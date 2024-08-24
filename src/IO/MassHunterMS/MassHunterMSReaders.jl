module MassHunterMSReaders

using Unitful

import ...JuChrom: GCMS
import ..InputOutput: FileFormat, Path, IOError, buildxic, importdata

export MassHunterMS


struct MassHunterMS{T1<:Integer, T2<:Integer, T3<:Integer, T4<:Integer} <: FileFormat
    scanmethodid::T1
    scantimedigits::T2
    mzdigits::T3
    intensitydigits::T4
    function MassHunterMS{T1, T2, T3, T4}(scanmethodid::T1, scantimedigits::T2, 
        mzdigits::T3, intensitydigits::T4) where {T1<:Integer, T2<:Integer, T3<:Integer, 
        T4<:Integer}
        new(scanmethodid, scantimedigits, mzdigits, intensitydigits)
    end
end


"""
    
    MassHunterMS(; scanmethodid::Integer=1, scantimedigits::Integer=1, 
    iondigits::Integer=1, intensitydigits::Integer=0)

Return an object representing the `MassHunterMS` file format. The optional `scanmethodid` 
keyword argument allows you to specify which scan method ID's data to read, especially 
when multiple scan methods have been applied simultaneously. The optional `scantimedigits` 
keyword argument allows you to round the scan time to a specified number of decimal digits. 
Similarly, the optional `iondigits` keyword argument rounds the mass-to-charge ratio values 
to the desired number of decimal places. Note that this is different from ion binning. The 
optional `intensitydigits` keyword argument rounds the intensity values to the specified 
number of decimals. If zero decimal places are specified, the values are returned as 
integers. The keyword arguments `scantimedigits`, `mzdigits`, and `intensitydigits` address 
the inconvenience that MassHunter files may store data in types different than the original 
data collected. A simple type conversion would be insufficient due to inaccuracies 
introduced during the data conversion before storing the data in MassHunter files. Note 
that when using the function importdata, the `MassHunterMS` data reader expects the Agilent 
.D folder as the source location.

See also [`FileFormat`](@ref), [`importdata`](@ref), [`JuChrom.binions`](@ref).

# Examples
```jldoctest
julia> dfolder = joinpath(JuChrom.agilent, "C7-C40_MassHunterMS.D");

julia> gcms = importdata(dfolder, MassHunterMS())
GCMS {scan times: Float64, ions: Float64, intensities: Int64}
2405 scans; scan time range: 3.199 minute - 31.651 minute
5174 ions; range: m/z 29.0 - 562.9
intensity range: 0 - 1187248
metadata: 0 entries

julia> gcms = importdata(dfolder, MassHunterMS(intensitydigits=1))
GCMS {scan times: Float64, ions: Float64, intensities: Float64}
2405 scans; scan time range: 3.199 minute - 31.651 minute
5174 ions; range: m/z 29.0 - 562.9
intensity range: 0.0 - 1.1872475e6
metadata: 0 entries
```
"""
function MassHunterMS(; scanmethodid::T1=1, scantimedigits::T2=3, mzdigits::T3=1, 
    intensitydigits::T4=0) where {T1<:Integer, T2<:Integer, T3<:Integer, T4<:Integer}
    MassHunterMS{T1, T2, T3, T4}(scanmethodid, scantimedigits, mzdigits, intensitydigits)
end


scanmethodid(fileformat::MassHunterMS) = fileformat.scanmethodid
scantimedigits(fileformat::MassHunterMS) = fileformat.scantimedigits
mzdigits(fileformat::MassHunterMS) = fileformat.mzdigits
intensitydigits(fileformat::MassHunterMS) = fileformat.intensitydigits


struct MSScanBin end
struct MSScanBinV1 end
struct MSPeakBin end
struct MSPeakBinV1 end


struct MSScanDataV1
    ScanID::Int32
    ScanMethodID::Int32
    TimeSegmentID::Int32
    ScanTime::Float64
    MSLevel::Int32
    ScanType::Int32
    TIC::Float64
    BasePeakMZ::Float64
    BasePeakValue::Float64
    CycleNumber::Int32
    Status::Int32
    IonMode::Int32
    IonPolarity::Int32
    Fragmentor::Float64
    CollisionEnergy::Float64
    MzOfInterest::Float64
    SamplingPeriod::Float64
    MeasuredMassRangeMin::Float64
    MeasuredMassRangeMax::Float64
    Threshold::Float64
    IsFragmentorDynamic::Int32
    IsCollisionEnergyDynamic::Int32
    SpectrumFormatID::Int32
    SpectrumOffset::Int64
    ByteCount::Int32
    PointCount::Int32
    MinY::Float64
    MaxY::Float64
    MinX::Float64
    MaxX::Float64
end


for type in (MSScanDataV1, )
    expr = Expr(:call, type, map(ft -> Expr(:call, :ltoh, Expr(:call, :read, :io, ft)), 
        fieldtypes(type))...)
    @eval Base.read(io::IO, ::Type{$type}) = $expr
end


function msscandata(::MSScanBinV1, file::AbstractString)
    size = sum(sizeof(fieldtype) for fieldtype in fieldtypes(MSScanDataV1))
    open(file, "r") do io
        seekend(io)
        endpos = position(io)
        seek(io, 88)
        startpos = ltoh(read(io, UInt16))
        seek(io, startpos)
        (endpos - startpos) % size == 0 || throw(
            IOError("MSScan.bin has unexpected data content"))
        scancount = div(endpos - startpos, size)
        data = Vector{MSScanDataV1}(undef, scancount)
        for i in 1:scancount
            data[i] = read(io, MSScanDataV1)
        end
        data
    end
end


pointcountsum(scandata::Vector{MSScanDataV1}) = sum(s.PointCount for s in scandata)
bytecountsum(scandata::Vector{MSScanDataV1}) = sum(s.ByteCount for s in scandata)
scantimes(scandata::Vector{MSScanDataV1}) = [s.ScanTime for s in scandata]
bytesperpoint(scandata::Vector{MSScanDataV1}) = (bytecountsum(scandata) 
    / pointcountsum(scandata))
scanmethodids(scandata::Vector{MSScanDataV1}) = [s.ScanMethodID for s in scandata]


function mspeakdata(::MSPeakBinV1, file::AbstractString, fileformat::FileFormat, 
    scandata::Vector{MSScanDataV1}, ::Type{T}) where {T}
    mzvec = Vector{T}(undef, pointcountsum(scandata))
    intsvec = Vector{T}(undef, pointcountsum(scandata))
    open(file, "r") do f
        i = 1
        for s in scandata
            seek(f, s.SpectrumOffset)
            mzvec[i:i+s.PointCount-1] = ltoh.(read!(f, Vector{T}(undef, s.PointCount)))
            intsvec[i:i+s.PointCount-1] = ltoh.(read!(f, Vector{T}(undef, s.PointCount)))
            i += s.PointCount
        end
    end
    mzvec_rounded = round.(mzvec, digits=mzdigits(fileformat))
    intsvec_rounded = round.(intsvec, digits=intensitydigits(fileformat))
    if intensitydigits(fileformat) == 0
        T2 = T == Float32 ? Int32 : Int64
        intsvec_rounded = convert(Vector{T2}, intsvec_rounded)
    end
    ions, xic = buildxic([scandata[i].PointCount for i in 1:length(scandata)], 
        mzvec_rounded, intsvec_rounded)
    scantimes_rounded_unitful = round.(scantimes(scandata), 
        digits=scantimedigits(fileformat)) * 1u"minute"
    GCMS(scantimes_rounded_unitful, ions, xic, Dict())
end


function magicnumber(file::AbstractString)
    open(file, "r") do io
        return ltoh(read(io, UInt16))
    end
end


function readfile(::MSPeakBinV1, path::AbstractString, fileformat::MassHunterMS, 
        scandata::Vector{MSScanDataV1})
    file = joinpath(path, "MSPeak.bin")
    signature = magicnumber(file)
    signature == 259 || throw(
        IOError("\"$file\" has unexpected file signature: \"$signature\""))
    methodscandata = filter(s -> s.ScanMethodID .== scanmethodid(fileformat), scandata)
    bpp = bytesperpoint(methodscandata)
    bpp == 8 || bpp == 16 || throw(
        IOError("\"$file\" has unexpected ByteCount: \"$bpp\""))
    type = bpp == 8 ? Float32 : Float64
    mspeakdata(MSPeakBinV1(), file, fileformat, methodscandata, type)
end


function readfile(::MSScanBinV1, path::AbstractString, 
    fileformat::MassHunterMS)::Vector{MSScanDataV1}
    file = joinpath(path, "MSScan.bin")
    signature = magicnumber(file)
    signature == 257 || throw(
        IOError("\"$file\" has unexpected file signature: \"$signature\""))
    scandata = msscandata(MSScanBinV1(), file)
    methodid = scanmethodid(fileformat)
    methodid in scanmethodids(scandata) || throw(
        IOError("\"$file\" contains no data for ScanMethodID \"$methodid\""))
    scandata
end


function readfiles(::MSScanBin, ::MSPeakBin, path::AbstractString, 
    fileformat::MassHunterMS)
    scandata = readfile(MSScanBinV1(), path, fileformat)
    readfile(MSPeakBinV1(), path, fileformat, scandata)
end


function importdata(::Path, path::AbstractString, fileformat::MassHunterMS)
    splitext(path)[end] == ".D" || throw(
        IOError("not an Agilent .D folder: \"$path\""))
    acqdatapath = joinpath(path, "AcqData")
    isdir(acqdatapath) || throw(
        IOError("cannot find folder AcqData in Agilent .D folder"))
    files = readdir(acqdatapath)
    if "MSScan.bin" ∉ files
        throw(IOError("cannot find \"MSScan.bin\" file in path: \"$path\""))
    elseif "MSPeak.bin" ∉ files
        throw(IOError("cannot find \"MSPeak.bin\" file in path: \"$path\""))
    else
        readfiles(MSScanBin(), MSPeakBin(), acqdatapath, fileformat)
    end
end


function importdata(fileformat::MassHunterMS, source::AbstractString)
    if isdir(source)
        importdata(Path(), source, fileformat)
    else
        throw(IOError("unsupported source for MassHunterMS data import"))
    end
end

end # module