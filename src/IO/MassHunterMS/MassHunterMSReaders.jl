module MassHunterMSReaders

using Unitful

import ...JuChrom: GCMS
import ..InputOutput: FileFormat, Path, IOError, buildxic, importdata

export MassHunterMS


struct MassHunterMS{T1<:Integer, T2<:Type, T3<:Type, T4<:Type} <: FileFormat
    scanmethodid::T1
    scantimetype::T2
    iontype::T3
    intensitytype::T4
    function MassHunterMS{T1, T2, T3, T4}(scanmethodid::T1, scantimetype::T2, 
        iontype::T3, intensitytype::T4) where {T1<:Integer, T2<:Type, T3<:Type, 
        T4<:Type}
        new(scanmethodid, scantimetype, iontype, intensitytype)
    end
end


"""
    
    MassHunterMS(; scanmethodid::Integer=1, scantimetype::Union{Float32, Float64}=Float32, 
    iontype::Union{Float32, Float64}=Float32, intensitytype::Union{Float32, Float64}=Float32)

Return an object representing the `MassHunterMS` file format. The optional `scanmethodid` 
keyword argument allows you to specify which scan method ID data to read, which is 
particularly useful when multiple scan methods have been applied simultaneously. The 
:ScanMethodID and :ScanMethodIDs metadata entries indicate which scan method ID data was 
extracted and which scan method IDs are available for a given run. The optional 
`scantimetype` keyword argument allows the scan time to be converted to the specified float 
type. Similarly, the optional `iontype` keyword argument allows the mass-to-charge ratio 
values to be converted to the desired float type. The optional `intensitytype` keyword 
argument allows the intensity values to be converted to the specified float type. The 
keyword arguments `scantimetype`, `iontype`, and `intensitytype` are provided to address 
the issue that MassHunter files may store data as Float64 when it was originally collected 
as Float32. The storage format (Float64 or Float32) of the ion and intensity values can be 
determined by checking the bytes per point value in the metadata: a value of 16 corresponds 
to Float64, while a value of 8 corresponds to Float32. Scan times are always stored as 
64-bit floats. Note that when using the importdata function, the MassHunterMS data reader 
expects the source location to be the Agilent .D folder.

See also [`FileFormat`](@ref), [`importdata`](@ref), [`JuChrom.binions`](@ref).

# Examples
```jldoctest
julia> dfolder = joinpath(JuChrom.agilent, "C7-C40_MassHunterMS.D");

julia> gcms = importdata(dfolder, MassHunterMS())
GCMS {scan times: Float32, ions: Float32, intensities: Float32}
2405 scans; scan time range: 3.1990166f0 minute - 31.650784f0 minute
50275 ions; range: m/z 29.02 - 562.89
intensity range: 0.0 - 1.1872475e6
metadata: 3 entries

julia> ions(gcms)[1:6]  # ions are converted per default to 32-bit floats
6-element Vector{Float32}:
 29.02
 29.03
 29.04
 29.05
 29.06
 29.07

julia> metadata(gcms)  # ion and intensity values are stored as 64-bit floats
Dict{Any, Any} with 3 entries:
  :ScanMethodID  => 1
  :BytesPerPoint => 16.0
  :ScanMethodIDs => [1]

julia> gcms = importdata(dfolder, MassHunterMS(iontype=Float64))  # convert ions to Float64
GCMS {scan times: Float32, ions: Float64, intensities: Float32}
2405 scans; scan time range: 3.1990166f0 minute - 31.650784f0 minute
50275 ions; range: m/z 29.020000457763672 - 562.8900146484375
intensity range: 0.0 - 1.1872475e6
metadata: 3 entries

julia> ions(gcms)[1:6]  # 64-bit floats do not display the ion values well
6-element Vector{Float64}:
 29.020000457763672
 29.030000686645508
 29.040000915527344
 29.049999237060547
 29.059999465942383
 29.06999969482422
```
"""
function MassHunterMS(; scanmethodid::T1=1, scantimetype::T2=Float32, iontype::T3=Float32, 
    intensitytype::T4=Float32) where {T1<:Integer, T2<:Type, T3<:Type, T4<:Type}
    MassHunterMS{T1, T2, T3, T4}(scanmethodid, scantimetype, iontype, intensitytype)
end


scanmethodid(fileformat::MassHunterMS) = fileformat.scanmethodid
scantimetype(fileformat::MassHunterMS) = fileformat.scantimetype
iontype(fileformat::MassHunterMS) = fileformat.iontype
intensitytype(fileformat::MassHunterMS) = fileformat.intensitytype


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
    scandata::Vector{MSScanDataV1}, ::Type{T}, ::Type{T2}, ::Type{T3}, ::Type{T4}, 
    bpp, scanmethodids) where {T, T2, T3, T4}
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
    mzvec_converted = convert(Vector{iontype(fileformat)}, mzvec)
    intsvec_converted = convert(Vector{intensitytype(fileformat)}, intsvec)
    ions, xic = buildxic([scandata[i].PointCount for i in 1:length(scandata)], 
        mzvec_converted, intsvec_converted)
    scantimes_converted_unitful = convert(Vector{scantimetype(fileformat)}, 
        scantimes(scandata)) * 1u"minute"
    GCMS(scantimes_converted_unitful, ions, xic, Dict(:BytesPerPoint => bpp, 
    :ScanMethodID => scanmethodid(fileformat), :ScanMethodIDs => scanmethodids))
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
    scanmethodids = sort(collect(Set(scanmethodid(fileformat))))
    methodscandata = length(scanmethodids) == 1 ? scandata : filter(
        s -> s.ScanMethodID .== scanmethodid(fileformat), scandata)
    bpp = bytesperpoint(methodscandata)
    bpp == 8 || bpp == 16 || throw(
        IOError("\"$file\" has unexpected ByteCount: \"$bpp\""))
    storagetype = bpp == 8 ? Float32 : Float64
    mspeakdata(MSPeakBinV1(), file, fileformat, methodscandata, storagetype, 
        scantimetype(fileformat), iontype(fileformat), intensitytype(fileformat), bpp, 
        scanmethodids)
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