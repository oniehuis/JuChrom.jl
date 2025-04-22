module AgilentFIDReaders

using StringEncodings
using Unitful

import ...JuChrom: Chrom
import ..InputOutput: FileFormat, File, Path, IOError, importdata

export AgilentFID

const encoding = "Windows-1252"

struct AgilentFID{T<:AbstractString} <: FileFormat
    datafilename::T
    AgilentFID{T}(datafilename::T) where {T<:AbstractString} = new(datafilename)
end


"""
    
    AgilentFID(; datafilename::AbstractString="FID1A.ch")

Returns an object representing the `AgilentFID` file format. The optional keyword argument, 
`datafilename`, allows you to specify a different name for the data file to be read, 
instead of the default `FID1A.ch`. Upper and lower case are not considered when specifying 
the data file name. Note that when using the function importdata, the `AgilentFID` data 
reader expects the Agilent .D folder as the source location, with the specified data file 
located at the top level of the .D folder. However, you can also provide the full pathname 
of the data file directly when using the function importdata, in which case the 
`datafilename` entry in the `AgilentFID` object is ignored.

See also [`FileFormat`](@ref), [`importdata`](@ref).

# Examples
```jldoctest
julia> AgilentFID()
AgilentFID{String}("FID1A.ch")

julia> AgilentFID(datafilename="fid_data_20240611.ch")
AgilentFID{String}("fid_data_20240611.ch")

julia> dfolder = joinpath(JuChrom.agilent, "ZK_ONUBE_Mix1_11.D");

julia> chrom = importdata(dfolder, AgilentFID())
Chrom {scan times: Float32, intensities: Float64}
4151 scans; scan time range: 0.437f0 ms - 830000.44f0 ms
intensity range: 0.0 - 1.0738316309895832e6
metadata: 10 entries
```
"""
function AgilentFID(; datafilename::T="FID1A.ch") where {T<:AbstractString}
    AgilentFID{T}(datafilename)
end



datafilename(fileformat::AgilentFID) = fileformat.datafilename


struct AgilentFIDv179 end


function readmetadata(f::IOStream, positionof::Dict{Symbol, Int}, stepwidth::Integer)
    metadata = Dict{Symbol, Union{String, Nothing}}()
    for (feature, pos) in pairs(positionof)
        seek(f, pos)
        bytes = ltoh.(read(f, ltoh(read(f, UInt8)) * stepwidth))[begin:stepwidth:end]
        value = decode(bytes, encoding)
        metadata[feature] = value ≠ "" ? value : nothing
    end
    metadata
end


function readfile(::AgilentFIDv179, file::AbstractString)

    # Location of string type metadata
    positionof = Dict{Symbol, Int}(
        :sample      =>  858,
        :description => 1369,
        :operator    => 1880,
        :datetime    => 2391,
        :inlet       => 2492,
        :type        => 2533,
        :method      => 2574,
        :software    => 3089,
        :signal      => 4213)

    open(file, "r") do f

        # Scan count
        seek(f, 278)  # 0x116
        scancount = ntoh(read(f, UInt32))
        scancount ≠ 0 || throw(FileFormatError("file has a scan count of zero: \"$file\""))

        # Start and stop time
        seek(f, 282)  # 0x11A
        start, stop = ntoh.(read!(f, Vector{Float32}(undef, 2)))
        scantimes = LinRange(start, stop, scancount) * 1u"ms"

        seek(f, 347)
        filetype = decode(ltoh.(read(f, ltoh(read(f, UInt8)) * 2))[begin:2:end], encoding)
        (filetype == "GC DATA FILE") || throw(
            FileFormatError(string("cannot read file with the file type signature ",
                "\"$filetype\": \"$file\"")))

        # Metadata
        metadata = readmetadata(f, positionof, 2)

        # Signal unit
        seek(f, 4172)
        signalunit = decode(ltoh.(read(f, ltoh(read(f, UInt8)) * 2))[begin:2:end], encoding)
        metadata[:signal_unit] = signalunit

        # Signal scaling factor
        seek(f, 4732)
        scalingfactor = ntoh(read(f, Float64))

        # Signals
        seek(f, 6144)
        signals = ltoh.(read!(f, Vector{Float64}(undef, scancount)))

        intensities = [abs(signal * scalingfactor) for signal in signals]

        Chrom(scantimes, intensities, metadata)
    end
end


function version(::AgilentFID, file::AbstractString)
    open(file, "r") do f
        return decode(ltoh.(read(f, ltoh(read(f, UInt8)))), encoding)
    end
end


function readfile(::AgilentFID, file::AbstractString)
    ver = version(AgilentFID(), file)
    ver == "179" || throw(FileFormatError(
        "cannot read file with a version \"$ver\": \"$file\""))
    readfile(AgilentFIDv179(), file)
end


function importdata(::Path, fileformat::AgilentFID, path::AbstractString)
    splitext(path)[end] == ".D" || throw(IOError("not an Agilent .D folder: \"$path\""))
    files = readdir(path)
    datafile = datafilename(fileformat)
    idcs = findall(file -> lowercase(file) == lowercase.(datafile), files)
    length(idcs) ≤ 1 || throw(IOError(string("found multiple data files in path ",
        "\"$path\": ", files[idcs])))
    if length(idcs) == 1
        filename = files[idcs[1]]
        importdata(File(), fileformat, joinpath(path, filename))
    else
        throw(IOError("cannot find specified data file \"$datafile\" in path \"$path\""))       
    end
end


function importdata(::File, fileformat::AgilentFID, file::AbstractString)
    readfile(fileformat, file)
end


function importdata(fileformat::AgilentFID, source::AbstractString)
    if isdir(source)
        importdata(Path(), fileformat, source)
    elseif isfile(source)
        importdata(File(), fileformat, source)
    else
        throw(IOError("unsupported source for AgilentFID data import: \"$source\""))
    end
end

end  # module
