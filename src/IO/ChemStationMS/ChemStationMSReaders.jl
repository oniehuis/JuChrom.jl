module ChemStationMSReaders

using StringEncodings
using Unitful

import ...JuChrom: ChromMS
import ..InputOutput: FileFormat, File, Path, IOError, buildxic, importdata

export ChemStationMS

const encoding = "Windows-1252"

struct ChemStationMS{T<:AbstractString} <: FileFormat
    datafilename::T
    ChemStationMS{T}(datafilename::T) where {T<:AbstractString} = new(datafilename)
end


"""
    
    ChemStationMS(; datafilename::AbstractString="data.ms")

Returns an object representing the `ChemStationMS` file format. The optional keyword 
argument, `datafilename`, allows you to specify a different name for the data file to be 
read, instead of the default `data.ms` (e.g., `datasim.ms`). Upper and lower case are not 
considered when specifying the data file name. Note that when using the function 
importdata, the `ChemStationMS` data reader expects the Agilent .D folder as the source 
location, with the specified data file located at the top level of the .D folder. 
However, you can also provide the full pathname of the data file directly when using 
the function importdata, in which case the `datafilename` entry in the `ChemStationMS` 
object is ignored.

See also [`FileFormat`](@ref), [`importdata`](@ref).

# Examples
```jldoctest
julia> ChemStationMS()
ChemStationMS{String}("data.ms")

julia> ChemStationMS(datafilename="DATASIM.MS")
ChemStationMS{String}("DATASIM.MS")

julia> dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D");

julia> chrom = importdata(dfolder, ChemStationMS())
ChromMS {scan times: Float32, ions: Float32, intensities: Int32}
2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms
5176 ions; range: m/z 29.0 - 562.9
intensity range: 0 - 1186816
metadata: 10 entries

julia> datafile = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D/data.ms");

julia> chrom = importdata(datafile, ChemStationMS())
ChromMS {scan times: Float32, ions: Float32, intensities: Int32}
2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms
5176 ions; range: m/z 29.0 - 562.9
intensity range: 0 - 1186816
metadata: 10 entries
```
"""
function ChemStationMS(; datafilename::T="data.ms") where {T<:AbstractString}
    ChemStationMS{T}(datafilename)
end


datafilename(fileformat::ChemStationMS) = fileformat.datafilename


function readmetadata(f::IOStream, positionof::Dict{Symbol, Int}, stepwidth::Integer)
    metadata = Dict{Any, Any}()
    for (feature, pos) in pairs(positionof)
        seek(f, pos)
        bytes = ltoh.(read(f, ltoh(read(f, UInt8)) * stepwidth))[begin:stepwidth:end]
        value = decode(bytes, encoding)
        metadata[feature] = value ≠ "" ? value : nothing
    end
    metadata
end


struct ChemStationMSV1 end


# Tested with data file from MSD ChemStation E.02.02.1431 1989-2011
# Tested with data file from MSD ChemStation F.01.03.2357 1989-2015
function readfile(::ChemStationMSV1, file::AbstractString)

    # Location of string type metadata
    positionof = Dict{Symbol, Int}(
        :sample      =>  24,
        :description =>  86,
        :operator    => 148,
        :datetime    => 178,
        :type        => 208,
        :inlet       => 218,
        :method      => 228)

    open(file, "r") do f

        # File version
        
        ver = decode(ltoh.(read(f, ltoh(read(f, UInt8)))), encoding)
        ver == "2" || throw(
            FileFormatError("cannot read file with a version \"$ver\": \"$file\""))

        # File type
        seek(f, 4)
        
        filetype = decode(ltoh.(read(f, ltoh(read(f, UInt8)))), encoding)
        (filetype == "GC / MS Data File" || filetype == "GC / MS DATA FILE") || throw(
            FileFormatError(string("cannot read file with the file type signature ",
                "\"$filetype\": \"$file\"")))

        # Extract string type metadata
        metadata = readmetadata(f, positionof, 1)

        # Extract integer type metadata
        seek(f, 252)
        metadata[:sequence] = convert(Int, ntoh(read(f, Int16)))
        metadata[:vial] = convert(Int, ntoh(read(f, Int16)))
        metadata[:replicate] = convert(Int, ntoh(read(f, Int16)))

        # Scan count
        seek(f, 278)
        scancount = ntoh(read(f, UInt32))

        # Move to scan data
        seek(f, 266)
        seek(f, 2 * ntoh(read(f, UInt16)) - 2)

        # Reading scan data
        scantimes = Vector{Float32}(undef, scancount)
        pointcounts = Vector{Int32}(undef, scancount) # Was Int! better concret? Int32 or Int64
        ionvec = Vector{Float32}()
        intsvec = Vector{Int32}()  # Was Int! Better concret? Int32 or Int64?

        for i = 1:scancount

            # Store next scan start position
            nextscanpos = position(f) + 2 * ntoh(read(f, UInt16))

            # Get scan time
            scantimes[i] = ntoh(read(f, Int32))

            # Get number of ion-intensity raw data pairs (= point counts) of scan
            pointcounts[i] = fld((ntoh(read(f, UInt16)) - 6), 2)
            
            # Get ion-intensity raw data pairs of scan
            skip(f, 10)
            rawdatapairs = ntoh.(read!(f, Vector{UInt16}(undef, 2 * pointcounts[i])))

            # Store raw ion values
            append!(ionvec, @view rawdatapairs[1:2:end])

            # Store raw intensity values
            append!(intsvec, @view rawdatapairs[2:2:end])

            # move to next scan position
            seek(f, nextscanpos)
        end

        # Calculate mass values
        @. ionvec /= 20

        # Calculate intensity values
        @. intsvec = (intsvec & 16383) * 8^(intsvec >> 14)
        
        ChromMS(scantimes * 1u"ms", buildxic(pointcounts, ionvec, intsvec)..., metadata)
    end
end


function version(::ChemStationMS, file::AbstractString)
    open(file, "r") do f
        ver = String(ltoh.(read(f, ltoh(read(f, UInt8)))))  # position 1 (= 0x1)
        seek(f, 4)
        ftype = String(ltoh.(read(f, ltoh(read(f, UInt8)))))
        return ver, ftype
    end
end


function readfile(::ChemStationMS, file::AbstractString)
    ver, ftype = version(ChemStationMS(), file)
    ver == "2" || throw("cannot read .ms files with a version \"$ver\": \"$file\"")
    readfile(ChemStationMSV1(), file)
end


function importdata(::Path, fileformat::ChemStationMS, path::AbstractString)
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


function importdata(::File, fileformat::ChemStationMS, file::AbstractString)
    readfile(fileformat, file)
end


function importdata(fileformat::ChemStationMS, source::AbstractString)
    if isdir(source)
        importdata(Path(), fileformat, source)
    elseif isfile(source)
        importdata(File(), fileformat, source)
    else
        throw(IOError("unsupported source for ChemStationMS data import"))
    end
end

end # module