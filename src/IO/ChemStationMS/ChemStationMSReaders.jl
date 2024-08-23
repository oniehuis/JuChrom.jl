module ChemStationMSReaders

using Unitful

import ...JuChrom: GCMS
import ..InputOutput: FileFormat, importdata 
import ..InputOutput: Path, IOError, FileExistsError, buildxic, File 

export ChemStationMS


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

# Example
```jldoctest
julia> ChemStationMS()
ChemStationMS{String}("data.ms")

julia> ChemStationMS(datafilename="DATASIM.MS")
ChemStationMS{String}("DATASIM.MS")

julia> dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D");

julia> gcms = importdata(dfolder, ChemStationMS())
GCMS {scan times: Float32, ions: Float32, intensities: Int64}
2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms
5176 ions; range: m/z 29.0 - 562.9
intensity range: 0 - 1186816
metadata: 0 entries

julia> datafile = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D/data.ms");

julia> gcms = importdata(datafile, ChemStationMS())
GCMS {scan times: Float32, ions: Float32, intensities: Int64}
2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms
5176 ions; range: m/z 29.0 - 562.9
intensity range: 0 - 1186816
metadata: 0 entries
```
"""
function ChemStationMS(; datafilename::T="data.ms") where {T<:AbstractString}
    ChemStationMS{T}(datafilename)
end


datafilename(fileformat::ChemStationMS) = fileformat.datafilename


struct ChemStationMSV1 end

# Reads data from MSD ChemStation F.01.03.2357 1989-2015 file provided by Wolf Haberer (Freiburg).
# Reads data from MSD ChemStation E.02.02.1431 1989-2011 file provided by Florian Menzel (Mainz).
# Reads data from Chemstation VERSION? file provided by Thomas Schmitt (Würzburg).
function readfile(::ChemStationMSV1, file::AbstractString)
    open(file, "r") do f

        # File version
        ver = String(ltoh.(read(f, ltoh(read(f, UInt8)))))  # position 1 (= 0x1)
        ver == "2" || throw(
                FileFormatError("cannot read file with a version \"$ver\": \"$file\""))

        # File type
        seek(f, 4)
        filetype = String(ltoh.(read(f, ltoh(read(f, UInt8)))))
        (filetype == "GC / MS Data File" || filetype == "GC / MS DATA FILE") || throw(
            FileFormatError(string("cannot read file with the file type signature ",
                "\"$filetype\": \"$file\"")))

        # Scan count
        seek(f, 278)
        scancount = ntoh(read(f, UInt32))

        # Move to scan data
        seek(f, 266)
        seek(f, 2 * ntoh(read(f, UInt16)) - 2)

        # Reading scan data
        scantimes = Vector{Float32}(undef, scancount)
        pointcounts = Vector{Int}(undef, scancount)
        ionvec = Vector{Float32}()
        intsvec = Vector{Int}()

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
        
        ions, xic = buildxic(pointcounts, ionvec, intsvec)

        GCMS(scantimes * 1u"ms", ions, xic, Dict())
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