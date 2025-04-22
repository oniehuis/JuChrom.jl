__precompile__()
module ShimadzuMSReaders

using PyCall
using Unitful

import ...JuChrom: ChromMS
import ..InputOutput: File, FileFormat, IOError, buildxic, importdata

const olefile = PyNULL()

function __init__()
    copy!(olefile, pyimport_conda("olefile", "olefile"))
end

export ShimadzuMS


"""
    
    ShimadzuMS()

Returns an object representing the `ShimadzuMS` file format. When using the importdata
function, the `ShimadzuMS` data reader expects the source to be a .qgd file.

See also [`FileFormat`](@ref), [`importdata`](@ref).

# Example
```jldoctest
julia> ShimadzuMS()
ShimadzuMS()

julia> qgdfile = joinpath(JuChrom.shimadzu, "AR190311.qgd");

julia> chrom = importdata(qgdfile, ShimadzuMS())
ChromMS {scan times: Int32, ions: Float64, intensities: UInt32}
9210 scans; scan time range: 183000 ms - 2024800 ms
3620 ions; range: m/z 34.9 - 600.4
intensity range: 0 - 924954
metadata: 0 entries
```
"""
struct ShimadzuMS <: FileFormat end 

function extractdata(bytesting)
    i = 1
    scantimes = Vector{Int32}()
    pointcounts = Vector{Int16}()
    mzvec = Vector{Float64}()
    intsvec = Vector{UInt32}()
    while i < length(bytesting)
        push!(scantimes, ltoh.(reinterpret(Int32, bytesting[i+4:i+7]))[1])
        bytecount = ltoh.(reinterpret(Int16, bytesting[i+20:i+21]))[1]
        @assert 1 ≤ bytecount ≤ 4 "unexpected byte count: $bytecount"
        pointcount = ltoh.(reinterpret(Int16, bytesting[i+22:i+23]))[1]
        push!(pointcounts, pointcount)
        start = i + 32
        stepwidth = 2 + bytecount
        stop = start + stepwidth * (pointcount - 1)
        scanmzs = Vector{Float64}(undef, pointcount)
        scanints = Vector{UInt32}(undef, pointcount)
        for (s, z) in enumerate(start:stepwidth:stop)
            scanmzs[s] = ltoh.(reinterpret(Int16, bytesting[z:z+1]))[1] / 20
            intbytes = zeros(UInt8, 4) 
            intbytes[1:bytecount] = bytesting[z+2:z+1+bytecount]
            int = ltoh.(reinterpret(UInt32, intbytes))[1]
            scanints[s] = int
        end
        append!(mzvec, scanmzs)
        append!(intsvec, scanints)
        i += 32 + (2 + bytecount) * pointcount
    end
    scantimes * u"ms", pointcounts, mzvec, intsvec
end

function readfile(::ShimadzuMS, file::AbstractString)
    # Read MS Raw Data from OLE file
    olefile.isOleFile(file) || throw(IOError("not an OLE file: $file"))
    ole::PyObject = olefile.OleFileIO(file)
    msrawdata::String = try
        stream = ["GCMS Raw Data", "MS Raw Data"] 
        ole.exists(stream)::Bool || throw(IOError("file does not contain stream $stream"))
        msrawdatastream::PyObject = ole.openstream(stream)
        msrawdatastream.read()
    finally
        ole.close()
    end

    # Transcode the data into a byte string
    bytesting::Base.CodeUnits{UInt8, String} = transcode(UInt8, msrawdata)

    # Extract data of interest
    scantimes, pointcounts, mzvec, intsvec = extractdata(bytesting)
    ions, xic = buildxic(pointcounts, mzvec, intsvec)
    ChromMS(scantimes, ions, xic)
end


function importdata(::File, fileformat::ShimadzuMS, file::AbstractString)
    readfile(ShimadzuMS(), file)
end


function importdata(fileformat::ShimadzuMS, source::AbstractString)
    if isfile(source)
        importdata(File(), fileformat, source)
    else
        throw(IOError("unsupported source for ShimadzuMS .qgd file data import"))
    end
end

end