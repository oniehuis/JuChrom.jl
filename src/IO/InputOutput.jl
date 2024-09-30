module InputOutput

using Unitful
export buildxic
export Source
export File
export Path
export IOError
export FileExistsError

import ..JuChrom: AbstractChrom, scantimes

export FileFormat
export AgilentFID
export ANDI
export ChemStationMS
export DelimitedText
export Excel
export MassHunterMS
export ShimadzuMS
export exportdata
export importdata

include("./Utilities.jl")


"""
    FileFormat

Supertype of all FileFormat implementations.

See also [`AgilentFID`](@ref), [`ANDI`](@ref), [`ChemStationMS`](@ref), 
[`MassHunterMS`](@ref).
"""
abstract type FileFormat end


"""
    importdata(source::AbstractString, fileformat::FileFormat)

Read data from the specified source using the appropriate file format reader. The source 
may be a file or a directory path, depending on the requirements of the chosen file format 
reader.

See also [`AgilentFID`](@ref), [`ANDI`](@ref), [`ChemStationMS`](@ref), 
[`MassHunterMS`](@ref).

# Examples
```jldoctest
julia> dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D");

julia> chrom = importdata(dfolder, ChemStationMS())
ChromMS {scan times: Float32, ions: Float32, intensities: Int32}
2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms
5176 ions; range: m/z 29.0 - 562.9
intensity range: 0 - 1186816
metadata: 10 entries
```
"""
function importdata(source::AbstractString, fileformat::FileFormat)
    importdata(fileformat, source)
end


include(joinpath(".", "AgilentFID", "AgilentFIDReaders.jl"))
using .AgilentFIDReaders

include(joinpath(".", "ANDI", "ANDIReaders.jl"))
using .ANDIReaders

include(joinpath(".", "ChemStationMS", "ChemStationMSReaders.jl"))
using .ChemStationMSReaders

include(joinpath(".", "MassHunterMS", "MassHunterMSReaders.jl"))
using .MassHunterMSReaders

include(joinpath(".", "ShimadzuMS", "ShimadzuMSReaders.jl"))
using .ShimadzuMSReaders


"""
    exportdata(chrom::AbstractChrom, file::AbstractString, fileformat::FileFormat; 
    timeunit::Unitful.TimeUnits, overwrite::Bool=false)


Write the scan times and associated intensity values of the AbstractChrom object to a  
file in the specified format. The optional parameter `timeunit` allows you 
to specify the unit for the exported scan times. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `overwrite` allows to specify whether an existing 
target file should be overwritten.

See also [`AbstractChrom`](@ref), [`DelimitedText`](@ref), [`Excel`](@ref).

# Examples
```julia-repl
julia> chrom = Chrom((1:3)u"minute", [123, 224, 103])
Chrom {scan times: Int64, intensities: Int64}
3 scans; scan times: 1 minute, 2 minute, 3 minute
intensity range: 103 - 224
metadata: 0 entries

julia> exportdata(chrom, "./delimtest", DelimitedText());

julia> exportdata(chrom, "./delimtest", Excel(sheetname="TIC"));

```
"""
function exportdata(chrom::AbstractChrom, file::AbstractString, fileformat::FileFormat; 
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))), overwrite::Bool=false)
    exportdata(fileformat, chrom, file; timeunit, overwrite)
end


include(joinpath(".", "DelimitedText", "DelimitedTextWriter.jl"))
using .DelimitedTextWriter

include(joinpath(".", "Excel", "ExcelWriter.jl"))
using .ExcelWriter

end  # module
