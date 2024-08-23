module InputOutput

using Unitful

# import ..JuChrom: AbstractGC, scantimes

export FileFormat
#export AgilentFID
#export ANDI
export ChemStationMS
#export DelimitedText
#export Excel
#export MassHunterMS
#export exportdata
export importdata


"""
    FileFormat

Supertype of all FileFormat implementations.

See also [`ChemStationMS`](@ref).
"""
abstract type FileFormat end
# [`AgilentFID`](@ref), [`ANDI`](@ref), , [`MassHunterMS`](@ref)


"""
    importdata(source::AbstractString, fileformat::FileFormat)

Read data from the specified source using the appropriate file format reader. The source 
may be a file or a directory path, depending on the requirements of the chosen file format 
reader.

# Example
```jldoctest
julia> examplefile = joinpath(JuChrom.exampledata, "C7-C40_ChemStationMS.D");

julia> gcms = importdata(examplefile, ChemStationMS())
GCMS {scan times: Float32, ions: Float32, intensities: Int64}
2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms
5176 ions; range: m/z 29.0 - 562.9
intensity range: 0 - 1186816
metadata: 0 entries
```
"""
function importdata(source::AbstractString, fileformat::FileFormat)
    importdata(fileformat, source)
end
# See also [`AgilentFID`](@ref), [`ANDI`](@ref), [`ChemStationMS`](@ref), [`MassHunterMS`](@ref), [`GCMS`](@ref), [`FID`](@ref), .
# julia> gcms = importdata(examplefile, ChemStationMS());
# julia> examplefile = joinpath(exampledata, "C7-C40_ChemStationMS.D");
#include(joinpath(".", "AgilentFID", "AgilentFIDReaders.jl"))
#using .AgilentFIDReaders

include(joinpath(".", "ChemStationMS", "ChemStationMSReaders.jl"))
using .ChemStationMSReaders

#include(joinpath(".", "MassHunterMS", "MassHunterMSReaders.jl"))
#using .MassHunterMSReaders

#include(joinpath(".", "ANDI", "ANDIReaders.jl"))
#using .ANDIReaders


# """
#     exportdata(gc::AbstractGC, file::AbstractString, fileformat::FileFormat; timeunit::Unitful.TimeUnits, overwrite::Bool=false)

# Write the scan times and associated intensity values of the AbstractGC object to the file in the specified file format. The optional keyword argument
# `timeunit` allows the scan times to be written in a specific time unit. All time units defined in the package [Unitful.jl](https://painterqubits.github.io/Unitful.jl)
# (e.g., u"s", u"minute") are supported. The optional keyword argument `overwrite` allows to specify that an existing target file should be overwritten.

# # Example
# ```julia-repl
# julia> tic = TIC([1u"s", 2u"s", 3u"s"], [123, 224, 103])

# julia> exportdata(tic, "/home/JuliaUser/tic", Excel())
# /home/JuliaUser/tic.xlsx

# julia> exportdata(tic, "/home/JuliaUser/tic", DelimitedText())
# /home/JuliaUser/tic.tsv"

# julia> exportdata(tic, "/home/JuliaUser/tic", DelimitedText(), timeunit=u"minute", overwrite=true)
# /home/JuliaUser/tic.tsv"
# ```

# See also [`AbstractGC`](@ref), [`DelimitedText`](@ref), [`Excel`](@ref).
# """
# function exportdata(gc::AbstractGC, file::AbstractString, fileformat::FileFormat; timeunit::Unitful.TimeUnits=unit(eltype(scantimes(gc))), overwrite::Bool=false)
#     exportdata(fileformat, gc, file; timeunit, overwrite)
# end


# include(joinpath(".", "DelimitedText", "DelimitedTextWriter.jl"))
# using .DelimitedTextWriter

# include(joinpath(".", "Excel", "ExcelWriter.jl"))
# using .ExcelWriter

end  # module
