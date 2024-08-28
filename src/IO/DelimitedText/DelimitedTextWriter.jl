module DelimitedTextWriter

using CSV
using Unitful

import ...JuMS: AbstractGC, scantimes, intensity
import ..InputOutput: FileFormat, exportdata

export DelimitedText

include("../Utilities.jl")

struct DelimitedText{T<:Union{Char, AbstractString}} <: FileFormat
    delim::T
    function DelimitedText{T}(delim::T) where T<:Union{Char, AbstractString}
        new(delim)
    end
end
# TYPE STABLE June 11, 2024


delim(fileformat::DelimitedText) = fileformat.delim
# TYPE STABLE June 11, 2024


"""
    
    DelimitedText(;delim::Union{Char, AbstractString}="\t") <: Fileformat

Return an DelimitedText file format object. The optional keyword argument `delim` allows you to specify the column delimiter as a character or a string.

# Example
```julia-repl
julia> DelimitedText()
DelimitedText{String}("\t")

julia> DelimitedText(delim=';')
DelimitedText{Char}(';')
```

See also [`Fileformat`](@ref), [`DelimitedText`](@ref), [`exportdata`](@ref).
"""
DelimitedText(;delim::T="\t") where T<:Union{Char, AbstractString} = DelimitedText{T}(delim)
# TYPE STABLE June 11, 2024


function exportdata(fileformat::DelimitedText, gc::AbstractGC, file::AbstractString, timeunit::Unitful.TimeUnits, overwrite::Bool)
    if string(delim(fileformat)) in (",", ";")
        file = string(first(splitext(file)), ".csv")
    elseif string(delim(fileformat)) == "\t"
        file = string(first(splitext(file)), ".tsv")
    else
        file = string(first(splitext(file)), ".txt")
    end
    overwrite || isfile(file) && throw(FileExistsError("a file with the same name already exists: \"$file\""))
    itr = ((scantimes=time, intensities=intensity(gc, i)) for (i, time) in enumerate(scantimes(gc, timeunit=timeunit, ustripped=true)))
    CSV.write(file, itr, delim=fileformat.delim, header=[string("scan time [", timeunit, "]"), "total intensity"])
end
# TYPE STABLE June 11, 2024


end  # module
