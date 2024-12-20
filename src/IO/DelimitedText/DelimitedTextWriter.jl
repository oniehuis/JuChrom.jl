module DelimitedTextWriter

using CSV
using Unitful

import ...JuChrom: AbstractChrom, AbstractChromMS, scantimes, intensities, intensity, ions, ioncount, scancount
import ..InputOutput: FileExistsError, FileFormat, exportdata

export DelimitedText


struct DelimitedText{T<:Union{Char, AbstractString}} <: FileFormat
    delim::T
    function DelimitedText{T}(delim::T) where T<:Union{Char, AbstractString}
        new(delim)
    end
end


delim(fileformat::DelimitedText) = fileformat.delim


"""
    
    DelimitedText(; delim::Union{AbstractChar, AbstractString}="\\t") <: Fileformat

Return a `DelimitedText` file format object. The optional `delim` keyword argument allows 
you to specify the column delimiter, which can be either a single character or a string. 
If the provided filename lacks a suffix, an appropriate one will be automatically appended.

See also [`FileFormat`](@ref), [`exportdata`](@ref).

# Examples
```julia-repl
julia> DelimitedText()
DelimitedText{String}("\t")

julia> DelimitedText(delim=';')
DelimitedText{Char}(';')
```
"""
function DelimitedText(; delim::T="\t") where T<:Union{AbstractChar, AbstractString}
    DelimitedText{T}(delim)
end


function exportdata(fileformat::DelimitedText, chrom::AbstractChrom, file::AbstractString; 
    timeunit::Unitful.TimeUnits, overwrite::Bool)
    pathstem, suffix = splitext(file)
    if suffix == ""
        if string(delim(fileformat)) in (",", ";")
            file = string(pathstem, ".csv")
        elseif string(delim(fileformat)) == "\t"
            file = string(pathstem, ".tsv")
        else
            file = string(pathstem, ".txt")
        end
    end
    overwrite || isfile(file) && throw(FileExistsError(
        "a file with the same name already exists: \"$file\""))
    itr = ((scantimes=time, intensities=intensity(chrom, i)) 
        for (i, time) in enumerate(scantimes(chrom, timeunit=timeunit, ustripped=true)))
    CSV.write(file, itr, delim=fileformat.delim, 
        header=[string("scan time [", timeunit, "]"), "intensity"])
end


function exportdata(fileformat::DelimitedText, chrom::AbstractChromMS, file::AbstractString; 
    timeunit::Unitful.TimeUnits, overwrite::Bool)
    pathstem, suffix = splitext(file)
    if suffix == ""
        if string(delim(fileformat)) in (",", ";")
            file = string(pathstem, ".csv")
        elseif string(delim(fileformat)) == "\t"
            file = string(pathstem, ".tsv")
        else
            file = string(pathstem, ".txt")
        end
    end
    overwrite || isfile(file) && throw(FileExistsError(
        "a file with the same name already exists: \"$file\""))

    header = [string("m/z ", i) for i in ions(chrom)]
    pushfirst!(header, string("scan times [", timeunit, "]"))

    itr = Vector{NamedTuple}(undef, scancount(chrom))
    for (idxₛ, time) in enumerate(scantimes(chrom, timeunit=timeunit, ustripped=true))
        nt = (scantime=time, zip(Symbol.(ions(chrom)), intensities(chrom)[idxₛ, :])...)
        itr[idxₛ] = nt
    end

    CSV.write(file, itr, delim=fileformat.delim, header=header)
end


end  # module
