module DelimitedTextWriter

using CSV
using Unitful

import ...JuChrom: AbstractChrom, AbstractChromMS, scantimes, intensities, intensity, ions, 
    ioncount, rimapper, retentionindexname, minretentiontime, maxretentiontime, 
    extrapolationmethod, retentionindex, scancount
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

    if isnothing(rimapper(chrom))
        header = [string("scan time [", timeunit, "]"), "intensity"]
        itr = ((scantimes=time, intensities=intensity(chrom, i)) 
            for (i, time) in enumerate(scantimes(chrom, timeunit=timeunit, ustripped=true)))
    else
        ris = Vector{Union{Float64, Missing}}(missing, scancount(chrom))
        ld = rimapper(chrom)
        for (i, rt) in enumerate(scantimes(chrom))
            if (minretentiontime(ld) ≤ rt ≤ maxretentiontime(ld) 
                || !isnothing(extrapolationmethod(ld)))
                ris[i] = retentionindex(ld, rt)
            end
        end
        header = [string("scan time [", timeunit, "]"), retentionindexname(ld), 
            "intensity"]
        itr = ((scantimes=time, retentionindex=ris[i], intensities=intensity(chrom, i)) 
            for (i, time) in enumerate(scantimes(chrom, timeunit=timeunit, ustripped=true)))
    end
    
    CSV.write(file, itr, delim=fileformat.delim, header=header)
    
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

    itr = Vector{NamedTuple}(undef, scancount(chrom))

    if isnothing(rimapper(chrom))
        header = vcat(string("scan times [", timeunit, "]"), 
            [string("m/z ", i) for i in ions(chrom)])

        for (idxₛ, time) in enumerate(scantimes(chrom, timeunit=timeunit, ustripped=true))
            nt = (scantime=time, zip(Symbol.(ions(chrom)), intensities(chrom)[idxₛ, :])...)
            itr[idxₛ] = nt
        end
    else
        header = vcat(string("scan times [", timeunit, "]"), 
            retentionindexname(rimapper(chrom)),
            [string("m/z ", i) for i in ions(chrom)])

        ris = Vector{Union{Float64, Missing}}(missing, scancount(chrom))
        ld = rimapper(chrom)
        for (i, rt) in enumerate(scantimes(chrom))
            if (minretentiontime(ld) ≤ rt ≤ maxretentiontime(ld) 
                || !isnothing(extrapolationmethod(ld)))
                ris[i] = retentionindex(ld, rt)
            end
        end

        for (idxₛ, time) in enumerate(scantimes(chrom, timeunit=timeunit, ustripped=true))
            nt = (scantime=time, retentionindex=ris[idxₛ], zip(Symbol.(ions(chrom)), 
                intensities(chrom)[idxₛ, :])...)
            itr[idxₛ] = nt
        end
    end

    CSV.write(file, itr, delim=fileformat.delim, header=header)
end


end  # module
