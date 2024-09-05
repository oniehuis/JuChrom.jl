module ExcelWriter

using XLSX
using Unitful

import ...JuChrom: AbstractChrom, scantimes, intensities
import ..InputOutput: FileExistsError, FileFormat, exportdata

export Excel


struct Excel{T<:Union{AbstractString, Nothing}} <: FileFormat
    sheetname::T
    function Excel{T}(sheetname::T) where T<:Union{AbstractString, Nothing}
        new(sheetname)
    end
end


sheetname(fileformat::FileFormat) = fileformat.sheetname


"""
    
    Excel(; sheetname::AbstractString) <: FileFormat

Returns an `Excel` file format object. The optional keyword argument `sheetname` allows you 
to specify the name of the sheet to operate on.

See also [`FileFormat`](@ref), [`exportdata`](@ref).

# Examples
```julia-repl
julia> Excel()
Excel{Nothing}(nothing)

julia> Excel(sheetname="TIC")
Excel{String}("TIC")
```
"""
Excel(; sheetname::T=nothing) where T<:Union{AbstractString, Nothing}= Excel{T}(sheetname)


function exportdata(fileformat::Excel, chrom::AbstractChrom, file::AbstractString;
    timeunit::Unitful.TimeUnits, overwrite::Bool)
    pathstem, suffix = splitext(file)
    if suffix ≠ ".xlsx"
        file = string(pathstem, ".xlsx")
    end
    overwrite || isfile(file) && throw(
        FileExistsError("a file with the same name already exists: \"$file\""))
    columns = [scantimes(chrom, timeunit=timeunit, ustripped=true), intensities(chrom)]
    XLSX.openxlsx(file, mode="w") do xf
        !isnothing(sheetname(fileformat)) && XLSX.rename!(xf[1], sheetname(fileformat))
        XLSX.writetable!(xf[1], columns, 
            [string("scan time [", timeunit, "]"), "intensity"])
    end
end

end  # module
