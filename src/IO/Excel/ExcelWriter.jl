module ExcelWriter

using XLSX
using Unitful

import ...JuMS: AbstractGC, scantimes, intensities
import ..InputOutput: FileFormat, exportdata

export Excel

include("../Utilities.jl")

struct Excel{T<:Union{AbstractString, Nothing}} <: FileFormat
    sheetname::T
    function Excel{T}(sheetname::T) where T<:Union{AbstractString, Nothing}
        new(sheetname)
    end
end
# TYPE STABLE June 11, 2024


sheetname(fileformat::FileFormat) = fileformat.sheetname
# TYPE STABLE June 11, 2024


"""
    
    Excel(;sheetname::AbstractString) <: Fileformat

Return an Excel file format object. The optional keyword argument `sheetname` allows you to specify the name of the sheet being operated on.

# Example
```julia-repl
julia> Excel()
Excel{Nothing}(nothing)

julia> Excel(sheetname="TIC")
Excel{String}("TIC")
```

See also [`Fileformat`](@ref), [`DelimitedText`](@ref), [`exportdata`](@ref).
"""
Excel(;sheetname::T=nothing) where T<:Union{AbstractString, Nothing}= Excel{T}(sheetname)
# TYPE STABLE June 11, 2024


function exportdata(fileformat::Excel, gc::AbstractGC, file::AbstractString, timeunit::Unitful.TimeUnits, overwrite::Bool)
    file = string(first(splitext(file)), ".xlsx")
    overwrite || isfile(file) && throw(FileExistsError("a file with the same name already exists: \"$file\""))
    columns = [scantimes(gc, timeunit=timeunit, ustripped=true), intensities(gc)]
    XLSX.openxlsx(file, mode="w") do xf
        !isnothing(sheetname(fileformat)) && XLSX.rename!(xf[1], sheetname(fileformat))
        XLSX.writetable!(xf[1], columns, [string("scan time [", timeunit, "]"), "total intensity"])
    end
end
# TYPE STABLE June 11, 2024

end  # module
