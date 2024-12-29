module ExcelWriter

using XLSX
using Unitful

import ...JuChrom: AbstractChrom, AbstractChromMS, scantimes, intensities, ions, rimapper, 
    retentionindexname, minretentiontime, maxretentiontime, extrapolationmethod,
    retentionindex, scancount
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
    if isnothing(rimapper(chrom))
        headers = [string("scan time [", timeunit, "]"), "intensity"]
        columns = [scantimes(chrom, timeunit=timeunit, ustripped=true), intensities(chrom)]
    else
        ris = Vector{Union{Float64, Missing}}(missing, scancount(chrom))
        ld = rimapper(chrom)
        for (i, rt) in enumerate(scantimes(chrom))
            if (minretentiontime(ld) ≤ rt ≤ maxretentiontime(ld) 
               || !isnothing(extrapolationmethod(ld)))
                ris[i] = retentionindex(ld, rt)
            end
        end
        headers = [string("scan time [", timeunit, "]"), retentionindexname(ld), 
            "intensity"]
        columns = [scantimes(chrom, timeunit=timeunit, ustripped=true), ris, 
            intensities(chrom)]
    end
    XLSX.openxlsx(file, mode="w") do xf
        !isnothing(sheetname(fileformat)) && XLSX.rename!(xf[1], sheetname(fileformat))
        XLSX.writetable!(xf[1], columns, headers)
    end
end


function exportdata(fileformat::Excel, chrom::AbstractChromMS, file::AbstractString;
    timeunit::Unitful.TimeUnits, overwrite::Bool)
    pathstem, suffix = splitext(file)
    if suffix ≠ ".xlsx"
        file = string(pathstem, ".xlsx")
    end
    overwrite || isfile(file) && throw(
        FileExistsError("a file with the same name already exists: \"$file\""))

    XLSX.openxlsx(file, mode="w") do xf
        !isnothing(sheetname(fileformat)) && XLSX.rename!(xf[1], sheetname(fileformat))
        sheet = xf[1]

        # Write scan time header
        sheet["A1"] = string("scan times [", timeunit, "]")

        # Write scan times
        sheet["A2", dim=1] = collect(scantimes(chrom, timeunit=timeunit, ustripped=true))

        if isnothing(rimapper(chrom))
            # Write ions
            sheet["B1", dim=2] = [string("m/z ", i) for i in ions(chrom)]

            # Write intensities
            sheet["B2"] = collect(intensities(chrom))
        else
            ris = Vector{Union{Float64, Missing}}(missing, scancount(chrom))
            ld = rimapper(chrom)
            for (i, rt) in enumerate(scantimes(chrom))
                if (minretentiontime(ld) ≤ rt ≤ maxretentiontime(ld) 
                   || !isnothing(extrapolationmethod(ld)))
                    ris[i] = retentionindex(ld, rt)
                end
            end

            # Write retention index header
            sheet["B1"] = retentionindexname(ld)

            # Write retention indices
            sheet["B2", dim=1] = ris

            # Write ions
            sheet["C1", dim=2] = [string("m/z ", i) for i in ions(chrom)]

            # Write intensities
            sheet["C2"] = collect(intensities(chrom))
        end
    end
end


end  # module
