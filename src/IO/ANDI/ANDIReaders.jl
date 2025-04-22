module ANDIReaders

using NetCDF
using Unitful

import ...JuChrom: ChromMS
import ..InputOutput: File, FileFormat, IOError, buildxic, importdata

export ANDI


"""
    
    ANDI(; scantimetype::Union{AbstractFloat, Nothing}=Nothing, iontype::Union{
    AbstractFloat, Nothing}=Nothing, intensitytype::Union{AbstractFloat, Nothing}=Nothing)

Returns an object representing the `ANDI` file format. When using the importdata function, 
the `ANDI` data reader expects the source location to be a netCDF file, which typically has 
a .CDF file extension. The keyword arguments `scantimetype`, `iontype`, and `intensitytype` 
are provided to address cases where `ANDI` files may store data as Float64 even though it 
was originally collected as Float32. The default type Nothing indicates that no conversion 
takes place and the type in which the data is stored is returned. Background: When decimal 
numbers are converted to binary floats, they are rounded to the nearest binary fraction 
rather than a decimal fraction. This rounding can cause slight changes in the values. If 
the data were initially recorded as Float32 and are then converted to Float64, these 
changes can become visible. However, the changes can usually be ignored. If one wishes to 
convert the float types, the optional `scantimetype` keyword argument allows the scan times 
to be converted to a specified float type. Similarly, the optional keyword arguments 
`iontype` and `intensitytype` allow the conversion of mass-to-charge ratio values and 
intensity values, respectively, to a desired float type. 

See also [`FileFormat`](@ref), [`importdata`](@ref).

# Example
```jldoctest
julia> ANDI()
ANDI{Nothing, Nothing, Nothing}(Nothing, Nothing, Nothing)

julia> cdffile = joinpath(JuChrom.andi, "C7-C40_13_Nov_1.CDF");

julia> chrom = importdata(cdffile, ANDI())
ChromMS {scan times: Float64, ions: Float32, intensities: Float32}
5221 scans; scan time range: 191.942 s - 3898.719 s
5248 ions; range: m/z 29.0 - 562.9
intensity range: 0.0 - 1.051136e6
metadata: 2 entries

julia> scantimes(chrom)[9:11]
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 197.623 s
 198.33300000000003 s
 199.043 s

julia> chrom = importdata(cdffile, ANDI(scantimetype=Float32))
ChromMS {scan times: Float32, ions: Float32, intensities: Float32}
5221 scans; scan time range: 191.942f0 s - 3898.719f0 s
5248 ions; range: m/z 29.0 - 562.9
intensity range: 0.0 - 1.051136e6
metadata: 2 entries

julia> scantimes(chrom)[9:11]
3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 197.623f0 s
 198.333f0 s
 199.043f0 s
```
"""
struct ANDI{T1<:Union{AbstractFloat, Nothing}, T2<:Union{AbstractFloat, Nothing}, 
    T3<:Union{AbstractFloat, Nothing}} <: FileFormat
    scantimetype::Type{T1}
    iontype::Type{T2}
    intensitytype::Type{T3}
    function ANDI{T1, T2, T3}(scantimetype::Type{T1}, iontype::Type{T2}, 
        intensitytype::Type{T3}) where {T1<:Union{AbstractFloat, Nothing}, 
        T2<:Union{AbstractFloat, Nothing}, T3<:Union{AbstractFloat, Nothing}}
        new(scantimetype, iontype, intensitytype)
    end
end

function ANDI(; scantimetype::Type{T1}=Nothing, iontype::Type{T2}=Nothing, 
    intensitytype::Type{T3}=Nothing) where {T1<:Union{AbstractFloat, Nothing}, 
    T2<:Union{AbstractFloat, Nothing}, T3<:Union{AbstractFloat, Nothing}}
    ANDI{T1, T2, T3}(scantimetype, iontype, intensitytype)
end

scantimetype(fileformat::FileFormat) = fileformat.scantimetype
iontype(fileformat::FileFormat) = fileformat.iontype
intensitytype(fileformat::FileFormat) = fileformat.intensitytype

struct ANDIV1 end

function readfile(::ANDIV1, file::AbstractString, ::Type{T1}, ::Type{T2}, ::Type{T3}
    ) where {T1, T2, T3}

    scantimes = ncread(file, "scan_acquisition_time")
    scantimeunit = ncgetatt(file, "time_values", "units")
    pointcounts = ncread(file, "point_count")
    mzvec = ncread(file, "mass_values")
    intsvec = ncread(file, "intensity_values")

    # Currently not implemented entries:
    ionunit = ncgetatt(file, "mass_values", "units")
    intensityunit = ncgetatt(file, "intensity_values", "units")

    # Sanity checks
    cond1::Bool = length(scantimes) == 0
    cond2::Bool = length(pointcounts) == 0
    cond3::Bool = length(mzvec) == 0
    cond4::Bool = length(intsvec) == 0

    (cond1 || cond2 || cond3 || cond4) && throw(IOError("cannot read file \"$file\""))

    # This may need to be updated in the future! So far I have only seen "Seconds" as a 
    # scantime unit
    @assert scantimeunit == "Seconds" string("scantimes in \"$file\" in an unexpected ", 
        "unit: \"$scantimeunit\"")

    scantimes_unitful = T1 <: Nothing ? (scantimes * 1u"s") : (convert(Vector{T1}, 
        scantimes) * 1u"s")

    ions, xic = buildxic(pointcounts, mzvec, intsvec)

    mzs = T2 <: Nothing ? ions : convert(Vector{T2}, ions)
    im = T3 <: Nothing ? xic : convert(Matrix{T3}, xic)

    ChromMS(scantimes_unitful, mzs, im, Dict(:ionunit => ionunit, 
        :intensityunit => intensityunit))
end

function importdata(::File, fileformat::ANDI, file::AbstractString)
    readfile(ANDIV1(), file, scantimetype(fileformat), iontype(fileformat), 
        intensitytype(fileformat))
end

function importdata(fileformat::ANDI, source::AbstractString)
    if isfile(source)
        importdata(File(), fileformat, source)
    else
        throw(IOError("unsupported source for ANDI data import"))
    end
end

end  # module
