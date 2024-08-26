module ANDIReaders

using NetCDF
using Unitful

import ...JuChrom: GCMS
import ..InputOutput: File, FileFormat, IOError, buildxic, importdata

export ANDI


"""
    
    ANDI(; scantimetype::Union{AbstractFloat, Nothing}=Nothing, iontype::Union{
    AbstractFloat, Nothing}=Nothing, intensitytype::Union{AbstractFloat, Nothing}=Nothing)

Returns an object representing the `ANDI` file format. `ANDI` files typically have a .CDF 
file extension. The optional keyword argument `scantimetype` allows the scan times to be 
converted to a specified float type. Similarly, the optional keyword arguments `iontype` 
and `intensitytype` enable the conversion of mass-to-charge ratio values and intensity 
values, respectively, to a desired float type. These keyword arguments — `scantimetype`, 
`iontype`, and `intensitytype` — are provided to address cases where `ANDI` files store 
data as Float64 even though it was originally collected as Float32. When using the 
importdata function, note that the `ANDI` data reader expects the source location to 
be a netCDF file.

See also [`FileFormat`](@ref), [`importdata`](@ref).

# Example
```jldoctest
julia> ANDI()
ANDI{Nothing, Nothing, Nothing}(Nothing, Nothing, Nothing)

julia> cdffile = joinpath(JuChrom.andi, "C7-C40_13_Nov_1.CDF");

julia> gcms = importdata(cdffile, ANDI())
GCMS {scan times: Float64, ions: Float32, intensities: Float32}
5221 scans; scan time range: 191.942 s - 3898.719 s
5248 ions; range: m/z 29.0 - 562.9
intensity range: 0.0 - 1.051136e6
metadata: 2 entries

julia> scantimes(gcms)[6:11]
6-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
            195.493 s
            196.203 s
            196.913 s
            197.623 s
 198.33300000000003 s
            199.043 s

julia> gcms = importdata(cdffile, ANDI(scantimetype=Float32))
GCMS {scan times: Float32, ions: Float32, intensities: Float32}
5221 scans; scan time range: 191.942f0 s - 3898.719f0 s
5248 ions; range: m/z 29.0 - 562.9
intensity range: 0.0 - 1.051136e6
metadata: 2 entries

julia> scantimes(gcms)[6:11]
6-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 195.493f0 s
 196.203f0 s
 196.913f0 s
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

    if cond1 || cond2 || cond3 || cond4
        error("cannot read file '$file'")
    end

    if scantimeunit == "Seconds"
        scantimes_unitful = T1 <: Nothing ? (scantimes * 1u"s") : (
            convert(Vector{T1}, scantimes) * 1u"s")
    else
        throw(ErrorException(string("scantimes in \"$file\" in an unexpected unit: ", 
            "\"$scantimeunit\"")))
    end

    ions, xic = buildxic(pointcounts, mzvec, intsvec)

    mzs = T2 <: Nothing ? ions : convert(Vector{T2}, ions)
    im = T3 <: Nothing ? xic : convert(Matrix{T2}, xic)

    GCMS(scantimes_unitful, mzs, im, Dict(:ionunit => ionunit, 
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
