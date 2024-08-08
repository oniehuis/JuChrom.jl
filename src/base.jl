using Reexport
@reexport using Unitful

"""
    AbstractChromatogram

The `AbstractChromatogram` type is the supertype of all chromatogram implementations in 
JuChrom. All subtypes of `AbstractChromatogram` (e.g., `FID`, `GCMS`, `TIC`) contain 
`scantimes`, `intensities`, and `source` information.

See also [`AbstractGC`](@ref), [`AbstractGCMS`](@ref), [`AbstractFID`](@ref), 
[`AbstractTIC`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), [`intensities`](@ref), 
[`scantimes`](@ref), [`source`](@ref).
"""
abstract type AbstractChromatogram end


"""
    AbstractGC <: AbstractChromatogram

The `AbstractGC` type is the supertype of all chromatogram implementations that have no 
mass-charge ratio (m/z) data (= `ions`) and therefore have a single intensity value 
associated with a given `scantime` in JuChrom (e.g., `FID`, `TIC`). The `intensities` are 
stored in a vector whose index corresponds to the index of the associated `scantime`.

See also [`AbstractChromatogram`](@ref), [`AbstractFID`](@ref), [`AbstractTIC`](@ref), 
[`FID`](@ref), [`TIC`](@ref), [`intensities`](@ref), [`scantimes`](@ref), [`source`](@ref).
"""
abstract type AbstractGC <: AbstractChromatogram end


"""
    AbstractFID <: AbstractGC

The `AbstractFID` type is the supertype of all flame ionization detector chromatogram 
implementations in JuMS (e.g., `FID`).

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`FID`](@ref), 
[`intensities`](@ref), [`scantimes`](@ref), [`source`](@ref).
"""
abstract type AbstractFID <: AbstractGC end


struct FID{
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real}, 
    T3<:Union{Nothing, AbstractString}} <: AbstractFID
    scantimes::T1
    intensities::T2
    source::T3
    function FID{T1, T2, T3}(scantimes::T1, intensities::T2, source::T3) where {
        T1<:AbstractVector{<:Unitful.Time},
        T2<:AbstractVector{<:Real},
        T3<:Union{Nothing, AbstractString}}
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan count"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain values less than zero"))
        new(scantimes, intensities, source)
    end
end


Base.broadcastable(fid::FID) = Ref(fid)


"""
    FID(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, 
    source::Union{Nothing, AbstractString}=nothing)

Return a FID object consisting of scan times, intensities, and optionally data source 
information.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractFID`](@ref), 
[`scantimes`](@ref), [`intensities`](@ref), [`source`](@ref).

In the following examples, the number types of the arrays passed to the object constructor
are explicitly annotated to illustrate that the `FID` object preserves the types.

# Example
```jldoctest
julia> FID(Int64[1, 2, 3]u"s", Int32[12, 956, 1])
FID {scantimes: Int64, intensities: Int32}
Source: nothing
3 scans; time range: 1 s - 3 s
intensity range: 1 - 956

julia> FID(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], "example data")
FID {scantimes: Int32, intensities: Float64}
Source: example data
3 scans; time range: 1 s - 3 s
intensity range: 1.0 - 956.0
```
"""
function FID(scantimes::T1, intensities::T2, source::T3=nothing) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real},
    T3<:Union{Nothing, AbstractString}}
    FID{T1, T2, T3}(scantimes, intensities, source)
end


"""
    AbstractGCMS <: AbstractChromatogram

The `AbstractGCMS` type is the supertype of all chromatogram implementations that contain 
mass-charge ratio (m/z) data (= `ions`) and associated abundance values (= `intensities`) 
and thus can have one or more ion intensity values associated with a given scan time in 
JuChrom (e.g., `GCMS`). The `intensities` are stored in a matrix where the row index 
corresponds to that of the associated `scantime` and where the column index corresponds
to that of the associated `ion`.

See also [`AbstractChromatogram`](@ref), [`GCMS`](@ref), [`intensities`](@ref), 
[`ions`](@ref), [`scantimes`](@ref), [`source`](@ref).
"""
abstract type AbstractGCMS <: AbstractChromatogram end


struct GCMS{
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real},
    T3<:AbstractMatrix{<:Real},
    T4<:Union{Nothing, AbstractString}} <: AbstractGCMS
    scantimes::T1
    ions::T2
    intensities::T3
    source::T4
    function GCMS{T1, T2, T3, T4}(scantimes::T1, ions::T2, intensities::T3, source::T4
        ) where {
            T1<:AbstractVector{<:Unitful.Time},
            T2<:AbstractVector{<:Real},
            T3<:AbstractMatrix{<:Real},
            T4<:Union{Nothing, AbstractString}}
        issorted(scantimes) || throw(ArgumentError("scan times not in ascending order"))
        issorted(ions) || throw(ArgumentError("ions not in ascending order"))
        size(intensities, 1) == length(scantimes) || throw(DimensionMismatch(
            "intensity matrix row count does not match scan count"))
        size(intensities, 2) == length(ions) || throw(DimensionMismatch(
            "intensity matrix column count does not match ion count"))
        count(i -> i < 0, intensities) == 0 || throw(ArgumentError(
            "intensity matrix contains at least one value less than zero"))
        new(scantimes, ions, intensities, source)
    end
end

Base.broadcastable(gcms::GCMS) = Ref(gcms)


"""
    GCMS(scantimes::AbstractVector{<:Unitful.Time}, ions::AbstractVector{<:Real}, 
    intensities::AbstractMatrix{<:Real}; source::Union{Nothing, AbstractString}) 

Return a GCMS object consisting of `scantimes`, `ions`, `intensities`, and optionally data 
`source` information.

See also [`AbstractChromatogram`](@ref), [`AbstractGCMS`](@ref), [`intensities`](@ref), 
[`ions`](@ref), [`scantimes`](@ref), [`source`](@ref).

In the following examples, the number types of the arrays passed to the object constructor
are explicitly annotated to illustrate that the `GCMS` object preserves the types.

# Example
```jldoctest
julia> GCMS(Int32[1, 2, 3]u"s", Int64[85, 100], Int32[0 12; 34 956; 23 1])
GCMS {scantimes: Int32, ions: Int64, intensities: Int32}
Source: nothing
3 scans; time range: 1 s - 3 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956

julia> GCMS([1.1f0, 2.1f0]u"s", Float32[35.1, 76.2], Int64[0 12; 34 956], "example data")
GCMS {scantimes: Float32, ions: Float32, intensities: Int64}
Source: example data
2 scans; time range: 1.1f0 s - 2.1f0 s
2 ions; range: m/z 35.1 - 76.2
intensity range: 0 - 956
```
"""
function GCMS(scantimes::T1, ions::T2, intensities::T3, source::T4=nothing) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real},
    T3<:AbstractMatrix{<:Real},
    T4<:Union{Nothing, AbstractString}}
    GCMS{T1, T2, T3, T4}(scantimes, ions, intensities, source)
end


"""
    AbstractTIC <: AbstractGC

The `AbstractTIC` type is the supertype of all total ion chromatogram implementations in 
JuChrom (e.g., `TIC`).

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`TIC`](@ref), 
[`intensities`](@ref), [`scantimes`](@ref), [`source`](@ref).
"""
abstract type AbstractTIC <: AbstractGC end


struct TIC{
    T1<:AbstractVector{<:Unitful.Time}, 
    T2<:AbstractVector{<:Real}, 
    T3<:Union{Nothing, AbstractString}} <: AbstractTIC
    scantimes::T1
    intensities::T2
    source::T3
    function TIC{T1, T2, T3}(scantimes::T1, intensities::T2, source::T3) where {
        T1<:AbstractVector{<:Unitful.Time},
        T2<:AbstractVector{<:Real},
        T3<:Union{Nothing, AbstractString}}
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan count"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain values less than zero"))
        new(scantimes, intensities, source)
    end
end


Base.broadcastable(tic::TIC) = Ref(tic)


"""
    TIC(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, 
    source::Union{Nothing, AbstractString}=nothing)

Return a TIC object consisting of scan times, intensities, and optionally data source 
information.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractTIC`](@ref), 
[`scantimes`](@ref), [`intensities`](@ref), [`source`](@ref).

In the following examples, the number types of the arrays passed to the object constructor
are explicitly annotated to illustrate that the `TIC` object preserves the types.

# Example
```jldoctest
julia> TIC(Int64[1, 2, 3]u"s", Int32[12, 956, 1])
TIC {scantimes: Int64, intensities: Int32}
Source: nothing
3 scans; time range: 1 s - 3 s
intensity range: 1 - 956

julia> TIC(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], "example data")
TIC {scantimes: Int32, intensities: Float64}
Source: example data
3 scans; time range: 1 s - 3 s
intensity range: 1.0 - 956.0
```
"""
function TIC(scantimes::T1, intensities::T2, source::T3=nothing) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real},
    T3<:Union{Nothing, AbstractString}}
    TIC{T1, T2, T3}(scantimes, intensities, source)
end


"""
    ions(gcms::AbstractGCMS)

Return the `ions` of an `AbstractGCMS` subtype object (e.g., `GCMS`).

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref).

In the following example, the number type of `ions` passed to the object constructor 
is explicitly annotated to illustrate that the GCMS object preserves the type.

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", Int64[85, 100], [0 12; 34 956; 23 1]);

julia> ions(gcms)
2-element Vector{Int64}:
  85
 100
```
"""
ions(gcms::AbstractGCMS) = gcms.ions


"""
    scantimes(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the `scantimes` of an `AbstractChromatogram` subtype object (e.g., `FID`, `GCMS`, 
`TIC`). The optional keyword argument `timeunit` allows you to change the unit of the 
returned `scantimes`. All time units defined in the package [Unitful.jl]
(https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are supported. 
The optional keyword argument `ustripped` allows you to specify whether the unit is 
stripped from the returned values. 

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref).

In the following example, the element type of the `scantime` vector passed to the object 
constructor is explicitly annotated to illustrate that the GCMS object preserves the type.

# Example
```jldoctest
julia> gcms = GCMS(Float32[1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> scantimes(gcms)
3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 1.0f0 s
 2.0f0 s
 3.0f0 s

julia> scantimes(gcms, timeunit=u"minute")
3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:
 0.016666668f0 minute
 0.033333335f0 minute
 0.050000004f0 minute

julia> scantimes(gcms, timeunit=u"minute", ustripped=true)
3-element Vector{Float32}:
 0.016666668
 0.033333335
 0.050000004
```
"""
function scantimes(chrom::AbstractChromatogram; 
    timeunit::Unitful.TimeUnits=unit(eltype(chrom.scantimes)), ustripped::Bool=false)
    ustripped ? ustrip.(timeunit, chrom.scantimes) : uconvert.(timeunit, chrom.scantimes)
end


"""
    intensities(chrom::AbstractChromatogram)

Return the intensities of an `AbstractChromatogram` subtype object.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractGCMS`](@ref),
[`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], Int64[0 12; 34 956; 23 1]);

julia> intensities(gcms)
3×2 Matrix{Int64}:
  0   12
 34  956
 23    1

julia> fid = FID([1, 2, 3]u"s", Float32[12.0, 956.0, 23.0]);

julia> intensities(fid)
3-element Vector{Float32}:
  12.0
 956.0
  23.0
```
"""
intensities(chrom::AbstractChromatogram) = chrom.intensities


"""
    source(chrom::AbstractChromatogram) -> Type{<:Union{AbstractString, Nothing}}

Return the source information of the AbstractChromatogram subtype object.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref).

# Example
```jldoctest
julia> gcms₁ = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1], "example data");

julia> source(gcms₁)
"example data"

julia> gcms₂ = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> source(gcms₂)

julia> isnothing(source(gcms₂))
true
```
"""
source(chrom::AbstractChromatogram) = chrom.source


"""
    scancount(chrom::AbstractChromatogram) -> Int

Return the number of scans from the AbstractChromatogram subtype object.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`scantimes`](@ref).

# Example
```jldoctest
julia> fid = FID([1, 2, 3]u"s", [12, 956, 23]);

julia> scancount(fid)
3
```
"""
scancount(chrom::AbstractChromatogram) = length(scantimes(chrom))


function minscantime(chrom::AbstractChromatogram;
    timeunit::Unitful.TimeUnits=unit(eltype(chrom.scantimes)), ustripped::Bool=false)
    ustripped ? ustrip(timeunit, first(chrom.scantimes)) : uconvert(timeunit, 
        first(chrom.scantimes))
end

function maxscantime(chrom::AbstractChromatogram; 
    timeunit::Unitful.TimeUnits=unit(eltype(chrom.scantimes)), ustripped::Bool=false)
    ustripped ? ustrip(timeunit, last(chrom.scantimes)) : uconvert(timeunit, 
    last(chrom.scantimes))
end

ioncount(gcms::AbstractGCMS) = length(ions(gcms))

minion(gcms::AbstractGCMS) = first(ions(gcms))
maxion(gcms::AbstractGCMS) = last(ions(gcms))
minintensity(chrom::AbstractChromatogram) = minimum(intensities(chrom))
maxintensity(chrom::AbstractChromatogram) = maximum(intensities(chrom))


function Base.show(io::IO, fid::FID)
    println(io, "FID {scantimes: ", eltype(ustrip.(scantimes(fid))), ", intensities: ", 
        eltype(intensities(fid)), "}")
    println(io, "Source: ", source(fid))
    println(io, scancount(fid), " scans; time range: ", minscantime(fid), " - ", 
        maxscantime(fid))
    print(io, "intensity range: ", minintensity(fid), " - ", maxintensity(fid))
end


function Base.show(io::IO, tic::TIC)
    println(io, "TIC {scantimes: ", eltype(ustrip.(scantimes(tic))), ", intensities: ", 
        eltype(intensities(tic)), "}")
    println(io, "Source: ", source(tic))
    println(io, scancount(tic), " scans; time range: ", minscantime(tic), " - ", 
        maxscantime(tic))
    print(io, "intensity range: ", minintensity(tic), " - ", maxintensity(tic))
end


function Base.show(io::IO, gcms::GCMS)
    println(io, "GCMS {scantimes: ", eltype(ustrip.(scantimes(gcms))), ", ions: ", 
        eltype(ions(gcms)), ", intensities: ", eltype(intensities(gcms)), "}")
    println(io, "Source: ", source(gcms))
    println(io, scancount(gcms), " scans; time range: ", minscantime(gcms), " - ", 
        maxscantime(gcms))
    println(io, ioncount(gcms), " ions; range: m/z ", minion(gcms), " - ", maxion(gcms))
    print(io, "intensity range: ", minintensity(gcms), " - ", maxintensity(gcms))
end
