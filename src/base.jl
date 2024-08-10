using Reexport
@reexport using Unitful

"""
    AbstractChromatogram

The `AbstractChromatogram` type is the supertype of all chromatogram implementations in 
JuChrom. All subtypes of `AbstractChromatogram` (e.g., `FID`, `GCMS`, `TIC`) contain 
`scantimes`, `intensities`, and `metadata`.

See also [`AbstractGC`](@ref), [`AbstractGCMS`](@ref), [`AbstractFID`](@ref), 
[`AbstractTIC`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), [`intensities`](@ref), 
[`scantimes`](@ref), [`metadata`](@ref).
"""
abstract type AbstractChromatogram end


"""
    AbstractGC <: AbstractChromatogram

The `AbstractGC` type is the supertype of all chromatogram implementations that have no 
mass-charge ratio (*m*/*z*) data (= `ions`) and therefore have a single intensity value 
associated with a given `scantime` in JuChrom (e.g., `FID`, `TIC`). The `intensities` are 
stored in a vector whose index corresponds to the index of the associated `scantime`.

See also [`AbstractChromatogram`](@ref), [`AbstractFID`](@ref), [`AbstractTIC`](@ref), 
[`FID`](@ref), [`TIC`](@ref), [`intensities`](@ref), [`scantimes`](@ref), 
[`metadata`](@ref).
"""
abstract type AbstractGC <: AbstractChromatogram end


"""
    AbstractFID <: AbstractGC

The `AbstractFID` type is the supertype of all flame ionization detector chromatogram 
implementations in JuMS (e.g., `FID`).

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`FID`](@ref), 
[`intensities`](@ref), [`scantimes`](@ref), [`metadata`](@ref).
"""
abstract type AbstractFID <: AbstractGC end


struct FID{
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real}} <: AbstractFID
    scantimes::T1
    intensities::T2
    metadata::Dict{Any, Any}
    function FID{T1, T2}(scantimes::T1, intensities::T2, metadata::Dict) where {
        T1<:AbstractVector{<:Unitful.Time},
        T2<:AbstractVector{<:Real}}
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan count"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain values less than zero"))
        new(scantimes, intensities, metadata)
    end
end


Base.broadcastable(fid::FID) = Ref(fid)


"""
    FID(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, 
    metadata::Dict=Dict{Any, Any})

Return a FID object consisting of scan times, intensities, and metadata.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractFID`](@ref), 
[`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref).

In the following examples, the number types of the arrays passed to the object constructor 
are explicitly annotated to illustrate that the `FID` object preserves the types.

# Example
```jldoctest
julia> FID(Int64[1, 2, 3]u"s", Int32[12, 956, 1])
FID {scantimes: Int64, intensities: Int32}
3 scans; time range: 1 s - 3 s
intensity range: 1 - 956
metadata: 

julia> FID(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], Dict(:id => 1, :name => "sample"))
FID {scantimes: Int32, intensities: Float64}
3 scans; time range: 1 s - 3 s
intensity range: 1.0 - 956.0
metadata: :id, :name
```
"""
function FID(scantimes::T1, intensities::T2, metadata::Dict=Dict{Any, Any}()) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real}}
    FID{T1, T2}(scantimes, intensities, metadata)
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
[`ions`](@ref), [`scantimes`](@ref), [`metadata`](@ref).
"""
abstract type AbstractGCMS <: AbstractChromatogram end


struct GCMS{
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real},
    T3<:AbstractMatrix{<:Real}} <: AbstractGCMS
    scantimes::T1
    ions::T2
    intensities::T3
    metadata::Dict{Any, Any}
    function GCMS{T1, T2, T3}(scantimes::T1, ions::T2, intensities::T3, metadata::Dict
        ) where {
            T1<:AbstractVector{<:Unitful.Time},
            T2<:AbstractVector{<:Real},
            T3<:AbstractMatrix{<:Real}}
        issorted(scantimes) || throw(ArgumentError("scan times not in ascending order"))
        issorted(ions) || throw(ArgumentError("ions not in ascending order"))
        size(intensities, 1) == length(scantimes) || throw(DimensionMismatch(
            "intensity matrix row count does not match scan count"))
        size(intensities, 2) == length(ions) || throw(DimensionMismatch(
            "intensity matrix column count does not match ion count"))
        count(i -> i < 0, intensities) == 0 || throw(ArgumentError(
            "intensity matrix contains at least one value less than zero"))
        new(scantimes, ions, intensities, metadata)
    end
end

Base.broadcastable(gcms::GCMS) = Ref(gcms)


"""
    GCMS(scantimes::AbstractVector{<:Unitful.Time}, ions::AbstractVector{<:Real}, 
    intensities::AbstractMatrix{<:Real}; metadata::Dict=Dict()) 

Return a GCMS object consisting of `scantimes`, `ions`, `intensities`, and  `metadata`.

See also [`AbstractChromatogram`](@ref), [`AbstractGCMS`](@ref), [`intensities`](@ref), 
[`ions`](@ref), [`scantimes`](@ref), [`metadata`](@ref).

In the following examples, the number types of the arrays passed to the object constructor 
are explicitly annotated to illustrate that the `GCMS` object preserves the types.

# Example
```jldoctest
julia> GCMS(Int32[1, 2, 3]u"s", Int64[85, 100], Int32[0 12; 34 956; 23 1])
GCMS {scantimes: Int32, ions: Int64, intensities: Int32}
3 scans; time range: 1 s - 3 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 

julia> GCMS([1.1f0, 2.1f0]u"s", Float32[35.1, 76.2], Int64[0 12; 34 956], Dict(:id => 1, :name => "sample"))
GCMS {scantimes: Float32, ions: Float32, intensities: Int64}
2 scans; time range: 1.1f0 s - 2.1f0 s
2 ions; range: m/z 35.1 - 76.2
intensity range: 0 - 956
metadata: :id, :name
```
"""
function GCMS(scantimes::T1, ions::T2, intensities::T3, metadata::Dict=Dict{Any, Any}()
    ) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real},
    T3<:AbstractMatrix{<:Real}}
    GCMS{T1, T2, T3}(scantimes, ions, intensities, metadata)
end


"""
    AbstractTIC <: AbstractGC

The `AbstractTIC` type is the supertype of all total ion chromatogram implementations in 
JuChrom (e.g., `TIC`).

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`TIC`](@ref), 
[`intensities`](@ref), [`scantimes`](@ref), [`metadata`](@ref).
"""
abstract type AbstractTIC <: AbstractGC end


struct TIC{
    T1<:AbstractVector{<:Unitful.Time}, 
    T2<:AbstractVector{<:Real}} <: AbstractTIC
    scantimes::T1
    intensities::T2
    metadata::Dict{Any, Any}
    function TIC{T1, T2}(scantimes::T1, intensities::T2, metadata::Dict) where {
        T1<:AbstractVector{<:Unitful.Time},
        T2<:AbstractVector{<:Real}}
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan count"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain values less than zero"))
        new(scantimes, intensities, metadata)
    end
end


Base.broadcastable(tic::TIC) = Ref(tic)


"""
    TIC(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, 
    metadata::Dict=Dict{Any, Any}())

Return a TIC object consisting of scan times, intensities, and metadata.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractTIC`](@ref), 
[`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref).

In the following examples, the number types of the arrays passed to the object constructor 
are explicitly annotated to illustrate that the `TIC` object preserves the types.

# Example
```jldoctest
julia> TIC(Int64[1, 2, 3]u"s", Int32[12, 956, 1])
TIC {scantimes: Int64, intensities: Int32}
3 scans; time range: 1 s - 3 s
intensity range: 1 - 956
metadata: 

julia> TIC(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], Dict(:id => 1, :name => "sample"))
TIC {scantimes: Int32, intensities: Float64}
3 scans; time range: 1 s - 3 s
intensity range: 1.0 - 956.0
metadata: :id, :name
```
"""
function TIC(scantimes::T1, intensities::T2, metadata::Dict=Dict{Any, Any}()) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real}}
    TIC{T1, T2}(scantimes, intensities, metadata)
end


"""
    ions(gcms::AbstractGCMS)

Return the `ions` of the `AbstractGCMS` subtype object.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`minion`](@ref), [`maxion`](@ref), 
[`ioncount`](@ref).

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

Return the `scantimes` of the `AbstractChromatogram` subtype object. The optional keyword 
argument `timeunit` allows you to change the unit of the returned `scantimes`. All time 
units defined in the package [Unitful.jl](https://painterqubits.github.io/Unitful.jl) 
(e.g., `u"s"`, `u"minute"`) are supported. The optional keyword argument `ustripped` 
allows you to specify whether the unit is stripped from the returned values.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`minscantime`](@ref), [`maxscantime`](@ref).

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

Return the intensities of the `AbstractChromatogram` subtype object.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`minintensity`](@ref), [`maxintensity`](@ref).

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
    metadata(chrom::AbstractChromatogram) -> Dict{Any, Any}

Return the metadata.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref).

# Example
```jldoctest
julia> gcms₁ = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1], Dict(:id => 1, :name => "sample"))
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
3 scans; time range: 1 s - 3 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: :id, :name

julia> metadata(gcms₁)
Dict{Any, Any} with 2 entries:
  :id   => 1
  :name => "sample"

julia> gcms₂ = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
3 scans; time range: 1 s - 3 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 

julia> metadata(gcms₂)
Dict{Any, Any}()
```
"""
metadata(chrom::AbstractChromatogram) = chrom.metadata


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


"""
    minscantime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the time of the first scan from the AbstractChromatogram subtype object. The 
optional keyword argument `timeunit` allows you to change the unit of the returned scan 
time. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` allows you to specify whether the unit 
is stripped from the returned value.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`maxscantime`](@ref), [`scantimes`](@ref), [`scancount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> minscantime(gcms)
1.0 s

julia> minscantime(gcms, timeunit=u"minute")
0.016666666666666666 minute

julia> minscantime(gcms, timeunit=u"minute", ustripped=true)
0.016666666666666666
```
"""
function minscantime(chrom::AbstractChromatogram;
    timeunit::Unitful.TimeUnits=unit(eltype(chrom.scantimes)), ustripped::Bool=false)
    ustripped ? ustrip(timeunit, first(chrom.scantimes)) : uconvert(timeunit, 
        first(chrom.scantimes))
end


"""
    maxscantime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
    ustripped::Bool)

Return the time of the last scan from the AbstractChromatogram subtype object. The optional 
keyword argument `timeunit` allows you to change the unit of the returned scan time. All 
time units defined in the package [Unitful.jl](https://painterqubits.github.io/Unitful.jl) 
(e.g., `u"s"`, `u"minute"`) are supported. The optional keyword argument `ustripped` allows 
you to specify whether the unit is stripped from the returned value.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`minscantime`](@ref), [`scantimes`](@ref), [`scancount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> maxscantime(gcms)
3.0 s

julia> maxscantime(gcms, timeunit=u"minute")
0.05 minute

julia> maxscantime(gcms, timeunit=u"minute", ustripped=true)
0.05
```
"""
function maxscantime(chrom::AbstractChromatogram; 
    timeunit::Unitful.TimeUnits=unit(eltype(chrom.scantimes)), ustripped::Bool=false)
    ustripped ? ustrip(timeunit, last(chrom.scantimes)) : uconvert(timeunit, 
    last(chrom.scantimes))
end


"""
    ioncount(gcms::AbstractGCMS) -> Int

Return the number of ions in the AbstractGCMS subtype object.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`ions`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> ioncount(gcms)
2
```
"""
ioncount(gcms::AbstractGCMS) = length(ions(gcms))


"""
    minion(gcms::AbstractGCMS)

Return the smallest ion in the AbstractGCMS subtype object.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`ions`](@ref), [`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", Int64[85, 100], [0 12; 34 956; 23 1]);

julia> minion(gcms)
85
```

"""
minion(gcms::AbstractGCMS) = first(ions(gcms))


"""
    maxion(gcms::AbstractGCMS)

Return the largest ion in the AbstractGCMS subtype object.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`minion`](@ref), [`ions`](@ref), 
[`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", Int64[85, 100], [0 12; 34 956; 23 1]);

julia> maxion(gcms)
100
```
"""
maxion(gcms::AbstractGCMS) = last(ions(gcms))


"""
    minintensity(chrom::AbstractChromatogram)

Return the minimum intensity in the AbstractChromatogram subtype object.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`maxintensity`](@ref), [`intensities`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], Int64[0 12; 34 956; 23 1]);

julia> minintensity(gcms)
0
```
"""
minintensity(chrom::AbstractChromatogram) = minimum(intensities(chrom))


"""
    maxintensity(chrom::AbstractChromatogram)

Return the maximum intensity in the AbstractChromatogram subtype object.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`minintensity`](@ref), [`intensities`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], Int64[0 12; 34 956; 23 1]);

julia> maxintensity(gcms)
956
```
"""
maxintensity(chrom::AbstractChromatogram) = maximum(intensities(chrom))


function Base.show(io::IO, fid::FID)
    println(io, "FID {scantimes: ", eltype(ustrip.(scantimes(fid))), ", intensities: ", 
        eltype(intensities(fid)), "}")
    println(io, scancount(fid), " scans; time range: ", minscantime(fid), " - ", 
        maxscantime(fid))
    println(io, "intensity range: ", minintensity(fid), " - ", maxintensity(fid))
    print(io, "metadata: ")
    first = true
    n = 0
    for key in sort(collect(keys(metadata(fid))))
        first ? (first = false) : print(io, ", ")
        show(io, key)
        n += 1
        n ≥ 10 && (print(io, " …"); break)
    end
end


function Base.show(io::IO, tic::TIC)
    println(io, "TIC {scantimes: ", eltype(ustrip.(scantimes(tic))), ", intensities: ", 
        eltype(intensities(tic)), "}")
    println(io, scancount(tic), " scans; time range: ", minscantime(tic), " - ", 
        maxscantime(tic))
    println(io, "intensity range: ", minintensity(tic), " - ", maxintensity(tic))
    print(io, "metadata: ")
    first = true
    n = 0
    for key in sort(collect(keys(metadata(tic))))
        first ? (first = false) : print(io, ", ")
        show(io, key)
        n += 1
        n ≥ 10 && (print(io, " …"); break)
    end
end


function Base.show(io::IO, gcms::GCMS)
    println(io, "GCMS {scantimes: ", eltype(ustrip.(scantimes(gcms))), ", ions: ", 
        eltype(ions(gcms)), ", intensities: ", eltype(intensities(gcms)), "}")
    println(io, scancount(gcms), " scans; time range: ", minscantime(gcms), " - ", 
        maxscantime(gcms))
    println(io, ioncount(gcms), " ions; range: m/z ", minion(gcms), " - ", maxion(gcms))
    println(io, "intensity range: ", minintensity(gcms), " - ", maxintensity(gcms))
    print(io, "metadata: ")
    first = true
    n = 0
    for key in sort(collect(keys(metadata(gcms))))
        first ? (first = false) : print(io, ", ")
        show(io, key)
        n += 1
        n ≥ 10 && (print(io, " …"); break)
    end
end
