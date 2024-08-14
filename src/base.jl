using Reexport
@reexport using Unitful
include("utilities.jl")

"""
    AbstractChromatogram

Supertype of all chromatogram implementations. All subtypes (e.g., `FID`, `GCMS`, `TIC`) 
contain `scantimes`, `intensities`, and `metadata`.

See also [`AbstractGC`](@ref), [`AbstractFID`](@ref), [`AbstractGCMS`](@ref), 
[`AbstractTIC`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), [`scantimes`](@ref), 
[`intensities`](@ref), [`metadata`](@ref).
"""
abstract type AbstractChromatogram end


"""
    AbstractGC <: AbstractChromatogram

Supertype of all chromatogram implementations that have no mass-charge ratio (*m*/*z*) data 
(= `ions`) and therefore have a single intensity value associated with a given `scantime` 
(e.g., `FID`, `TIC`). The `intensities` are stored in a vector whose index corresponds to 
the index of the associated `scantime`.

See also [`AbstractChromatogram`](@ref), [`AbstractFID`](@ref), [`AbstractTIC`](@ref), 
[`FID`](@ref), [`TIC`](@ref), [`scantimes`](@ref), [`intensities`](@ref), 
[`metadata`](@ref).
"""
abstract type AbstractGC <: AbstractChromatogram end


"""
    AbstractFID <: AbstractGC

Supertype of all flame ionization detector chromatogram implementations in JuMS (e.g., 
`FID`).

See also [`AbstractChromatogram`](@ref), [`AbstractFID`](@ref), [`AbstractGC`](@ref), 
[`FID`](@ref), [`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref).
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
            ArgumentError("scantimes not in ascending order"))
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

Construct a `FID` object consisting of `scantimes`, `intensities`, and `metadata`. Note that 
the `scantimes` must be in ascending order and the `intensities` must not contain values 
less than zero.

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
metadata: 0 entries

julia> FID(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], Dict("name" => "sample"))
FID {scantimes: Int32, intensities: Float64}
3 scans; time range: 1 s - 3 s
intensity range: 1.0 - 956.0
metadata: 1 entry

julia> FID([2, 1, 3]u"s", [12.0, 956.0, 1.0])
ERROR: ArgumentError: scantimes not in ascending order
[...]

julia> FID([1, 2, 3]u"s", [-12.0, 956.0, 1.0])
ERROR: ArgumentError: intensity values contain values less than zero
[...]
```
"""
function FID(scantimes::T1, intensities::T2, metadata::Dict=Dict{Any, Any}()) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real}}
    FID{T1, T2}(scantimes, intensities, metadata)
end


"""
    AbstractGCMS <: AbstractChromatogram

Supertype of all chromatogram implementations that contain mass-charge ratio (*m*/*z*) data 
(= `ions`) and associated abundance values (= `intensities`) and thus can have one or more 
ion intensity values associated with a given scan time (e.g., `GCMS`). The `intensities` 
are stored in a matrix where the row index corresponds to that of the associated `scantime` 
and where the column index corresponds to that of the associated `ion`.

See also [`AbstractChromatogram`](@ref), [`AbstractGCMS`](@ref), [`GCMS`](@ref), 
[`scantimes`](@ref), [`ions`](@ref), [`intensities`](@ref), [`metadata`](@ref).
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
        issorted(scantimes) || throw(ArgumentError("scantimes not in ascending order"))
        issorted(ions) || throw(ArgumentError("ions not in ascending order"))
        size(intensities, 1) == length(scantimes) || throw(DimensionMismatch(
            "intensity matrix row count does not match scan count"))
        size(intensities, 2) == length(ions) || throw(DimensionMismatch(
            "intensity matrix column count does not match ion count"))
        count(i -> i < 0, intensities) == 0 || throw(ArgumentError(
            "intensity values contain at least one value less than zero"))
        new(scantimes, ions, intensities, metadata)
    end
end

Base.broadcastable(gcms::GCMS) = Ref(gcms)


"""
    GCMS(scantimes::AbstractVector{<:Unitful.Time}, ions::AbstractVector{<:Real}, 
    intensities::AbstractMatrix{<:Real}; metadata::Dict=Dict{Any, Any}()) 

Construct a `GCMS` object consisting of `scantimes`, `ions`, `intensities`, and `metadata`. 
Note that the `scantimes` and the `ions` must be in ascending order and the `intensities` 
must not contain values less than zero.

See also [`AbstractChromatogram`](@ref), [`AbstractGCMS`](@ref), [`scantimes`](@ref), 
[`ions`](@ref), [`intensities`](@ref), [`metadata`](@ref), [`totalionchromatogram`](@ref).

In the following examples, the number types of the arrays passed to the object constructor 
are explicitly annotated to illustrate that the `GCMS` object preserves the types.

# Example
```jldoctest
julia> GCMS(Int32[1, 2, 3]u"s", Int64[85, 100], Int32[0 12; 34 956; 23 1])
GCMS {scantimes: Int32, ions: Int64, intensities: Int32}
3 scans; time range: 1 s - 3 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 0 entries

julia> GCMS([1.1f0, 2.1f0]u"s", [35.1f0, 76.2f0], Int64[0 12; 34 956], Dict(:id => 4))
GCMS {scantimes: Float32, ions: Float32, intensities: Int64}
2 scans; time range: 1.1f0 s - 2.1f0 s
2 ions; range: m/z 35.1 - 76.2
intensity range: 0 - 956
metadata: 1 entry

julia> GCMS([2, 1, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
ERROR: ArgumentError: scantimes not in ascending order
[...]

julia> GCMS([1, 2, 3]u"s", [100, 85], [0 12; 34 956; 23 1])
ERROR: ArgumentError: ions not in ascending order
[...]

julia> GCMS([1, 2, 3]u"s", [85, 100], [0 -12; 34 956; 23 1])
ERROR: ArgumentError: intensity values contain at least one value less than zero
[...]
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

Supertype of all total ion chromatogram implementations (e.g., `TIC`).

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractTIC`](@ref), 
[`TIC`](@ref), [`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref).
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
            ArgumentError("scantimes not in ascending order"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan count"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain at least one value less than zero"))
        new(scantimes, intensities, metadata)
    end
end


Base.broadcastable(tic::TIC) = Ref(tic)


"""
    TIC(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, 
    metadata::Dict=Dict{Any, Any}())

Construct a `TIC` object consisting of `scantimes`, `intensities`, and `metadata`. Note 
that the `scantimes` must be in ascending order and the `intensities` must not contain 
values less than zero.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractTIC`](@ref), 
[`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref), 
[`totalionchromatogram`](@ref).

In the following examples, the number types of the arrays passed to the object constructor 
are explicitly annotated to illustrate that the `TIC` object preserves the types.

# Example
```jldoctest
julia> TIC(Int64[1, 2, 3]u"s", Int32[12, 956, 1])
TIC {scantimes: Int64, intensities: Int32}
3 scans; time range: 1 s - 3 s
intensity range: 1 - 956
metadata: 0 entries

julia> TIC(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], Dict("name" => "sample"))
TIC {scantimes: Int32, intensities: Float64}
3 scans; time range: 1 s - 3 s
intensity range: 1.0 - 956.0
metadata: 1 entry

julia> TIC([2, 1, 3]u"s", [12.0, 956.0, 1.0])
ERROR: ArgumentError: scantimes not in ascending order
[...]

julia> TIC([1, 2, 3]u"s", [-12.0, 956.0, 1.0])
ERROR: ArgumentError: intensity values contain at least one value less than zero
[...]
```
"""
function TIC(scantimes::T1, intensities::T2, metadata::Dict=Dict{Any, Any}()) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real}}
    TIC{T1, T2}(scantimes, intensities, metadata)
end


"""
    binions(gcms::AbstractGCMS; ionbin::Function=integerion)

Return a `GCMS` object in which the `ions` are binned according to the `ionbin` function 
(default function is `integer`) and the `intensities` of the binned `ions` are summed.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`integer`](@ref), [`intensities`](@ref), 
[`ions`](@ref), [`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS((1:3)u"s", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])
GCMS {scantimes: Int64, ions: Float64, intensities: Int64}
3 scans; time range: 1 s - 3 s
3 ions; range: m/z 84.8 - 100.9
intensity range: 0 - 956
metadata: 0 entries

julia> intensities(gcms)
3×3 Matrix{Int64}:
  0  24   12
  0   0  956
 23   0    1

julia> gcmsᵢ = binions(gcms)
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
3 scans; time range: 1 s - 3 s
2 ions; range: m/z 85 - 101
intensity range: 0 - 956
metadata: 0 entries

julia> ions(gcmsᵢ)
2-element Vector{Int64}:
  85
 101

julia> intensities(gcmsᵢ)
3×2 Matrix{Int64}:
 24   12
  0  956
 23    1


julia> custom_ionbin(ion) = integer(ion, start=0.9);

julia> gcmsᵢ = binions(gcms, ionbin=custom_ionbin)
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
3 scans; time range: 1 s - 3 s
3 ions; range: m/z 84 - 101
intensity range: 0 - 956
metadata: 0 entries

julia> ions(gcmsᵢ)
3-element Vector{Int64}:
  84
  85
 101
```
"""
function binions(gcms::AbstractGCMS; ionbin::Function=integer)
    binnedions = sort(unique(ionbin.(ions(gcms))))
    imₙ = zeros(eltype(intensities(gcms)), (size(intensities(gcms), 1), length(binnedions)))
    for (i, ion) in enumerate(ions(gcms))
        binnedion = ionbin(ion)
        imₙ[:, searchsortedlast(binnedions, binnedion)] .+= @view intensities(gcms)[:, i]
    end
    imₙ
    GCMS(scantimes(gcms), binnedions, imₙ, metadata(gcms))
end


"""
    ions(gcms::AbstractGCMS)

Return the `ions`.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`ion`](@ref), [`minion`](@ref), 
[`maxion`](@ref), [`ionindex`](@ref), [`ioncount`](@ref).

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
ions(gcms::GCMS) = gcms.ions


"""
    scantime(chrom::AbstractChromatogram, index::Integer; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the `scantime` by specifying the scan `index`. The optional parameter 
`timeunit` allows you to change the unit of the returned `scantime`. All time units defined 
in the package [Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, 
`u"minute"`) are supported. The optional keyword argument `ustripped` allows you to specify 
whether the unit is stripped from the returned value.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`scantimes`](@ref), [`minscantime`](@ref), [`maxscantime`](@ref), [`scancount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> scantime(gcms, 2)
2 s

julia> scantime(gcms, 2, timeunit=u"minute")
1//30 minute

julia> scantime(gcms, 2, timeunit=u"minute", ustripped=true)
1//30

julia> scantime(gcms, 5, timeunit=u"minute", ustripped=true)
ERROR: BoundsError: attempt to access 3-element Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}} at index [5]
```
"""
function scantime(chrom::AbstractChromatogram, index::Integer; 
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))), ustripped::Bool=false)
    firstindex(scantimes(chrom)) ≤ index ≤ lastindex(scantimes(chrom)) || throw(
        BoundsError(scantimes(chrom)[index]))
    ustripped ? ustrip(timeunit, scantimes(chrom)[index]) : uconvert(timeunit, 
        scantimes(chrom)[index])
end


"""
    scantimeindex(gcms::AbstractGCMS, time::Unitful.Time; precisetime::Bool=false) -> Int

Return the index of the `scantime` closest to `time` in the `scantimes`. 
All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl)
(e.g., `u"s"`, `u"minute"`) are supported. If there is a tie, the larger `scantime` is 
returned. If the optional parameter `precisetime` is set to `true`, the specified time 
must exist in the `scantimes`, otherwise an error is thrown.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`scantimes`](@ref), [`scantime`](@ref), [`minscantime`](@ref), [`maxscantime`](@ref), 
[`scancount`](@ref), [`ionindex`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1.1f0, 2.1f0, 3.1f0]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scantimes: Float32, ions: Int64, intensities: Int64}
3 scans; time range: 1.1f0 s - 3.1f0 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 0 entries

julia> scantimeindex(gcms, 1.1f0u"s", precisetime=true)
1

julia> scantimeindex(gcms, 2.1u"s", precisetime=true)
2

julia> scantimeindex(gcms, 2.2u"s", precisetime=true)
ERROR: ArgumentError: scantime 2.2 s does not exist
[...]

julia> scantimeindex(gcms, 2.2u"s")
2
```
"""
function scantimeindex(chrom::AbstractChromatogram, time::Unitful.Time; 
    precisetime::Bool=false)
    precisetime || return findclosest(scantimes(chrom), time)
    t = ustrip(uconvert(unit(eltype(scantimes(chrom))), time))
    for (index, element) in enumerate(scantimes(chrom, ustripped=true))
        element ≈ t && return index
    end
    throw(ArgumentError("scantime $time does not exist"))
end


"""
    scantimes(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the `scantimes`. The optional keyword argument `timeunit` allows you to change the 
unit of the returned `scantimes`. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` allows you to specify whether the unit 
is stripped from the returned values.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`scantime`](@ref), [`minscantime`](@ref), [`maxscantime`](@ref), [`scancount`](@ref).

In the following example, the element type of the `scantimes` vector passed to the object 
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

Return the `intensities`.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractGCMS`](@ref), 
[`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), [`minintensity`](@ref), 
[`maxintensity`](@ref), [`scantimes`](@ref), [`scancount`](@ref), [`ions`](@ref), 
[`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], Int64[0 12; 34 956; 23 1]);

julia> intensities(gcms)
3×2 Matrix{Int64}:
  0   12
 34  956
 23    1


julia> intensities(gcms)[3, :]  # all intensites of 3rd scan
2-element Vector{Int64}:
 23
  1

julia> intensities(gcms)[:, 1]  # all intensites of 1st ion
3-element Vector{Int64}:
  0
 34
 23

julia> intensities(gcms)[2, 1]  # intensity of 2nd scan, 1st ion
34

julia> fid = FID([1, 2, 3]u"s", Float32[12.0, 956.0, 23.0]);

julia> intensities(fid)
3-element Vector{Float32}:
  12.0
 956.0
  23.0

julia> intensities(fid)[2]  # intensity of 2nd scan
956.0f0
```
"""
intensities(chrom::AbstractChromatogram) = chrom.intensities


"""
    metadata(chrom::AbstractChromatogram) -> Dict{Any, Any}

Return the `metadata`.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref).

# Example
```jldoctest
julia> gcms₁ = GCMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956], Dict(:id => 4))
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
2 scans; time range: 1 s - 2 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 1 entry

julia> metadata(gcms₁)
Dict{Any, Any} with 1 entry:
  :id => 4

julia> gcms₂ = GCMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956])
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
2 scans; time range: 1 s - 2 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 0 entries

julia> metadata(gcms₂)
Dict{Any, Any}()

julia> metadata(gcms₂)["name"] = "sample"
"sample"

julia> metadata(gcms₂)[:id] = 123
123

julia> metadata(gcms₂)
Dict{Any, Any} with 2 entries:
  "name" => "sample"
  :id    => 123

julia> gcms₂
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
2 scans; time range: 1 s - 2 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 2 entries

julia> delete!(metadata(gcms₂), "name")
Dict{Any, Any} with 1 entry:
  :id => 123

julia> gcms₂
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
2 scans; time range: 1 s - 2 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 1 entry
```
"""
metadata(chrom::AbstractChromatogram) = chrom.metadata


"""
    scancount(chrom::AbstractChromatogram) -> Int

Return the number of scans.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`scantimes`](@ref), [`minscantime`](@ref), [`maxscantime`](@ref).

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

Return the time of the first scan. The optional keyword argument `timeunit` allows you to 
change the unit of the returned `scantime`. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` allows you to specify whether the unit 
is stripped from the returned value.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`maxscantime`](@ref), [`scantimes`](@ref), [`scancount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

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

Return the time of the last scan. The optional keyword argument `timeunit` allows you to 
change the unit of the returned `scantime`. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` allows you to specify whether the unit 
is stripped from the returned value.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), 
[`minscantime`](@ref), [`scantimes`](@ref), [`scancount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

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
    ion(gcms::AbstractGCMS, index::Integer)

Return the `ion` at the specified `index`.

See also [`AbstractGCMS`](@ref), [`ions`](@ref), [`ionindex`](@ref), [`minion`](@ref), 
[`maxion`](@ref), [`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS((1:3)u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> ion(gcms, 1)
85

julia> ion(gcms, 2)
100

julia> ion(gcms, 3)
ERROR: BoundsError: attempt to access 2-element Vector{Int64} at index [3]
[...]
```
"""
function ion(gcms::AbstractGCMS, index::Integer)
    firstindex(ions(gcms)) ≤ index ≤ lastindex(ions(gcms)) || throw(
        BoundsError(ions(gcms)[index]))
    ions(gcms)[index]
end


"""
    ioncount(gcms::AbstractGCMS) -> Int

Return the number of `ions`.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`ions`](@ref), [`ion`](@ref), 
[`minion`](@ref), [`maxion`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> ioncount(gcms)
2
```
"""
ioncount(gcms::AbstractGCMS) = length(ions(gcms))


"""
    ionindex(gcms::AbstractGCMS, ion::Real) -> Int

Return the index of the `ion`. If the `ion` does not exist, an error is thrown.

See also [`AbstractGCMS`](@ref), [`ions`](@ref), [`ion`](@ref), [`minion`](@ref), 
[`maxion`](@ref), [`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS((1:3)u"s", [85.2f0, 100.1f0], [0 12; 34 956; 23 1]);

julia> ionindex(gcms, 100.1)
2

julia> ionindex(gcms, 201.1)
ERROR: ArgumentError: ion 201.1 does not exist
[...]
```
"""
function ionindex(gcms::AbstractGCMS, ion::Real)
    for (index, element) in enumerate(ions(gcms))
        element ≈ ion && return index
    end
    throw(ArgumentError("ion $ion does not exist"))
end


"""
    minion(gcms::AbstractGCMS)

Return the smallest `ion`.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`maxion`](@ref), [`ions`](@ref), 
[`ion`](@ref), [`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> minion(gcms)
85
```

"""
minion(gcms::AbstractGCMS) = first(ions(gcms))


"""
    maxion(gcms::AbstractGCMS)

Return the largest `ion`.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`minion`](@ref), [`ions`](@ref), 
[`ion`](@ref), [`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> maxion(gcms)
100
```
"""
maxion(gcms::AbstractGCMS) = last(ions(gcms))


"""
    minintensity(chrom::AbstractChromatogram)

Return the minimum `intensity`.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractGCMS`](@ref), 
[`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), [`maxintensity`](@ref), 
[`intensities`](@ref), [`scantimes`](@ref), [`scancount`](@ref), [`ions`](@ref), 
[`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> minintensity(gcms)
0
```
"""
minintensity(chrom::AbstractChromatogram) = minimum(intensities(chrom))


"""
    maxintensity(chrom::AbstractChromatogram)

Return the maximum `intensity`.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractGCMS`](@ref), 
[`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), [`minintensity`](@ref), 
[`intensities`](@ref), [`scantimes`](@ref), [`scancount`](@ref), [`ions`](@ref), 
[`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> maxintensity(gcms)
956
```
"""
maxintensity(chrom::AbstractChromatogram) = maximum(intensities(chrom))


"""
    runduration(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the `runduration`. The optional keyword argument `timeunit` allows you to change the 
unit of the returned time interval. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl)(e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` allows you to specify whether the unit 
is stripped from the returned value. 

See also [`AbstractChromatogram`](@ref), [`scantimes`](@ref), [`minscantime`](@ref), 
[`maxscantime`](@ref), [`scancount`](@ref).

# Example
```jldoctest
julia> fid = FID([30.1u"minute", 40.8u"minute", 51.5u"minute"], [12, 956, 23])
FID {scantimes: Float64, intensities: Int64}
3 scans; time range: 30.1 minute - 51.5 minute
intensity range: 12 - 956
metadata: 0 entries

julia> runduration(fid)
21.4 minute

julia> runduration(fid, timeunit=u"s")
1284.0 s

julia> runduration(fid, timeunit=u"s", ustripped=true)
1284.0
```
"""
function runduration(chrom::AbstractChromatogram; 
    timeunit::Unitful.TimeUnits=unit(eltype(chrom.scantimes)), ustripped::Bool=false)
    Δt = maxscantime(chrom) - minscantime(chrom)
    ustripped ? ustrip(timeunit, Δt) : uconvert(timeunit, Δt)
end


"""
    totalionchromatogram(gcms::GCMS)

Compute the total ion chromatrogram.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`GCMS`](@ref), 
[`AbstractTIC`](@ref), [`TIC`](@ref), [`scantimes`](@ref), [`intensities`](@ref), 
[`metadata`](@ref).

In the following example, the element type of `intensities` passed to the object 
constructor  is explicitly annotated to show that the `TIC` object has the same type.

# Example
```jldoctest
julia> gcms = GCMS((1.1:3.1)u"s", [85, 100], Int64[0 12; 34 956; 23 1]);

julia> tic = totalionchromatogram(gcms)
TIC {scantimes: Float64, intensities: Int64}
3 scans; time range: 1.1 s - 3.1 s
intensity range: 12 - 990
metadata: 0 entries

julia> intensities(tic)
3-element Vector{Int64}:
  12
 990
  24

julia> gcms = GCMS((1.1:3.1)u"s", [85, 100], Float64[0 12; 34 956; 23 1]);

julia> tic = totalionchromatogram(gcms)
TIC {scantimes: Float64, intensities: Float64}
3 scans; time range: 1.1 s - 3.1 s
intensity range: 12.0 - 990.0
metadata: 0 entries

julia> intensities(tic)
3-element Vector{Float64}:
  12.0
 990.0
  24.0
```
"""
totalionchromatogram(gcms::GCMS) = TIC(scantimes(gcms), vec(sum(intensities(gcms), 
    dims=2)), metadata(gcms))


function Base.show(io::IO, fid::FID)
    println(io, "FID {scantimes: ", eltype(ustrip.(scantimes(fid))), ", intensities: ", 
        eltype(intensities(fid)), "}")
    println(io, scancount(fid), " scans; time range: ", minscantime(fid), " - ", 
        maxscantime(fid))
    println(io, "intensity range: ", minintensity(fid), " - ", maxintensity(fid))
    n = length(metadata(fid))
    print(io, "metadata: ", n, (n == 0 || n > 1) ? " entries" : " entry" )
end


function Base.show(io::IO, gcms::GCMS)
    println(io, "GCMS {scantimes: ", eltype(ustrip.(scantimes(gcms))), ", ions: ", 
        eltype(ions(gcms)), ", intensities: ", eltype(intensities(gcms)), "}")
    println(io, scancount(gcms), " scans; time range: ", minscantime(gcms), " - ", 
        maxscantime(gcms))
    println(io, ioncount(gcms), " ions; range: m/z ", minion(gcms), " - ", maxion(gcms))
    println(io, "intensity range: ", minintensity(gcms), " - ", maxintensity(gcms))
    n = length(metadata(gcms))
    print(io, "metadata: ", n, (n == 0 || n > 1) ? " entries" : " entry" )
end


function Base.show(io::IO, tic::TIC)
    println(io, "TIC {scantimes: ", eltype(ustrip.(scantimes(tic))), ", intensities: ", 
        eltype(intensities(tic)), "}")
    println(io, scancount(tic), " scans; time range: ", minscantime(tic), " - ", 
        maxscantime(tic))
    println(io, "intensity range: ", minintensity(tic), " - ", maxintensity(tic))
    n = length(metadata(tic))
    print(io, "metadata: ", n, (n == 0 || n > 1) ? " entries" : " entry" )
end


############################################################################################
# Functions related to ion scan time shifts
############################################################################################

"""
    IonScanOrder

Supertype of all IonScanOrder implementations in JuMS.

See also [`LinearAscending`](@ref), [`LinearDescending`](@ref).
"""
abstract type IonScanOrder end
# , [`timeshift`](@ref)

struct LinearAscending{T1<:Real, T2<:Real} <: IonScanOrder
    start::T1
    stop::T2
    function LinearAscending{T1, T2}(start::T1, stop::T2) where {T1<:Real, T2<:Real}
        0 ≤ start < stop ≤ 1 || throw(ArgumentError(string("start=$start and stop=$stop do", 
            " not satisfy condition 0 ≤ start < stop ≤ 1")))
        new(start, stop)
    end
end


"""
    LinearAscending(; start::Real=0, stop::Real=1) <: IonScanOrder

Construct a `LinearAscending` ion scan order object. It specifies that the ions were scanned 
in linear ascending order (i.e., smallest ion first, last ion last) during each scan. The 
time to scan each ion is assumed to be equal and is the result of dividing the total scan 
interval time equally among the ions. The optional `start` and `stop` parameters allow you 
to limit the interval time during which the ions were scanned in each scan. They specify 
relative points in the scan interval: 0 ≤ start < stop ≤ 1. The default values are 
`start`=0 and `stop`=1, which means that the scan of the smallest ion started at the 
beginning of the scan interval and the scan of the largest ion ended at the end of the 
scan interval. In contrast, setting the start value to 0.5 would indicate that the ions 
were only scanned during the second half of the scan interval (e.g. because the instrument 
switched between SIM mode and Scan mode during each scan interval).

# Example
```jldoctest
julia> LinearAscending()
LinearAscending{Int64, Int64}(0, 1)

julia> LinearAscending(start=0.5)
LinearAscending{Float64, Int64}(0.5, 1)

julia> LinearAscending(start=0.1, stop=0.5)
LinearAscending{Float64, Float64}(0.1, 0.5)

julia> LinearAscending(start=0.5, stop=0.5)
ERROR: ArgumentError: start=0.5 and stop=0.5 do not satisfy condition 0 ≤ start < stop ≤ 1
[...]
```
"""
function LinearAscending(; start::T1=0, stop::T2=1) where {T1<:Real, T2<:Real}
    LinearAscending{T1, T2}(start, stop)
end

# See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`LinearDescending`](@ref), 
# [`timeshift`](@ref), [`ions`](@ref), [`minion`](@ref), [`maxion`](@ref), 
# [`ioncount`](@ref), [`meanscantime`](@ref), [`scantimes`](@ref), [`minscantime`](@ref), 
# [`maxscantime`](@ref).


struct LinearDescending{T1<:Real, T2<:Real} <: IonScanOrder
    start::T1
    stop::T2
    function LinearDescending{T1, T2}(start::T1, stop::T2) where {T1<:Real, T2<:Real}
        0 ≤ start < stop ≤ 1 || throw(ArgumentError(string("start=$start and stop=$stop do", 
            " not satisfy condition 0 ≤ start < stop ≤ 1")))
        new(start, stop)
    end
end


"""
    LinearDescending(; start::Real=0, stop::Real=1) <: IonScanOrder

Construct a `LinearDescending` ion scan order object. It specifies that the ions were 
scanned in linear descending order (i.e., largest ion first, smallest ion last) during each 
scan. The time to scan each ion is assumed to be equal and is the result of dividing the 
total scan interval time equally among the ions. The optional `start` and `stop` parameters 
allow you to limit the interval time during which the ions were scanned in each scan. They 
specify relative points in the scan interval: 0 ≤ start < stop ≤ 1. The default values are 
`start`=0 and `stop`=1, which means that the scan of the largest ion started at the 
beginning of the scan interval and the scan of the smallest ion ended at the end of the 
scan interval. In contrast, setting the start value to 0.5 would indicate that the ions 
were only scanned during the second half of the scan interval (e.g. because the instrument 
switched between SIM mode and Scan mode during each scan interval).

# Example
```jldoctest
julia> LinearDescending()
LinearDescending{Int64, Int64}(0, 1)

julia> LinearDescending(start=0.5)
LinearDescending{Float64, Int64}(0.5, 1)

julia> LinearDescending(start=0.1, stop=0.5)
LinearDescending{Float64, Float64}(0.1, 0.5)

julia> LinearDescending(start=0.5, stop=0.5)
ERROR: ArgumentError: start=0.5 and stop=0.5 do not satisfy condition 0 ≤ start < stop ≤ 1
[...]
```
"""
function LinearDescending(; start::T1=0, stop::T2=1) where {T1<:Real, T2<:Real}
    LinearDescending{T1, T2}(start, stop)
end

# See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`LinearAscending`](@ref), 
# [`timeshift`](@ref), [`ions`](@ref), [`minion`](@ref), [`maxion`](@ref), 
# [`ioncount`](@ref), [`meanscantime`](@ref), [`scantimes`](@ref), [`minscantime`](@ref), 
# [`maxscantime`](@ref).


# function timeshift(::LinearAscending, gcms::AbstractGCMS)
#     index::Integer -> begin
#         firstindex(ions(gcms)) ≤ index ≤ lastindex(ions(gcms)) || error("index $index outside of ion index range")
#         (index - ioncount(gcms)) * meanscantime(gcms) / ioncount(gcms)
#     end
# end
# # TYPE STABLE August 3, 2024


# function timeshift(::LinearDescending, gcms::AbstractGCMS)
#     index::Integer -> begin
#         firstindex(ions(gcms)) ≤ index ≤ lastindex(ions(gcms)) || error("index $index outside of ion index range")
#         (1 - index) * meanscantime(gcms) / ioncount(gcms)
#     end
# end
# # TYPE STABLE August 3, 2024


# """
#     timeshift(gcms::AbstractGCMS, ionscanorder::IonScanOrder)

# Return a function that returns the time shift of an ion relative to the scan timestamp if the ion index is provided
# and the `scanindex` and `ionindex` are specified. The function returned is determined by the ionscanorder type. If the
# ionscanorder is LinearAscending(), it is assumed that the ions are scanned in ascending order over the entire scan time
# interval. This means that the time shift relative to the timestamp is greatest for the smallest ion and zero for the largest
# ion. If the ionscanorder is LinearDescending(), it is assumed that the ions are scanned in descending order over the entire
# scan time interval. This means that the time shift relative to the scan timestamp is greatest for the largest ion and zero
# for the smallest ion. Both functions assume that the time for scanning each ion is evenly distributed over the average scan
# time of a scan in the AbstractGCMS object and that the timestamp is the time at which the scanning of the last ion was completed..

# # Example
# ```julia-repl
# julia> gcms = GCMS([1.0u"s", 2.0u"s", 3.0u"s"], [85, 100], [0 12; 34 956; 23 1])
# GCMS {scantimes: Quantity{Float64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}, ions: Int64, intensities: Int64}
# Source: nothing
# 3 scans; time range: 1.0 s - 3.0 s
# 2 ions; range: m/z 85 - 100
# intensity range: 0 - 956

# julia> δt = timeshift(gcms, LinearAscending())
# #60 (generic function with 1 method)

# julia> δt(1)
# -0.5 s

# julia> δt(2)
# 0.0 s

# julia> δt = timeshift(gcms, LinearDescending())
# #62 (generic function with 1 method)

# julia> δt(1)
# 0.0 s

# julia> δt(2)
# -0.5 s
# ```

# See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`scantimeindex`](@ref), [`ionindex`](@ref), [`timeshift`](@ref), [`IonScanOrder`](@ref), [`AscendingIons`](@ref), [`DescendingIons`](@ref).
# """
# timeshift(gcms::AbstractGCMS, ionscanorder::IonScanOrder) = timeshift(ionscanorder, gcms)
# # TYPE STABLE August 3, 2024