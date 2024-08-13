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

Construct a FID object consisting of `scantimes`, `intensities`, and `metadata`. Note that 
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
metadata: 0

julia> FID(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], Dict("name" => "sample"))
FID {scantimes: Int32, intensities: Float64}
3 scans; time range: 1 s - 3 s
intensity range: 1.0 - 956.0
metadata: 1

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

Construct a GCMS object consisting of `scantimes`, `ions`, `intensities`, and `metadata`. 
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
metadata: 0

julia> GCMS([1.1f0, 2.1f0]u"s", [35.1f0, 76.2f0], Int64[0 12; 34 956], Dict(:id => 1))
GCMS {scantimes: Float32, ions: Float32, intensities: Int64}
2 scans; time range: 1.1f0 s - 2.1f0 s
2 ions; range: m/z 35.1 - 76.2
intensity range: 0 - 956
metadata: 1

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

Construct a TIC object consisting of `scantimes`, `intensities`, and `metadata`. Note that 
the `scantimes` must be in ascending order and the `intensities` must not contain values 
less than zero.

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
metadata: 0

julia> TIC(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], Dict("name" => "sample"))
TIC {scantimes: Int32, intensities: Float64}
3 scans; time range: 1 s - 3 s
intensity range: 1.0 - 956.0
metadata: 1

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
ions(gcms::AbstractGCMS) = gcms.ions


"""
    scantime(chrom::AbstractChromatogram, index::Integer; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the `scantime` by specifying the scan `index`. The optional parameter 
`timeunit` allows you to change the unit of the returned time value. All time units defined 
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

Return the index of the scantime in the `scantimes` that is closest to the given `time`. 
All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl)
(e.g., `u"s"`, `u"minute") are supported. If there is a tie, the larger scan time is used. 
If the optional parameter `precisetime` is set to `true`, the specified time must exist in 
`scantimes`, otherwise an error is thrown.

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
metadata: 0

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

Return the intensities.

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

Return the metadata.

See also [`AbstractChromatogram`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref).

# Example
```jldoctest
julia> gcms₁ = GCMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956], Dict(:id => 1))
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
2 scans; time range: 1 s - 2 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 1

julia> metadata(gcms₁)
Dict{Any, Any} with 1 entry:
  :id => 1

julia> gcms₂ = GCMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956])
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
2 scans; time range: 1 s - 2 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 0

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
metadata: 2

julia> delete!(metadata(gcms₂), "name")
Dict{Any, Any} with 1 entry:
  :id => 123

julia> gcms₂
GCMS {scantimes: Int64, ions: Int64, intensities: Int64}
2 scans; time range: 1 s - 2 s
2 ions; range: m/z 85 - 100
intensity range: 0 - 956
metadata: 1
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
change the unit of the returned scan time. All time units defined in the package 
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
change the unit of the returned scan time. All time units defined in the package 
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

Return the number of ions.

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

Return the smallest ion.

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

Return the largest ion.

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

Return the minimum intensity.

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

Return the maximum intensity.

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
    totalionchromatogram(gcms::GCMS)

Compute the total ion chromatrogram.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`GCMS`](@ref), 
[`AbstractTIC`](@ref), [`TIC`](@ref), [`scantimes`](@ref), [`intensities`](@ref), 
[`metadata`](@ref).

In the following example, the element type of `intensities` passed to the object 
constructor  is explicitly annotated to show that the TIC object has the same type.

# Example
```jldoctest
julia> gcms = GCMS((1.1:3.1)u"s", [85, 100], Int64[0 12; 34 956; 23 1]);

julia> tic = totalionchromatogram(gcms)
TIC {scantimes: Float64, intensities: Int64}
3 scans; time range: 1.1 s - 3.1 s
intensity range: 12 - 990
metadata: 0

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
metadata: 0

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
    print(io, "metadata: ", length(metadata(fid)))
end


function Base.show(io::IO, gcms::GCMS)
    println(io, "GCMS {scantimes: ", eltype(ustrip.(scantimes(gcms))), ", ions: ", 
        eltype(ions(gcms)), ", intensities: ", eltype(intensities(gcms)), "}")
    println(io, scancount(gcms), " scans; time range: ", minscantime(gcms), " - ", 
        maxscantime(gcms))
    println(io, ioncount(gcms), " ions; range: m/z ", minion(gcms), " - ", maxion(gcms))
    println(io, "intensity range: ", minintensity(gcms), " - ", maxintensity(gcms))
    print(io, "metadata: ", length(metadata(gcms)))
end


function Base.show(io::IO, tic::TIC)
    println(io, "TIC {scantimes: ", eltype(ustrip.(scantimes(tic))), ", intensities: ", 
        eltype(intensities(tic)), "}")
    println(io, scancount(tic), " scans; time range: ", minscantime(tic), " - ", 
        maxscantime(tic))
    println(io, "intensity range: ", minintensity(tic), " - ", maxintensity(tic))
    print(io, "metadata: ", length(metadata(tic)))
end
