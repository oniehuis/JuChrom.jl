include("utilities.jl")

"""
    AbstractRiMapper

Supertype for all retention index mapper implementations. All subtypes (e.g., RiMapper) 
have fields to store the retention index name and metadata.

See also [`RiMapper`](@ref), [`retentionindexname`](@ref), [`metadata`](@ref).
"""
abstract type AbstractRiMapper end


"""
    AbstractChromatogram

Supertype for all chromatogram implementations. All subtypes (e.g., Chrom, ChromMS) 
include scan times and intensities and, optionally, metadata and a retention index mapper.

See also [`AbstractChrom`](@ref), [`AbstractChromMS`](@ref), [`ChromMS`](@ref), [`Chrom`](@ref), 
[`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref), [`rimapper`](@ref),.
"""
abstract type AbstractChromatogram end


"""
    AbstractChrom <: AbstractChromatogram

Supertype for all chromatogram implementations that have a single intensity value 
associated with each scan time (e.g., Chrom). The intensities are stored in a vector, with 
the index corresponding to the scan time index.

See also [`AbstractChromatogram`](@ref), [`Chrom`](@ref), [`scantimes`](@ref), 
[`intensities`](@ref), [`metadata`](@ref), [`rimapper`](@ref), .
"""
abstract type AbstractChrom <: AbstractChromatogram end


"""
    AbstractChromMS <: AbstractChromatogram

Supertype for all chromatogram implementations that include mass-charge ratio (*m*/*z*) 
data (ions) and associated abundance values (intensities) (e.g., ChromMS). This type can 
have one or more ion intensity values associated with a given scan time. Intensities are 
stored in a matrix where the row index represents the scan time and the column index 
represents the ion.

See also [`AbstractChromatogram`](@ref), [`AbstractChromMS`](@ref), [`ChromMS`](@ref), 
[`scantimes`](@ref), [`ions`](@ref), [`intensities`](@ref), [`metadata`](@ref), 
[`rimapper`](@ref).
"""
abstract type AbstractChromMS <: AbstractChromatogram end


mutable struct ChromMS{
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real},
    T3<:AbstractMatrix{<:Real},
    } <: AbstractChromMS
    const scantimes::T1
    const ions::T2
    const intensities::T3
    const metadata::Dict{Any, Any}
    rimapper::Union{AbstractRiMapper, Nothing}
    function ChromMS{T1, T2, T3}(
        scantimes::T1,
        ions::T2,
        intensities::T3,
        metadata::Dict, 
        rimapper::Union{AbstractRiMapper, Nothing}
        ) where {
            T1<:AbstractVector{<:Unitful.Time},
            T2<:AbstractVector{<:Real},
            T3<:AbstractMatrix{<:Real}}
        length(scantimes) > 0 || throw(ArgumentError("no scan time(s) provided"))
        length(intensities) > 0 || throw(ArgumentError("no intensity value(s) provided"))
        length(ions) > 0 || throw(ArgumentError("no ion(s) provided"))
        size(intensities, 1) == length(scantimes) || throw(DimensionMismatch(
            "intensity matrix row count does not match scan time count"))
        size(intensities, 2) == length(ions) || throw(DimensionMismatch(
            "intensity matrix column count does not match ion count"))
        length(Set(scantimes)) == length(scantimes) || throw(
            ArgumentError("scan times contain identical values"))
        issorted(scantimes) || throw(ArgumentError("scan times not in ascending order"))
        length(Set(ions)) == length(ions) || throw(
            ArgumentError("ions contain identical values"))
        issorted(ions) || throw(ArgumentError("ions not in ascending order"))
        new(scantimes, ions, intensities, metadata, rimapper)
    end
end


Base.broadcastable(chrom::ChromMS) = Ref(chrom)


"""
    ChromMS(scantimes::AbstractVector{<:Unitful.Time}, ions::AbstractVector{<:Real}, 
    intensities::AbstractMatrix{<:Real}; metadata::Dict=Dict{Any, Any}(), 
    rimapper::Union{AbstractRiMapper, Nothing}=nothing) <: AbstractChromMS

Construct a `ChromMS` object that includes `scantimes`, `ions`, `intensities`, and 
`metadata`. Ensure that both `scantimes` and `ions` are in ascending order. The scan times 
must include a time unit. All time units supported by the 
[Unitful.jl package](https://painterqubits.github.io/Unitful.jl) 
(e.g., `u"s"`, `u"minute"`) are accepted. You can optionally use the keyword argument 
`rimapper` to include a retention index mapper.

See also [`AbstractChromatogram`](@ref), [`AbstractChromMS`](@ref), [`scantimes`](@ref), 
[`ions`](@ref), [`intensities`](@ref), [`metadata`](@ref).

In the following examples, the types of `scantimes`, `ions`, and `intensities` are 
explicitely annotated to demonstrate that the `ChromMS` object preserves these types.

# Examples
```jldoctest
julia> ChromMS(Int32[1, 2, 3]u"s", Int64[85, 100], Int32[0 12; 34 956; 23 1])
ChromMS {scan times: Int32, ions: Int64, intensities: Int32}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> ChromMS([1.1f0, 2.1f0]u"s", [35.1f0, 76.2f0], Int64[0 12; 34 956], Dict(:id => 4))
ChromMS {scan times: Float32, ions: Float32, intensities: Int64}
2 scans; scan times: 1.1f0 s, 2.1f0 s
2 ions: m/z 35.1, 76.2
intensity range: 0 - 956
metadata: 1 entry

julia> ChromMS([2, 1, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
ERROR: ArgumentError: scan times not in ascending order
[...]

julia> ChromMS([1, 2, 3]u"s", [100, 85], [0 12; 34 956; 23 1])
ERROR: ArgumentError: ions not in ascending order
[...]
```
"""
function ChromMS(scantimes::T1, ions::T2, intensities::T3, metadata::Dict=Dict{Any, Any}();
    rimapper::Union{AbstractRiMapper, Nothing}=nothing
    ) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real},
    T3<:AbstractMatrix{<:Real}}
    ChromMS{T1, T2, T3}(scantimes, ions, intensities, metadata, rimapper)
end


mutable struct Chrom{
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real}
    } <: AbstractChrom
    const scantimes::T1
    const intensities::T2
    const metadata::Dict{Any, Any}
    rimapper::Union{AbstractRiMapper, Nothing}
    function Chrom{T1, T2}(
        scantimes::T1,
        intensities::T2, 
        metadata::Dict,
        rimapper::Union{AbstractRiMapper, Nothing}
        ) where { 
            T1<:AbstractVector{<:Unitful.Time}, 
            T2<:AbstractVector{<:Real}}
        length(scantimes) > 0 || throw(ArgumentError("no scan time(s) provided"))
        length(intensities) > 0 || throw(ArgumentError("no intensity value(s) provided"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan time count"))
        length(Set(scantimes)) == length(scantimes) || throw(
            ArgumentError("scan times contain identical values"))
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        new(scantimes, intensities, metadata, rimapper)
    end
end


Base.broadcastable(chrom::Chrom) = Ref(chrom)


"""
    Chrom(scantimes::AbstractVector{<:Unitful.Time}, 
    intensities::AbstractVector{<:Real}, metadata::Dict=Dict{Any, Any}; 
    rimapper::Union{AbstractRiMapper, Nothing}=nothing) <: AbstractChrom

Construct a `Chrom` object that includes `scantimes`, `intensities`, and `metadata`. 
Ensure that both `scantimes` and `ions` are in ascending order. The scan times must 
include a time unit. All time units supported by the 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) package (e.g., `u"s"`, 
`u"minute"`) are accepted. You can optionally use the keyword argument `rimapper` to 
include a retention index mapper.

See also [`AbstractChromatogram`](@ref), [`AbstractChrom`](@ref), [`scantimes`](@ref), 
[`intensities`](@ref), [`metadata`](@ref), [`rimapper`](@ref).

In the following examples, the types of `scantimes` and `intensities` are explicitely 
annotated to demonstrate that the `Chrom` object preserves these types.

# Examples
```jldoctest
julia> Chrom(Int64[1, 2, 3]u"s", Int32[12, 956, 1])
Chrom {scan times: Int64, intensities: Int32}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 1 - 956
metadata: 0 entries

julia> Chrom(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], Dict("name" => "sample"))
Chrom {scan times: Int32, intensities: Float64}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 1.0 - 956.0
metadata: 1 entry

julia> Chrom([2, 1, 3]u"s", [12.0, 956.0, 1.0])
ERROR: ArgumentError: scan times not in ascending order
[...]
```
"""
function Chrom(scantimes::T1, intensities::T2, metadata::Dict=Dict{Any, Any}();
    rimapper::Union{AbstractRiMapper, Nothing}=nothing) where {
        T1<:AbstractVector{<:Unitful.Time}, T2<:AbstractVector{<:Real}}
    Chrom{T1, T2}(scantimes, intensities, metadata, rimapper)
end


"""
    IonScanOrder

Supertype of all ion scan order implementations.

See also [`LinearAscending`](@ref), [`LinearDescending`](@ref), [`ionscantimeshift`](@ref), 
[`ionscantimes`](@ref), [`ionscantime`](@ref), [`ionscantimeindex`](@ref).
"""
abstract type IonScanOrder end


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

Construct a `LinearAscending` ion scan order object. This object specifies that ions were 
scanned in linear ascending order (i.e., smallest ion first, largest ion last) during each 
scan. The time allocated to scan each ion is assumed to be equal, resulting from dividing 
the total scan interval time equally among the ions. The optional `start` and `stop` 
parameters allow you to limit the time interval during which ions were scanned in each 
scan. These parameters specify relative points within the scan interval: `0 ≤ start < stop 
≤ 1`. The default values are `start=0` and `stop=1`, meaning the scan of the smallest ion 
began at the start of the scan interval and the scan of the largest ion concluded at the 
end of the interval. In contrast, setting the start value to 0.5 indicates that ions were 
scanned only during the second half of the scan interval. For instance, this could occur if 
the instrument switched between SIM mode and Scan mode during each scan interval, operating 
in Scan mode only during the second half, which generated the data in question.

See also [`AbstractChromMS`](@ref), [`LinearDescending`](@ref), [`ionscantimeshift`](@ref), 
[`ionscantimes`](@ref), [`ionscantime`](@ref), [`ionscantimeindex`](@ref), [`ions`](@ref), 
[`minion`](@ref), [`maxion`](@ref), [`ioncount`](@ref), [`scanduration`](@ref), 
[`scantimes`](@ref), [`minscantime`](@ref), [`maxscantime`](@ref).

# Examples
```jldoctest
julia> LinearAscending()
LinearAscending{Int64, Int64}(0, 1)

julia> LinearAscending(; start=0.5)
LinearAscending{Float64, Int64}(0.5, 1)

julia> LinearAscending(; start=0.1, stop=0.5)
LinearAscending{Float64, Float64}(0.1, 0.5)

julia> LinearAscending(; start=0.5, stop=0.5)
ERROR: ArgumentError: start=0.5 and stop=0.5 do not satisfy condition 0 ≤ start < stop ≤ 1
[...]
```
"""
function LinearAscending(; start::T1=0, stop::T2=1) where {T1<:Real, T2<:Real}
    LinearAscending{T1, T2}(start, stop)
end


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

Construct a `LinearDescending` ion scan order object. This object specifies that ions were 
scanned in linear descending order (i.e., largest ion first, smallest ion last) during each 
scan. The time allocated to scan each ion is assumed to be equal, resulting from dividing 
the total scan interval time equally among the ions. The optional `start` and `stop` 
parameters allow you to limit the time interval during which ions were scanned in each 
scan. These parameters specify relative points within the scan interval: `0 ≤ start < stop 
≤ 1`. The default values are `start=0` and `stop=1`, meaning the scan of the largest ion 
began at the start of the scan interval and the scan of the smallest ion concluded at the 
end of the interval. In contrast, setting the start value to 0.5 indicates that ions were 
scanned only during the second half of the scan interval. For instance, this could occur if 
the instrument switched between SIM mode and Scan mode during each scan interval, operating 
in Scan mode only during the second half, which generated the data in question.

See also [`AbstractChromMS`](@ref), [`LinearAscending`](@ref), [`ionscantimeshift`](@ref), 
[`ionscantimes`](@ref), [`ionscantime`](@ref), [`ionscantimeindex`](@ref), [`ions`](@ref), 
[`minion`](@ref), [`maxion`](@ref), [`ioncount`](@ref), [`scanduration`](@ref), 
[`scantimes`](@ref), [`minscantime`](@ref), [`maxscantime`](@ref).

# Examples
```jldoctest
julia> LinearDescending()
LinearDescending{Int64, Int64}(0, 1)

julia> LinearDescending(; start=0.5)
LinearDescending{Float64, Int64}(0.5, 1)

julia> LinearDescending(; start=0.1, stop=0.5)
LinearDescending{Float64, Float64}(0.1, 0.5)

julia> LinearDescending(; start=0.5, stop=0.5)
ERROR: ArgumentError: start=0.5 and stop=0.5 do not satisfy condition 0 ≤ start < stop ≤ 1
[...]
```
"""
function LinearDescending(; start::T1=0, stop::T2=1) where {T1<:Real, T2<:Real}
    LinearDescending{T1, T2}(start, stop)
end


"""
    binions(chrom::ChromMS; ionbin::Function=integer)

Return a ChromMS object in which the ions are binned according to the `ionbin` function (the 
default is the integer function), and the intensities of the binned ions are summed.

See also [`AbstractChromMS`](@ref), [`integer`](@ref), [`intensities`](@ref), [`ions`](@ref), 
[`ioncount`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS((1:3)u"s", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])
ChromMS {scan times: Int64, ions: Float64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
3 ions: m/z 84.8, 85.2, 100.9
intensity range: 0 - 956
metadata: 0 entries

julia> intensities(chrom)
3×3 Matrix{Int64}:
  0  24   12
  0   0  956
 23   0    1

julia> chrom_binnedions = binions(chrom)
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 101
intensity range: 0 - 956
metadata: 0 entries

julia> ions(chrom_binnedions)
2-element Vector{Int64}:
  85
 101

julia> intensities(chrom_binnedions)
3×2 Matrix{Int64}:
 24   12
  0  956
 23    1


julia> custom_ionbin(ion) = integer(ion, start=0.9);

julia> chrom_binnedions = binions(chrom, ionbin=custom_ionbin)
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
3 ions: m/z 84, 85, 101
intensity range: 0 - 956
metadata: 0 entries

julia> ions(chrom_binnedions)
3-element Vector{Int64}:
  84
  85
 101
```
"""
function binions(chrom::ChromMS; ionbin::Function=integer)
    binnedions = sort(unique(ionbin.(ions(chrom))))
    imₙ = zeros(eltype(intensities(chrom)), (size(intensities(chrom), 1), length(binnedions)))
    for (i, ion) in enumerate(ions(chrom))
        binnedion = ionbin(ion)
        imₙ[:, searchsortedlast(binnedions, binnedion)] .+= @view intensities(chrom)[:, i]
    end
    ChromMS(scantimes(chrom), binnedions, imₙ, metadata(chrom), rimapper=rimapper(chrom))
end


"""
    integer(value::Real; start::Real=0.7) -> Int

Return the integer for the given `value` that satisfies the following condition: 
`integer - 1 + start ≤ value < integer + start`, where `0 ≤ start < 1`.

See also [`AbstractChromMS`](@ref), [`ChromMS`](@ref), [`binions`](@ref), [`ions`](@ref).

# Examples
```jldoctest
julia> integer(29.7)
30

julia> integer(30.0)
30

julia> integer(30.69)
30

julia> integer(29.7, start=0.8)
29
```
"""
function integer(value::Real; start::Real=0.7)
    0 ≤ start < 1 || throw(ArgumentError(string("fractional digits of binning interval ",
        "start is outside the interval [0,1) of allowed values: $start")))
    start == 0 && return trunc(Int, value)
    trunc(Int, value + (1 - start))
end


"""
    intensities(chrom::AbstractChrom; scanindexrange::OrdinalRange{T, S}) where {T<:Integer, 
    S<:Integer}

Return the intensities. The optional keyword argument `scanindexrange` lets you select a 
subset of scans for which the intensities will be returned. Note that the function will 
return either a reference to the intensity vector or a view into the intensity vector, 
depending on whether a subset of scans is selected.

See also [`AbstractChrom`](@ref), [`intensity`](@ref), [`minintensity`](@ref), 
[`maxintensity`](@ref), [`scancount`](@ref).

# Examples
```jldoctest
julia> chrom = Chrom([1, 2, 3]u"s", [123, 224, 103])
Chrom {scan times: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 103 - 224
metadata: 0 entries

julia> intensities(chrom)  # reference to the data structure
3-element Vector{Int64}:
 123
 224
 103

julia> intensities(chrom)[:]  # a copy of these values
3-element Vector{Int64}:
 123
 224
 103

julia> intensities(chrom; scanindexrange=2:3)  # view into the data structure
2-element view(::Vector{Int64}, 2:3) with eltype Int64:
 224
 103

julia> intensities(chrom; scanindexrange=2:3)[:]  # a copy of these values
2-element Vector{Int64}:
 224
 103
```
"""
function intensities(
    chrom::AbstractChrom; scanindexrange::OrdinalRange{T, S}=firstindex(scantimes(
        chrom)):lastindex(scantimes(chrom))) where {T<:Integer, S<:Integer}

    if scanindexrange == firstindex(scantimes(chrom)):lastindex(scantimes(chrom))
        chrom.intensities
    else
        @view chrom.intensities[scanindexrange]
    end
end


"""
    intensities(chrom::AbstractChromMS; scanindexrange::OrdinalRange{T1, S1}, 
    ionindexrange::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, 
    S2<:Integer}

Return the intensities. The optional keyword arguments `scanindexrange` and `ionindexrange` 
allow you to select specific parts of the intensity matrix to be returned. Note that the 
function returns either a reference to the matrix or a view into it, depending on whether 
the keyword arguments specify subranges of the matrix.

See also [`AbstractChromMS`](@ref), [`scancount`](@ref), [`ioncount`](@ref), 
[`minintensity`](@ref), [`maxintensity`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> intensities(chrom)  # reference to the data structure
3×2 Matrix{Int64}:
  0   12
 34  956
 23    1

julia> intensities(chrom)[:, :]  # a copy of these values
3×2 Matrix{Int64}:
  0   12
 34  956
 23    1

julia> intensities(chrom, ionindexrange=1:1)  # all intensities of the ion at index 1
3×1 view(::Matrix{Int64}, 1:3, 1:1) with eltype Int64:
  0
 34
 23

julia> intensities(chrom, ionindexrange=1:1)[:]  # a copy of these values
3-element Vector{Int64}:
  0
 34
 23

julia> intensities(chrom, scanindexrange=1:1)  # intensities of all ions from scan 1
1×2 view(::Matrix{Int64}, 1:1, 1:2) with eltype Int64:
 0  12

julia> intensities(chrom, scanindexrange=1:1)[:]  # a copy of these values
2-element Vector{Int64}:
  0
 12

julia> intensities(chrom, scanindexrange=1:2, ionindexrange=1:2)
2×2 view(::Matrix{Int64}, 1:2, 1:2) with eltype Int64:
  0   12
 34  956

julia> intensities(chrom, scanindexrange=1:2, ionindexrange=1:2)[:, :]
2×2 Matrix{Int64}:
  0   12
 34  956
```
"""
function intensities(
    chrom::AbstractChromMS;
    scanindexrange::OrdinalRange{T1, S1}=firstindex(scantimes(chrom)):lastindex(
        scantimes(chrom)), 
    ionindexrange::OrdinalRange{T2, S2}=firstindex(ions(chrom)):lastindex(ions(chrom))
    ) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}

    if (scanindexrange == firstindex(scantimes(chrom)):lastindex(scantimes(chrom)) 
        && ionindexrange == firstindex(ions(chrom)):lastindex(ions(chrom)))
        chrom.intensities
    else
        @view intensities(chrom)[scanindexrange, ionindexrange]
    end
end


"""
    intensity(chrom::AbstractChrom, scanindex::Integer)

Return the intensity for a scan by specifying its `scanindex`.

See also [`AbstractChrom`](@ref), [`intensities`](@ref), [`minintensity`](@ref), 
[`maxintensity`](@ref).

# Examples
```jldoctest
julia> chrom = Chrom([1, 2, 3]u"s", [123, 224, 103])
Chrom {scan times: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 103 - 224
metadata: 0 entries

julia> intensity(chrom, 1)
123

julia> intensity(chrom, 2)
224
```
"""
intensity(chrom::AbstractChrom, scanindex::Integer) = intensities(chrom)[scanindex]


"""
    intensity(chrom::AbstractChrom, time::Unitful.Time; precisetime::Bool=false)

Return the intensity at a given `time`. All time units defined in the package
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. By default, the intensity associated with the scan whose timestamp is closest 
to the given `time` is returned. In case of a tie, the intensity of the scan with the 
later timestamp is used. If the optional parameter `precisetime` is set to true, the 
specified `time` must exactly match a timestamp in the vector; otherwise, an error is 
thrown.

See also [`AbstractChrom`](@ref), [`intensities`](@ref), [`scantimes`](@ref), 
[`minscantime`](@ref), [`maxscantime`](@ref).

# Examples
```jldoctest
julia> chrom = Chrom([1.0, 2.0, 3.0]u"s", [123, 224, 103])
Chrom {scan times: Float64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
intensity range: 103 - 224
metadata: 0 entries

julia> intensity(chrom, 1.5u"s")
224

julia> intensity(chrom, 1u"s", precisetime=true)
123

julia> intensity(chrom, 1.5u"s", precisetime=true)
ERROR: ArgumentError: scantime 1.5 s does not exist
[...]
```
"""
function intensity(chrom::AbstractChrom, time::Unitful.Time; precisetime::Bool=false)
    intensity(chrom, scantimeindex(chrom, time, precisetime=precisetime))
end


"""
    intensity(chrom::AbstractChromMS, scanindex::Integer, ionindex::Integer)

Return the intensity of an ion in a scan, given the `scanindex` of the scan and the 
`ionindex` of the ion.

See also [`AbstractChromMS`](@ref), [`scancount`](@ref), [`ions`](@ref), [`ioncount`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> intensity(chrom, 2, 1)
34

julia> intensity(chrom, 1, 2)
12
```
"""
intensity(chrom::AbstractChromMS, scanindex::Integer, ionindex::Integer) = intensities(
    chrom)[scanindex, ionindex]


"""
    intensity(chrom::AbstractChromMS, time::Unitful.Time, ion::Real; precisetime::Bool=false)

Return the intensity of an `ion` at a given `time`. All time units defined in the package
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. By default, the intensity associated with the scan whose timestamp is closest 
to the given `time` is returned. In case of a tie, the intensity of the scan with the 
later timestamp is used. If the optional parameter `precisetime` is set to true, the 
specified `time` must exactly match a timestamp in the vector; otherwise, an error is 
thrown.

See also [`AbstractChromMS`](@ref), [`intensities`](@ref), [`ions`](@ref), [`scantimes`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Float64, ions: Int64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> intensity(chrom, 1.9u"s", 85)
34

julia> intensity(chrom, 2.9u"s", 85, precisetime=true)
ERROR: ArgumentError: scantime 2.9 s does not exist
[...]

julia> intensity(chrom, 3u"s", 85, precisetime=true)
23
```
"""
function intensity(chrom::AbstractChromMS, time::Unitful.Time, ion::Real; 
    precisetime::Bool=false)
    intensity(chrom, scantimeindex(chrom, time; precisetime=precisetime), 
        ionindex(chrom, ion))
end


"""
    ion(chrom::AbstractChromMS, ionindex::Integer)

Return the ion at the specified `ionindex`.

See also [`AbstractChromMS`](@ref), [`ions`](@ref), [`ionindex`](@ref), [`minion`](@ref), 
[`maxion`](@ref), [`ioncount`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS((1:3)u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> ion(chrom, 1)
85

julia> ion(chrom, 2)
100

julia> ion(chrom, 3)
ERROR: BoundsError: attempt to access 2-element Vector{Int64} at index [3]
[...]
```
"""
ion(chrom::AbstractChromMS, ionindex::Integer) = ions(chrom)[ionindex]


"""
    ioncount(chrom::AbstractChromMS) -> Int

Return the number of ions.

See also [`AbstractChromMS`](@ref), [`ions`](@ref), [`ion`](@ref), [`minion`](@ref), 
[`maxion`](@ref).

# Example
```jldoctest
julia> chrom = ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> ioncount(chrom)
2
```
"""
ioncount(chrom::AbstractChromMS) = length(ions(chrom))


"""
    ionindex(chrom::AbstractChromMS, ion::Real) -> Int

Return the index of the specified `ion`. An error is thrown if the `ion` does not exist.

See also [`AbstractChromMS`](@ref), [`ions`](@ref), [`ioncount`](@ref), [`ion`](@ref), 
[`minion`](@ref), [`maxion`](@ref), [`ionscantime`](@ref), [`ionscantimeshift`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS((1:3)u"s", [85.2f0, 100.1f0], [0 12; 34 956; 23 1]);

julia> ionindex(chrom, 100.1)
2

julia> ionindex(chrom, 201.1)
ERROR: ArgumentError: ion 201.1 does not exist
[...]
```
"""
function ionindex(chrom::AbstractChromMS, ion::Real)
    for (index, element) in enumerate(ions(chrom))
        element ≈ ion && return index
    end
    throw(ArgumentError("ion $ion does not exist"))
end


"""
    ionscantime(δtᵢ::Function, chrom::AbstractChromMS, scanindex::Integer, ionindex::Integer; 
    timeunit::Unitful.TimeUnits, ustripped::Bool=false)

Return the time at which an ion was actually scanned, given the `scanindex`, `ionindex`, 
and a function `δtᵢ` that computes the time difference between the timestamp of a scan and 
the scan time of the ion from the `ionindex`. The optional parameter `timeunit` allows you 
to specify the unit of the returned scan time. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether the unit is 
included in the returned value. Note that the timestamp of a scan is assumed to be the time 
when the scanning of ion intensities associated with that scan was completed.

See also [`AbstractChromMS`](@ref), [`scantimes`](@ref), [`scantime`](@ref), 
[`scantimeindex`](@ref), [`ions`](@ref), [`ionindex`](@ref), [`ionscantimeshift`](@ref), 
[`IonScanOrder`](@ref), [`LinearAscending`](@ref), [`LinearDescending`](@ref).

# Examples
```julia-repl
julia> chrom = ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Float64, ions: Int64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> δtᵢ = ionscantimeshift(chrom, LinearDescending());

julia> ionscantime(δtᵢ, chrom, 2, 1)
2.0 s

julia> ionscantime(δtᵢ, chrom, 2, 2)
1.5 s

julia> ionscantime(δtᵢ, chrom, 2, 2; timeunit=u"minute")
0.025 minute

julia> ionscantime(δtᵢ, chrom, 2, 2; timeunit=u"minute", ustripped=true)
0.025
```
"""
function ionscantime(δtᵢ::Function, chrom::AbstractChromMS, scanindex::Integer, 
    ionindex::Integer; timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))),
    ustripped::Bool=false)
    t = scantime(chrom, scanindex) + δtᵢ(ionindex)
    ustripped ? ustrip(timeunit, t) : uconvert(timeunit, t)
end


"""
    ionscantimeindex(δtᵢ::Function, chrom::AbstractChromMS, ionindex::Integer, 
    time::Unitful.Time; precisetime::Bool=false) -> Int

Return the index of the scan where the scan time for the ion is closest to the specified 
`time`, given the `ionindex` and a function `δtᵢ` that computes the time difference between 
the timestamp of a scan and the scan time of the ion from the `ionindex`. All time units 
defined in the package [Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., 
`u"s"`, `u"minute"`) are supported. In case of a tie, the larger scan index is returned. 
If the optional parameter `precisetime` is set to true, the ion must have been scanned 
exactly at the specified `time`; otherwise, an error is thrown.

See also [`AbstractChromMS`](@ref), [`scantimeindex`](@ref), [`ionscantime`](@ref), 
[`ionscantimeshift`](@ref), [`IonScanOrder`](@ref), [`LinearAscending`](@ref), 
[`LinearDescending`](@ref), [`scantimes`](@ref), [`scantime`](@ref), [`ions`](@ref), 
[`ion`](@ref), [`ionindex`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> δtᵢ = ionscantimeshift(chrom, LinearDescending());

julia> ionscantime(δtᵢ, chrom, 2, 2)
1.5 s

julia> ionscantimeindex(δtᵢ, chrom, 2, 1.5u"s")
2

julia> ionscantimeindex(δtᵢ, chrom, 2, 1.6u"s")
2

julia> ionscantimeindex(δtᵢ, chrom, 2, 1.6u"s", precisetime=true)
ERROR: ArgumentError: ion has not been scanned at the time 1.6 s
[...]

julia> ionscantimeindex(δtᵢ, chrom, 2, 1.5u"s", precisetime=true)
2

julia> ionscantimeindex(δtᵢ, chrom, 1, 1.4u"s")
1

julia> ionscantime(δtᵢ, chrom, 1, 1)
1.0 s
```
"""
function ionscantimeindex(δtᵢ::Function, chrom::AbstractChromMS, ionindex::Integer, 
    time::Unitful.Time; precisetime::Bool=false)
    precisetime || return findclosest(ionscantimes(δtᵢ, chrom, ionindex), time)
    t = ustrip(uconvert(unit(eltype(scantimes(chrom))), time))
    for (index, element) in enumerate(ionscantimes(δtᵢ, chrom, ionindex, ustripped=true))
        element ≈ t && return index
    end
    throw(ArgumentError("ion has not been scanned at the time $time"))
end


"""
    ionscantimes(δtᵢ::Function, chrom::AbstractChromMS, ionindex::Integer; 
    timeunit::Unitful.TimeUnits, ustripped::Bool=false)

Return the times at which an ion was actually scanned, given the `ionindex` and a function 
`δtᵢ` that computes the time difference between the timestamp of a scan and the scan time 
of the ion from the `ionindex`. The optional parameter `timeunit` allows you to specify the 
unit for the returned scan times. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether to include the 
unit in the returned value. Note that the timestamp of a scan is assumed to be the time 
when the scanning of ion intensities associated with that scan was completed.

See also [`AbstractChromMS`](@ref), [`ionscantime`](@ref), [`ionscantimeshift`](@ref), 
[`IonScanOrder`](@ref), [`LinearAscending`](@ref), [`LinearDescending`](@ref), 
[`scantimes`](@ref), [`scantimeindex`](@ref), [`ions`](@ref), [`ionindex`](@ref).

# Examples
```julia-repl
julia> chrom = ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Float64, ions: Int64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> δtᵢ = ionscantimeshift(chrom, LinearDescending());

julia> ionscantimes(δtᵢ, chrom, 1)
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 1.0 s
 2.0 s
 3.0 s

julia> ionscantimes(δtᵢ, chrom, 2)
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 0.5 s
 1.5 s
 2.5 s

julia> ionscantimes(δtᵢ, chrom, 2; timeunit=u"minute")
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:
 0.008333333333333333 minute
                0.025 minute
 0.041666666666666664 minute

julia> ionscantimes(δtᵢ, chrom, 2; timeunit=u"minute", ustripped=true)
3-element Vector{Float64}:
 0.008333333333333333
 0.025
 0.041666666666666664
```
"""
function ionscantimes(δtᵢ::Function, chrom::AbstractChromMS, ionindex::Integer; 
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))), ustripped::Bool=false)
    ts = scantimes(chrom) .+ δtᵢ(ionindex)
    ustripped ? ustrip.(timeunit, ts) : uconvert.(timeunit, ts)
end


function ionscantimeshift(ionscanorder::LinearAscending, chrom::AbstractChromMS, error::Real)
    intervalsize = scanduration(chrom, error=error)
    index::Integer -> begin
        firstindex(ions(chrom)) ≤ index ≤ lastindex(ions(chrom)) || throw(
            BoundsError(ions(chrom), index))
        intervalsize * ((ionscanorder.stop - 1) + ((index - ioncount(chrom))
            * (ionscanorder.stop - ionscanorder.start)) / ioncount(chrom))
    end
end


function ionscantimeshift(ionscanorder::LinearDescending, chrom::AbstractChromMS, error::Real)
    intervalsize = scanduration(chrom, error=error)
    index::Integer -> begin
        firstindex(ions(chrom)) ≤ index ≤ lastindex(ions(chrom)) || throw(
            BoundsError(ions(chrom), index))
        intervalsize * ((ionscanorder.stop - 1) + ((1 - index) 
            * (ionscanorder.stop - ionscanorder.start)) / ioncount(chrom))
    end
end


"""
    ionscantimeshift(chrom::AbstractChromMS, ionscanorder::IonScanOrder; error::Real=0.001)

Return a function, based on the `ionscanorder`, that calculates the time difference between 
the timestamp of a scan and the time when an ion was actually scanned, given the ion index 
as an argument. The time difference will be zero or negative, since the timestamp of a scan 
is considered to be when the scanning of the last ion was completed. The returned function 
assumes that the duration of each scan is consistent throughout the run. The optional 
keyword argument `error` lets you specify the maximum allowed deviation of the scan 
duration, as a fraction of the average scan time, between the timestamps of two consecutive 
scans.

See also [`AbstractChromMS`](@ref), [`IonScanOrder`](@ref), [`LinearAscending`](@ref), 
[`LinearDescending`](@ref), [`ionscantime`](@ref), [`ions`](@ref), [`minion`](@ref), 
[`maxion`](@ref), [`ioncount`](@ref), [`scanduration`](@ref), [`scantimes`](@ref), 
[`minscantime`](@ref), [`maxscantime`](@ref).

# Examples
```julia-repl
julia> chrom = ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Float64, ions: Int64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> δtᵢ = ionscantimeshift(chrom, LinearAscending());

julia> δtᵢ(1)
-0.5 s

julia> δtᵢ(2)
0.0 s

julia> δtᵢ = ionscantimeshift(chrom, LinearDescending());

julia> δtᵢ(1)
0.0 s

julia> δtᵢ(2)
-0.5 s

julia> δtᵢ = ionscantimeshift(chrom, LinearDescending(start=0.5));

julia> δtᵢ(1)
0.0 s

julia> δtᵢ(2)
-0.25 s

julia> δtᵢ = ionscantimeshift(chrom, LinearDescending(stop=0.5));

julia> δtᵢ(1)
-0.5 s

julia> δtᵢ(2)
-0.75 s

julia> δtᵢ = ionscantimeshift(chrom, LinearDescending(start=0.25, stop=0.75));

julia> δtᵢ(1)
-0.25 s

julia> δtᵢ(2)
-0.5 s
```
"""
ionscantimeshift(chrom::AbstractChromMS, ionscanorder::IonScanOrder, error::Real=0.01
    ) = ionscantimeshift(ionscanorder, chrom, error)


"""
    ions(chrom::AbstractChromMS)

Return the ions.

See also [`AbstractChromMS`](@ref), [`ioncount`](@ref), [`ion`](@ref), [`minion`](@ref), 
[`maxion`](@ref), [`ionindex`](@ref), [`ionscantime`](@ref), [`ionscantimeshift`](@ref). 

In the following examples, the type of `ions` is explicitely annotated to demonstrate that 
the ChromMS object preserves this type.

# Examples
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", Int64[85, 100], [0 12; 34 956; 23 1]);

julia> ions(chrom)
2-element Vector{Int64}:
  85
 100

julia> chrom = ChromMS([1, 2, 3]u"s", Float32[85, 100], [0 12; 34 956; 23 1]);

julia> ions(chrom)
2-element Vector{Float32}:
  85.0
 100.0
```
"""
ions(chrom::AbstractChromMS) = chrom.ions


"""
    maxintensity(chrom::AbstractChrom[; scanindexrange::OrdinalRange{T, S}]) 
    where {T<:Integer, S<:Integer})

Return the minimum intensity. The optional keyword argument `scanindexrange` allows you to 
specify a range of scans from which the minimum intensity is determined.

See also [`AbstractChrom`](@ref), [`minintensity`](@ref), [`intensities`](@ref), 
[`intensity`](@ref).

# Examples
```jldoctest
julia> chrom = Chrom([1, 2, 3]u"s", [12, 1, 956]);

julia> maxintensity(chrom)
956

julia> maxintensity(chrom, scanindexrange=1:2)
12
```
"""
function maxintensity(chrom::AbstractChrom; scanindexrange::OrdinalRange{T, S}=(
    firstindex(scantimes(chrom)):lastindex(scantimes(chrom)))) where {
    T<:Integer, S<:Integer}
    maximum(intensities(chrom, scanindexrange=scanindexrange))
end


"""
    maxintensity(chrom::AbstractChromMS[; scanindexrange::OrdinalRange{T1, S1}, 
    ionindexrange::OrdinalRange{T2, S2}]) where {T1<:Integer, S1<:Integer, 
    T2<:Integer, S2<:Integer})

Return the maximum intensity. The optional keyword arguments `scanindexrange` and 
`ionindexrange` allow you to select a range of scans and ions from which the minimum 
intensity is determined.

See also [`AbstractChromMS`](@ref), [`minintensity`](@ref), [`intensities`](@ref), 
[`intensity`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 1; 23 956]);

julia> maxintensity(chrom)
956

julia> maxintensity(chrom, scanindexrange=1:2)
34

julia> maxintensity(chrom, ionindexrange=1:1)
34

julia> maxintensity(chrom, scanindexrange=1:2, ionindexrange=2:2)
12
```
"""
function maxintensity(chrom::AbstractChromMS; scanindexrange::OrdinalRange{T1, S1}=(
    firstindex(scantimes(chrom)):lastindex(scantimes(chrom))), 
    ionindexrange::OrdinalRange{T2, S2}=(firstindex(ions(chrom)):lastindex(ions(chrom)))
    ) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
    maximum(intensities(chrom, scanindexrange=scanindexrange, ionindexrange=ionindexrange))
end


"""
    maxion(chrom::AbstractChromMS)

Return the largest ion.

See also [`AbstractChromMS`](@ref), [`minion`](@ref), [`ions`](@ref), [`ion`](@ref), 
[`ioncount`](@ref).

# Example
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> maxion(chrom)
100
```
"""
maxion(chrom::AbstractChromMS) = last(ions(chrom))


"""
    maxscantime(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{<:Integer, 
    <:Integer}]; timeunit::Unitful.TimeUnits, ustripped::Bool)

Return the maximum scan time. The optional second positional argument allows you to 
specify the scan range for which the maximum scan time is returned. The optional keyword 
argument `timeunit` lets you change the unit of the returned scan time. All time units 
defined in the package [Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., 
u"s", u"minute") are supported. The optional keyword argument `ustripped` lets you choose 
whether to include the unit in the returned value.

See also [`AbstractChromatogram`](@ref), [`minscantime`](@ref), [`scantimes`](@ref), 
[`scantime`](@ref),[`scancount`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> maxscantime(chrom)
3.0 s

julia> maxscantime(chrom, timeunit=u"minute")
0.05 minute

julia> maxscantime(chrom, timeunit=u"minute", ustripped=true)
0.05

julia> maxscantime(chrom, 1:2)
2.0 s

julia> maxscantime(chrom, 1:2, timeunit=u"minute", ustripped=true)
0.03333333333333333
```
"""
function maxscantime(chrom::AbstractChromatogram,
    scanindexrange::OrdinalRange{T, S}=(firstindex(scantimes(chrom)):lastindex(
        scantimes(chrom))); 
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))), 
    ustripped::Bool=false
    ) where {T<:Integer, S<:Integer}

    t = last(@view scantimes(chrom)[scanindexrange])
    ustripped ? ustrip(timeunit, t) : uconvert(timeunit, t)
end


"""
    metadata(chrom::AbstractChromatogram) -> Dict{Any, Any}

Return the metadata.

See also [`AbstractChromatogram`](@ref).

# Examples
```jldoctest
julia> chrom₁ = ChromMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956], Dict(:id => 4))
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
2 scans; scan times: 1 s, 2 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 1 entry

julia> metadata(chrom₁)
Dict{Any, Any} with 1 entry:
  :id => 4

julia> chrom₂ = ChromMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956])
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
2 scans; scan times: 1 s, 2 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> metadata(chrom₂)
Dict{Any, Any}()

julia> metadata(chrom₂)["name"] = "sample"
"sample"

julia> metadata(chrom₂)[:id] = 123
123

julia> metadata(chrom₂)
Dict{Any, Any} with 2 entries:
  "name" => "sample"
  :id    => 123

julia> chrom₂
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
2 scans; scan times: 1 s, 2 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 2 entries

julia> delete!(metadata(chrom₂), "name")
Dict{Any, Any} with 1 entry:
  :id => 123

julia> chrom₂
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
2 scans; scan times: 1 s, 2 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 1 entry
```
"""
metadata(chrom::AbstractChromatogram) = chrom.metadata


"""
    minintensity(chrom::AbstractChrom[, greaterthan::Real; 
    scanindexrange::OrdinalRange{T, S}]) where {T<:Integer, S<:Integer})

Return the minimum intensity. The optional positional argument `greaterthan` allows you to 
specify a threshold value; the returned minimum intensity will be greater than this 
threshold. If no intensity exceeds the specified value, nothing is returned. The optional 
keyword argument `scanindexrange` allows you to specify a range of scans from which the 
minimum intensity is determined.

See also [`AbstractChrom`](@ref), [`maxintensity`](@ref), [`intensities`](@ref), 
[`intensity`](@ref).

# Examples
```jldoctest
julia> chrom = Chrom([1, 2, 3]u"s", [12, 956, 1]);

julia> minintensity(chrom)
1

julia> minintensity(chrom, 1)
12

julia> minintensity(chrom, scanindexrange=1:2)
12

julia> minintensity(chrom, 1, scanindexrange=2:3)
956
```
"""
function minintensity(chrom::AbstractChrom; scanindexrange::OrdinalRange{T, S}=(
    firstindex(scantimes(chrom)):lastindex(scantimes(chrom)))) where {
    T<:Integer, S<:Integer}
    minimum(intensities(chrom, scanindexrange=scanindexrange))
end


function minintensity(chrom::AbstractChrom, greaterthan::Real; 
    scanindexrange::OrdinalRange{T, S}=(firstindex(scantimes(chrom)):lastindex(
    scantimes(chrom)))) where {T<:Integer, S<:Integer}
    filtered_intensities = filter(num -> num > greaterthan, intensities(chrom, 
        scanindexrange=scanindexrange))
    length(filtered_intensities) > 0 ? minimum(filtered_intensities) : nothing
end


"""
    minintensity(chrom::AbstractChromMS[, greaterthan::Real; 
    scanindexrange::OrdinalRange{T1, S1}, ionindexrange::OrdinalRange{T2, S2}]) 
    where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer})

Return the minimum intensity. The optional positional argument `greaterthan` allows you to 
specify a threshold value; the returned minimum intensity will be greater than this 
threshold. If no intensity exceeds the specified value, nothing is returned. The optional 
keyword arguments `scanindexrange` and `ionindexrange` allow you to select a range of scans
and ions from which the minimum intensity is determined.

See also [`AbstractChromMS`](@ref), [`maxintensity`](@ref), [`intensities`](@ref), 
[`intensity`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> minintensity(chrom)
0

julia> minintensity(chrom, 0)
1

julia> minintensity(chrom, scanindexrange=2:3)
1

julia> minintensity(chrom, ionindexrange=1:1)
0

julia> minintensity(chrom, scanindexrange=2:3, ionindexrange=1:1)
23

julia> minintensity(chrom, 25, scanindexrange=2:3, ionindexrange=1:1)
34
```
"""
function minintensity(chrom::AbstractChromMS; scanindexrange::OrdinalRange{T1, S1}=(
    firstindex(scantimes(chrom)):lastindex(scantimes(chrom))), 
    ionindexrange::OrdinalRange{T2, S2}=(firstindex(ions(chrom)):lastindex(ions(chrom)))
    ) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
    minimum(intensities(chrom, scanindexrange=scanindexrange, ionindexrange=ionindexrange))
end


function minintensity(chrom::AbstractChromMS, greaterthan::Real; 
    scanindexrange::OrdinalRange{T1, S1}=(firstindex(scantimes(chrom)):lastindex(
    scantimes(chrom))), ionindexrange::OrdinalRange{T2, S2}=(firstindex(
    ions(chrom)):lastindex(ions(chrom)))) where {T1<:Integer, S1<:Integer, T2<:Integer, 
    S2<:Integer}
    filtered_intensities = filter(num -> num > greaterthan, intensities(chrom, 
        scanindexrange=scanindexrange, ionindexrange=ionindexrange))
    length(filtered_intensities) > 0 ? minimum(filtered_intensities) : nothing
end


"""
    minion(chrom::AbstractChromMS)

Return the smallest ion.

See also [`AbstractChromMS`](@ref), [`maxion`](@ref), [`ions`](@ref), [`ion`](@ref), 
[`ioncount`](@ref).

# Example
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> minion(chrom)
85
```

"""
minion(chrom::AbstractChromMS) = first(ions(chrom))


"""
    minscantime(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{<:Integer, 
    <:Integer}}]; timeunit::Unitful.TimeUnits, ustripped::Bool=false)

Return the minimum scan time. The optional second positional argument `scanindexrange` 
allows you to specify a scan range for which the minimum scan time is returned. The 
optional keyword argument `timeunit` lets you change the unit of the returned scan time. 
All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether to include the 
unit in the returned value.

See also [`AbstractChromatogram`](@ref), [`maxscantime`](@ref), [`scantimes`](@ref), 
[`scantime`](@ref), [`scancount`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> minscantime(chrom)
1.0 s

julia> minscantime(chrom, timeunit=u"minute")
0.016666666666666666 minute

julia> minscantime(chrom, timeunit=u"minute", ustripped=true)
0.016666666666666666

julia> minscantime(chrom, 2:3)
2.0 s

julia> minscantime(chrom, 2:3, timeunit=u"minute", ustripped=true)
0.03333333333333333
```
"""
function minscantime(chrom::AbstractChromatogram, 
    scanindexrange::OrdinalRange{T, S}=(firstindex(scantimes(chrom)):lastindex(
        scantimes(chrom)));
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))), 
    ustripped::Bool=false
    ) where {T<:Integer, S<:Integer}

    t = first(@view scantimes(chrom)[scanindexrange])
    ustripped ? ustrip(timeunit, t) : uconvert(timeunit, t)
end


"""
    rimapper(chrom::AbstractChromatogram)

Return the retention index mapper object. If no object is stored, the function returns 
nothing.

See also [`AbstractChromatogram`](@ref), [`AbstractRiMapper`](@ref), [`RiMapper`](@ref).

# Example
```jldoctest
julia> chrom = Chrom([1, 2, 3, 4, 5]u"s", [12, 956, 23, 45, 25]);

julia> rimapper(chrom) === nothing
true

julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000);

julia> chrom = Chrom([1, 2, 3, 4, 5]u"s", [12, 956, 23, 45, 25], rimapper=ld);

julia> rimapper(chrom)
RiMapper {index name: Kovats, calibration points: 5}
retention times: 1 minute, 2 minute, 3 minute, 4 minute, 5 minute
retention indices: 1000, 2000, 3000, 4000, 5000
interpolation method: NaturalCubicBSpline(false)
extrapolation method: nothing
metadata: 0 entries
```
"""
rimapper(chrom::AbstractChromatogram) = chrom.rimapper


"""
    rimapper!(chrom::AbstractChromatogram, rim::AbstractRiMapper)

Assign an retention index mapper to the AbstractChromatogram object.

See also [`AbstractChromatogram`](@ref), [`AbstractRiMapper`](@ref), [`RiMapper`](@ref).

# Example
```jldoctest
julia> chrom = Chrom([1, 2, 3, 4, 5]u"s", [12, 956, 23, 45, 25]);

julia> rimapper(chrom) === nothing
true

julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000);

julia> rimapper!(chrom, ld);

julia> rimapper(chrom)
RiMapper {index name: Kovats, calibration points: 5}
retention times: 1 minute, 2 minute, 3 minute, 4 minute, 5 minute
retention indices: 1000, 2000, 3000, 4000, 5000
interpolation method: NaturalCubicBSpline(false)
extrapolation method: nothing
metadata: 0 entries

julia> retentionindex(chrom, 2.2u"minute") ≈ 2200.0
true
```
"""
rimapper!(chrom::AbstractChromatogram, rim::AbstractRiMapper) = chrom.rimapper = rim


"""
    runduration(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the duration of the run. The optional keyword argument `timeunit` lets you specify 
the unit for the returned time interval. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether to include 
the unit in the returned value.

See also [`AbstractChromatogram`](@ref), [`scantimes`](@ref), [`minscantime`](@ref), 
[`maxscantime`](@ref), [`scancount`](@ref).

# Examples
```jldoctest
julia> chrom = Chrom([30.1u"minute", 40.8u"minute", 51.5u"minute"], [12, 956, 23])
Chrom {scan times: Float64, intensities: Int64}
3 scans; scan times: 30.1 minute, 40.8 minute, 51.5 minute
intensity range: 12 - 956
metadata: 0 entries

julia> runduration(chrom)
21.4 minute

julia> runduration(chrom, timeunit=u"s")
1284.0 s

julia> runduration(chrom, timeunit=u"s", ustripped=true)
1284.0
```
"""
function runduration(chrom::AbstractChromatogram; 
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))), ustripped::Bool=false)
    Δt = maxscantime(chrom) - minscantime(chrom)
    ustripped ? ustrip(timeunit, Δt) : uconvert(timeunit, Δt)
end


"""
    scancount(chrom::AbstractChromatogram) -> Int

Return the number of scans.

See also [`AbstractChromatogram`](@ref), [`scantimes`](@ref).

# Example
```jldoctest
julia> chrom = Chrom([1, 2, 3]u"s", [12, 956, 23]);

julia> scancount(chrom)
3
```
"""
scancount(chrom::AbstractChromatogram) = length(scantimes(chrom))


"""
    scanduration(chrom::AbstractChromatogram; error::Real=0.01, 
    timeunit::Unitful.TimeUnits, ustripped::Bool=false)

Calculates the periodicity with which the scans were recorded over time. The optional 
keyword argument `error` allows you to specify the maximum allowable deviation of the 
time interval between consecutive scans from the average scan time, as a fraction of 
the average scan time. The optional keyword argument `timeunit` lets you specify the unit 
for the returned value. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether the unit is 
include in the returned value.

See also [`AbstractChromatogram`](@ref), [`scantimes`](@ref), [`minscantime`](@ref), 
[`maxscantime`](@ref), [`scancount`](@ref), [`runduration`](@ref).

# Examples
```jldoctest
julia> chrom = Chrom([1.0, 2.0, 3.0]u"s", [12, 956, 1]);

julia> scanduration(chrom)
1.0 s

julia> scanduration(chrom, timeunit=u"minute")
0.016666666666666666 minute

julia> scanduration(chrom, timeunit=u"minute", ustripped=true)
0.016666666666666666

julia> scanduration(Chrom([1.0, 1.99, 3.0]u"s", [12, 956, 1]))
ERROR: ArgumentError: maximum scan duration variation above threshold: 0.010000000000000009 > 0.01
[...]

julia> scanduration(Chrom([1.0, 1.99, 3.0]u"s", [12, 956, 1]), error=0.02)
1.0 s
```
"""
function scanduration(chrom::AbstractChromatogram; error::Real=0.01,
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))), ustripped::Bool=false)
    scancount(chrom) > 1 || throw(
        ArgumentError("cannot calculate the scan duration from a single scan"))
    Δts = Set{eltype(scantimes(chrom))}()
    first = true
    for i in eachindex(scantimes(chrom))
        first && (first = false; continue)
        push!(Δts, scantime(chrom, i) - scantime(chrom, i - 1))
    end
    mean = runduration(chrom) / (scancount(chrom) - 1)
    error_obs = max(abs(maximum(Δts) - mean), abs(minimum(Δts) - mean)) / mean
    error_obs > error && throw(ArgumentError(string("maximum scan duration ", 
        "variation above threshold: $error_obs > $error")))
    ustripped ? ustrip.(timeunit, mean) : uconvert.(timeunit, mean)
end


"""
    scantime(chrom::AbstractChromatogram, scanindex::Integer; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the scan time for a given `scanindex`. The optional parameter `timeunit` lets you 
specify the unit for the returned scan time. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether the unit is 
included in the returned value.

See also [`AbstractChromatogram`](@ref), [`scantimes`](@ref), [`minscantime`](@ref), 
[`maxscantime`](@ref), [`scancount`](@ref), [`ionscantime`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> scantime(chrom, 2)
2 s

julia> scantime(chrom, 2, timeunit=u"minute")
1//30 minute

julia> scantime(chrom, 2, timeunit=u"minute", ustripped=true)
1//30
```
"""
function scantime(chrom::AbstractChromatogram, scanindex::Integer; 
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))), ustripped::Bool=false)
    ustripped ? ustrip(timeunit, scantimes(chrom)[scanindex]) : uconvert(timeunit, 
        scantimes(chrom)[scanindex])
end


"""
    scantimeindex(chrom::AbstractChromatogram, time::Unitful.Time; 
    precisetime::Bool=false) -> Int

Return the index of the scan time closest to `time` in the scan times. All time units 
defined in the package [Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., 
`u"s"`, `u"minute"`) are supported. In case of a tie, the larger scan time is returned. If 
the optional parameter `precisetime` is set to true, the specified time must exactly match 
a scan time value; otherwise, an error is thrown.

See also [`AbstractChromatogram`](@ref), [`scantimes`](@ref), [`scantime`](@ref), 
[`minscantime`](@ref), [`maxscantime`](@ref), [`scancount`](@ref), [`ionindex`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1.1f0, 2.1f0, 3.1f0]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Float32, ions: Int64, intensities: Int64}
3 scans; scan times: 1.1f0 s, 2.1f0 s, 3.1f0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> scantimeindex(chrom, 1.1f0u"s", precisetime=true)
1

julia> scantimeindex(chrom, 2.1u"s", precisetime=true)
2

julia> scantimeindex(chrom, 2.2u"s", precisetime=true)
ERROR: ArgumentError: scantime 2.2 s does not exist
[...]

julia> scantimeindex(chrom, 2.2u"s")
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
    scantimes(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{<:Integer, 
    <:Integer}]; timeunit::Unitful.TimeUnits, ustripped::Bool=false)

Return the scan times. The optional second positional argument allows you to specify a 
range of scan indices. The optional keyword argument `timeunit` lets you change the unit of 
the returned scan times. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether to include the 
unit in the returned values. Note: If no time unit conversion is applied and the unit is 
not stripped, the function returns a reference to or a view into the data structure.

See also [`AbstractChromatogram`](@ref), [`scantime`](@ref), [`minscantime`](@ref), 
[`maxscantime`](@ref), [`scancount`](@ref).

In the following example, the type of `scantimes` is explicitely annotated to demonstrate 
that the `ChromMS` object preserves this type.

# Example
```jldoctest
julia> chrom = ChromMS(Float32[1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> scantimes(chrom)  # reference to the data structure
3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 1.0f0 s
 2.0f0 s
 3.0f0 s

julia> scantimes(chrom, timeunit=u"minute")
3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:
 0.016666668f0 minute
 0.033333335f0 minute
 0.050000004f0 minute

julia> scantimes(chrom, timeunit=u"minute", ustripped=true)
3-element Vector{Float32}:
 0.016666668
 0.033333335
 0.050000004

julia> scantimes(chrom, 2:3)  # view into the data structure
2-element view(::Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}, 2:3) with eltype Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}:
 2.0f0 s
 3.0f0 s

julia> scantimes(chrom, 2:3)[:]  # a copy of these values
2-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 2.0f0 s
 3.0f0 s

julia> scantimes(chrom, 2:3, timeunit=u"minute", ustripped=true)
2-element Vector{Float32}:
 0.033333335
 0.050000004
```
"""
function scantimes(chrom::AbstractChromatogram;
    timeunit::Unitful.TimeUnits=unit(eltype(chrom.scantimes)), ustripped::Bool=false)
    convert = unit(eltype(chrom.scantimes)) ≠ timeunit
    if !convert && !ustripped
        chrom.scantimes
    elseif !convert && ustripped
        ustrip.(chrom.scantimes)
    elseif convert && !ustripped
        uconvert.(timeunit, chrom.scantimes)
    else
        ustrip.(timeunit, chrom.scantimes)
    end
end


function scantimes(chrom::AbstractChromatogram, scanindexrange::OrdinalRange{T, S}; 
    timeunit::Unitful.TimeUnits=unit(eltype(chrom.scantimes)), ustripped::Bool=false
    ) where {T<:Integer, S<:Integer}
    convert = unit(eltype(chrom.scantimes)) ≠ timeunit
    if !convert && !ustripped
        @view chrom.scantimes[scanindexrange]
    elseif !convert && ustripped
        ustrip.(@view chrom.scantimes[scanindexrange])
    elseif convert && !ustripped
        uconvert.(timeunit, @view chrom.scantimes[scanindexrange])
    else
        ustrip.(timeunit, @view chrom.scantimes[scanindexrange])
    end
end


"""
    totalionchromatogram(chrom::ChromMS)

Compute the total ion chromatrogram.

See also [`ChromMS`](@ref), [`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref).

In the following example, the type of the intensities is explicitely annotated to 
demonstrate that the returned Chrom object preserves this type.

# Example
```jldoctest
julia> chrom = ChromMS((1.1f0:3.1f0)u"s", [85, 100], Int64[0 12; 34 956; 23 1])
ChromMS {scan times: Float32, ions: Int64, intensities: Int64}
3 scans; scan times: 1.1f0 s, 2.1f0 s, 3.1f0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> tic = totalionchromatogram(chrom)
Chrom {scan times: Float32, intensities: Int64}
3 scans; scan times: 1.1f0 s, 2.1f0 s, 3.1f0 s
intensity range: 12 - 990
metadata: 0 entries

julia> intensities(tic)
3-element Vector{Int64}:
  12
 990
  24
```
"""
totalionchromatogram(chrom::ChromMS) = Chrom(scantimes(chrom), vec(sum(intensities(chrom), 
    dims=2)), metadata(chrom), rimapper=rimapper(chrom))


function Base.show(io::IO, chrom::AbstractChromatogram)
    # Information about the element types
    typename = name(typeof(chrom))
    print(io, "$typename {scan times: ", eltype(ustrip.(scantimes(chrom)))) 
    isa(chrom, AbstractChromMS) && print(io, ", ions: ", eltype(ions(chrom)))
    println(io, ", intensities: ", eltype(intensities(chrom)), "}")

    # Information about the scan count and the scan times
    if scancount(chrom) > 10
        println(io, scancount(chrom), " scans; scan time range: ", minscantime(chrom), 
            " - ", maxscantime(chrom))
    elseif 1 < scancount(chrom) ≤ 10
        print(io, scancount(chrom), " scans; scan times: ")
        for i in eachindex(scantimes(chrom))
            print(io, scantime(chrom, i))
            i < scancount(chrom) ? print(io, ", ") : println(io)
        end
    else
        println(io, scancount(chrom), " scan; scan time: ", scantime(chrom, 1))
    end

    # Information about the ions
    if isa(chrom, AbstractChromMS)
        if ioncount(chrom) > 10
            println(io, ioncount(chrom), " ions; range: m/z ", minion(chrom), " - ", 
                maxion(chrom))
        elseif 1 < ioncount(chrom) ≤ 10
            print(io, ioncount(chrom), " ions: m/z ")
            for (i, ion) in enumerate(ions(chrom))
                print(io, ion)
                i < ioncount(chrom) ? print(io, ", ") : println(io)
            end
        else
            println(io, ioncount(chrom), " ion: m/z ", ion(chrom, 1))
        end
    end

    # Information about the intensities
    if length(intensities(chrom)) > 1
        println(io, "intensity range: ", minintensity(chrom), " - ", 
            maxintensity(chrom))
    else
        println(io, "intensity: ", first(intensities(chrom)))
    end

    # Information about the metadata
    n = length(metadata(chrom))
    print(io, "metadata: ", n, (n == 0 || n > 1) ? " entries" : " entry" )

    # ADD INFORMATION AN RIMAPPER
    if rimapper(chrom) !== nothing
        print(io, "\nretention index mapper: ", retentionindexname(rimapper(chrom)))
    end
end
