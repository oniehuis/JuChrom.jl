include("utilities.jl")

"""
    AbstractChromatogram

Supertype for all chromatogram implementations. All subtypes (e.g., FID, GCMS, TIC) 
include scan times, intensities, and metadata.

See also [`AbstractGC`](@ref), [`AbstractFID`](@ref), [`AbstractGCMS`](@ref), 
[`AbstractTIC`](@ref), [`FID`](@ref), [`GCMS`](@ref), [`TIC`](@ref), [`scantimes`](@ref), 
[`intensities`](@ref), [`metadata`](@ref).
"""
abstract type AbstractChromatogram end



"""
    JuChrom.RetentionIndexStyle

Supertype for traits related to the availability of retention index data.

See also [`HasRetentionIndexData`](@ref), [`HasNoRetentionIndexData`](@ref).
"""
abstract type RetentionIndexStyle end


"""
    JuChrom.HasRetentionIndexData <: JuChrom.RetentionIndexStyle

Data type that functions as a trait indicating the availability of retention index data.

See also [`RetentionIndexStyle`](@ref), [`HasNoRetentionIndexData`](@ref), 
[`RiFID`](@ref), [`retentionindexname`](@ref), [`retentionindices`](@ref), 
[`maxretentionindex`](@ref), [`minretentionindex`](@ref).
"""
struct HasRetentionIndexData <: RetentionIndexStyle end


"""
    JuChrom.HasNoRetentionIndexData <: JuChrom.RetentionIndexStyle

Data type that functions as a trait indicating the absence of retention index data.

See also [`RetentionIndexStyle`](@ref), [`HasNoRetentionIndexData`](@ref).
"""
struct HasNoRetentionIndexData <: RetentionIndexStyle end


"""
    JuChrom.RetentionIndexStyle(::Type)

Return the RetentionIndexStyle trait assigned to the specified data type.

See also [`AbstractChromatogram`](@ref), [`HasRetentionIndexData`](@ref), 
[`HasNoRetentionIndexData`](@ref).

# Examples
```jldoctest
julia> tic = TIC([1, 2, 3]u"s", [12, 956, 1]);

julia> JuChrom.RetentionIndexStyle(typeof(tic))
JuChrom.HasNoRetentionIndexData()

julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", Int64[100, 200, 300], [12, 956, 1]);

julia> JuChrom.RetentionIndexStyle(typeof(tic))
JuChrom.HasRetentionIndexData()
```
"""
RetentionIndexStyle(::Type) = HasNoRetentionIndexData()


"""
    AbstractGC <: AbstractChromatogram

Supertype for all chromatogram implementations that lack mass-charge ratio (*m*/*z*) data 
(ions) and thus have a single intensity value associated with each scan time (e.g., FID, 
TIC). The intensities are stored in a vector, with the index corresponding to the scan time 
index.

See also [`AbstractChromatogram`](@ref), [`AbstractFID`](@ref), [`AbstractTIC`](@ref), 
[`FID`](@ref), [`TIC`](@ref), [`scantimes`](@ref), [`intensities`](@ref), 
[`metadata`](@ref).
"""
abstract type AbstractGC <: AbstractChromatogram end


"""
    AbstractFID <: AbstractGC

Supertype for all flame ionization detector chromatogram implementations (e.g., FID).

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractFID`](@ref),  
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
        length(scantimes) > 0 || throw(ArgumentError("no scan time(s) provided"))
        length(intensities) > 0 || throw(ArgumentError("no intensity value(s) provided"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan time count"))
        length(Set(scantimes)) == length(scantimes) || throw(
            ArgumentError("scan times contain identical values"))
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain at least one value less than zero"))
        new(scantimes, intensities, metadata)
    end
end


Base.broadcastable(fid::FID) = Ref(fid)


"""
    FID(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, 
    metadata::Dict=Dict{Any, Any})

Construct a `FID` object that includes `scantimes`, `intensities`, and `metadata`. Note 
that `scantimes` must be in ascending order, and `intensities` must not contain any values 
less than zero.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractFID`](@ref), 
[`RiFID`](@ref), [`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref).

In the following examples, the types of `scantimes` and `intensities` are explicitely 
annotated to demonstrate that the `FID` object preserves these types.

# Examples
```jldoctest
julia> FID(Int64[1, 2, 3]u"s", Int32[12, 956, 1])
FID {scan times: Int64, intensities: Int32}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 1 - 956
metadata: 0 entries

julia> FID(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], Dict("name" => "sample"))
FID {scan times: Int32, intensities: Float64}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 1.0 - 956.0
metadata: 1 entry

julia> FID([2, 1, 3]u"s", [12.0, 956.0, 1.0])
ERROR: ArgumentError: scan times not in ascending order
[...]

julia> FID([1, 2, 3]u"s", [-12.0, 956.0, 1.0])
ERROR: ArgumentError: intensity values contain at least one value less than zero
[...]
```
"""
function FID(scantimes::T1, intensities::T2, metadata::Dict=Dict{Any, Any}()) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractVector{<:Real}}
    FID{T1, T2}(scantimes, intensities, metadata)
end


struct RiFID{
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractString,
    T3<:AbstractVector{<:Union{Real, Missing}},
    T4<:AbstractVector{<:Real}} <: AbstractFID
    scantimes::T1
    retentionindexname::T2
    retentionindices::T3
    intensities::T4
    metadata::Dict{Any, Any}
    function RiFID{T1, T2, T3, T4}(scantimes::T1, retentionindexname::T2, 
        retentionindices::T3, intensities::T4, metadata::Dict) where {
        T1<:AbstractVector{<:Unitful.Time},
        T2<:AbstractString,
        T3<:AbstractVector{<:Union{Real, Missing}},
        T4<:AbstractVector{<:Real}}
        length(scantimes) > 0 || throw(ArgumentError("no scan time(s) provided"))
        length(retentionindexname) > 0 || throw(ArgumentError(
            "no retention index name provided"))
        numerical_ris = collect(skipmissing(retentionindices))
        length(numerical_ris) > 0 || throw(ArgumentError(
            "no retention index value(s) provided"))
        length(Set(numerical_ris)) == length(numerical_ris) || throw(ArgumentError(
            "retention indices contain identical values"))
        issorted(numerical_ris) || throw(
            ArgumentError("retention indices not in ascending order"))
        minval = minimum(numerical_ris)
        minidx = findfirst(x -> !ismissing(x) && x == minval, retentionindices)
        maxval = maximum(numerical_ris)
        maxidx = findfirst(x -> !ismissing(x) && x == maxval, retentionindices)
        count(x -> ismissing(x), retentionindices[minidx:maxidx]) == 0 ||
            throw(ArgumentError("numerical retention index values are discontinuous"))
        length(intensities) > 0 || throw(ArgumentError("no intensity value(s) provided"))
        length(retentionindices) == length(scantimes) || throw(
            DimensionMismatch("retention index count does not match scan time count"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan time count"))
        length(Set(scantimes)) == length(scantimes) || throw(
            ArgumentError("scan times contain identical values"))
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain at least one value less than zero"))
        new(scantimes, retentionindexname, retentionindices, intensities, metadata)
    end
end


Base.broadcastable(fid::RiFID) = Ref(fid)


RetentionIndexStyle(::Type{<:RiFID}) = HasRetentionIndexData()


"""
    RiFID(scantimes::AbstractVector{<:Unitful.Time}, retentionindexname::AbstractString, 
    retentionindices::{<:Union{Real, Missing}}, intensities::AbstractVector{<:Real}, 
    metadata::Dict=Dict{Any, Any})

Create an `RiFID` object that includes `scantimes`, `retentionindexname`, 
`retentionindices`, `intensities`, and `metadata`. The `retentionindices` array may have 
missing values at either end, but it must contain at least one numerical entry, and all 
numerical entries must be continuous.  Ensure that `scantimes` and `retentionindices` are 
both in ascending order and have the same length. Additionally, ensure that all values in 
`intensities` are non-negative.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractFID`](@ref), 
[`scantimes`](@ref), [`retentionindexname`](@ref), [`retentionindices`](@ref), 
[`intensities`](@ref), [`metadata`](@ref).

# Examples
```jldoctest
julia> RiFID(Int32[1, 2, 3]u"s", "Kovats", Int64[100, 200, 300], Int64[12, 956, 1])
RiFID {scan times: Int32, retention indices: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
3 retention indices: 100 (≘ 1 s), 200 (≘ 2 s), 300 (≘ 3 s)
intensity range: 1 - 956
metadata: 0 entries

julia> RiFID(Int32[1, 2, 3]u"s", "Kovats", [missing, 200.0, 300.0], Int64[12, 956, 1])
RiFID {scan times: Int32, retention indices: Union{Missing, Float64}, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 retention indices: 200.0 (≘ 2 s), 300.0 (≘ 3 s)
intensity range: 1 - 956
metadata: 0 entries

julia> RiFID([1, 2, 3]u"s", "Kovats", [missing, 300.0, 200.0], [12, 956, 1])
ERROR: ArgumentError: retention indices not in ascending order
[...]

julia> RiFID([1, 2, 3]u"s", "Kovats", [100.0, missing, 300.0], [12, 956, -1])
ERROR: ArgumentError: numerical retention index values are discontinuous
[...]

julia> RiFID([1, 2, 3]u"s", "Kovats", [100.0, 200.0, 300.0], [12, 956, -1])
ERROR: ArgumentError: intensity values contain at least one value less than zero
[...]
```
"""
function RiFID(scantimes::T1, retentionindexname::T2, retentionindices::T3,  
    intensities::T4, metadata::Dict=Dict{Any, Any}()) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractString,
    T3<:AbstractVector{<:Union{Real, Missing}},
    T4<:AbstractVector{<:Real}}
    RiFID{T1, T2, T3, T4}(scantimes, retentionindexname, retentionindices, intensities, 
        metadata)
end


"""
    AbstractTIC <: AbstractGC

Supertype for all total ion chromatogram implementations (e.g., `TIC`).

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
        length(scantimes) > 0 || throw(ArgumentError("no scan time(s) provided"))
        length(intensities) > 0 || throw(ArgumentError("no intensity value(s) provided"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan time count"))
        length(Set(scantimes)) == length(scantimes) || throw(
            ArgumentError("scan times contain identical values"))
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain at least one value less than zero"))
        new(scantimes, intensities, metadata)
    end
end


Base.broadcastable(tic::TIC) = Ref(tic)


"""
    TIC(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, 
    metadata::Dict=Dict{Any, Any}())

Construct a `TIC` object that includes `scantimes`, `intensities`, and `metadata`. Note 
that `scantimes` must be in ascending order, and `intensities` must not contain any values 
less than zero.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractTIC`](@ref), 
[`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref), 
[`totalionchromatogram`](@ref).

In the following examples, the types of `scantimes` and `intensities` are explicitely 
annotated to demonstrate that the `TIC` object preserves these types.

# Examples
```jldoctest
julia> TIC(Int64[1, 2, 3]u"s", Int32[12, 956, 1])
TIC {scan times: Int64, intensities: Int32}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 1 - 956
metadata: 0 entries

julia> TIC(Int32[1, 2, 3]u"s", Float64[12.0, 956.0, 1.0], Dict("name" => "sample"))
TIC {scan times: Int32, intensities: Float64}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 1.0 - 956.0
metadata: 1 entry

julia> TIC([2, 1, 3]u"s", [12.0, 956.0, 1.0])
ERROR: ArgumentError: scan times not in ascending order
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


struct RiTIC{
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractString,
    T3<:AbstractVector{<:Union{Real, Missing}},
    T4<:AbstractVector{<:Real}} <: AbstractTIC
    scantimes::T1
    retentionindexname::T2
    retentionindices::T3
    intensities::T4
    metadata::Dict{Any, Any}
    function RiTIC{T1, T2, T3, T4}(scantimes::T1, retentionindexname::T2, 
        retentionindices::T3, intensities::T4, metadata::Dict) where {
        T1<:AbstractVector{<:Unitful.Time},
        T2<:AbstractString,
        T3<:AbstractVector{<:Union{Real, Missing}},
        T4<:AbstractVector{<:Real}}
        length(scantimes) > 0 || throw(ArgumentError("no scan time(s) provided"))
        length(retentionindexname) > 0 || throw(ArgumentError(
            "no retention index name provided"))
        numerical_ris = collect(skipmissing(retentionindices))
        length(numerical_ris) > 0 || throw(ArgumentError(
            "no retention index value(s) provided"))
        length(Set(numerical_ris)) == length(numerical_ris) || throw(
            ArgumentError("retention indices contain identical values"))
        issorted(numerical_ris) || throw(ArgumentError(
            "retention indices not in ascending order"))
        minval = minimum(numerical_ris)
        minidx = findfirst(x -> !ismissing(x) && x == minval, retentionindices)
        maxval = maximum(numerical_ris)
        maxidx = findfirst(x -> !ismissing(x) && x == maxval, retentionindices)
        count(x -> ismissing(x), retentionindices[minidx:maxidx]) == 0 ||
            throw(ArgumentError("numerical retention index values are discontinuous"))
        length(intensities) > 0 || throw(ArgumentError("no intensity value(s) provided"))
        length(retentionindices) == length(scantimes) || throw(
            DimensionMismatch("retention index count does not match scan time count"))
        length(intensities) == length(scantimes) || throw(
            DimensionMismatch("intensity count does not match scan time count"))
        length(scantimes) == length(scantimes) || throw(
            ArgumentError("scan times contain identical values"))
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain at least one value less than zero"))
        new(scantimes, retentionindexname, retentionindices, intensities, metadata)
    end
end


Base.broadcastable(fid::RiTIC) = Ref(fid)


RetentionIndexStyle(::Type{<:RiTIC}) = HasRetentionIndexData()


"""
    RiTIC(scantimes::AbstractVector{<:Unitful.Time}, retentionindexname::AbstractString, 
    retentionindices::{<:Union{Real, Missing}}, intensities::AbstractVector{<:Real}, 
    metadata::Dict=Dict{Any, Any})

Create an `RiTIC` object that includes `scantimes`, `retentionindexname`, 
`retentionindices`, `intensities`, and `metadata`. The `retentionindices` array may have 
missing values at either end, but it must contain at least one numerical entry, and all 
numerical entries must be continuous.  Ensure that `scantimes` and `retentionindices` are 
both in ascending order and have the same length. Additionally, ensure that all values in 
`intensities` are non-negative.

See also [`AbstractChromatogram`](@ref), [`AbstractGC`](@ref), [`AbstractTIC`](@ref), 
[`scantimes`](@ref), [`retentionindexname`](@ref), [`retentionindices`](@ref), 
[`intensities`](@ref), [`metadata`](@ref).

# Examples
```jldoctest
julia> RiTIC(Int32[1, 2, 3]u"s", "Kovats", Int64[100, 200, 300], Int64[12, 956, 1])
RiTIC {scan times: Int32, retention indices: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
3 retention indices: 100 (≘ 1 s), 200 (≘ 2 s), 300 (≘ 3 s)
intensity range: 1 - 956
metadata: 0 entries

julia> RiTIC(Int32[1, 2, 3]u"s", "Kovats", [missing, 200.0, 300.0], Int64[12, 956, 1])
RiTIC {scan times: Int32, retention indices: Union{Missing, Float64}, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 retention indices: 200.0 (≘ 2 s), 300.0 (≘ 3 s)
intensity range: 1 - 956
metadata: 0 entries

julia> RiTIC([1, 2, 3]u"s", "Kovats", [missing, 300.0, 200.0], [12, 956, 1])
ERROR: ArgumentError: retention indices not in ascending order
[...]

julia> RiTIC([1, 2, 3]u"s", "Kovats", [100.0, missing, 300.0], [12, 956, 1])
ERROR: ArgumentError: numerical retention index values are discontinuous
[...]

julia> RiTIC([1, 2, 3]u"s", "Kovats", [100.0, 200.0, 300.0], [12, 956, -1])
ERROR: ArgumentError: intensity values contain at least one value less than zero
[...]
```
"""
function RiTIC(scantimes::T1, retentionindexname::T2, retentionindices::T3,  
    intensities::T4, metadata::Dict=Dict{Any, Any}()) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractString,
    T3<:AbstractVector{<:Union{Real, Missing}},
    T4<:AbstractVector{<:Real}}
    RiTIC{T1, T2, T3, T4}(scantimes, retentionindexname, retentionindices, intensities, 
        metadata)
end


"""
    AbstractGCMS <: AbstractChromatogram

Supertype for all chromatogram implementations that include mass-charge ratio (*m*/*z*) 
data (ions) and associated abundance values (intensities) (e.g., GCMS). This type can have 
one or more ion intensity values associated with a given scan time. Intensities are stored 
in a matrix where the row index represents the scan time and the column index represents 
the ion.

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
        count(i -> i < 0, intensities) == 0 || throw(ArgumentError(
            "intensity values contain at least one value less than zero"))
        new(scantimes, ions, intensities, metadata)
    end
end


Base.broadcastable(gcms::GCMS) = Ref(gcms)


"""
    GCMS(scantimes::AbstractVector{<:Unitful.Time}, ions::AbstractVector{<:Real}, 
    intensities::AbstractMatrix{<:Real}; metadata::Dict=Dict{Any, Any}()) 

Construct a `GCMS` object that includes `scantimes`, `ions`, `intensities`, and `metadata`. 
Note that `scantimes` and `ions` must be in ascending order, and `intensities` must not 
contain any values less than zero.

See also [`AbstractChromatogram`](@ref), [`AbstractGCMS`](@ref), [`scantimes`](@ref), 
[`ions`](@ref), [`intensities`](@ref), [`metadata`](@ref).

In the following examples, the types of `scantimes`, `ions`, and `intensities` are 
explicitely annotated to demonstrate that the `GCMS` object preserves these types.

# Examples
```jldoctest
julia> GCMS(Int32[1, 2, 3]u"s", Int64[85, 100], Int32[0 12; 34 956; 23 1])
GCMS {scan times: Int32, ions: Int64, intensities: Int32}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> GCMS([1.1f0, 2.1f0]u"s", [35.1f0, 76.2f0], Int64[0 12; 34 956], Dict(:id => 4))
GCMS {scan times: Float32, ions: Float32, intensities: Int64}
2 scans; scan times: 1.1f0 s, 2.1f0 s
2 ions: m/z 35.1, 76.2
intensity range: 0 - 956
metadata: 1 entry

julia> GCMS([2, 1, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
ERROR: ArgumentError: scan times not in ascending order
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


struct RiGCMS{
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractString,
    T3<:AbstractVector{<:Union{Real, Missing}},
    T4<:AbstractVector{<:Real},
    T5<:AbstractMatrix{<:Real}} <: AbstractGCMS
    scantimes::T1
    retentionindexname::T2
    retentionindices::T3
    ions::T4
    intensities::T5
    metadata::Dict{Any, Any}
    function RiGCMS{T1, T2, T3, T4, T5}(scantimes::T1, retentionindexname::T2, 
        retentionindices::T3, ions::T4, intensities::T5, metadata::Dict) where {
        T1<:AbstractVector{<:Unitful.Time},
        T2<:AbstractString,
        T3<:AbstractVector{<:Union{Real, Missing}},
        T4<:AbstractVector{<:Real},
        T5<:AbstractMatrix{<:Real}}
        length(scantimes) > 0 || throw(ArgumentError("no scan time(s) provided"))
        length(retentionindexname) > 0 || throw(ArgumentError(
            "no retention index name provided"))
        numerical_ris = collect(skipmissing(retentionindices))
        length(numerical_ris) > 0 || throw(ArgumentError(
            "no retention index value(s) provided"))
        length(Set(numerical_ris)) == length(numerical_ris) || throw(
            ArgumentError("retention indices contain identical values"))
        issorted(numerical_ris) || throw(ArgumentError(
            "retention indices not in ascending order"))
        minval = minimum(numerical_ris)
        minidx = findfirst(x -> !ismissing(x) && x == minval, retentionindices)
        maxval = maximum(numerical_ris)
        maxidx = findfirst(x -> !ismissing(x) && x == maxval, retentionindices)
        count(x -> ismissing(x), retentionindices[minidx:maxidx]) == 0 ||
            throw(ArgumentError("numerical retention index values are discontinuous"))
        length(ions) > 0 || throw(ArgumentError("no ion(s) provided"))
        length(intensities) > 0 || throw(ArgumentError("no intensity value(s) provided"))
        length(retentionindices) == length(scantimes) || throw(
            DimensionMismatch("retention index count does not match scan time count"))
        size(intensities, 1) == length(scantimes) || throw(DimensionMismatch(
            "intensity matrix row count does not match scan time count"))
        size(intensities, 2) == length(ions) || throw(DimensionMismatch(
            "intensity matrix column count does not match ion count"))
        length(Set(scantimes)) == length(scantimes) || throw(
            ArgumentError("scan times contain identical values"))
        issorted(scantimes) || throw(
            ArgumentError("scan times not in ascending order"))
        length(Set(ions)) == length(ions) || throw(
            ArgumentError("ions contain identical values"))
        issorted(ions) || throw(ArgumentError("ions not in ascending order"))
        count(i -> i < 0, intensities) == 0 || throw(
            ArgumentError("intensity values contain at least one value less than zero"))
        new(scantimes, retentionindexname, retentionindices, ions, intensities, metadata)
    end
end


Base.broadcastable(gcms::RiGCMS) = Ref(gcms)


RetentionIndexStyle(::Type{<:RiGCMS}) = HasRetentionIndexData()


"""
    RiGCMS(scantimes::AbstractVector{<:Unitful.Time}, retentionindexname::AbstractString, 
    retentionindices::{<:Union{Real, Missing}}, ions::AbstractVector{<:Real}, 
    intensities::AbstractVector{<:Real}, metadata::Dict=Dict{Any, Any})

Create an `RiTIC` object that includes `scantimes`, `retentionindexname`, 
`retentionindices`, `ions`, `intensities`, and `metadata`. The `retentionindices` array may 
have missing values at either end, but it must contain at least one numerical entry, and 
all numerical entries must be continuous.  Ensure that `scantimes`, 
`retentionindices`, and `ions` are all in ascending order. Ensure that `scantimes` and 
`retentionindices` have the same length. Additionally, ensure that all values in 
`intensities` are non-negative.

See also [`AbstractChromatogram`](@ref), [`AbstractGCMS`](@ref), [`scantimes`](@ref), 
[`retentionindexname`](@ref), [`retentionindices`](@ref), [`ions`](@ref), 
[`intensities`](@ref), [`metadata`](@ref).

# Examples
```jldoctest
julia> RiGCMS([1.0, 2.0]u"s", "Kovats", [10.0, 20.0], Int64[85, 100], Int64[0 12; 34 6])
RiGCMS {scan times: Float64, retention indices: Float64, ions: Int64, intensities: Int64}
2 scans; scan times: 1.0 s, 2.0 s
2 retention indices: 10.0 (≘ 1.0 s), 20.0 (≘ 2.0 s)
2 ions: m/z 85, 100
intensity range: 0 - 34
metadata: 0 entries

julia> RiGCMS([1.0, 2.0]u"s", "Kovats", [missing, 10.0], Int64[85, 100], Int64[0 12; 34 6])
RiGCMS {scan times: Float64, retention indices: Union{Missing, Float64}, ions: Int64, intensities: Int64}
2 scans; scan times: 1.0 s, 2.0 s
1 retention index: 10.0 (≘ 2.0 s)
2 ions: m/z 85, 100
intensity range: 0 - 34
metadata: 0 entries

julia> RiGCMS([1, 2, 3]u"s", "Kovats", [200, 100, 300], [85, 100], [0 12; 34 6; 23 1]);
ERROR: ArgumentError: retention indices not in ascending order
[...]

julia> RiGCMS([1, 2, 3]u"s", "Kovats", [100, missing, 300], [85, 100], [0 12; 34 6; 23 1]);
ERROR: ArgumentError: numerical retention index values are discontinuous
[...]

julia> RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], [0 12; 34 6; 23 -1]);
ERROR: ArgumentError: intensity values contain at least one value less than zero
[...]
```
"""
function RiGCMS(scantimes::T1, retentionindexname::T2, retentionindices::T3,  
    ions::T4, intensities::T5, metadata::Dict=Dict{Any, Any}()) where {
    T1<:AbstractVector{<:Unitful.Time},
    T2<:AbstractString,
    T3<:AbstractVector{<:Union{Real, Missing}},
    T4<:AbstractVector{<:Real},
    T5<:AbstractMatrix{<:Real}}
    RiGCMS{T1, T2, T3, T4, T5}(scantimes, retentionindexname, retentionindices, ions, 
        intensities, metadata)
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

See also [`AbstractGCMS`](@ref), [`LinearDescending`](@ref), [`ionscantimeshift`](@ref), 
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

See also [`AbstractGCMS`](@ref), [`LinearAscending`](@ref), [`ionscantimeshift`](@ref), 
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
    binions(gcms::Union{GCMS, RiGCMS}; ionbin::Function=integer)

Return a GCMS object in which the ions are binned according to the `ionbin` function (the 
default is the integer function), and the intensities of the binned ions are summed.

See also [`AbstractGCMS`](@ref), [`integer`](@ref), [`intensities`](@ref), [`ions`](@ref), 
[`ioncount`](@ref).

# Examples
```jldoctest
julia> gcms = GCMS((1:3)u"s", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])
GCMS {scan times: Int64, ions: Float64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
3 ions: m/z 84.8, 85.2, 100.9
intensity range: 0 - 956
metadata: 0 entries

julia> intensities(gcms)
3×3 Matrix{Int64}:
  0  24   12
  0   0  956
 23   0    1

julia> gcms_binnedions = binions(gcms)
GCMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 101
intensity range: 0 - 956
metadata: 0 entries

julia> ions(gcms_binnedions)
2-element Vector{Int64}:
  85
 101

julia> intensities(gcms_binnedions)
3×2 Matrix{Int64}:
 24   12
  0  956
 23    1


julia> custom_ionbin(ion) = integer(ion, start=0.9);

julia> gcms_binnedions = binions(gcms, ionbin=custom_ionbin)
GCMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
3 ions: m/z 84, 85, 101
intensity range: 0 - 956
metadata: 0 entries

julia> ions(gcms_binnedions)
3-element Vector{Int64}:
  84
  85
 101
```
"""
function binions(gcms::GCMS; ionbin::Function=integer)
    binnedions = sort(unique(ionbin.(ions(gcms))))
    imₙ = zeros(eltype(intensities(gcms)), (size(intensities(gcms), 1), length(binnedions)))
    for (i, ion) in enumerate(ions(gcms))
        binnedion = ionbin(ion)
        imₙ[:, searchsortedlast(binnedions, binnedion)] .+= @view intensities(gcms)[:, i]
    end
    GCMS(scantimes(gcms), binnedions, imₙ, metadata(gcms))
end


function binions(gcms::RiGCMS; ionbin::Function=integer)
    binnedions = sort(unique(ionbin.(ions(gcms))))
    imₙ = zeros(eltype(intensities(gcms)), (size(intensities(gcms), 1), length(binnedions)))
    for (i, ion) in enumerate(ions(gcms))
        binnedion = ionbin(ion)
        imₙ[:, searchsortedlast(binnedions, binnedion)] .+= @view intensities(gcms)[:, i]
    end
    RiGCMS(scantimes(gcms), retentionindexname(gcms), retentionindices(gcms), binnedions, 
        imₙ, metadata(gcms))
end


"""
    cosine(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Return the angle between two non-zero vectors, which can be considered a measure of the
similarity (i.e., `cosine` similarity) between the two vectors.

# Examples
```jldoctest
julia> cosine([100, 500, 250], [200, 1000, 0])
0.8978872704229618

julia> cosine([100, 0, 50], [0, 20, 0])
0.0

julia> cosine([10, 50, 25], [100, 500, 250])
1.0
```
"""
function cosine(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || throw(
        DimensionMismatch("vectors x and y have different lengths"))
    length(x) == 0 && throw(ArgumentError("vectors x and y are empty"))
    iszero(x) && throw(ArgumentError("vector x contains only zeros"))
    iszero(y) && throw(ArgumentError("vector y contains only zeros"))
    s = sum(x .* y) / (sqrt(sum(x.^2)) * sqrt(sum(y.^2)))
    0 ≤ s ≤ 1 && return s
    (s < 0 || isnan(s)) ? zero(typeof(s)) : one(typeof(s))
end
## See also [`similarity`](@ref).


"""
    integer(value:::Real; start::Real=0.7) -> Int

Return the integer for the given `value` that satisfies the following condition: 
`integer - 1 + start ≤ value < integer + start`, where `0 ≤ start < 1`.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`binions`](@ref), [`ions`](@ref).

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
    intensities(chrom::AbstractGC; scanindexrange::OrdinalRange{T, S}) where {T<:Integer, 
    S<:Integer}

Return the intensities. The optional keyword argument `scanindexrange` lets you select a 
subset of scans for which the intensities will be returned. Note that the function will 
return either a reference to the intensity vector or a view into the intensity vector, 
depending on whether a subset of scans is selected.

See also [`AbstractGC`](@ref), [`intensity`](@ref), [`minintensity`](@ref), 
[`maxintensity`](@ref), [`scancount`](@ref).

# Examples
```jldoctest
julia> tic = TIC([1, 2, 3]u"s", [123, 224, 103])
TIC {scan times: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 103 - 224
metadata: 0 entries

julia> intensities(tic)  # reference to the data structure
3-element Vector{Int64}:
 123
 224
 103

julia> intensities(tic)[:]  # a copy of these values
3-element Vector{Int64}:
 123
 224
 103

julia> intensities(tic; scanindexrange=2:3)  # view into the data structure
2-element view(::Vector{Int64}, 2:3) with eltype Int64:
 224
 103

julia> intensities(tic; scanindexrange=2:3)[:]  # a copy of these values
2-element Vector{Int64}:
 224
 103
```
"""
function intensities(
    chrom::AbstractGC; scanindexrange::OrdinalRange{T, S}=firstindex(scantimes(
        chrom)):lastindex(scantimes(chrom))) where {T<:Integer, S<:Integer}

    if scanindexrange == firstindex(scantimes(chrom)):lastindex(scantimes(chrom))
        chrom.intensities
    else
        @view chrom.intensities[scanindexrange]
    end
end


"""
    intensities(gcms::AbstractGCMS; scanindexrange::OrdinalRange{T1, S1}, 
    ionindexrange::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, 
    S2<:Integer}

Return the intensities. The optional keyword arguments `scanindexrange` and 
`ionindexrange` allow you to select specific parts of the intensity matrix to be returned. 
Note that the function returns either a reference to the matrix or a view into it, 
depending on whether the keyword arguments specify subranges of the matrix.

See also [`AbstractGCMS`](@ref), [`scancount`](@ref), [`ioncount`](@ref), 
[`minintensity`](@ref), [`maxintensity`](@ref).

# Examples
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> intensities(gcms)  # reference to the data structure
3×2 Matrix{Int64}:
  0   12
 34  956
 23    1

julia> intensities(gcms)[:, :]  # a copy of these values
3×2 Matrix{Int64}:
  0   12
 34  956
 23    1

julia> intensities(gcms, ionindexrange=1:1)  # all intensities of the ion at index 1
3×1 view(::Matrix{Int64}, 1:3, 1:1) with eltype Int64:
  0
 34
 23

julia> intensities(gcms, ionindexrange=1:1)[:]  # a copy of these values
3-element Vector{Int64}:
  0
 34
 23

julia> intensities(gcms, scanindexrange=1:1)  # intensities of all ions from scan 1
1×2 view(::Matrix{Int64}, 1:1, 1:2) with eltype Int64:
 0  12

julia> intensities(gcms, scanindexrange=1:1)[:]  # a copy of these values
2-element Vector{Int64}:
  0
 12

julia> intensities(gcms, scanindexrange=1:2, ionindexrange=1:2)
2×2 view(::Matrix{Int64}, 1:2, 1:2) with eltype Int64:
  0   12
 34  956

julia> intensities(gcms, scanindexrange=1:2, ionindexrange=1:2)[:, :]
2×2 Matrix{Int64}:
  0   12
 34  956
```
"""
function intensities(
    gcms::AbstractGCMS;
    scanindexrange::OrdinalRange{T1, S1}=firstindex(scantimes(gcms)):lastindex(scantimes(gcms)), 
    ionindexrange::OrdinalRange{T2, S2}=firstindex(ions(gcms)):lastindex(ions(gcms))
    ) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}

    if (scanindexrange == firstindex(scantimes(gcms)):lastindex(scantimes(gcms)) 
        && ionindexrange == firstindex(ions(gcms)):lastindex(ions(gcms)))
        gcms.intensities
    else
        @view intensities(gcms)[scanindexrange, ionindexrange]
    end
end


"""
    intensity(chrom::AbstractGC, scanindex::Integer)

Return the intensity for a scan by specifying its `scanindex`.

See also [`AbstractGC`](@ref), [`intensities`](@ref), [`minintensity`](@ref), 
[`maxintensity`](@ref).

# Examples
```jldoctest
julia> tic = TIC([1, 2, 3]u"s", [123, 224, 103])
TIC {scan times: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
intensity range: 103 - 224
metadata: 0 entries

julia> intensity(tic, 1)
123

julia> intensity(tic, 2)
224
```
"""
intensity(chrom::AbstractGC, scanindex::Integer) = intensities(chrom)[scanindex]


"""
    intensity(chrom::AbstractGC, time::Unitful.Time; precisetime::Bool=false)

Return the intensity at a given `time`. All time units defined in the package
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. By default, the intensity associated with the scan whose timestamp is closest 
to the given `time` is returned. In case of a tie, the intensity of the scan with the 
later timestamp is used. If the optional parameter `precisetime` is set to true, the 
specified `time` must exactly match a timestamp in the vector; otherwise, an error is 
thrown.

See also [`AbstractGC`](@ref), [`intensities`](@ref), [`scantimes`](@ref), 
[`minscantime`](@ref), [`maxscantime`](@ref).

# Examples
```jldoctest
julia> tic = TIC([1.0, 2.0, 3.0]u"s", [123, 224, 103])
TIC {scan times: Float64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
intensity range: 103 - 224
metadata: 0 entries

julia> intensity(tic, 1.5u"s")
224

julia> intensity(tic, 1u"s", precisetime=true)
123

julia> intensity(tic, 1.5u"s", precisetime=true)
ERROR: ArgumentError: scantime 1.5 s does not exist
[...]
```
"""
function intensity(chrom::AbstractGC, time::Unitful.Time; precisetime::Bool=false)
    intensity(chrom, scantimeindex(chrom, time, precisetime=precisetime))
end


"""
    intensity(gcms::AbstractGCMS, scanindex::Integer, ionindex::Integer)

Return the intensity of an ion in a scan, given the `scanindex` of the scan and the 
`ionindex` of the ion.

See also [`AbstractGCMS`](@ref), [`scancount`](@ref), [`ions`](@ref), [`ioncount`](@ref).

# Examples
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> intensity(gcms, 2, 1)
34

julia> intensity(gcms, 1, 2)
12
```
"""
intensity(gcms::AbstractGCMS, scanindex::Integer, ionindex::Integer) = intensities(
    gcms)[scanindex, ionindex]


"""
    intensity(gcms::AbstractGCMS, time::Unitful.Time, ion::Real; precisetime::Bool=false)

Return the intensity of an `ion` at a given `time`. All time units defined in the package
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. By default, the intensity associated with the scan whose timestamp is closest 
to the given `time` is returned. In case of a tie, the intensity of the scan with the 
later timestamp is used. If the optional parameter `precisetime` is set to true, the 
specified `time` must exactly match a timestamp in the vector; otherwise, an error is 
thrown.

See also [`AbstractGCMS`](@ref), [`intensities`](@ref), [`ions`](@ref), [`scantimes`](@ref).

# Examples
```jldoctest
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scan times: Float64, ions: Int64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> intensity(gcms, 1.9u"s", 85)
34

julia> intensity(gcms, 2.9u"s", 85, precisetime=true)
ERROR: ArgumentError: scantime 2.9 s does not exist
[...]

julia> intensity(gcms, 3u"s", 85, precisetime=true)
23
```
"""
function intensity(gcms::AbstractGCMS, time::Unitful.Time, ion::Real; 
    precisetime::Bool=false)
    intensity(gcms, scantimeindex(gcms, time; precisetime=precisetime), 
        ionindex(gcms, ion))
end


"""
    ion(gcms::AbstractGCMS, ionindex::Integer)

Return the ion at the specified `ionindex`.

See also [`AbstractGCMS`](@ref), [`ions`](@ref), [`ionindex`](@ref), [`minion`](@ref), 
[`maxion`](@ref), [`ioncount`](@ref).

# Examples
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
ion(gcms::AbstractGCMS, ionindex::Integer) = ions(gcms)[ionindex]


"""
    ioncount(gcms::AbstractGCMS) -> Int

Return the number of ions.

See also [`AbstractGCMS`](@ref), [`ions`](@ref), [`ion`](@ref), [`minion`](@ref), 
[`maxion`](@ref).

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

Return the index of the specified `ion`. An error is thrown if the `ion` does not exist.

See also [`AbstractGCMS`](@ref), [`ions`](@ref), [`ioncount`](@ref), [`ion`](@ref), 
[`minion`](@ref), [`maxion`](@ref), [`ionscantime`](@ref), [`ionscantimeshift`](@ref).

# Examples
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
    ionscantime(δtᵢ::Function, gcms::AbstractGCMS, scanindex::Integer, ionindex::Integer; 
    timeunit::Unitful.TimeUnits, ustripped::Bool=false)

Return the time at which an ion was actually scanned, given the `scanindex`, `ionindex`, 
and a function `δtᵢ` that computes the time difference between the timestamp of a scan and 
the scan time of the ion from the `ionindex`. The optional parameter `timeunit` allows you 
to specify the unit of the returned scan time. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether the unit is 
included in the returned value. Note that the timestamp of a scan is assumed to be the time 
when the scanning of ion intensities associated with that scan was completed.

See also [`AbstractGCMS`](@ref), [`scantimes`](@ref), [`scantime`](@ref), 
[`scantimeindex`](@ref), [`ions`](@ref), [`ionindex`](@ref), [`ionscantimeshift`](@ref), 
[`IonScanOrder`](@ref), [`LinearAscending`](@ref), [`LinearDescending`](@ref).

# Examples
```julia-repl
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scan times: Float64, ions: Int64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> δtᵢ = ionscantimeshift(gcms, LinearDescending());

julia> ionscantime(δtᵢ, gcms, 2, 1)
2.0 s

julia> ionscantime(δtᵢ, gcms, 2, 2)
1.5 s

julia> ionscantime(δtᵢ, gcms, 2, 2; timeunit=u"minute")
0.025 minute

julia> ionscantime(δtᵢ, gcms, 2, 2; timeunit=u"minute", ustripped=true)
0.025
```
"""
function ionscantime(δtᵢ::Function, gcms::AbstractGCMS, scanindex::Integer, 
    ionindex::Integer; timeunit::Unitful.TimeUnits=unit(eltype(scantimes(gcms))),
    ustripped::Bool=false)
    t = scantime(gcms, scanindex) + δtᵢ(ionindex)
    ustripped ? ustrip(timeunit, t) : uconvert(timeunit, t)
end


"""
    ionscantimeindex(δtᵢ::Function, gcms::AbstractGCMS, ionindex::Integer, 
    time::Unitful.Time; precisetime::Bool=false) -> Int

Return the index of the scan where the scan time for the ion is closest to the specified 
`time`, given the `ionindex` and a function `δtᵢ` that computes the time difference between 
the timestamp of a scan and the scan time of the ion from the `ionindex`. All time units 
defined in the package [Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., 
`u"s"`, `u"minute"`) are supported. In case of a tie, the larger scan index is returned. 
If the optional parameter `precisetime` is set to true, the ion must have been scanned 
exactly at the specified `time`; otherwise, an error is thrown.

See also [`AbstractGCMS`](@ref), [`scantimeindex`](@ref), [`ionscantime`](@ref), 
[`ionscantimeshift`](@ref), [`IonScanOrder`](@ref), [`LinearAscending`](@ref), 
[`LinearDescending`](@ref), [`scantimes`](@ref), [`scantime`](@ref), [`ions`](@ref), 
[`ion`](@ref), [`ionindex`](@ref).

# Examples
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> δtᵢ = ionscantimeshift(gcms, LinearDescending());

julia> ionscantime(δtᵢ, gcms, 2, 2)
1.5 s

julia> ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s")
2

julia> ionscantimeindex(δtᵢ, gcms, 2, 1.6u"s")
2

julia> ionscantimeindex(δtᵢ, gcms, 2, 1.6u"s", precisetime=true)
ERROR: ArgumentError: ion has not been scanned at the time 1.6 s
[...]

julia> ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s", precisetime=true)
2

julia> ionscantimeindex(δtᵢ, gcms, 1, 1.4u"s")
1

julia> ionscantime(δtᵢ, gcms, 1, 1)
1.0 s
```
"""
function ionscantimeindex(δtᵢ::Function, gcms::AbstractGCMS, ionindex::Integer, 
    time::Unitful.Time; precisetime::Bool=false)
    precisetime || return findclosest(ionscantimes(δtᵢ, gcms, ionindex), time)
    t = ustrip(uconvert(unit(eltype(scantimes(gcms))), time))
    for (index, element) in enumerate(ionscantimes(δtᵢ, gcms, ionindex, ustripped=true))
        element ≈ t && return index
    end
    throw(ArgumentError("ion has not been scanned at the time $time"))
end


"""
    ionscantimes(δtᵢ::Function, gcms::AbstractGCMS, ionindex::Integer; 
    timeunit::Unitful.TimeUnits, ustripped::Bool=false)

Return the times at which an ion was actually scanned, given the `ionindex` and a function 
`δtᵢ` that computes the time difference between the timestamp of a scan and the scan time 
of the ion from the `ionindex`. The optional parameter `timeunit` allows you to specify the 
unit for the returned scan times. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether to include the 
unit in the returned value. Note that the timestamp of a scan is assumed to be the time 
when the scanning of ion intensities associated with that scan was completed.

See also [`AbstractGCMS`](@ref), [`ionscantime`](@ref), [`ionscantimeshift`](@ref), 
[`IonScanOrder`](@ref), [`LinearAscending`](@ref), [`LinearDescending`](@ref), 
[`scantimes`](@ref), [`scantimeindex`](@ref), [`ions`](@ref), [`ionindex`](@ref).

# Examples
```julia-repl
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scan times: Float64, ions: Int64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> δtᵢ = ionscantimeshift(gcms, LinearDescending());

julia> ionscantimes(δtᵢ, gcms, 1)
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 1.0 s
 2.0 s
 3.0 s

julia> ionscantimes(δtᵢ, gcms, 2)
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 0.5 s
 1.5 s
 2.5 s

julia> ionscantimes(δtᵢ, gcms, 2; timeunit=u"minute")
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:
 0.008333333333333333 minute
                0.025 minute
 0.041666666666666664 minute

julia> ionscantimes(δtᵢ, gcms, 2; timeunit=u"minute", ustripped=true)
3-element Vector{Float64}:
 0.008333333333333333
 0.025
 0.041666666666666664
```
"""
function ionscantimes(δtᵢ::Function, gcms::AbstractGCMS, ionindex::Integer; 
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(gcms))), ustripped::Bool=false)
    ts = scantimes(gcms) .+ δtᵢ(ionindex)
    ustripped ? ustrip.(timeunit, ts) : uconvert.(timeunit, ts)
end


function ionscantimeshift(ionscanorder::LinearAscending, gcms::AbstractGCMS, error::Real)
    intervalsize = scanduration(gcms, error=error)
    index::Integer -> begin
        firstindex(ions(gcms)) ≤ index ≤ lastindex(ions(gcms)) || throw(
            BoundsError(ions(gcms), index))
        intervalsize * ((ionscanorder.stop - 1) + ((index - ioncount(gcms))
            * (ionscanorder.stop - ionscanorder.start)) / ioncount(gcms))
    end
end


function ionscantimeshift(ionscanorder::LinearDescending, gcms::AbstractGCMS, error::Real)
    intervalsize = scanduration(gcms, error=error)
    index::Integer -> begin
        firstindex(ions(gcms)) ≤ index ≤ lastindex(ions(gcms)) || throw(
            BoundsError(ions(gcms), index))
        intervalsize * ((ionscanorder.stop - 1) + ((1 - index) 
            * (ionscanorder.stop - ionscanorder.start)) / ioncount(gcms))
    end
end


"""
    ionscantimeshift(gcms::AbstractGCMS, ionscanorder::IonScanOrder; error::Real=0.001)

Return a function, based on the `ionscanorder`, that calculates the time difference between 
the timestamp of a scan and the time when an ion was actually scanned, given the ion index 
as an argument. The time difference will be zero or negative, since the timestamp of a scan 
is considered to be when the scanning of the last ion was completed. The returned function 
assumes that the duration of each scan is consistent throughout the run. The optional 
keyword argument `error` lets you specify the maximum allowed deviation of the scan 
duration, as a fraction of the average scan time, between the timestamps of two consecutive 
scans.

See also [`AbstractGCMS`](@ref), [`IonScanOrder`](@ref), [`LinearAscending`](@ref), 
[`LinearDescending`](@ref), [`ionscantime`](@ref), [`ions`](@ref), [`minion`](@ref), 
[`maxion`](@ref), [`ioncount`](@ref), [`scanduration`](@ref), [`scantimes`](@ref), 
[`minscantime`](@ref), [`maxscantime`](@ref).

# Examples
```julia-repl
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scan times: Float64, ions: Int64, intensities: Int64}
3 scans; scan times: 1.0 s, 2.0 s, 3.0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> δtᵢ = ionscantimeshift(gcms, LinearAscending());

julia> δtᵢ(1)
-0.5 s

julia> δtᵢ(2)
0.0 s

julia> δtᵢ = ionscantimeshift(gcms, LinearDescending());

julia> δtᵢ(1)
0.0 s

julia> δtᵢ(2)
-0.5 s

julia> δtᵢ = ionscantimeshift(gcms, LinearDescending(start=0.5));

julia> δtᵢ(1)
0.0 s

julia> δtᵢ(2)
-0.25 s

julia> δtᵢ = ionscantimeshift(gcms, LinearDescending(stop=0.5));

julia> δtᵢ(1)
-0.5 s

julia> δtᵢ(2)
-0.75 s

julia> δtᵢ = ionscantimeshift(gcms, LinearDescending(start=0.25, stop=0.75));

julia> δtᵢ(1)
-0.25 s

julia> δtᵢ(2)
-0.5 s
```
"""
ionscantimeshift(gcms::AbstractGCMS, ionscanorder::IonScanOrder, error::Real=0.001
    ) = ionscantimeshift(ionscanorder, gcms, error)


"""
    ions(gcms::AbstractGCMS)

Return the ions.

See also [`AbstractGCMS`](@ref), [`ioncount`](@ref), [`ion`](@ref), [`minion`](@ref), 
[`maxion`](@ref), [`ionindex`](@ref), [`ionscantime`](@ref), [`ionscantimeshift`](@ref). 

In the following examples, the type of `ions` is explicitely annotated to demonstrate that 
the GCMS object preserves this type.

# Examples
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", Int64[85, 100], [0 12; 34 956; 23 1]);

julia> ions(gcms)
2-element Vector{Int64}:
  85
 100

julia> gcms = GCMS([1, 2, 3]u"s", Float32[85, 100], [0 12; 34 956; 23 1]);

julia> ions(gcms)
2-element Vector{Float32}:
  85.0
 100.0
```
"""
ions(gcms::AbstractGCMS) = gcms.ions


"""
    maxintensity(chrom::AbstractGC[; scanindexrange::OrdinalRange{T, S}]) 
    where {T<:Integer, S<:Integer})

Return the minimum intensity. The optional keyword argument `scanindexrange` allows you to 
specify a range of scans from which the minimum intensity is determined.

See also [`AbstractGC`](@ref), [`minintensity`](@ref), [`intensities`](@ref), 
[`intensity`](@ref).

# Examples
```jldoctest
julia> tic = TIC([1, 2, 3]u"s", [12, 1, 956]);

julia> maxintensity(tic)
956

julia> maxintensity(tic, scanindexrange=1:2)
12
```
"""
function maxintensity(chrom::AbstractGC; scanindexrange::OrdinalRange{T, S}=(
    firstindex(scantimes(chrom)):lastindex(scantimes(chrom)))) where {
    T<:Integer, S<:Integer}
    maximum(intensities(chrom, scanindexrange=scanindexrange))
end


"""
    maxintensity(gcms::AbstractGCMS[; scanindexrange::OrdinalRange{T1, S1}, 
    ionindexrange::OrdinalRange{T2, S2}]) where {T1<:Integer, S1<:Integer, 
    T2<:Integer, S2<:Integer})

Return the maximum intensity. The optional keyword arguments `scanindexrange` and 
`ionindexrange` allow you to select a range of scans and ions from which the minimum 
intensity is determined.

See also [`AbstractGCMS`](@ref), [`minintensity`](@ref), [`intensities`](@ref), 
[`intensity`](@ref).

# Examples
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 1; 23 956]);

julia> maxintensity(gcms)
956

julia> maxintensity(gcms, scanindexrange=1:2)
34

julia> maxintensity(gcms, ionindexrange=1:1)
34

julia> maxintensity(gcms, scanindexrange=1:2, ionindexrange=2:2)
12
```
"""
function maxintensity(gcms::AbstractGCMS; scanindexrange::OrdinalRange{T1, S1}=(
    firstindex(scantimes(gcms)):lastindex(scantimes(gcms))), 
    ionindexrange::OrdinalRange{T2, S2}=(firstindex(ions(gcms)):lastindex(ions(gcms)))
    ) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
    maximum(intensities(gcms, scanindexrange=scanindexrange, ionindexrange=ionindexrange))
end


"""
    maxion(gcms::AbstractGCMS)

Return the largest ion.

See also [`AbstractGCMS`](@ref), [`minion`](@ref), [`ions`](@ref), [`ion`](@ref), 
[`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> maxion(gcms)
100
```
"""
maxion(gcms::AbstractGCMS) = last(ions(gcms))


"""
    maxretentionindex(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{T<:Integer, 
    <:Integer}]; scanindex::Bool=false)

Return the maximum retention index, ignoring any missing values. An optional second 
argument lets you specify a range of scans for which the maximum retention index will be 
returned. The optional keyword argument `scanindex` lets you choose whether to return the
corresponding scan index as well. Note that this function only applies to certain 
AbstractChromatogram subtypes that store retention index data, such as RiFID.

See also [`AbstractChromatogram`](@ref), [`RiFID`](@ref), [`maxretentionindex`](@ref), 
[`retentionindices`](@ref), [`retentionindexname`](@ref).

# Examples
```jldoctest
julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", Int64[100, 200, 300], [12, 956, 1]);

julia> maxretentionindex(tic)
300

julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", [100.0, 200.0, missing], [12, 956, 1]);

julia> maxretentionindex(tic, 2:3)
200.0

julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", [missing, 200.0, 300.0], [12, 956, 1]);

julia> maxretentionindex(tic, scanindex=true)
(300.0, 3)
```
"""
function maxretentionindex(chrom::T1, scanindexrange::OrdinalRange{T2, T3}=firstindex(
    scantimes(chrom)):lastindex(scantimes(chrom)); scanindex::Bool=false) where {T1<:AbstractChromatogram, 
        T2<:Integer, T3<:Integer}
    maxretentionindex(RetentionIndexStyle(T1), chrom, scanindexrange, scanindex)
end


function maxretentionindex(::HasRetentionIndexData, chrom, scanindexrange, scanindex)
    indices = collect(skipmissing(chrom.retentionindices[scanindexrange]))
    length(indices) > 0 || throw(
        ArgumentError("scan index range contains no numerical values"))
    maxval = maximum(skipmissing(chrom.retentionindices[scanindexrange]))
    if scanindex
        maxindex = findfirst(x -> !ismissing(x) && x == maxval, chrom.retentionindices)
        maxval, maxindex
    else
        maxval
    end
end


maxretentionindex(::HasNoRetentionIndexData, chrom, scanindexrange, scanindex) = throw(
    MethodError(minretentionindex, (chrom,)))


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
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> maxscantime(gcms)
3.0 s

julia> maxscantime(gcms, timeunit=u"minute")
0.05 minute

julia> maxscantime(gcms, timeunit=u"minute", ustripped=true)
0.05

julia> maxscantime(gcms, 1:2)
2.0 s

julia> maxscantime(gcms, 1:2, timeunit=u"minute", ustripped=true)
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
julia> gcms₁ = GCMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956], Dict(:id => 4))
GCMS {scan times: Int64, ions: Int64, intensities: Int64}
2 scans; scan times: 1 s, 2 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 1 entry

julia> metadata(gcms₁)
Dict{Any, Any} with 1 entry:
  :id => 4

julia> gcms₂ = GCMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956])
GCMS {scan times: Int64, ions: Int64, intensities: Int64}
2 scans; scan times: 1 s, 2 s
2 ions: m/z 85, 100
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
GCMS {scan times: Int64, ions: Int64, intensities: Int64}
2 scans; scan times: 1 s, 2 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 2 entries

julia> delete!(metadata(gcms₂), "name")
Dict{Any, Any} with 1 entry:
  :id => 123

julia> gcms₂
GCMS {scan times: Int64, ions: Int64, intensities: Int64}
2 scans; scan times: 1 s, 2 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 1 entry
```
"""
metadata(chrom::AbstractChromatogram) = chrom.metadata



"""
    minintensity(chrom::AbstractGC[, greaterthan::Real; 
    scanindexrange::OrdinalRange{T, S}]) where {T<:Integer, S<:Integer})

Return the minimum intensity. The optional positional argument `greaterthan` allows you to 
specify a threshold value; the returned minimum intensity will be greater than this 
threshold. If no intensity exceeds the specified value, nothing is returned. The optional 
keyword argument `scanindexrange` allows you to specify a range of scans from which the 
minimum intensity is determined.

See also [`AbstractGC`](@ref), [`maxintensity`](@ref), [`intensities`](@ref), 
[`intensity`](@ref).

# Examples
```jldoctest
julia> tic = TIC([1, 2, 3]u"s", [12, 956, 1]);

julia> minintensity(tic)
1

julia> minintensity(tic, 1)
12

julia> minintensity(tic, scanindexrange=1:2)
12

julia> minintensity(tic, 1, scanindexrange=2:3)
956
```
"""
function minintensity(chrom::AbstractGC; scanindexrange::OrdinalRange{T, S}=(
    firstindex(scantimes(chrom)):lastindex(scantimes(chrom)))) where {
    T<:Integer, S<:Integer}
    minimum(intensities(chrom, scanindexrange=scanindexrange))
end


function minintensity(chrom::AbstractGC, greaterthan::Real; 
    scanindexrange::OrdinalRange{T, S}=(firstindex(scantimes(chrom)):lastindex(
    scantimes(chrom)))) where {T<:Integer, S<:Integer}
    filtered_intensities = filter(num -> num > greaterthan, intensities(chrom, 
        scanindexrange=scanindexrange))
    length(filtered_intensities) > 0 ? minimum(filtered_intensities) : nothing
end


"""
    minintensity(gcms::AbstractGCMS[, greaterthan::Real; 
    scanindexrange::OrdinalRange{T1, S1}, ionindexrange::OrdinalRange{T2, S2}]) 
    where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer})

Return the minimum intensity. The optional positional argument `greaterthan` allows you to 
specify a threshold value; the returned minimum intensity will be greater than this 
threshold. If no intensity exceeds the specified value, nothing is returned. The optional 
keyword arguments `scanindexrange` and `ionindexrange` allow you to select a range of scans
and ions from which the minimum intensity is determined.

See also [`AbstractGCMS`](@ref), [`maxintensity`](@ref), [`intensities`](@ref), 
[`intensity`](@ref).

# Examples
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> minintensity(gcms)
0

julia> minintensity(gcms, 0)
1

julia> minintensity(gcms, scanindexrange=2:3)
1

julia> minintensity(gcms, ionindexrange=1:1)
0

julia> minintensity(gcms, scanindexrange=2:3, ionindexrange=1:1)
23

julia> minintensity(gcms, 25, scanindexrange=2:3, ionindexrange=1:1)
34
```
"""
function minintensity(gcms::AbstractGCMS; scanindexrange::OrdinalRange{T1, S1}=(
    firstindex(scantimes(gcms)):lastindex(scantimes(gcms))), 
    ionindexrange::OrdinalRange{T2, S2}=(firstindex(ions(gcms)):lastindex(ions(gcms)))
    ) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
    minimum(intensities(gcms, scanindexrange=scanindexrange, ionindexrange=ionindexrange))
end


function minintensity(gcms::AbstractGCMS, greaterthan::Real; 
    scanindexrange::OrdinalRange{T1, S1}=(firstindex(scantimes(gcms)):lastindex(
    scantimes(gcms))), ionindexrange::OrdinalRange{T2, S2}=(firstindex(
    ions(gcms)):lastindex(ions(gcms)))) where {T1<:Integer, S1<:Integer, T2<:Integer, 
    S2<:Integer}
    filtered_intensities = filter(num -> num > greaterthan, intensities(gcms, 
        scanindexrange=scanindexrange, ionindexrange=ionindexrange))
    length(filtered_intensities) > 0 ? minimum(filtered_intensities) : nothing
end


"""
    minion(gcms::AbstractGCMS)

Return the smallest ion.

See also [`AbstractGCMS`](@ref), [`maxion`](@ref), [`ions`](@ref), [`ion`](@ref), 
[`ioncount`](@ref).

# Example
```jldoctest
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> minion(gcms)
85
```

"""
minion(gcms::AbstractGCMS) = first(ions(gcms))


"""
    minretentionindex(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{T<:Integer, 
    <:Integer}]; scanindex::Bool=false)

Return the minimum retention index, ignoring any missing values. An optional second 
argument lets you specify a range of scans for which the minimum retention index will be 
returned. The optional keyword argument `scanindex` lets you choose whether to return the
corresponding scan index as well. Note that this function only applies to certain 
AbstractChromatogram subtypes that store retention index data, such as RiFID.

See also [`AbstractChromatogram`](@ref), [`RiFID`](@ref), [`maxretentionindex`](@ref), 
[`retentionindices`](@ref), [`retentionindexname`](@ref).

# Examples
```jldoctest
julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", Int64[100, 200, 300], [12, 956, 1]);

julia> minretentionindex(tic)
100

julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", [100.0, 200.0, missing], [12, 956, 1]);

julia> minretentionindex(tic, 2:3)
200.0

julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", [missing, 200.0, 300.0], [12, 956, 1]);

julia> minretentionindex(tic, scanindex=true)
(200.0, 2)
```
"""
function minretentionindex(chrom::T1, scanindexrange::OrdinalRange{T2, T3}=firstindex(
    scantimes(chrom)):lastindex(scantimes(chrom)); scanindex::Bool=false) where {
    T1<:AbstractChromatogram, T2<:Integer, T3<:Integer}
    minretentionindex(RetentionIndexStyle(T1), chrom, scanindexrange, scanindex)
end


function minretentionindex(::HasRetentionIndexData, chrom, scanindexrange, scanindex)
    indices = collect(skipmissing(chrom.retentionindices[scanindexrange]))
    length(indices) > 0 || throw(
        ArgumentError("scan index range contains no numerical values"))
    minval = minimum(skipmissing(chrom.retentionindices[scanindexrange]))
    if scanindex
        minindex = findfirst(x -> !ismissing(x) && x == minval, chrom.retentionindices)
        minval, minindex
    else
        minval
    end
end


minretentionindex(::HasNoRetentionIndexData, chrom, scanindexrange, scanindex) = throw(
    MethodError(minretentionindex, (chrom,)))


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
julia> gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> minscantime(gcms)
1.0 s

julia> minscantime(gcms, timeunit=u"minute")
0.016666666666666666 minute

julia> minscantime(gcms, timeunit=u"minute", ustripped=true)
0.016666666666666666

julia> minscantime(gcms, 2:3)
2.0 s

julia> minscantime(gcms, 2:3, timeunit=u"minute", ustripped=true)
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
    retentionindex(chrom::AbstractChromatogram, scanindex::Integer)

Return the retention index for a given `scanindex`. Note that this function only applies to 
certain AbstractChromatogram subtypes that store retention index data, such as RiFID. Also 
note that the retention index returned may be missing a value.

See also [`AbstractChromatogram`](@ref), [`RiFID`](@ref), [`retentionindices`](@ref), 
[`maxretentionindex`](@ref), [`minretentionindex`](@ref) [`retentionindexname`](@ref).

# Examples
```jldoctest
julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, missing], [12, 956, 1]);

julia> retentionindex(tic, 1)
100

julia> retentionindex(tic, 3)
missing
```
"""
function retentionindex(chrom::T, scanindex::Integer) where {T<:AbstractChromatogram}
    retentionindex(RetentionIndexStyle(T), chrom, scanindex)
end


retentionindex(::HasRetentionIndexData, chrom, scanindex) = chrom.retentionindices[
    scanindex]


retentionindex(::HasNoRetentionIndexData, chrom, scanindex) = throw(
    MethodError(retentionindex, (chrom,)))


"""
    retentionindexname(chrom::AbstractChromatogram)

Return the retention index name. Note that this function only applies to certain 
AbstractChromatogram subtypes that store retention index data, such as RiFID.

See also [`AbstractChromatogram`](@ref), [`RiFID`](@ref), [`retentionindices`](@ref), 
[`maxretentionindex`](@ref), [`minretentionindex`](@ref).

# Example
```jldoctest
julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", Int64[100, 200, 300], [12, 956, 1]);

julia> retentionindexname(tic)
"Kovats"
```
"""
retentionindexname(chrom::T) where {T<:AbstractChromatogram} = retentionindexname(
    RetentionIndexStyle(T), chrom)
retentionindexname(::HasRetentionIndexData, chrom) = chrom.retentionindexname
retentionindexname(::HasNoRetentionIndexData, chrom) = throw(
    MethodError(retentionindexname, (chrom,)))


"""
    retentionindices(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{<:Integer, 
    <:Integer}])

Return the retention indices. An optional second argument lets you specify a range of scans 
for which the retention indices will be returned. Note that this function only 
applies to certain AbstractChromatogram subtypes that store retention index data, such as 
RiFID. Additionally, be aware that the returned retention indices may include missing 
values. Depending on whether a subset of scans is selected, the function will return either 
a reference to the data structure holding the retention indices or a view into it.

See also [`AbstractChromatogram`](@ref), [`RiFID`](@ref), [`retentionindexname`](@ref), 
[`maxretentionindex`](@ref), [`minretentionindex`](@ref).

# Examples
```jldoctest
julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", Int64[100, 200, 300], [12, 956, 1]);

julia> retentionindices(tic)
3-element Vector{Int64}:
 100
 200
 300

julia> tic = RiTIC([1, 2, 3]u"s", "Kovats", [missing, 200.0, 300.0], [12, 956, 1]);

julia> retentionindices(tic)
3-element Vector{Union{Missing, Float64}}:
    missing
 200.0
 300.0

julia> retentionindices(tic, 2:3)  # view into the data structure
2-element view(::Vector{Union{Missing, Float64}}, 2:3) with eltype Union{Missing, Float64}:
 200.0
 300.0

julia> retentionindices(tic, 2:3)[:]  # a copy of these values
2-element Vector{Union{Missing, Float64}}:
 200.0
 300.0
```
"""
function retentionindices(chrom::T1, scanindexrange::OrdinalRange{T2, T3}=firstindex(
    scantimes(chrom)):lastindex(scantimes(chrom))) where {
    T1<:AbstractChromatogram, T2<:Integer, T3<:Integer}
        retentionindices(RetentionIndexStyle(T1), chrom, scanindexrange)
end


function retentionindices(::HasRetentionIndexData, chrom, scanindexrange)
    if scanindexrange == firstindex(scantimes(chrom)):lastindex(scantimes(chrom))
        chrom.retentionindices
    else
        @view chrom.retentionindices[scanindexrange]
    end
end


retentionindices(::HasNoRetentionIndexData, chrom, scanindexrange) = throw(
    MethodError(retentionindices, (chrom,)))


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
julia> fid = FID([30.1u"minute", 40.8u"minute", 51.5u"minute"], [12, 956, 23])
FID {scan times: Float64, intensities: Int64}
3 scans; scan times: 30.1 minute, 40.8 minute, 51.5 minute
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
julia> fid = FID([1, 2, 3]u"s", [12, 956, 23]);

julia> scancount(fid)
3
```
"""
scancount(chrom::AbstractChromatogram) = length(scantimes(chrom))


"""
    scanduration(chrom::AbstractChromatogram; error::Real=0.001, 
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
julia> fid = FID([1.0, 2.0, 3.0]u"s", [12, 956, 1]);

julia> scanduration(fid)
1.0 s

julia> scanduration(fid, timeunit=u"minute")
0.016666666666666666 minute

julia> scanduration(fid, timeunit=u"minute", ustripped=true)
0.016666666666666666

julia> scanduration(FID([1.0, 1.99, 3.0]u"s", [12, 956, 1]))
ERROR: ArgumentError: maximum scan duration variation above threshold: 0.010000000000000009 > 0.001
[...]

julia> scanduration(FID([1.0, 1.99, 3.0]u"s", [12, 956, 1]), error=0.02)
1.0 s
```
"""
function scanduration(chrom::AbstractChromatogram; error::Real=0.001,
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
julia> gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> scantime(gcms, 2)
2 s

julia> scantime(gcms, 2, timeunit=u"minute")
1//30 minute

julia> scantime(gcms, 2, timeunit=u"minute", ustripped=true)
1//30
```
"""
function scantime(chrom::AbstractChromatogram, scanindex::Integer; 
    timeunit::Unitful.TimeUnits=unit(eltype(scantimes(chrom))), ustripped::Bool=false)
    ustripped ? ustrip(timeunit, scantimes(chrom)[scanindex]) : uconvert(timeunit, 
    scantimes(chrom)[scanindex])
end


"""
    scantimeindex(gcms::AbstractChromatogram, time::Unitful.Time; 
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
julia> gcms = GCMS([1.1f0, 2.1f0, 3.1f0]u"s", [85, 100], [0 12; 34 956; 23 1])
GCMS {scan times: Float32, ions: Int64, intensities: Int64}
3 scans; scan times: 1.1f0 s, 2.1f0 s, 3.1f0 s
2 ions: m/z 85, 100
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
that the `GCMS` object preserves this type.

# Example
```jldoctest
julia> gcms = GCMS(Float32[1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]);

julia> scantimes(gcms)  # reference to the data structure
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

julia> scantimes(gcms, 2:3)  # view into the data structure
2-element view(::Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}, 2:3) with eltype Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}:
 2.0f0 s
 3.0f0 s

julia> scantimes(gcms, 2:3)[:]  # a copy of these values
2-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
 2.0f0 s
 3.0f0 s

julia> scantimes(gcms, 2:3, timeunit=u"minute", ustripped=true)
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
    totalionchromatogram(gcms::Union{GCMS, RiGCMS})

Compute the total ion chromatrogram.

See also [`GCMS`](@ref), [`scantimes`](@ref), [`intensities`](@ref), [`metadata`](@ref).

In the following example, the type of the intensities is explicitely annotated to 
demonstrate that the TIC object preserves this type.

# Example
```jldoctest
julia> gcms = GCMS((1.1f0:3.1f0)u"s", [85, 100], Int64[0 12; 34 956; 23 1])
GCMS {scan times: Float32, ions: Int64, intensities: Int64}
3 scans; scan times: 1.1f0 s, 2.1f0 s, 3.1f0 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> tic = totalionchromatogram(gcms)
TIC {scan times: Float32, intensities: Int64}
3 scans; scan times: 1.1f0 s, 2.1f0 s, 3.1f0 s
intensity range: 12 - 990
metadata: 0 entries

julia> intensities(tic)
3-element Vector{Int64}:
  12
 990
  24

julia> gcms = RiGCMS((1:3)u"s", "Kovats", 1:1.0:3, [85, 100], Float64[0 12; 34 956; 23 1])
RiGCMS {scan times: Int64, retention indices: Float64, ions: Int64, intensities: Float64}
3 scans; scan times: 1 s, 2 s, 3 s
3 retention indices: 1.0 (≘ 1 s), 2.0 (≘ 2 s), 3.0 (≘ 3 s)
2 ions: m/z 85, 100
intensity range: 0.0 - 956.0
metadata: 0 entries

julia> tic = totalionchromatogram(gcms)
RiTIC {scan times: Int64, retention indices: Float64, intensities: Float64}
3 scans; scan times: 1 s, 2 s, 3 s
3 retention indices: 1.0 (≘ 1 s), 2.0 (≘ 2 s), 3.0 (≘ 3 s)
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


totalionchromatogram(gcms::RiGCMS) = RiTIC(scantimes(gcms), retentionindexname(gcms), 
    retentionindices(gcms), vec(sum(intensities(gcms), dims=2)), metadata(gcms))


function Base.show(io::IO, chrom::AbstractChromatogram)
    # Information about the element types
    type = typeof(chrom)
    typename = name(type)
    hasridata = HasRetentionIndexData() == RetentionIndexStyle(type)
    hasiondata = isa(chrom, AbstractGCMS)
    print(io, "$typename {scan times: ", eltype(ustrip.(scantimes(chrom)))) 
    hasridata && print(io, ", retention indices: ", eltype(retentionindices(chrom)))
    hasiondata && print(io, ", ions: ", eltype(ions(chrom)))
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

    # Information about the retention indices
    if hasridata
        ricount = length(collect(skipmissing(retentionindices(chrom))))
        minri, minri_idx = minretentionindex(chrom, scanindex=true)
        maxri, maxri_idx = maxretentionindex(chrom, scanindex=true)
        if ricount > 10
            println(io, ricount, " retention indices; retention index range: ", minri, 
                " (≘ ", scantime(chrom, minri_idx), ") – ", maxri, " (≘ ", scantime(chrom, 
                maxri_idx), ")")
        elseif 1 < ricount ≤ 10
            print(io, ricount, " retention indices: ")
            for (count, i) in enumerate(minri_idx:maxri_idx)
                print(io, retentionindex(chrom, i), " (≘ ", scantime(chrom, i), ")")
                count < ricount ? print(io, ", ") : println(io)
            end
        else
            println(io, ricount, " retention index: ", minri, " (≘ ", 
            scantime(chrom, minri_idx), ")")
        end
    end

    # Information about the ions
    if hasiondata
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
end
