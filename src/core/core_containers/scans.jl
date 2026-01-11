"""
    AbstractScan{R, I}

Abstract supertype for all scan objects, such as individual chromatographic or mass
spectrometric scan points. 

`R` is the separation unit type (`Unitful.Units` subtype or `Nothing`), and `I` is the 
signal intensity unit type (`Unitful.Units` subtype or `Nothing`).

Concrete subtypes are expected to define `retention::Real`, `retentionunit::Union{
Unitful.Units, Nothing}`, `intensityunit::Union{Unitful.Units, Nothing}`, and 
`attrs::NamedTuple`. The intensity itself is not prescribed: subtypes may store it under 
`intensity`, `intensities`, or another field, and it may be a scalar, vector, matrix, or 
other form.

The type parameter `I` indicates only the unit of the intensity data, not its structure
or storage location. Subtypes should document how to access their intensity values.

See also: [`AbstractChromScan`](@ref), [`AbstractMassScan`](@ref), [`JuChrom.ChromScan`](@ref),
[`attrs`](@ref), [`JuChrom.MassScan`](@ref), [`intensityunit(::AbstractScan)`](@ref),
[`rawretention(::AbstractScan)`](@ref), [`retention(::AbstractScan)`](@ref), 
[`retentionunit(::AbstractScan)`](@ref).
"""
abstract type AbstractScan{R, I} end

"""
    AbstractChromScan{R, I} <: AbstractScan{R, I}

Abstract supertype for individual chromatographic scan points. 

`R` is the separation unit type (`Unitful.Units` subtype or `Nothing`), and `I` is the 
signal intensity unit type (`Unitful.Units` subtype or `Nothing`).

Concrete subtypes are expected to define `retention::Real`, `retentionunit::Union{
Unitful.Units, Nothing}`, `intensity::Real`, `intensityunit::Union{Unitful.Units, Nothing}`, 
and `attrs::NamedTuple`. Subtypes may define additional fields as needed.

See also: [`AbstractScan`](@ref), [`ChromScan`](@ref), [`attrs`](@ref), 
[`intensity`](@ref JuChrom.intensity(::AbstractChromScan{<:Any, Nothing})), 
[`intensityunit(::AbstractScan)`](@ref), 
[`rawintensity(::AbstractChromScan{<:Any, Nothing})`](@ref), 
[`rawretention(::AbstractScan)`](@ref), [`retention(::AbstractScan)`](@ref), 
[`retentionunit(::AbstractScan)`](@ref).
"""
abstract type AbstractChromScan{R, I} <: AbstractScan{R, I} end

"""
    ChromScan{R, I} <: AbstractChromScan{R, I}

Concrete chromatographic scan point with a single intensity value at a retention
coordinate. `ChromScan` is a subtype of `AbstractChromScan`.

`R` is the retention unit type (`Unitful.Units` subtype or `Nothing`), and `I` is the
intensity unit type (`Unitful.Units` subtype or `Nothing`).

Fields include `retention::Real`, `retentionunit::Union{Unitful.Units, Nothing}`,
`intensity::Real`, `intensityunit::Union{Unitful.Units, Nothing}`, and
`attrs::NamedTuple`.

See also: [`AbstractChromScan`](@ref), [`AbstractScan`](@ref), [`attrs`](@ref),
[`intensity(::AbstractChromScan{<:Any, Nothing})`](@ref),
[`intensityunit(::AbstractScan)`](@ref),
[`rawintensity(::AbstractChromScan{<:Any, Nothing})`](@ref),
[`rawretention(::AbstractScan)`](@ref), [`retention(::AbstractScan)`](@ref),
[`retentionunit(::AbstractScan)`](@ref).
"""
struct ChromScan{
    T1<:Real,
    T2<:Union{Nothing, Unitful.Units},
    T3<:Real,
    T4<:Union{Nothing, Unitful.Units},
    T5<:NamedTuple
    } <: AbstractChromScan{T2, T4}

    retention::T1
    retentionunit::T2
    intensity::T3
    intensityunit::T4
    attrs::T5

    function ChromScan{T1, T2, T3, T4, T5}(
        retention::T1,
        retentionunit::T2,
        intensity::T3,
        intensityunit::T4,
        attrs::T5
        ) where {
            T1<:Real,
            T2<:Union{Nothing, Unitful.Units},
            T3<:Real,
            T4<:Union{Nothing, Unitful.Units},
            T5<:NamedTuple
        }

        isfinite(retention) || throw(ArgumentError("retention must be finite"))
        isfinite(intensity) || throw(ArgumentError("intensity must be finite"))

        new{T1, T2, T3, T4, T5}(
            retention, retentionunit, 
            intensity, intensityunit, 
            attrs)
    end
end

"""
    ChromScan(
        retention::Union{Real, AbstractQuantity{<:Real}},
        intensity::Union{Real, AbstractQuantity{<:Real}};
        attrs::NamedTuple=NamedTuple()
    )

Construct a `ChromScan` object representing a single chromatographic scan point acquired at
a given value along the separation axis (e.g. retention time, index, or position).

Both `retention` and `intensity` may be given as plain numbers or `Unitful.AbstractQuantity` 
values. If units are provided, they are stripped from the values and stored separately in 
the `retentionunit` and `intensityunit` fields. If no unit is provided, `nothing` is 
stored.

`retention` is the separation coordinate (for example `2u"s"` or `5000`), `intensity` is
the signal intensity (for example `1000` or `10u"pA"`), and `attrs` provides optional
scan-level metadata as a `NamedTuple`. Returns a `ChromScan` instance where numeric values
and their units are stored separately. Throws `ArgumentError` if `retention` or `intensity`
is not finite.

# Examples
```jldoctest
julia> using Unitful

julia> csc = ChromScan(5u"minute", 1000u"pA");

julia> retention(csc)
5 minute

julia> rawretention(csc)
5

julia> intensity(csc)
1000 pA

julia> rawintensity(csc)
1000

julia> attrs(csc)
NamedTuple()

julia> csc2 = ChromScan(5000, 100, attrs=(scan_id=1,));

julia> attrs(csc2)
(scan_id = 1,)
```
"""
@inline ChromScan(rt::AbstractQuantity, I::AbstractQuantity; attrs::NamedTuple=NamedTuple()) = begin
    rt_v = ustrip(rt); rt_u = unit(rt)
    I_v  = ustrip(I);  I_u  = unit(I)
    ChromScan{typeof(rt_v), typeof(rt_u), typeof(I_v), typeof(I_u), typeof(attrs)}(
        rt_v, rt_u, I_v, I_u, attrs)
end

@inline ChromScan(rt::AbstractQuantity, I::Real; attrs::NamedTuple=NamedTuple()) = begin
    rt_v = ustrip(rt); rt_u = unit(rt)
    ChromScan{typeof(rt_v), typeof(rt_u), typeof(I), Nothing, typeof(attrs)}(
        rt_v, rt_u, I, nothing, attrs)
end

@inline ChromScan(rt::Real, I::AbstractQuantity; attrs::NamedTuple=NamedTuple()) = begin
    I_v = ustrip(I); I_u = unit(I)
    ChromScan{typeof(rt), Nothing, typeof(I_v), typeof(I_u), typeof(attrs)}(
        rt, nothing, I_v, I_u, attrs)
end

@inline ChromScan(rt::Real, I::Real; attrs::NamedTuple=NamedTuple()) =
    ChromScan{typeof(rt), Nothing, typeof(I), Nothing, typeof(attrs)}(
        rt, nothing, I, nothing, attrs)

"""
    ChromScan(
        retention::Real,
        retentionunit::Union{Unitful.Units, Nothing},
        intensity::Real,
        intensityunit::Union{Unitful.Units, Nothing};
        attrs::NamedTuple=NamedTuple()
    )

Construct a `ChromScan` object by explicitly providing the unit-stripped numeric values and
their associated units.

This constructor bypasses automatic unit inference. It is intended for advanced use cases
such as programmatic data loading, conversion, or optimization, where unit parsing has
already been performed externally.

`retention` is the numeric separation coordinate, `retentionunit` is its unit or
`nothing`, `intensity` is the numeric signal intensity, `intensityunit` is its unit or
`nothing`, and `attrs` provides optional scan-level metadata. Throws `ArgumentError` if
`retention` or `intensity` is not finite. No unit consistency checks are performed; values
are assumed to be prevalidated.
"""
@inline ChromScan(rt_val::Real, rt_unit::Union{Nothing, Unitful.Units},
    I_val::Real, I_unit::Union{Nothing, Unitful.Units}; attrs::NamedTuple=NamedTuple()) =
    ChromScan{typeof(rt_val), typeof(rt_unit), typeof(I_val), typeof(I_unit), typeof(attrs)}(
        rt_val, rt_unit, I_val, I_unit, attrs)

"""
    Base.:(==)(a::ChromScan, b::ChromScan)

Return `true` if two `ChromScan` objects have approximately equal retention and intensity
values (`isapprox` defaults), identical retention and intensity units, and identical
`attrs`.

See also: [`AbstractChromScan`](@ref), [`ChromScan`](@ref), [`attrs`](@ref),
[`intensity(::AbstractChromScan{<:Any, Nothing})`](@ref),
[`intensityunit(::AbstractScan)`](@ref),
[`rawintensity(::AbstractChromScan{<:Any, Nothing})`](@ref),
[`rawretention(::AbstractScan)`](@ref), [`retention(::AbstractScan)`](@ref),
[`retentionunit(::AbstractScan)`](@ref).
"""
Base.:(==)(a::ChromScan, b::ChromScan) = 
    retention(a) ≈ retention(b) &&
    retentionunit(a) == retentionunit(b) &&
    intensity(a) ≈ intensity(b) &&
    intensityunit(a) == intensityunit(b) &&
    attrs(a) == attrs(b)

"""
    AbstractMassScan{R, M, I} <: AbstractScan{R, I}

Abstract supertype for individual mass spectrometric scan points. 

`R` is the separation unit type (`Unitful.Units` subtype or `Nothing`), `M` is the m/z unit 
type (`Unitful.Units` subtype or `Nothing`), and `I` is the signal intensity unit type 
(`Unitful.Units` subtype or `Nothing`).

Concrete subtypes are expected to define `retention::Real`, `retentionunit::Union{
Unitful.Units, Nothing}`, `mzvalues::AbstractVector{<:Real}`, `mzunit::Union{
Unitful.Units, Nothing}`, `intensities::AbstractVector{<:Real}`, `intensityunit::Union{
Unitful.Units, Nothing}`, `level::Integer`, and `attrs::NamedTuple`. Subtypes may define 
additional fields as needed.

See also: [`AbstractScan`](@ref), [`MassScan`](@ref), [`attrs`](@ref), 
[`intensities(::AbstractMassScan)`](@ref), [`intensityunit(::AbstractScan)`](@ref), 
[`mzcount(::AbstractMassScan)`](@ref), [`mzvalues(::AbstractMassScan)`](@ref), 
[`mzunit(::AbstractMassScan)`](@ref), [`rawintensities(::AbstractMassScan)`](@ref), 
[`rawmzvalues(::AbstractMassScan)`](@ref), [`rawretention(::AbstractScan)`](@ref), 
[`retention(::AbstractScan)`](@ref), [`retentionunit(::AbstractScan)`](@ref).
"""
abstract type AbstractMassScan{R, M, I} <: AbstractScan{R, I} end

"""
    MassScan{R, M, I} <: AbstractMassScan{R, M, I}

Concrete mass spectrometric scan point with m/z and intensity vectors. `MassScan` is a
subtype of `AbstractMassScan`.

`R` is the retention unit type (`Unitful.Units` subtype or `Nothing`), `M` is the m/z unit
type (`Unitful.Units` subtype or `Nothing`), and `I` is the intensity unit type
(`Unitful.Units` subtype or `Nothing`).

Fields include `retention::Real`, `retentionunit::Union{Unitful.Units, Nothing}`,
`mzvalues::AbstractVector{<:Real}`, `mzunit::Union{Unitful.Units, Nothing}`,
`intensities::AbstractVector{<:Real}`, `intensityunit::Union{Unitful.Units, Nothing}`,
`level::Integer`, and `attrs::NamedTuple`.

See also: [`AbstractMassScan`](@ref), [`AbstractScan`](@ref), [`attrs`](@ref), 
[`intensities(::AbstractMassScan)`](@ref), [`intensityunit(::AbstractScan)`](@ref), 
[`mzcount(::AbstractMassScan)`](@ref), [`mzunit(::AbstractMassScan)`](@ref),
[`mzvalues(::AbstractMassScan)`](@ref), [`rawintensities(::AbstractMassScan)`](@ref), 
[`rawretention(::AbstractScan)`](@ref), [`rawmzvalues(::AbstractMassScan)`](@ref),
[`retention(::AbstractScan)`](@ref), [`retentionunit(::AbstractScan)`](@ref).
"""
struct MassScan{
    T1<:Real,
    T2<:Union{Nothing, Unitful.Units},
    T3<:AbstractVector{<:Real},
    T4<:Union{Nothing, Unitful.Units},
    T5<:AbstractVector{<:Real},
    T6<:Union{Nothing, Unitful.Units},
    T7<:Integer,
    T8<:NamedTuple
    } <: AbstractMassScan{T2, T4, T6}

    retention::T1
    retentionunit::T2
    mzvalues::T3
    mzunit::T4
    intensities::T5
    intensityunit::T6
    level::T7
    attrs::T8

    function MassScan{T1, T2, T3, T4, T5, T6, T7, T8}(
        retention::T1,
        retentionunit::T2,
        mzvalues::T3,
        mzunit::T4,
        intensities::T5,
        intensityunit::T6,
        level::T7,
        attrs::T8
        ) where {
            T1<:Real,
            T2<:Union{Nothing, Unitful.Units},
            T3<:AbstractVector{<:Real},
            T4<:Union{Nothing, Unitful.Units},
            T5<:AbstractVector{<:Real},
            T6<:Union{Nothing, Unitful.Units},
            T7<:Integer,
            T8<:NamedTuple
        }

        isfinite(retention) || throw(ArgumentError("retention must be finite"))
        length(mzvalues) ≥ 1 || throw(ArgumentError("no m/z value(s) provided"))
        all(isfinite, mzvalues) || throw(ArgumentError("all m/z values must be finite"))
        all(mz -> mz > 0, mzvalues) || throw(
            ArgumentError("all m/z values must be greater than zero"))
        all(diff(mzvalues) .> 0) || throw(ArgumentError(
            "m/z values must be strictly increasing (no duplicates allowed)"))
        length(intensities) ≥ 1 || throw(ArgumentError("no intensity value(s) provided"))
        all(isfinite, intensities) || throw(ArgumentError("all intensities must be finite"))
        length(mzvalues) == length(intensities) || throw(
            DimensionMismatch("m/z value count does not match intensity count"))
        level ≥ 1 || throw(ArgumentError("level must be ≥ 1"))

        new{T1, T2, T3, T4, T5, T6, T7, T8}(
            retention, retentionunit,
            mzvalues, mzunit,
            intensities, intensityunit,
            level,
            attrs)
    end
end

"""
    MassScan(
        retention::Union{Real, AbstractQuantity{<:Real}},
        mzvalues::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}},
        intensities::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}};
        level::Integer=1,
        attrs::NamedTuple=NamedTuple()
    )

Construct a `MassScan` object representing a single scan acquired at a specific point along
the separation axis (e.g. time, index, or physical position).

All three arguments — `retention`, `mzvalues`, and `intensities` — may be provided as
plain numbers or `Unitful.AbstractQuantity` values. If units are present, they are stripped 
from the values and stored separately. Each of mzvalues and intensities must be either 
entirely unitless or share a consistent unit within their respective vectors.

`retention` is the separation coordinate (for example `2.0u"s"` or `5000`), `mzvalues`
is a vector of m/z values, and `intensities` is the matching vector of signal intensities.
`level` sets the MS level (default `1`), and `attrs` provides optional scan-level metadata.
Returns a `MassScan` instance where numeric values and their units (if present) are stored
separately.

Throws `ArgumentError` if `retention` is not finite, if `mzvalues` are empty, non-positive,
non-finite, or not strictly increasing, if `mzvalues` or `intensities` mix unitful and
unitless values or have inconsistent units, if `intensities` are empty or contain non-finite
values, or if `level < 1`. Throws `DimensionMismatch` if 
`length(mzvalues) ≠ length(intensities)`.

# Examples
```jldoctest
julia> msc = MassScan(2.0u"s", [100.0, 150.0], [1000, 2000]);

julia> retention(msc)
2.0 s

julia> msc.mzvalues == [100.0, 150.0]
true

julia> intensities(msc) == [1000, 2000]
true

julia> attrs(msc)
NamedTuple()

julia> msc2 = MassScan(3.0u"s", [120.0, 160.0], [600, 1100], attrs=(scan_id=42,));

julia> attrs(msc2)
(scan_id = 42,)
```
"""
@inline MassScan(rt::AbstractQuantity,
                 mz::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}},
                 I::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}};
                 level::Integer=1,
                 attrs::NamedTuple=NamedTuple()) = begin
    rt_v = ustrip(rt); rt_u = unit(rt)
    mz_v, mz_u = strip_units_checked(mz, "m/z values")
    I_v,  I_u  = strip_units_checked(I,  "intensities")
    MassScan{typeof(rt_v), typeof(rt_u),
             typeof(mz_v), typeof(mz_u),
             typeof(I_v),  typeof(I_u),
             typeof(level), typeof(attrs)}(
        rt_v, rt_u, mz_v, mz_u, I_v, I_u, level, attrs)
end

@inline MassScan(rt::Real,
                 mz::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}},
                 I::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}};
                 level::Integer=1,
                 attrs::NamedTuple=NamedTuple()) = begin
    mz_v, mz_u = strip_units_checked(mz, "m/z values")
    I_v,  I_u  = strip_units_checked(I,  "intensities")
    MassScan{typeof(rt), Nothing,
             typeof(mz_v), typeof(mz_u),
             typeof(I_v),  typeof(I_u),
             typeof(level), typeof(attrs)}(
        rt, nothing, mz_v, mz_u, I_v, I_u, level, attrs)
end

"""
    MassScan(
        retention::Real,
        retentionunit::Union{Unitful.Units, Nothing},
        mzvalues::AbstractVector{<:Real},
        mzunit::Union{Unitful.Units, Nothing},
        intensities::AbstractVector{<:Real},
        intensityunit::Union{Unitful.Units, Nothing};
        level::Integer=1,
        attrs::NamedTuple=NamedTuple()
    )

Construct a `MassScan` object from pre-parsed and unit-stripped values. This constructor is 
intended for advanced use cases where unit inference and validation are handled externally 
(e.g. during deserialization or optimized data loading).

`retention` is the numeric retention value, `retentionunit` is its unit (or `nothing`),
`mzvalues` is the vector of m/z values, `mzunit` is its unit (or `nothing`),
`intensities` is the vector of numeric intensity values, and `intensityunit` is its unit
(or `nothing`). `level` sets the MS level (default `1`), and `attrs` provides optional
metadata. Throws `ArgumentError` or `DimensionMismatch` if values are invalid (see inner
constructor). This constructor performs no unit consistency checks; use it only when values
have already been validated.
"""
@inline MassScan(rt_val::Real, rt_unit,
                 mz::AbstractVector{<:Real}, mzunit,
                 I::AbstractVector{<:Real}, I_unit;
                 level::Integer=1,
                 attrs::NamedTuple=NamedTuple()) =
    MassScan{typeof(rt_val), typeof(rt_unit),
             typeof(mz), typeof(mzunit),
             typeof(I),  typeof(I_unit),
             typeof(level), typeof(attrs)}(
        rt_val, rt_unit, mz, mzunit, I, I_unit, level, attrs)

"""
    Base.:(==)(a::MassScan, b::MassScan)

Return `true` if two `MassScan` objects have approximately equal retention, m/z values,
and intensities (`isapprox` defaults), identical m/z, retention, and intensity units, the
same m/z length, the same `level`, and identical `attrs`.

See also: [`AbstractMassScan`](@ref), [`MassScan`](@ref),[`attrs`](@ref), 
[`intensities(::AbstractMassScan)`](@ref), [`intensityunit(::AbstractScan)`](@ref), 
[`mzunit(::AbstractMassScan)`](@ref), [`mzvalues(::AbstractMassScan)`](@ref), 
[`rawintensities(::AbstractMassScan)`](@ref), [`rawmzvalues(::AbstractMassScan)`](@ref), 
[`rawretention(::AbstractScan)`](@ref), [`retention(::AbstractScan)`](@ref), 
[`retentionunit(::AbstractScan)`](@ref).
"""
Base.:(==)(a::MassScan, b::MassScan) = 
    isapprox(retention(a), retention(b)) &&
    retentionunit(a) == retentionunit(b) &&
    length(mzvalues(a)) == length(mzvalues(b)) &&
    all(isapprox.(rawmzvalues(a), rawmzvalues(b))) &&
    mzunit(a) == mzunit(b) &&
    all(isapprox.(rawintensities(a), rawintensities(b)))  &&
    intensityunit(a) == intensityunit(b) &&
    level(a) == level(b) &&
    attrs(a) == attrs(b)
