"""
    AbstractScan{R, I}

Abstract supertype for all scan objects, such as individual chromatographic or mass
spectrometric scan points.

# Type Parameters
- `R`: Type of the separation unit (i.e. `Unitful.Units` subtype or `Nothing`)
- `I`: Type of the signal intensity unit (i.e. `Unitful.Units` subtype or `Nothing`)

# Required Fields for Subtypes
Concrete subtypes must define the following fields:

- `retention::Real` — separation coordinate (e.g. time, index, or position), stored without 
  units
- `retention_unit::Union{Unitful.Units, Nothing}` — unit of the separation axis, or 
  `nothing` if not provided
- `intensity_unit::Union{Unitful.Units, Nothing}` — unit of the signal intensity, or 
  `nothing` if not provided
- `attrs::NamedTuple` — additional scan-level attrs

**Note:** This abstract type does **not** enforce a specific field name or structure for the 
intensity data. Subtypes may define fields like `intensity`, `intensities`, or others, and 
the intensity may be a scalar, a vector, a matrix, or any other form.

The type parameter `I` only indicates the unit of the intensity data, not its structure or 
where it is stored. Subtypes are responsible for documenting how to access their intensity 
values.
"""
abstract type AbstractScan{R, I} end

"""
    AbstractChromScan{R, I} <: AbstractScan{R, I}

Abstract supertype for individual chromatographic scan points.

# Type Parameters
- `R`: Type of the separation unit (i.e. `Unitful.Units` subtype or `Nothing`)
- `I`: Type of the signal intensity unit (i.e. `Unitful.Units` subtype or `Nothing`)

# Required Fields for Subtypes
Concrete subtypes must define the following fields:

- `retention::Real` — separation coordinate (e.g. time, index, position), stored without 
  units
- `retention_unit::Union{Unitful.Units, Nothing}` — unit of retention, or `nothing` if not 
  provided
- `intensity::Real` — signal intensity value, stored without units
- `intensity_unit::Union{Unitful.Units, Nothing}` — unit of intensity, or `nothing` if not 
  provided
- `attrs::NamedTuple` — additional scan-level attrs

Subtypes may define additional fields as necessary.
"""
abstract type AbstractChromScan{R, I} <: AbstractScan{R, I} end

struct ChromScan{
    T1<:Real,
    T2<:Union{Nothing, Unitful.Units},
    T3<:Real,
    T4<:Union{Nothing, Unitful.Units},
    T5<:NamedTuple
    } <: AbstractChromScan{T2, T4}

    retention::T1
    retention_unit::T2
    intensity::T3
    intensity_unit::T4
    attrs::T5

    function ChromScan{T1, T2, T3, T4, T5}(
        retention::T1,
        retention_unit::T2,
        intensity::T3,
        intensity_unit::T4,
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
            retention, retention_unit, 
            intensity, intensity_unit, 
            attrs)
    end
end

"""
    ChromScan(retention::Union{Real, AbstractQuantity{<:Real}},
              intensity::Union{Real, AbstractQuantity{<:Real}};
              attrs::NamedTuple = NamedTuple())

Construct a `ChromScan` object representing a single chromatographic scan point acquired at
a given value along the separation axis (e.g. retention time, index, or position).

Both `retention` and `intensity` may be given as plain numbers or `Unitful.AbstractQuantity` 
values. If units are provided, they are stripped from the values and stored separately in 
the `retention_unit` and `intensity_unit` fields. If no unit is provided, `nothing` is 
stored.

# Arguments
- `retention`: Separation coordinate (e.g. `2u"s"` or `5000`)
- `intensity`: Signal intensity (e.g. `1000` or `10u"pA"`)
- `attrs`: Optional scan-level attrs as a `NamedTuple` (default: empty)

# Returns
A `ChromScan` instance where numeric values and their units are stored separately.

# Throws
- `ArgumentError` if `retention` or `intensity` is not finite

See also [`AbstractChromScan`](@ref), [`AbstractScan`](@ref), [`retention`](@ref), 
[`rawretention`](@ref), [`retentionunit`](@ref), [`intensity`](@ref), 
[`rawintensity`](@ref), [`intensityunit`](@ref), [`attrs`](@ref).

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
    ChromScan(retention::Real, retention_unit::Union{Unitful.Units, Nothing},
              intensity::Real, intensity_unit::Union{Unitful.Units, Nothing};
              attrs::NamedTuple = NamedTuple())

Construct a `ChromScan` object by explicitly providing the unit-stripped numeric values and
their associated units.

This constructor bypasses automatic unit inference. It is intended for advanced use cases
such as programmatic data loading, conversion, or optimization, where unit parsing has
already been performed externally.

# Arguments
- `retention`: Numeric value representing the separation coordinate
- `retention_unit`: Unit of retention (`Unitful.Units` subtype) or `nothing`
- `intensity`: Numeric intensity value
- `intensity_unit`: Unit of intensity (`Unitful.Units` subtype) or `nothing`
- `attrs`: Optional scan-level attrs (`NamedTuple`, default: empty)

# Throws
- `ArgumentError` if `retention` or `intensity` is not finite

# Note
No unit consistency checks are performed in this constructor. It assumes all values are 
prevalidated.
"""
@inline ChromScan(rt_val::Real, rt_unit::Union{Nothing, Unitful.Units},
    I_val::Real, I_unit::Union{Nothing, Unitful.Units}; attrs::NamedTuple=NamedTuple()) =
    ChromScan{typeof(rt_val), typeof(rt_unit), typeof(I_val), typeof(I_unit), typeof(attrs)}(
        rt_val, rt_unit, I_val, I_unit, attrs)

Base.:(==)(a::ChromScan, b::ChromScan) = 
    retention(a) ≈ retention(b) &&
    retentionunit(a) == retentionunit(b) &&
    intensity(a) ≈ intensity(b) &&
    intensityunit(a) == intensityunit(b) &&
    attrs(a) == attrs(b)

"""
    AbstractMassScan{R, M, I} <: AbstractScan{R, I}

Abstract supertype for individual mass spectrometric scan points.

# Type Parameters
- `R`: Type of the separation unit (i.e. `Unitful.Units` subtype or `Nothing`)
- `M`: Type of the m/z value unit (i.e. `Unitful.Units` subtype or `Nothing`)
- `I`: Type of the signal intensity unit (i.e. `Unitful.Units` subtype or `Nothing`)

# Required Fields for Subtypes
Concrete subtypes must define the following fields:

- `retention::Real` — separation coordinate (e.g. time, index, position), stored without 
  units
- `retention_unit::Union{Unitful.Units, Nothing}` — unit of separation, or `nothing` if not 
  provided
- `mz_values::AbstractVector{<:Real}` — vector of m/z values (must be non-empty, finite, 
  larger than zero, strictly increasing)
- `mz_unit::Union{Unitful.Units, Nothing}` — optional unit of m/z values (e.g. Da/e); 
  typically `nothing` since m/z is unitless by convention
- `intensities::AbstractVector{<:Real}` — vector of intensity values (same length as 
  `mz_values`, all finite)
- `intensity_unit::Union{Unitful.Units, Nothing}` — unit of intensity, or `nothing` if not 
  provided
- `level::Integer` — MS level (must be ≥ 1)
- `attrs::NamedTuple` — additional scan-level attrs

Subtypes may define additional fields as necessary.
"""
abstract type AbstractMassScan{R, M, I} <: AbstractScan{R, I} end

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
    retention_unit::T2
    mz_values::T3
    mz_unit::T4
    intensities::T5
    intensity_unit::T6
    level::T7
    attrs::T8

    function MassScan{T1, T2, T3, T4, T5, T6, T7, T8}(
        retention::T1,
        retention_unit::T2,
        mz_values::T3,
        mz_unit::T4,
        intensities::T5,
        intensity_unit::T6,
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
        length(mz_values) ≥ 1 || throw(ArgumentError("no m/z value(s) provided"))
        all(isfinite, mz_values) || throw(ArgumentError("all m/z values must be finite"))
        all(mz -> mz > 0, mz_values) || throw(
            ArgumentError("all m/z values must be greater than zero"))
        all(diff(mz_values) .> 0) || throw(ArgumentError(
            "m/z values must be strictly increasing (no duplicates allowed)"))  # consider removing this constraint to reflect different ms scan directions
        length(intensities) ≥ 1 || throw(ArgumentError("no intensity value(s) provided"))
        all(isfinite, intensities) || throw(ArgumentError("all intensities must be finite"))
        length(mz_values) == length(intensities) || throw(
            DimensionMismatch("m/z value count does not match intensity count"))
        level ≥ 1 || throw(ArgumentError("level must be ≥ 1"))

        new{T1, T2, T3, T4, T5, T6, T7, T8}(
            retention, retention_unit,
            mz_values, mz_unit,
            intensities, intensity_unit,
            level,
            attrs)
    end
end

"""
    MassScan(
        retention::Union{Real, AbstractQuantity{<:Real}},
        mz_values::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}},
        intensities::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}};
        level::Integer = 1,
        attrs::NamedTuple = NamedTuple()
    )

Construct a `MassScan` object representing a single scan acquired at a specific point along
the separation axis (e.g. time, index, or physical position).

All three arguments — `retention`, `mz_values`, and `intensities` — may be provided as
plain numbers or `Unitful.AbstractQuantity` values. If units are present, they are stripped 
from the values and stored separately. Each of mz_values and intensities must be either 
entirely unitless or share a consistent unit within their respective vectors.

# Arguments
- `retention`: Separation coordinate (e.g. `2.0u"s"` or `5000`)
- `mz_values`: Vector of mass-to-charge (m/z) values; must be non-empty, finite, strictly 
  increasing, and > 0
- `intensities`: Vector of signal intensities (same length as `mz_values`, all finite)
- `level`: MS level (default: `1`; must be ≥ 1)
- `attrs`: Optional scan-level attrs as a `NamedTuple`

# Returns
A `MassScan` instance where numeric values and their units (if present) are stored 
separately.

# Throws
- `ArgumentError` if:
  - `retention` is not finite
  - `mz_values` is empty, contains non-positive or non-finite values, or is not strictly 
    increasing
  - `mz_values` mix unitful and unitless values
  - `mz_values` have inconsistent units
  - `intensities` is empty
  - `intensities` contain non-finite values
  - `intensities` mix unitful and unitless values
  - `intensities` have inconsistent units
  - `level` is less than 1
- `DimensionMismatch` if `length(mz_values) ≠ length(intensities)`

See also [`AbstractMassScan`](@ref), [`AbstractScan`](@ref), [`retention`](@ref), 
[`rawretention`](@ref), [`retentionunit`](@ref), [`mzvalues`](@ref), [`rawmzvalues`](@ref), 
[`mzcount`](@ref), [`mzunit`](@ref) [`intensities`](@ref), [`rawintensities`](@ref), 
[`intensityunit`](@ref), [`level`](@ref), [`attrs`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(2.0u"s", [100.0, 150.0], [1000, 2000]);

julia> retention(msc)
2.0 s

julia> msc.mz_values == [100.0, 150.0]
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
    MassScan(retention::Real, retention_unit::Union{Unitful.Units, Nothing},
             mz_values::AbstractVector{<:Real},
             mz_unit::Union{Unitful.Units, Nothing},
             intensities::AbstractVector{<:Real},
             intensity_unit::Union{Unitful.Units, Nothing};
             level::Integer = 1,
             attrs::NamedTuple = NamedTuple())

Construct a `MassScan` object from pre-parsed and unit-stripped values.
This constructor is intended for advanced use cases where unit inference and validation
are handled externally (e.g. during deserialization or optimized data loading).

# Arguments
- `retention`: Numeric retention value (unit-stripped)
- `retention_unit`: Unit of retention (or `nothing`)
- `mz_values`: Vector of m/z values (must be validated externally)
- `mz_unit`: Unit of m/z values (or `nothing`)
- `intensities`: Vector of numeric intensity values (unit-stripped)
- `intensity_unit`: Unit of intensity (or `nothing`)
- `level`: MS level (default: 1)
- `attrs`: Optional attrs (`NamedTuple`)

# Throws
- `ArgumentError` or `DimensionMismatch` if values are invalid (see inner constructor)

# Note
This constructor performs no unit consistency checks. Use only when values
have already been validated.
"""
@inline MassScan(rt_val::Real, rt_unit,
                 mz::AbstractVector{<:Real}, mz_unit,
                 I::AbstractVector{<:Real}, I_unit;
                 level::Integer=1,
                 attrs::NamedTuple=NamedTuple()) =
    MassScan{typeof(rt_val), typeof(rt_unit),
             typeof(mz), typeof(mz_unit),
             typeof(I),  typeof(I_unit),
             typeof(level), typeof(attrs)}(
        rt_val, rt_unit, mz, mz_unit, I, I_unit, level, attrs)

Base.:(==)(a::MassScan, b::MassScan) = 
    isapprox(retention(a), retention(b)) &&
    retentionunit(a) == retentionunit(b) &&
    length(mzvalues(a)) == length(mzvalues(b)) &&
    all(isapprox.(rawmzvalues(a), rawmzvalues(b))) &&
    mzunit(a) == mzunit(b) &&
    all(isapprox.(rawintensities(a), rawintensities(b)))  &&
    intensityunit(a) == intensityunit(b) &&
    attrs(a) == attrs(b)
