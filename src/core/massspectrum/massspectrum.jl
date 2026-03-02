"""
    AbstractMassSpectrum{M, I} <: Any

Abstract supertype for mass spectrum containers holding paired m/z values and
intensity values.

`M` is the m/z unit type (`Unitful.Units` subtype or `Nothing`), and `I` is the
intensity unit type (`Unitful.Units` subtype or `Nothing`).

Concrete subtypes are expected to define `mzvalues::AbstractVector{<:Real}`,
`mzunit::Union{Unitful.Units, Nothing}`, `intensities::AbstractVector{<:Real}`,
`intensityunit::Union{Unitful.Units, Nothing}`, and `attrs::NamedTuple`.

See also
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`attrs`](@ref JuChrom.attrs(::AbstractMassSpectrum)),
[`intensities`](@ref JuChrom.intensities(::AbstractMassSpectrum{<:Any, Nothing})),
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassSpectrum)),
[`mzcount`](@ref JuChrom.mzcount(::AbstractMassSpectrum)),
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassSpectrum)),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassSpectrum{Nothing, <:Any})).
"""
abstract type AbstractMassSpectrum{M, I} <: Any end

"""
    MassSpectrum{T1, T2, T3, T4, T5} <: AbstractMassSpectrum{T2, T4}

Concrete container for a single mass spectrum.

`T2` and `T4` encode the m/z and intensity unit types (each a `Unitful.Units` subtype
or `Nothing`). The remaining type parameters capture the stored numeric vectors and
metadata.

Fields include `mzvalues::AbstractVector{<:Real}`, `mzunit::Union{Unitful.Units, Nothing}`,
`intensities::AbstractVector{<:Real}`, `intensityunit::Union{Unitful.Units, Nothing}`,
and `attrs::NamedTuple`.

The constructor validates that m/z values are strictly increasing, positive, finite,
and that intensity values are finite with matching length.

See also
[`AbstractMassSpectrum`](@ref JuChrom.AbstractMassSpectrum),
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`attrs`](@ref JuChrom.attrs(::AbstractMassSpectrum)),
[`intensities`](@ref JuChrom.intensities(::AbstractMassSpectrum{<:Any, Nothing})),
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassSpectrum)),
[`mzcount`](@ref JuChrom.mzcount(::AbstractMassSpectrum)),
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassSpectrum)),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassSpectrum{Nothing, <:Any})).
"""
struct MassSpectrum{
    T1<:AbstractVector{<:Real},
    T2<:Union{Nothing, Unitful.Units},
    T3<:AbstractVector{<:Real},
    T4<:Union{Nothing, Unitful.Units},
    T5<:NamedTuple
    } <: AbstractMassSpectrum{T2, T4}

    mzvalues::T1
    mzunit::T2
    intensities::T3
    intensityunit::T4
    attrs::T5

    function MassSpectrum{T1, T2, T3, T4, T5}(
        mzvalues::T1,
        mzunit::T2,
        intensities::T3,
        intensityunit::T4,
        attrs::T5
        ) where {
            T1<:AbstractVector{<:Real},
            T2<:Union{Nothing, Unitful.Units},
            T3<:AbstractVector{<:Real},
            T4<:Union{Nothing, Unitful.Units},
            T5<:NamedTuple
        }

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

        new{T1, T2, T3, T4, T5}(
            mzvalues, mzunit,
            intensities, intensityunit,
            attrs)
    end
end

"""
    MassSpectrum(
        mzvalues::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}},
        intensities::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}};
        attrs::NamedTuple=NamedTuple()
    )

Construct a `MassSpectrum` from m/z and intensity vectors that may be unitless or
unitful. If units are present, they are stored and the numeric values are stripped.

Inputs must have strictly increasing, positive, finite m/z values, finite intensities,
and matching lengths.

See also
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassSpectrum{Nothing, <:Any})),
[`intensities`](@ref JuChrom.intensities(::AbstractMassSpectrum{<:Any, Nothing})),
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassSpectrum)),
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassSpectrum)),
[`attrs`](@ref JuChrom.attrs(::AbstractMassSpectrum)).

# Examples
```jldoctest
julia> ms = MassSpectrum([100.0, 150.0, 200.0]u"Th", [10.0, 20.0, 5.0]u"pA");

julia> mzunit(ms)
Th

julia> intensityunit(ms)
pA

julia> mzvalues(ms) == [100.0, 150.0, 200.0]u"Th"
true
```
"""
@inline MassSpectrum(
        mzvalues::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}},
        intensities::AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}};
        attrs::NamedTuple=NamedTuple()) = begin
    mz_v, mz_u = JuChrom.strip_units_checked(mzvalues, "m/z values")
    ints_v,  ints_u  = JuChrom.strip_units_checked(intensities,  "intensities")
    MassSpectrum{typeof(mz_v), typeof(mz_u),
                 typeof(ints_v),  typeof(ints_u),
                 typeof(attrs)}(
        mz_v, mz_u, ints_v, ints_u, attrs)
end

"""
    MassSpectrum(
        mzvalues::AbstractVector{<:Real},
        mzunit::Union{Nothing, Unitful.Units},
        intensities::AbstractVector{<:Real},
        intensityunit::Union{Nothing, Unitful.Units};
        attrs::NamedTuple=NamedTuple()
    )

Construct a `MassSpectrum` from unitless numeric vectors and explicit unit metadata.

See also
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassSpectrum{Nothing, <:Any})),
[`intensities`](@ref JuChrom.intensities(::AbstractMassSpectrum{<:Any, Nothing})),
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassSpectrum)),
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassSpectrum)),
[`attrs`](@ref JuChrom.attrs(::AbstractMassSpectrum)).

# Examples
```jldoctest
julia> ms = MassSpectrum([100.0, 150.0], u"Th", [10.0, 20.0], u"pA"; attrs=(scan_id=7,));

julia> attrs(ms).scan_id == 7
true

julia> intensities(ms, unit=u"nA") == [0.01, 0.02]u"nA"
true
```
"""
@inline MassSpectrum(
    mzvalues::AbstractVector{<:Real},
    mzunit::Union{Nothing, Unitful.Units},
    intensities::AbstractVector{<:Real},
    intensityunit::Union{Nothing, Unitful.Units};
    attrs::NamedTuple=NamedTuple()
) = MassSpectrum{typeof(mzvalues), typeof(mzunit),
                 typeof(intensities), typeof(intensityunit),
                 typeof(attrs)}(
    mzvalues, mzunit, intensities, intensityunit, attrs)

"""
    attrs(ms::AbstractMassSpectrum)

Return the metadata `NamedTuple` associated with a mass spectrum.

See also
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassSpectrum{Nothing, <:Any})),
[`intensities`](@ref JuChrom.intensities(::AbstractMassSpectrum{<:Any, Nothing})).

# Examples
```jldoctest
julia> ms = MassSpectrum([100.0, 150.0], [10.0, 20.0]; attrs=(sample="QC",));

julia> attrs(ms).sample == "QC"
true
```
"""
attrs(ms::AbstractMassSpectrum) = ms.attrs

"""
    intensityunit(ms::AbstractMassSpectrum)

Return the intensity unit for a mass spectrum, or `nothing` if unspecified.

See also
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`intensities`](@ref JuChrom.intensities(::AbstractMassSpectrum{<:Any, Nothing})).

# Examples
```jldoctest
julia> ms = MassSpectrum([100.0, 150.0], [10.0, 20.0]);

julia> intensityunit(ms) === nothing
true
```
"""
intensityunit(ms::AbstractMassSpectrum) = ms.intensityunit

"""
    intensities(ms::AbstractMassSpectrum; unit=nothing)

Return the intensity values from a mass spectrum.

If `unit` is specified, intensities are converted with `Unitful.uconvert`. If the
spectrum is unitless and a unit is requested, an error is thrown.

See also
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassSpectrum)),
[`mzcount`](@ref JuChrom.mzcount(::AbstractMassSpectrum)),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassSpectrum{Nothing, <:Any})).

# Examples
```jldoctest
julia> ms = MassSpectrum([100.0, 150.0], [10.0, 20.0]u"pA");

julia> intensities(ms) == [10.0, 20.0]u"pA"
true

julia> intensities(ms, unit=u"nA") == [0.01, 0.02]u"nA"
true
```
"""
@inline function intensities(
    ms::AbstractMassSpectrum{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(ms.intensities, unit, "intensities")
end

@inline function intensities(
    ms::AbstractMassSpectrum{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(ms.intensities, ms.intensityunit, unit)
end

"""
    mzcount(ms::AbstractMassSpectrum) -> Integer

Return the number of m/z points in the spectrum.

See also
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassSpectrum{Nothing, <:Any})),
[`intensities`](@ref JuChrom.intensities(::AbstractMassSpectrum{<:Any, Nothing})).

# Examples
```jldoctest
julia> ms = MassSpectrum([100.0, 150.0, 200.0], [10.0, 20.0, 5.0]);

julia> mzcount(ms)
3
```
"""
mzcount(ms::AbstractMassSpectrum) = length(ms.mzvalues)

"""
    mzunit(ms::AbstractMassSpectrum)

Return the m/z unit for a mass spectrum, or `nothing` if unspecified.

See also
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassSpectrum{Nothing, <:Any})).

# Examples
```jldoctest
julia> ms = MassSpectrum([100.0, 150.0], [10.0, 20.0]);

julia> mzunit(ms) === nothing
true
```
"""
mzunit(ms::AbstractMassSpectrum) = ms.mzunit

"""
    mzvalues(ms::AbstractMassSpectrum; unit=nothing)

Return the m/z values from a mass spectrum.

If `unit` is specified, values are converted with `Unitful.uconvert`. If the spectrum
is unitless and a unit is requested, an error is thrown.

See also
[`MassSpectrum`](@ref JuChrom.MassSpectrum),
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassSpectrum)),
[`intensities`](@ref JuChrom.intensities(::AbstractMassSpectrum{<:Any, Nothing})),
[`mzcount`](@ref JuChrom.mzcount(::AbstractMassSpectrum)).

# Examples
```jldoctest
julia> ms = MassSpectrum([100.0, 150.0]u"Th", [10.0, 20.0]);

julia> mzvalues(ms) == [100.0, 150.0]u"Th"
true

julia> mzvalues(ms, unit=u"kTh") == [0.1, 0.15]u"kTh"
true
```
"""
@inline function mzvalues(
    ms::AbstractMassSpectrum{Nothing, <:Any};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(ms.mzvalues, unit, "m/z values")
end

@inline function mzvalues(
    ms::AbstractMassSpectrum{<:Unitful.Units, <:Any};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(ms.mzvalues, ms.mzunit, unit)
end
