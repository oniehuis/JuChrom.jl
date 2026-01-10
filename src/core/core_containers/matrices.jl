"""
    AbstractMassScanMatrix{R, M, I} <: Any

Abstract supertype for all matrix-based representations of mass spectrometry scan data.

`R` is the separation axis unit type (`Unitful.Units` subtype or `Nothing`), `M` is the
m/z unit type (`Unitful.Units` subtype or `Nothing`), and `I` is the intensity unit type
(`Unitful.Units` subtype or `Nothing`).

Concrete subtypes are expected to define `retentions`, `mz_values`, and `intensities`
arrays, plus `level` and metadata fields (`instrument`, `acquisition`, `sample`, `user`)
as `NamedTuple`s, with `extras::Dict{String, Any}` for unstructured metadata. Subtypes may
include additional fields as needed.

See also: [`MassScanMatrix`](@ref).
"""
abstract type AbstractMassScanMatrix{R, M, I} end

"""
    MassScanMatrix{R, M, I} <: AbstractMassScanMatrix{R, M, I}

Concrete matrix-based representation of mass spectrometry scan data, stored as aligned
retention and m/z axes with an intensity matrix. `MassScanMatrix` is a subtype of
`AbstractMassScanMatrix`.

`R` is the separation axis unit type (`Unitful.Units` subtype or `Nothing`), `M` is the
m/z unit type (`Unitful.Units` subtype or `Nothing`), and `I` is the intensity unit type
(`Unitful.Units` subtype or `Nothing`).

Fields include `retentions::AbstractVector{<:Real}`, `mz_values::AbstractVector{<:Real}`,
`intensities::AbstractMatrix{<:Real}`, optional units stored in 
`retention_unit::Union{Unitful.Units, Nothing}`, `mz_unit::Union{Unitful.Units, Nothing}`,
and `intensity_unit::Union{Unitful.Units, Nothing}`, an MS `level::Integer`, and metadata
in `instrument::NamedTuple`, `acquisition::NamedTuple`, `user::NamedTuple`,
`sample::NamedTuple`, and `extras::Dict{String, Any}`.
"""
struct MassScanMatrix{
    T1<:AbstractVector{<:Real},
    T2<:Union{Nothing, Unitful.Units},
    T3<:AbstractVector{<:Real},
    T4<:Union{Nothing, Unitful.Units},
    T5<:AbstractMatrix{<:Real},
    T6<:Union{Nothing, Unitful.Units},
    T7<:Integer,
    T8<:NamedTuple,
    T9<:NamedTuple,
    T10<:NamedTuple,
    T11<:NamedTuple,
    T12<:Dict{String, Any}
    } <: AbstractMassScanMatrix{T2, T4, T6}

    retentions::T1
    retention_unit::T2
    mz_values::T3
    mz_unit::T4
    intensities::T5
    intensity_unit::T6
    level::T7
    instrument::T8
    acquisition::T9
    user::T10
    sample::T11
    extras::T12

    # Inner constructor with validation to ensure consistent, valid data
    function MassScanMatrix{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12}(
        retentions::T1,
        retention_unit::T2,
        mz_values::T3,
        mz_unit::T4,
        intensities::T5,
        intensity_unit::T6,
        level::T7,
        instrument::T8,
        acquisition::T9,
        user::T10,
        sample::T11,
        extras::T12,
        ) where {
            T1<:AbstractVector{<:Real},
            T2<:Union{Nothing, Unitful.Units},
            T3<:AbstractVector{<:Real},
            T4<:Union{Nothing, Unitful.Units},
            T5<:AbstractMatrix{<:Real},
            T6<:Union{Nothing, Unitful.Units},
            T7<:Integer,
            T8<:NamedTuple,
            T9<:NamedTuple,
            T10<:NamedTuple,
            T11<:NamedTuple,
            T12<:Dict{String, Any}}

        # Check unitfree retentions
        !isempty(retentions) || throw(ArgumentError("No retention(s) provided."))
        all(isfinite, retentions) || throw(ArgumentError("All retentions must be finite."))

        # Check unitfree m/z values
        !isempty(mz_values) || throw(ArgumentError("No m/z value(s) provided."))
        all(isfinite, mz_values) || throw(ArgumentError("All m/z values must be finite."))
        all(mz -> mz > 0, mz_values) || throw(
            ArgumentError("All m/z values must be greater than zero."))
        all(diff(mz_values) .> 0) || throw(ArgumentError(
            "m/z values must be strictly increasing (no duplicates allowed)."))

        # Check unitfree intensities
        !isempty(intensities) || throw(ArgumentError("No intensity value(s) provided."))
        all(isfinite, intensities) || throw(ArgumentError(
            "All intensities must be finite."))

        # Confirm the number of retentions matches the number of rows in intensity matrix
        length(retentions) == size(intensities, 1) || throw(
           ArgumentError("Retention count does not match the row count of intensities."))

        # Confirm the number of m/z values matches the number of columns in intensity matrix
        length(mz_values) == size(intensities, 2) || throw(
           ArgumentError("m/z value count does not match the column count of intensities."))

        # Validate that scan level is positive
        level ≥ 1 || throw(ArgumentError("`level` must be ≥ 1"))

        # Construct the MassScanMatrix instance
        new{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12}(
            retentions, retention_unit,
            mz_values, mz_unit,
            intensities, intensity_unit,
            level,
            instrument, acquisition, user, sample,
            extras)
    end
end

"""
    MassScanMatrix(
        retention_unitfree::AbstractVector{<:Real},
        retention_unit::Union{Nothing, Unitful.Units},
        mz_values_unitfree::AbstractVector{<:Real},
        mz_unit::Union{Nothing, Unitful.Units},
        intensities_unitfree::AbstractMatrix{<:Real},
        intensity_unit::Union{Nothing, Unitful.Units};
        level::Integer=1,
        instrument::NamedTuple=NamedTuple(),
        acquisition::NamedTuple=NamedTuple(),
        user::NamedTuple=NamedTuple(),
        sample::NamedTuple=NamedTuple(),
        extras::AbstractDict=Dict{String, Any}()
    )

Construct a `MassScanMatrix` from already unitless numeric arrays and explicit unit fields.

This method is intended for advanced use, when you have already separated units from values 
(e.g. after pre-processing or deserialization). All arrays must be strictly numeric 
(`Real`), and unit arguments must be either a compatible `Unitful.Units` object or `nothing` 
if unitless.

`retention_unitfree` provides the unitless separation coordinates and `retention_unit`
their unit (or `nothing`). `mz_values_unitfree` provides unitless m/z values and `mz_unit`
their unit (or `nothing`). `intensities_unitfree` provides the unitless intensity matrix
and `intensity_unit` its unit (or `nothing`). `level` sets the MS level (default `1`), and
`instrument`, `acquisition`, `user`, `sample`, and `extras` carry optional metadata.
Returns a `MassScanMatrix` with unitless arrays and units stored separately.

Throws `ArgumentError` if retention or m/z values are empty or non-finite, if m/z values are
not strictly increasing or are non-positive, if intensities are empty or non-finite, or if
`level < 1`. Throws `DimensionMismatch` if the intensity matrix shape does not match
`length(retention_unitfree)` × `length(mz_values_unitfree)`.

# Examples
```jldoctest
julia> ret = [1.0, 2.0];
       mzs = [100.0, 200.0, 300.0];
       ints = [1.0 2.0 3.0
               4.0 5.0 6.0];

julia> msm = MassScanMatrix(ret, u"s", mzs, nothing, ints, nothing);

julia> msm.retention_unit
s
julia> isnothing(msm.mz_unit)
true
julia> isnothing(msm.intensity_unit)
true
julia> size(msm.intensities)
(2, 3)
```
"""
function MassScanMatrix(
    retention_unitfree::T1,
    retention_unit::T2,
    mz_values_unitfree::T3,
    mz_unit::T4,
    intensities_unitfree::T5,
    intensity_unit::T6;
    level::T7=1,
    instrument::T8=NamedTuple(),
    acquisition::T9=NamedTuple(),
    user::T10=NamedTuple(),
    sample::T11=NamedTuple(),
    extras::T12=Dict{String, Any}()) where {
        T1<:AbstractVector{<:Real},
        T2<:Union{Nothing, Unitful.Units},
        T3<:AbstractVector{<:Real},
        T4<:Union{Nothing, Unitful.Units},
        T5<:AbstractMatrix{<:Real},
        T6<:Union{Nothing, Unitful.Units},
        T7<:Integer,
        T8<:NamedTuple,
        T9<:NamedTuple,
        T10<:NamedTuple,
        T11<:NamedTuple,
        T12<:AbstractDict{<:AbstractString, <:Any}}

    # Convert metadata to Dict{String, Any}
    converted_extras = Dict{String, Any}(string(k) => v for (k, v) in extras)

    # Call the inner constructor with all provided and default arguments
    MassScanMatrix{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, Dict{String, Any}}(
        retention_unitfree, retention_unit,
        mz_values_unitfree, mz_unit,
        intensities_unitfree, intensity_unit,
        level,
        instrument, acquisition, user, sample,
        converted_extras)
end

"""
    MassScanMatrix(
        retentions::AbstractVector{<:Union{Real, AbstractQuantity}},
        mz_values::AbstractVector{<:Union{Real, AbstractQuantity}},
        intensities::AbstractMatrix{<:Union{Real, AbstractQuantity}};
        level::Integer=1,
        instrument::NamedTuple=NamedTuple(),
        acquisition::NamedTuple=NamedTuple(),
        user::NamedTuple=NamedTuple(),
        sample::NamedTuple=NamedTuple(),
        extras::AbstractDict=Dict{String, Any}()
    )

Construct a `MassScanMatrix` representing a collection of aligned mass spectrometry scans.

All inputs may contain raw numeric values or `Unitful.AbstractQuantity` values. If units are 
present, they are stripped and stored separately. Each input must either be entirely 
unitless or use consistent units across its values.

`retentions` provides separation coordinates, `mz_values` provides m/z values, and
`intensities` provides the intensity matrix; each may be unitless or unitful with
consistent units. `level` sets the MS level (default `1`), and `instrument`, `acquisition`,
`user`, `sample`, and `extras` carry optional metadata. Returns a `MassScanMatrix` with
unitless arrays and units stored separately.

Throws `ArgumentError` if retentions are empty or non-finite, if m/z values are empty,
non-strictly increasing, non-positive, or non-finite, if intensities are empty or contain
non-finite values, if units are inconsistent, or if `level < 1`. Throws
`DimensionMismatch` if the intensity matrix shape does not match
`length(retentions)` × `length(mz_values)`.

See also: [`AbstractMassScanMatrix`](@ref), [`retentions`](@ref), [`rawretentions`](@ref), 
[`retentionunit`](@ref), [`mzvalues`](@ref), [`rawmzvalues`](@ref), [`mzunit`](@ref), 
[`intensities`](@ref), [`rawintensities`](@ref), [`intensityunit`](@ref), [`level`](@ref),
[`instrument`](@ref), [`acquisition`](@ref), [`user`](@ref), [`sample`](@ref), 
[`extras`](@ref).

# Examples
```jldoctest
julia> ret = [1.0, 2.0]u"s"
       mzs = [100.0, 200.0, 300.0];
       ints = [1.0 2.0 3.0
               4.0 5.0 6.0];

julia> msm = MassScanMatrix(ret, mzs, ints);

julia> msm.retention_unit
s

julia> isnothing(msm.mz_unit)
true

julia> isnothing(msm.intensity_unit)
true

julia> size(msm.intensities)
(2, 3)

julia> msm.level
1

julia> msm.extras
Dict{String, Any}()
```
"""
function MassScanMatrix(
    retentions::T1,
    mz_values::T2,
    intensities::T3;
    level::T4=1,
    instrument::T5=NamedTuple(),
    acquisition::T6=NamedTuple(),
    user::T7=NamedTuple(),
    sample::T8=NamedTuple(),
    extras::T9=Dict{String, Any}()) where {
        T1<:AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}},
        T2<:AbstractVector{<:Union{Real, AbstractQuantity{<:Real}}},
        T3<:AbstractMatrix{<:Union{Real, AbstractQuantity{<:Real}}},
        T4<:Integer,
        T5<:NamedTuple,
        T6<:NamedTuple,
        T7<:NamedTuple,
        T8<:NamedTuple,
        T9<:AbstractDict{<:AbstractString, <:Any}}

    retention_unitfree, retention_unit = strip_units_checked(retentions, "retentions")
    mz_values_unitfree, mz_unit = strip_units_checked(mz_values, "m/z values")
    intensities_unitfree, intensity_unit = strip_units_checked(intensities, "intensities")

    # Convert metadata to Dict{String, Any}
    converted_extras = Dict{String, Any}(string(k) => v for (k, v) in extras)

    # Call the inner constructor with all provided and default arguments
    MassScanMatrix{
        typeof(retention_unitfree), typeof(retention_unit),
        typeof(mz_values_unitfree), typeof(mz_unit),
        typeof(intensities_unitfree), typeof(intensity_unit),
        typeof(level),
        typeof(instrument), typeof(acquisition), typeof(user), typeof(sample),
        Dict{String, Any}}(
        retention_unitfree, retention_unit,
        mz_values_unitfree, mz_unit,
        intensities_unitfree, intensity_unit,
        level,
        instrument, acquisition, user, sample,
        converted_extras)
end

"""
    Base.:(-)(msm₁::MassScanMatrix, msm₂::MassScanMatrix)

Return a `MassScanMatrix` whose intensities are the elementwise difference
`intensities(msm₁) .- intensities(msm₂)`, after verifying that retention coordinates, m/z
values, MS level, and metadata match.

Throws `DimensionMismatch` if any of `retentions`, `mz_values`, `level`, or metadata
(`instrument`, `acquisition`, `user`, `sample`, `extras`) differ between the inputs.
"""
function Base.:(-)(msm₁::MassScanMatrix, msm₂::MassScanMatrix)
    retentions(msm₁) == retentions(msm₂) || throw(DimensionMismatch(
        "MassScanMatrix differ in their retention coordinates."))
    mzvalues(msm₁) == mzvalues(msm₂) || throw(DimensionMismatch(
        "MassScanMatrix differ in their m/z values."))
    level(msm₁) == level(msm₂) || throw(DimensionMismatch(
        "MassScanMatrix differ in their MS levels."))
    instrument(msm₁) == instrument(msm₂) || throw(DimensionMismatch(
        "MassScanMatrix differ in the instrument data"))
    acquisition(msm₁) == acquisition(msm₂) || throw(DimensionMismatch(
        "MassScanMatrix differ in the acquisition data"))
    user(msm₁) == user(msm₂) || throw(DimensionMismatch(
        "MassScanMatrix differ in the user data"))
    sample(msm₁) == sample(msm₂) || throw(DimensionMismatch(
        "MassScanMatrix differ in the sample data"))
    isequal(extras(msm₁), extras(msm₂)) || throw(DimensionMismatch(
        "MassScanMatrix differ in their extras"))

    new_intensities = intensities(msm₁) .- intensities(msm₂)

    MassScanMatrix(
      deepcopy(retentions(msm₁)),
      deepcopy(mzvalues(msm₁)),
      new_intensities;
      level       = deepcopy(level(msm₁)),
      instrument  = deepcopy(instrument(msm₁)),
      acquisition = deepcopy(acquisition(msm₁)),
      user        = deepcopy(user(msm₁)),
      sample      = deepcopy(sample(msm₁)),
      extras      = deepcopy(extras(msm₁))
    )
end
