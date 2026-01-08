"""
    AbstractMassScanMatrix{R, M, I}

Abstract supertype for all matrix-based representations of mass spectrometry scan data.

This type is parameterized as follows:
- `R`: Type of the separation axis unit (a subtype of `Unitful.Units` or `Nothing`)
- `M`: Type of the m/z unit (a subtype of `Unitful.Units` or `Nothing`)
- `I`: Type of the intensity unit (a subtype of `Unitful.Units` or `Nothing`)

Concrete subtypes must provide, at a minimum:
- `retentions`: A vector of real values (with or without units) representing separation 
  coordinates (e.g. time, index, or position)
- `mz_values`: A strictly increasing vector of real-valued (with or without units) m/z 
  values (no duplicates)
- `intensities`: A matrix of real-valued (with or without units) signal intensities; 
  each row corresponds to a scan, each column to an m/z value
- `level`: The MS level of the scans (e.g., `1` for MS1, `2` for MS2)
- `instrument`, `acquisition`, `sample`, `user`: structured metadata as `NamedTuple`s
- `extras`: Optional unstructured metadata (`Dict{String, Any}`)

Subtypes may include additional fields if needed.
"""
abstract type AbstractMassScanMatrix{R, M, I} end

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
            "m/z values must be strictly increasing (no duplicates allowed)."))  # consider removing this constraint to reflect different ms scan directions

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
        level::Integer = 1,
        instrument::NamedTuple = NamedTuple(),
        acquisition::NamedTuple = NamedTuple(),
        user::NamedTuple = NamedTuple(),
        sample::NamedTuple = NamedTuple(),
        extras::AbstractDict = Dict{String, Any}()
    )

Construct a `MassScanMatrix` from already unitless numeric arrays and explicit unit fields.

This method is intended for advanced use, when you have already separated units from 
values (e.g., after pre-processing or deserialization). All arrays must be strictly 
numeric (`Real`), and unit arguments must be either a compatible `Unitful.Units` object 
or `nothing` if unitless.

# Arguments
- `retention_unitfree`: Vector of retention/separation coordinates (unitless, non-empty, 
  finite).
- `retention_unit`: Unit for retention values (`Unitful.Units` or `nothing`).
- `mz_values_unitfree`: Vector of mass-to-charge (m/z) values (unitless, strictly 
  increasing, positive, non-empty, finite).
- `mz_unit`: Unit for m/z values (`Unitful.Units` or `nothing`).
- `intensities_unitfree`: 2D matrix of intensities (unitless, finite, size: 
  `length(retention_unitfree)` × `length(mz_values_unitfree)`).
- `intensity_unit`: Unit for intensities (`Unitful.Units` or `nothing`).
- `level`: MS level (default: `1`; must be ≥ 1).
- `instrument`: Optional instrument metadata as a `NamedTuple`.
- `acquisition`: Optional acquisition metadata as a `NamedTuple`.
- `user`: Optional user metadata as a `NamedTuple`.
- `sample`: Optional sample metadata as a `NamedTuple`.
- `extras`: Optional metadata as a `Dict{String, Any}`.

# Returns
A `MassScanMatrix` instance with all values stored as unitless arrays and units preserved 
in separate fields.

# Throws
- `ArgumentError` if:
    - `retention_unitfree` or `mz_values_unitfree` is empty or contains non-finite values
    - `mz_values_unitfree` is not strictly increasing or contains non-positive values
    - `intensities_unitfree` is empty or contains non-finite values
    - `level` is less than 1
- `DimensionMismatch` if the shape of `intensities_unitfree` does not match
  `length(retention_unitfree)` × `length(mz_values_unitfree)`

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
        retentions::AbstractVector{<:Union{Real, Quantity}},
        mz_values::AbstractVector{<:Union{Real, Quantity}},
        intensities::AbstractMatrix{<:Union{Real, Quantity}};
        level::Integer = 1,
        instrument::NamedTuple = NamedTuple(),
        acquisition::NamedTuple = NamedTuple(),
        user::NamedTuple = NamedTuple(),
        sample::NamedTuple = NamedTuple(),
        extras::AbstractDict = Dict{String, Any}()
    )

Construct a `MassScanMatrix` representing a collection of aligned mass spectrometry scans.

All inputs may contain raw numeric values or `Unitful.Quantity` values. If units are 
present, they are stripped and stored separately. Each input must either be entirely 
unitless or use consistent units across its values.

# Arguments
- `retentions`: A vector of separation coordinates (e.g. times, positions, or indices);
  must be non-empty and finite; may include units.
- `mz_values`: A strictly increasing, positive, non-empty vector of mass-to-charge (m/z)
  values; may include units.
- `intensities`: A 2D matrix of intensities, where rows match `retentions` and columns
  match `mz_values`; all values must be finite and units (if any) must be consistent.
- `level`: MS level (default: `1`; must be ≥ 1).
- `instrument`: Optional instrument metadata as a `NamedTuple`.
- `acquisition`: Optional acquisition metadata as a `NamedTuple`.
- `user`: Optional user metadata as a `NamedTuple`.
- `sample`: Optional sample metadata as a `NamedTuple`.
- `extras`: Optional metadata as a `Dict{String, Any}`.

# Returns
A `MassScanMatrix` instance with values stored unitless, and units preserved in separate
fields.

# Throws
- `ArgumentError` if:
  - `retentions` is empty or contains non-finite values
  - `mz_values` is empty, non-strictly increasing, non-positive, or non-finite
  - `intensities` is empty or contains non-finite values
  - units (if present) are inconsistent across any of the inputs
  - `level` is less than 1
- `DimensionMismatch` if the shape of `intensities` does not match
  `length(retentions)` × `length(mz_values)`

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
        T1<:AbstractVector{<:Union{Real, Quantity{<:Real}}},
        T2<:AbstractVector{<:Union{Real, Quantity{<:Real}}},
        T3<:AbstractMatrix{<:Union{Real, Quantity{<:Real}}},
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

function Base.:(-)(msm₁::MassScanMatrix, msm₂::MassScanMatrix)
    retentions(msm₁) == retentions(msm₂)    || throw(DimensionMismatch(
        "MassScanMatrix differ in their retention coordinates."))
    mzvalues(msm₁)   == mzvalues(msm₂)      || throw(DimensionMismatch(
        "MassScanMatrix differ in their m/z values."))
    level(msm₁)      == level(msm₂)         || throw(DimensionMismatch(
        "MassScanMatrix differ in their MS levels."))
    instrument(msm₁) == instrument(msm₂)    || throw(DimensionMismatch(
        "MassScanMatrix differ in the instrument data"))
    acquisition(msm₁) == acquisition(msm₂)  || throw(DimensionMismatch(
        "MassScanMatrix differ in the acquisition data"))
    user(msm₁)       == user(msm₂)          || throw(DimensionMismatch(
        "MassScanMatrix differ in the user data"))
    sample(msm₁)     == sample(msm₂)        || throw(DimensionMismatch(
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
