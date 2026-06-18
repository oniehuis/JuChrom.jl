"""
    AbstractVarianceMassScanMatrix{R, M, I, V} <: AbstractMassScanMatrix{R, M, I}

Abstract supertype for mass scan matrices that carry per-cell intensity variances.

`R`, `M`, and `I` have the same meaning as for `AbstractMassScanMatrix`: retention,
m/z, and intensity unit types. `V` is the variance unit type (`Unitful.Units` subtype
or `Nothing`). Concrete subtypes are expected to support the ordinary mass-scan matrix
interface and additionally define `variances`, `rawvariances`, and `varianceunit`.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix),
[`VarianceMassScanMatrix`](@ref JuChrom.VarianceMassScanMatrix),
[`variances`](@ref JuChrom.variances),
[`rawvariances`](@ref JuChrom.rawvariances),
[`varianceunit`](@ref JuChrom.varianceunit).
"""
abstract type AbstractVarianceMassScanMatrix{R, M, I, V} <:
    AbstractMassScanMatrix{R, M, I} end

"""
    VarianceMassScanMatrix{R, M, I, T, A, V} <:
        AbstractVarianceMassScanMatrix{R, M, I, V}

Matrix-based mass spectrometry scan data with a per-cell variance matrix.

`VarianceMassScanMatrix` wraps an `AbstractMassScanMatrix` and stores a variance matrix with
the same shape as `rawintensities(msm)`. It preserves the ordinary mass-scan matrix interface
through forwarding methods while adding variance-specific accessors.

The stored variance values are unitless numeric values. Their unit is stored separately in
`varianceunit`, following the same representation used by `MassScanMatrix` for intensities.
"""
struct VarianceMassScanMatrix{
    R,
    M,
    I,
    T<:AbstractMassScanMatrix{R, M, I},
    A<:AbstractMatrix{<:Real},
    V<:Union{Nothing, Unitful.Units}
    } <: AbstractVarianceMassScanMatrix{R, M, I, V}

    msm::T
    variances::A
    varianceunit::V

    function VarianceMassScanMatrix{R, M, I, T, A, V}(
        msm::T,
        variances::A,
        varianceunit::V,
        ) where {
            R,
            M,
            I,
            T<:AbstractMassScanMatrix{R, M, I},
            A<:AbstractMatrix{<:Real},
            V<:Union{Nothing, Unitful.Units}}

        !isempty(variances) || throw(ArgumentError("No variance value(s) provided."))
        size(variances) == size(rawintensities(msm)) || throw(DimensionMismatch(
            "Variance matrix shape must match the intensity matrix shape."))
        all(isfinite, variances) || throw(ArgumentError("All variances must be finite."))
        all(v -> v ≥ zero(v), variances) || throw(ArgumentError(
            "All variances must be nonnegative."))
        validate_varianceunit(msm, varianceunit)

        new{R, M, I, T, A, V}(msm, variances, varianceunit)
    end
end

Base.broadcastable(vmsm::VarianceMassScanMatrix) = Base.RefValue(vmsm)
Base.parent(vmsm::VarianceMassScanMatrix) = vmsm.msm

function default_varianceunit(msm::AbstractMassScanMatrix)
    unit = intensityunit(msm)
    isnothing(unit) && return nothing
    unit^2
end

function validate_varianceunit(msm::AbstractMassScanMatrix, varianceunit::Nothing)
    isnothing(intensityunit(msm)) || throw(ArgumentError(
        "varianceunit cannot be `nothing` when intensities have a unit."))
    nothing
end

function validate_varianceunit(msm::AbstractMassScanMatrix, varianceunit::Unitful.Units)
    int_unit = intensityunit(msm)
    isnothing(int_unit) && throw(ArgumentError(
        "varianceunit must be `nothing` when intensities are unitless."))

    expected = (one(Float64) * int_unit)^2
    observed = one(Float64) * varianceunit
    Unitful.dimension(observed) == Unitful.dimension(expected) || throw(
        Unitful.DimensionError(observed, expected))
    nothing
end

function VarianceMassScanMatrix(
    msm::T,
    variances::A,
    varianceunit::V,
    ) where {
        R,
        M,
        I,
        T<:AbstractMassScanMatrix{R, M, I},
        A<:AbstractMatrix{<:Real},
        V<:Union{Nothing, Unitful.Units}}

    VarianceMassScanMatrix{R, M, I, T, A, V}(msm, variances, varianceunit)
end

function VarianceMassScanMatrix(
    msm::AbstractMassScanMatrix,
    variances::AbstractMatrix{<:Union{Real, AbstractQuantity{<:Real}}},
    )

    variances_unitfree, inferred_unit = strip_units_checked(variances, "variances")
    varianceunit = isnothing(inferred_unit) ? default_varianceunit(msm) : inferred_unit
    VarianceMassScanMatrix(msm, variances_unitfree, varianceunit)
end

# ── Display ───────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, vmsm::VarianceMassScanMatrix)
    n_scans = scancount(vmsm)

    retentions_unitfree = rawretentions(vmsm)
    retention_unit = retentionunit(vmsm)
    retention_range = (minimum(retentions_unitfree), maximum(retentions_unitfree))

    mz_unitfree = rawmzvalues(vmsm)
    mz_unit = mzunit(vmsm)
    mz_range = (minimum(mz_unitfree), maximum(mz_unitfree))

    intensities_unitfree = rawintensities(vmsm)
    intensity_unit = intensityunit(vmsm)
    intensity_range = (minimum(intensities_unitfree), maximum(intensities_unitfree))

    variances_unitfree = rawvariances(vmsm)
    variance_unit = varianceunit(vmsm)
    variance_range = (minimum(variances_unitfree), maximum(variances_unitfree))

    ret_unit_str = isnothing(retention_unit) ? "unitless" : string(retention_unit)
    mz_unit_str = isnothing(mz_unit) ? "unitless" : string(mz_unit)
    int_unit_str = isnothing(intensity_unit) ? "unitless" : string(intensity_unit)
    var_unit_str = isnothing(variance_unit) ? "unitless" : string(variance_unit)

    metadata_sections = [
        ("Instrument", instrument(vmsm)),
        ("Acquisition", acquisition(vmsm)),
        ("User", user(vmsm)),
        ("Sample", sample(vmsm))
    ]

    non_empty_sections = filter(section -> !isempty(section[2]), metadata_sections)
    has_metadata = !isempty(non_empty_sections) || !isempty(extras(vmsm))

    scan_word = n_scans == 1 ? "scan" : "scans"
    println(io, "VarianceMassScanMatrix with $n_scans $scan_word")

    println(io, "├─ Retention:")
    println(io, "│  ├─ Range: $(retention_range[1]) to $(retention_range[2]) ($ret_unit_str)")
    println(io, "│  └─ Type: $(eltype(retentions_unitfree))")

    println(io, "├─ M/Z values:")
    println(io, "│  ├─ Range: $(mz_range[1]) to $(mz_range[2]) ($mz_unit_str)")
    println(io, "│  ├─ Total data points: $(length(mz_unitfree))")
    println(io, "│  └─ Type: $(eltype(mz_unitfree))")

    println(io, "├─ Intensity:")
    println(io, "│  ├─ Range: $(intensity_range[1]) to $(intensity_range[2]) ($int_unit_str)")
    println(io, "│  ├─ Type: $(eltype(intensities_unitfree))")
    println(io, "│  └─ MS Level: $(level(vmsm))")

    variance_prefix = has_metadata ? "├─" : "└─"
    println(io, "$variance_prefix Variance:")
    variance_sub_prefix = has_metadata ? "│" : " "
    println(io, "$variance_sub_prefix  ├─ Range: $(variance_range[1]) to $(variance_range[2]) ($var_unit_str)")
    println(io, "$variance_sub_prefix  └─ Type: $(eltype(variances_unitfree))")

    ScanSeriesDisplay.print_annotations(io, metadata_sections, extras(vmsm))
end

function Base.show(io::IO, ::MIME"text/plain", vmsm::VarianceMassScanMatrix)
    show(io, vmsm)
end
