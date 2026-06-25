
# ── whiten ────────────────────────────────────────────────────────────────────────────────

function isdimensionlessvarianceunit(unit::Nothing)
    unit === nothing
end

function isdimensionlessvarianceunit(unit::Unitful.Units)
    Unitful.dimension(one(Float64) * unit) == Unitful.dimension(one(Float64))
end

"""
    whiten(vmsm::AbstractVarianceMassScanMatrix; sigmafloor=:auto, floorquantile=0.05)
    whiten(vmsm::AbstractVarianceMassScanMatrix, sigmafloor::Real)
    whiten(vmsm::AbstractVarianceMassScanMatrix, sigmafloor::Unitful.AbstractQuantity)
    whiten(vmsm::AbstractVarianceMassScanMatrix, :auto; floorquantile=0.05)

Whiten a variance mass scan matrix by dividing each intensity by its propagated standard
deviation, with a standard-deviation floor.

For each cell, `denominator_variance = max(variance, sigmafloor^2)`. The returned
intensity is `intensity / sqrt(denominator_variance)`, and the returned variance is
`variance / denominator_variance`. The returned `VarianceMassScanMatrix` keeps the
retention axis, m/z axis, level, metadata, and extras of the input. Both `intensityunit`
and `varianceunit` of the returned object are `nothing`.

By default, `sigmafloor = :auto` infers a positive floor from the `floorquantile` quantile
of the positive propagated standard deviations. Dimensionless data produce a dimensionless
floor. Physical data infer the floor in `intensityunit(vmsm)`.

A bare real `sigmafloor` is accepted only when `varianceunit(vmsm)` is `nothing` or
dimensionless, for example after [`clr`](@ref JuChrom.clr). For physical intensity data,
pass a unitful floor such as `0.1u"pA"`; intensities and variances are converted to that
floor unit before units are stripped.

# Throws
- `ArgumentError`: if `sigmafloor <= 0`, `floorquantile` is outside `(0, 1]`, `:auto`
  cannot infer a positive floor, or a bare real floor is used for physical variances.
- `DimensionMismatch`: if the variance matrix shape does not match the intensity matrix.
"""
function whiten(
    vmsm::AbstractVarianceMassScanMatrix;
    sigmafloor=:auto,
    floorquantile::Real=0.05)

    sigmafloor === :auto && return whiten(vmsm, :auto; floorquantile=floorquantile)
    whiten(vmsm, sigmafloor)
end

function whiten(
    vmsm::AbstractVarianceMassScanMatrix,
    sigmafloor::Symbol;
    floorquantile::Real=0.05)

    sigmafloor === :auto || throw(ArgumentError(
        "sigmafloor must be positive, a Unitful quantity, or :auto."))
    whiten(vmsm, inferwhitensigmafloor(vmsm, floorquantile))
end

function whiten(vmsm::AbstractVarianceMassScanMatrix, sigmafloor::Real)
    sigmafloor > 0 || throw(ArgumentError("sigmafloor must be positive."))
    isdimensionlessvarianceunit(varianceunit(vmsm)) || throw(ArgumentError(
        "whiten with a bare real sigmafloor requires varianceunit(vmsm) to be nothing " *
        "or dimensionless. Use a unitful sigmafloor for physical variances."))

    x = rawintensities(vmsm)
    variances = rawvariances(vmsm)
    whitenfromvalues(vmsm, x, variances, sigmafloor)
end

function whiten(vmsm::AbstractVarianceMassScanMatrix, sigmafloor::AbstractQuantity{<:Real})
    sigmafloor > zero(sigmafloor) || throw(ArgumentError("sigmafloor must be positive."))

    floorunit = Unitful.unit(sigmafloor)
    x = rawintensities(vmsm; unit=floorunit)
    variances = rawvariances(vmsm; unit=floorunit^2)
    whitenfromvalues(vmsm, x, variances, Unitful.ustrip(floorunit, sigmafloor))
end

function inferwhitensigmafloor(
    vmsm::AbstractVarianceMassScanMatrix,
    floorquantile::Real)

    q = Float64(floorquantile)
    isfinite(q) && 0 < q <= 1 || throw(ArgumentError(
        "floorquantile must be finite and in the interval (0, 1]."))

    if isdimensionlessvarianceunit(varianceunit(vmsm))
        return inferwhitensigmafloorvalue(rawvariances(vmsm), q)
    end

    intunit = intensityunit(vmsm)
    isnothing(intunit) && throw(ArgumentError(
        "cannot infer a unitful sigmafloor without intensityunit(vmsm)."))
    inferwhitensigmafloorvalue(rawvariances(vmsm; unit=intunit^2), q) * intunit
end

function inferwhitensigmafloorvalue(variances::AbstractMatrix{<:Real}, q::Real)
    sigmas = sqrt.(filter(v -> isfinite(v) && v > zero(v), vec(variances)))
    isempty(sigmas) && throw(ArgumentError(
        "cannot infer sigmafloor from variances without positive finite values."))

    sigmafloor = quantile(sigmas, q)
    isfinite(sigmafloor) && sigmafloor > 0 || throw(ArgumentError(
        "cannot infer a positive finite sigmafloor from the variance values."))
    sigmafloor
end

function whitenfromvalues(
    vmsm::AbstractVarianceMassScanMatrix,
    x::AbstractMatrix{<:Real},
    variances::AbstractMatrix{<:Real},
    sigmafloor::Real)

    size(x) == size(variances) || throw(
        DimensionMismatch("intensities and variances must have identical sizes"))

    denominator_variances = max.(variances, abs2(sigmafloor))
    whitened_intensities = x ./ sqrt.(denominator_variances)
    whitened_variances = variances ./ denominator_variances

    msm_out = MassScanMatrix(
        copy(rawretentions(vmsm)),
        retentionunit(vmsm),
        copy(rawmzvalues(vmsm)),
        mzunit(vmsm),
        whitened_intensities,
        nothing;
        level=level(vmsm),
        instrument=deepcopy(instrument(vmsm)),
        acquisition=deepcopy(acquisition(vmsm)),
        user=deepcopy(user(vmsm)),
        sample=deepcopy(sample(vmsm)),
        extras=deepcopy(extras(vmsm)),
    )

    VarianceMassScanMatrix(msm_out, whitened_variances, nothing)
end
