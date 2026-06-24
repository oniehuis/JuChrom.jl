
# ── whiten ────────────────────────────────────────────────────────────────────────────────

function isdimensionlessvarianceunit(unit::Nothing)
    unit === nothing
end

function isdimensionlessvarianceunit(unit::Unitful.Units)
    Unitful.dimension(one(Float64) * unit) == Unitful.dimension(one(Float64))
end

"""
    whiten(vmsm::AbstractVarianceMassScanMatrix, sigmafloor::Real) -> VarianceMassScanMatrix

Whiten a variance mass scan matrix by dividing each intensity by its propagated standard
deviation, with a standard-deviation floor.

This method is intended for dimensionless transformed data, typically the result of
[`clr`](@ref JuChrom.clr). The variance unit of `vmsm` must be `nothing` or dimensionless;
physical intensity variances such as `pA^2` are rejected so raw signal data are not
accidentally whitened before log-ratio transformation.

For each cell, `denominator_variance = max(variance, sigmafloor^2)`. The returned
intensity is `intensity / sqrt(denominator_variance)`, and the returned variance is
`variance / denominator_variance`. The returned `VarianceMassScanMatrix` keeps the
retention axis, m/z axis, level, metadata, and extras of the input. Both `intensityunit`
and `varianceunit` of the returned object are `nothing`.

# Throws
- `ArgumentError`: if `sigmafloor <= 0` or `varianceunit(vmsm)` is not dimensionless.
- `DimensionMismatch`: if the variance matrix shape does not match the intensity matrix.
"""
function whiten(vmsm::AbstractVarianceMassScanMatrix, sigmafloor::Real)
    sigmafloor > 0 || throw(ArgumentError("sigmafloor must be positive."))
    isdimensionlessvarianceunit(varianceunit(vmsm)) || throw(ArgumentError(
        "whiten requires varianceunit(vmsm) to be nothing or dimensionless."))

    x = rawintensities(vmsm)
    variances = rawvariances(vmsm)
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
