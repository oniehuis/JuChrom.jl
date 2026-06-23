# ── clr ───────────────────────────────────────────────────────────────────────────────────

"""
    clr(vmsm::AbstractVarianceMassScanMatrix) -> VarianceMassScanMatrix

Apply the centered log-ratio (CLR) transformation to the intensity matrix of a variance
mass scan matrix and propagate the associated per-cell variances.

The input intensities must be strictly positive. `clr` does not replace zeros or negative
values; use [`replacecensored`](@ref JuChrom.replacecensored) before `clr` when censored
or baseline-corrected data contain non-positive cells.

The returned `VarianceMassScanMatrix` keeps the retention axis, m/z axis, level, metadata,
and extras of the input. CLR values are log-ratio coordinates, so both `intensityunit` and
`varianceunit` of the returned object are `nothing`.

Variance propagation uses the delta-method variance of `log(x)` followed by the finite
CLR centering correction:
`sigma2_clr[i] = sigma2_log[i] * (1 - 2/N) + sum(sigma2_log) / N^2`, where
`sigma2_log[i] = variance[i] / intensity[i]^2` and `N` is the number of cells.

# Throws
- `DomainError`: if any intensity is zero or negative.
- `DimensionMismatch`: if the variance matrix shape does not match the intensity matrix.
"""
function clr(vmsm::AbstractVarianceMassScanMatrix)
    x = rawintensities(vmsm)
    variances = rawvariances(vmsm; unit=default_varianceunit(vmsm))

    size(x) == size(variances) || throw(
        DimensionMismatch("intensities and variances must have identical sizes"))
    all(>(0), x) || throw(DomainError(
        x,
        "all intensities must be > 0; call replacecensored first if zeros or negative values are present",
    ))

    logx = log.(x)
    clr_intensities = logx .- mean(logx)

    sigma2_log = variances ./ abs2.(x)
    N = length(x)
    total_sigma2_log = sum(sigma2_log)
    clr_variances = @. sigma2_log * (1 - 2 / N) + total_sigma2_log / N^2

    msm_out = MassScanMatrix(
        copy(rawretentions(vmsm)),
        retentionunit(vmsm),
        copy(rawmzvalues(vmsm)),
        mzunit(vmsm),
        clr_intensities,
        nothing;
        level=level(vmsm),
        instrument=deepcopy(instrument(vmsm)),
        acquisition=deepcopy(acquisition(vmsm)),
        user=deepcopy(user(vmsm)),
        sample=deepcopy(sample(vmsm)),
        extras=deepcopy(extras(vmsm)),
    )

    VarianceMassScanMatrix(msm_out, clr_variances, nothing)
end
