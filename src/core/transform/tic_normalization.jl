# ── ticnormalize ──────────────────────────────────────────────────────────────────────────

"""
    ticnormalize(msm::MassScanMatrix) -> MassScanMatrix
    ticnormalize(vmsm::AbstractVarianceMassScanMatrix) -> VarianceMassScanMatrix

Normalize a mass-scan matrix by the sum of positive total-ion-chromatogram (TIC) values.

For an intensity matrix `X`, the TIC is the row-wise sum of intensities. `ticnormalize`
computes

```julia
tic = vec(sum(X; dims=2))
scale = sum(t for t in tic if t > 0)
```

and returns intensities `X ./ scale`. All cells are scaled, including zero and negative
cells, but only positive TIC values contribute to `scale`. The returned intensities are
dimensionless, so `intensityunit` is `nothing`. Retention coordinates, m/z values, MS
level, metadata, and extras are preserved.

For `AbstractVarianceMassScanMatrix`, `ticnormalize` returns a `VarianceMassScanMatrix`.
Variances are propagated by treating the observed TIC sum as a fixed scaling factor:

```julia
V_normalized = V ./ scale^2
```

This deliberately does not propagate uncertainty in `scale`; doing so would induce
covariances between all normalized cells, while `VarianceMassScanMatrix` stores only
per-cell variances. With this fixed-scaling rule, `clr(ticnormalize(vmsm))` has the same
CLR intensities and propagated CLR variances as `clr(vmsm)`, provided all intensities are
strictly positive. If zeros or negative values are present and the result is intended for
CLR, call [`replacecensored`](@ref JuChrom.replacecensored) before `ticnormalize`.

Throws `ArgumentError` if the input has no finite positive TIC sum.
"""
function ticnormalize(msm::MassScanMatrix)
    scale = tic_normalization_scale(rawintensities(msm))
    normalized = rawintensities(msm) ./ scale

    MassScanMatrix(
        copy(rawretentions(msm)),
        retentionunit(msm),
        copy(rawmzvalues(msm)),
        mzunit(msm),
        normalized,
        nothing;
        level=level(msm),
        instrument=deepcopy(instrument(msm)),
        acquisition=deepcopy(acquisition(msm)),
        user=deepcopy(user(msm)),
        sample=deepcopy(sample(msm)),
        extras=deepcopy(extras(msm))
    )
end

function ticnormalize(vmsm::AbstractVarianceMassScanMatrix)
    scale = tic_normalization_scale(rawintensities(vmsm))
    normalized_msm = tic_normalize_msm(vmsm, scale)
    normalized_variances =
        rawvariances(vmsm; unit=default_varianceunit(vmsm)) ./ abs2(scale)

    VarianceMassScanMatrix(normalized_msm, normalized_variances, nothing)
end

function tic_normalize_msm(msm::AbstractMassScanMatrix, scale::Real)
    MassScanMatrix(
        copy(rawretentions(msm)),
        retentionunit(msm),
        copy(rawmzvalues(msm)),
        mzunit(msm),
        rawintensities(msm) ./ scale,
        nothing;
        level=level(msm),
        instrument=deepcopy(instrument(msm)),
        acquisition=deepcopy(acquisition(msm)),
        user=deepcopy(user(msm)),
        sample=deepcopy(sample(msm)),
        extras=deepcopy(extras(msm))
    )
end

function tic_normalization_scale(intensities::AbstractMatrix{<:Real})
    tic = vec(sum(intensities; dims=2))
    scale = sum(t for t in tic if isfinite(t) && t > zero(t); init=zero(eltype(tic)))
    isfinite(scale) && scale > zero(scale) || throw(ArgumentError(
        "ticnormalize requires at least one finite positive TIC value"))

    scale
end
