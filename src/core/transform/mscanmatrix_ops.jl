# ── subtractbaseline ─────────────────────────────────────────────────────────────────────

"""
    subtractbaseline(msm::MassScanMatrix, baseline::MassScanMatrix) -> MassScanMatrix
    subtractbaseline(
        vmsm::AbstractVarianceMassScanMatrix,
        baseline::MassScanMatrix
    ) -> VarianceMassScanMatrix

Subtract a same-grid baseline estimate from a mass-scan matrix.

The baseline must match the input retention coordinates, m/z values, retention unit, m/z
unit, intensity unit, intensity matrix shape, and MS level. Structured annotations and
`extras` are not compared. The returned `MassScanMatrix` preserves the input metadata and
`extras`, while the baseline object's metadata are used only for the compatibility checks
listed above.

For `AbstractVarianceMassScanMatrix` inputs, the baseline is subtracted from the parent
mass-scan matrix and the original variance matrix and variance unit are preserved.
Baseline-estimation uncertainty is not propagated.
"""
function subtractbaseline(msm::MassScanMatrix, baseline::MassScanMatrix)
    _subtractbaseline(msm, baseline)
end

function subtractbaseline(vmsm::AbstractVarianceMassScanMatrix, baseline::MassScanMatrix)
    corrected = _subtractbaseline(parent(vmsm), baseline)
    VarianceMassScanMatrix(
        corrected,
        copy(rawvariances(vmsm)),
        deepcopy(varianceunit(vmsm))
    )
end

function _subtractbaseline(msm::MassScanMatrix, baseline::MassScanMatrix)
    validate_baseline_compatible(msm, baseline)
    corrected_intensities = rawintensities(msm) .- rawintensities(baseline)

    MassScanMatrix(
        copy(rawretentions(msm)),
        retentionunit(msm),
        copy(rawmzvalues(msm)),
        mzunit(msm),
        corrected_intensities,
        intensityunit(msm);
        level=level(msm),
        instrument=deepcopy(instrument(msm)),
        acquisition=deepcopy(acquisition(msm)),
        user=deepcopy(user(msm)),
        sample=deepcopy(sample(msm)),
        extras=deepcopy(extras(msm))
    )
end

function validate_baseline_compatible(msm::MassScanMatrix, baseline::MassScanMatrix)
    retentionunit(msm) == retentionunit(baseline) || throw(DimensionMismatch(
        "baseline retention unit does not match input."))
    mzunit(msm) == mzunit(baseline) || throw(DimensionMismatch(
        "baseline m/z unit does not match input."))
    intensityunit(msm) == intensityunit(baseline) || throw(DimensionMismatch(
        "baseline intensity unit does not match input."))
    rawretentions(msm) == rawretentions(baseline) || throw(DimensionMismatch(
        "baseline retention coordinates do not match input."))
    rawmzvalues(msm) == rawmzvalues(baseline) || throw(DimensionMismatch(
        "baseline m/z values do not match input."))
    size(rawintensities(msm)) == size(rawintensities(baseline)) ||
        throw(DimensionMismatch("baseline intensity matrix shape does not match input."))
    level(msm) == level(baseline) || throw(DimensionMismatch(
        "baseline MS level does not match input."))

    nothing
end
