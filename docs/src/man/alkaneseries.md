# Alkane Series Finding

Use `findalkaneseries` to identify an n-alkane series in GC-MS data for retention-index
calibration. The production API is under active development.

The result reports whether the default alkane calibration points are mapper-ready:

```julia
asr = findalkanes(msm)
if asr.success
    mapper = fitmap(asr)
end
```

`asr.success` is `true` and `asr.status === :ok` only when the default calibration-point
selection yields at least three mapper anchors. A detected partial ladder can therefore
still return `asr.status === :too_few_mapper_steps`.

The alkane-ladder wrapper uses its own default smoothing strength (`λ=3e-9`), based on
leave-one-out prediction across ladder runs. Pass `λ=...` explicitly to override it.

To curate the calibration anchors without mutating the result, extract calibration points:

```julia
points = alkaneladdercalibrationpoints(asr; exclude=[27])
mapper = fitmap(points)
```

Manual points can be appended with `extra`; this is the preferred way to replace or extend
calibration anchors while keeping `AlkaneSeriesResult` as the unmodified audit trail of the
algorithm output.

```@docs
AlkaneStandard
AlkaneReferenceChannels
AlkaneChannelInfo
AlkaneAbundanceWindow
AlkaneAbundanceInfo
AlkaneMolecularIonContrast
AlkaneMolecularIonInfo
AlkaneSeriesResult
AlkaneLadderStep
AlkaneLadderCalibrationPoint
defaultalkanestandard
alkanereferencespectra
alkanereferencespectrum
findalkanes
findalkaneseries
alkaneladdersteps
alkaneladdercalibrationpoints
alkaneladderscanorder
alkaneladderapexes
alkaneladdermassspectra
```
