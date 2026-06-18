# Alkane Series Finding

Use `findalkaneseries` to identify an n-alkane series in GC-MS data for retention-index
calibration. The production API is under active development.

The result can be used directly to fit a retention mapper:

```julia
asr = findalkanes(msm)
mapper = fitmap(asr)
```

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
