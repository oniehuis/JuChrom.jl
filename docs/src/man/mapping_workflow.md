# Retention mapping workflow

This page shows a complete retention-mapping workflow using one fitted mapper throughout.
The reference pages describe the individual APIs in more detail.

## Fitting a mapper

Fitting a mapper requires two vectors of equal length, with at least three points. One
vector contains the potentially unitful retentions in domain A, for example retention
times. The second vector contains the retentions in domain B. For alkane RT -> RI
calibration, domain B is the Kováts retention index, which is dimensionless.

For GC-MS n-alkane standards, these pairs usually do not need to be entered manually:
JuChrom can infer the ladder with
[`findalkanes`](@ref) and return an [`AlkaneSeriesResult`](@ref). See
[Alkane Series Finding](alkaneseries.md) for the ladder-detection details.

```@example retention_mapping_ladder
using JuChrom
using JuChrom.ChemStationMSLoader

# Load and preprocess an example n-alkane standard run
file = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
mss = load(ChemStationMS(file; mode=:ms))
tmss = retentiontrim(mss; start=4.1u"minute", stop=29.7u"minute")
btmss = binmzvalues(tmss)
msm = mscanmatrix(btmss)
nmsm = dwellnormalize(msm)

asr = findalkanes(nmsm)

asr.success || error("Failed to find a suitable ladder: $(asr.status)")

mapper = fitmap(asr)
```

This direct path is the usual workflow. `asr.success` is `true` only when the default
calibration-point selection yields at least three mapper-ready anchors. If it is `false`,
inspect `asr.status` and the ladder diagnostics before fitting; common causes are too few
accepted ladder steps, failed or ambiguous ladder detection, or anchor points that need
manual curation.

Only extract the calibration list when you want to inspect or curate it, for example to
remove a suspect step or add a manually confirmed step while keeping the original
`AlkaneSeriesResult` as the audit trail.

```@example retention_mapping_ladder
detected_points = alkaneladdercalibrationpoints(asr)
detected_C27 = only(filter(point -> point.ladderstep == 27, detected_points))

manual_C27 = AlkaneLadderCalibrationPoint(
    ladderstep=27,
    retention=detected_C27.retention,           # replace by the reviewed apex value
    retentionunit=detected_C27.retentionunit,   # unit of domain A; may be `nothing`
    retentionindex=2700.0,                      # domain B: unitless Kováts RI
    source=:manual,
)

points = alkaneladdercalibrationpoints(asr; exclude=[27], extra=[manual_C27])
mapper = fitmap(points)
```

`AlkaneLadderCalibrationPoint` is specific to alkane RT -> RI calibration, so the target
coordinate is stored as the unitless field `retentionindex`; there is no separate unit for
domain B. For a different mapping where domain B has a physical unit, use the general
[`fitmap`](@ref) call with two vectors instead. Manual points default to
`goodforcalibration=true`, meaning that the point is intended as a trusted mapper anchor.
If a point should not be fitted, leave it out of the curated list passed to `fitmap`.

If you need the coordinates as ordinary vectors, extract them from the calibration points.
This is mostly useful for reporting, plotting, or custom edits outside the
`AlkaneLadderCalibrationPoint` interface; otherwise pass `asr` or `points` directly to
[`fitmap`](@ref).

```julia
ordered = sort(points; by = point -> point.retention)
retention_values = [point.retention for point in ordered]
retention_unit = first(ordered).retentionunit

retention_times = isnothing(retention_unit) ?
    retention_values :
    retention_values .* retention_unit

kovats_indices = [point.retentionindex for point in ordered]
mapper = fitmap(retention_times, kovats_indices)
```

The extracted vectors must still satisfy the usual mapper requirements: both coordinate
vectors need at least three values, the retention coordinates and retention indices must be
strictly increasing after sorting by retention, and all retention coordinates must share
the same unit. If a curated list violates these constraints, remove or replace the
problematic anchors before fitting.

To keep this workflow self-contained, the remaining examples use explicit coordinates. In
this example, domain A contains retention times and domain B contains Kováts retention
indices.

```@example retention_mapping_workflow
using JuChrom
using CairoMakie

retention_times = [
    27.972, 30.043, 32.011, 33.884, 35.679, 37.398, 39.049, 40.636, 42.164,
    43.636, 45.058, 46.431, 47.759, 49.046, 50.296, 51.669, 53.266, 55.149,
]u"minute"

kovats_indices = [
    1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0,
    2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0,
]

mapper = fitmap(retention_times, kovats_indices)
```

The printed [`RetentionMapper`](@ref RetentionMapper) summary provides a compact overview
of the fitted anchor range and units. The default `λ=3e-9` is tuned for normalized
RT -> RI calibration and is intended to avoid boundary-of-monotonicity wiggles.

## Diagnostic plotting

Plot the forward map, its derivative, and the residuals to check whether the fitted curve
tracks the anchors without adding implausible oscillation.

```@example retention_mapping_workflow
fig = plot(mapper; size=(900, 600))
save("retention_mapping_workflow.svg", fig)
nothing # hide
```

![](retention_mapping_workflow.svg)

The black and red numbers at some anchor points show how much the fitted mapping departs
from those anchors. Visible oscillation in the derivative can indicate that a larger `λ`
is appropriate.

Fit a smoother mapper with stronger regularization and compare the diagnostics.

```@example retention_mapping_workflow
mapper_smooth = fitmap(retention_times, kovats_indices, λ=1e-6)
fig = plot(mapper_smooth; size=(900, 600))
save("retention_mapping_workflow_smooth.svg", fig)
nothing # hide
```

![](retention_mapping_workflow_smooth.svg)

The smoother fit deviates more from the anchor points but suppresses derivative wiggles.
Whether that is desirable depends on whether the anchor deviations reflect measurement noise
or real chromatographic structure.

## Applying the mapper

Use [`applymap`](@ref) to map retention times to retention indices.

```@example retention_mapping_workflow
ri = applymap(mapper, 41.5u"minute")
```

Broadcast over vectors with dot syntax.

```@example retention_mapping_workflow
rts = [10, 29.3, 35.0]u"minute"
ris = applymap.(mapper, rts)
```

Values outside the anchor domain are linearly extrapolated using the slope at the nearest
boundary. Set `warn=true` during analysis when extrapolation should be surfaced.

```@example retention_mapping_workflow
ri_extrapolated = applymap(mapper, 10u"minute"; warn=true)
```

Use [`invmap`](@ref) to map from retention index back to retention time.

```@example retention_mapping_workflow
ris = [1853.2, 3137.3, 3501.0]
rts = invmap.(mapper, ris)
```

For direct intensity correction, evaluate the local derivative. If `ri = f(t)`, then
`d(ri)/dt` is the local axis-stretch factor; area-preserving intensity transforms divide
the intensity by that factor.

```@example retention_mapping_workflow
scantimes = [1802.5, 1803.0, 1803.5]u"s"
intensities = [1000, 4000, 3500] / 0.5u"s"
dridt = derivmap.(mapper, scantimes, rA_unit=u"s")
ints_transformed = intensities ./ dridt
```

For `MassScanMatrix` and `VarianceMassScanMatrix`, [`applymap`](@ref) applies the Jacobian
correction and updates intensity and variance units automatically.

## Persisting mappers

When `JLD2` is loaded, JuChrom's JLD2 extension can store and restore
[`RetentionMapper`](@ref) objects.

```@example retention_mapping_workflow
using JLD2

save_object("retention_mapper.jld2", mapper)
mapper_loaded = load_object("retention_mapper.jld2")

jldsave("retention_mappers.jld2"; mapper, mapper_smooth)
mapper_reloaded = JLD2.load("retention_mappers.jld2", "mapper")
mapper_smooth_loaded = JLD2.load("retention_mappers.jld2", "mapper_smooth")

(mapper_loaded isa RetentionMapper, mapper_reloaded isa RetentionMapper,
    mapper_smooth_loaded isa RetentionMapper)
```

For API details, see [Fitting maps](mapping_fitting.md),
[Diagnostic plots](mapping_plotting.md), [Applying maps](mapping_application.md), and
[Mapper containers](mapping_containers.md).
