# Overview

## Mapping rationale

Retention mapping converts an irregular axis (time or distance) into a stable, comparable 
index so chromatograms from different runs can be aligned, interpolated, and compared on a 
common grid. This makes retention behavior reproducible across batches, instruments, and 
methods while preserving ordering and monotonicity.

A continuous, differentiable mapping does more than align coordinates: it enables 
intensity-aware transforms via the 
[`Jacobian`](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant). When you warp 
the axis, the mapping derivative provides a principled way to rescale intensities so areas 
and peak shapes remain physically meaningful. This yields consistent peak integrals across 
transformed domains, supports smooth interpolation, and avoids artifacts from piecewise or 
discontinuous mappings.

Mapping points are typically collected from standards with known reference positions. In 
gas chromatography, a common example is an n-alkane ladder, which yields paired arrays of 
retention times and 
[Kováts retention indices](https://en.wikipedia.org/wiki/Kovats_retention_index) 
(Kováts 1958). JuChrom can infer these points from GC-MS alkane ladder data with
[`findalkanes`](@ref) and expose the accepted mapper anchors with
[`alkaneladdercalibrationpoints`](@ref). Another approach uses a curated set of internal
standards (e.g., a few stable compounds spiked into every run), producing matched
retention pairs that anchor the mapping across batches (Skoog et al. 2007).

## JuChrom retention mapping

JuChrom provides [`fitmap`](@ref JuChrom.fitmap) to infer an empirical, smooth mapping
function from paired retention values in domain A (e.g., retention times) and reference
values in domain B (e.g., Kováts indices). The fit constructs a cubic B-spline with a
fixed smoothing penalty `λ` and validates that the resulting curve is strictly increasing.
Concretely, the objective minimizes squared residuals plus a curvature penalty based on
the spline’s second derivative: large changes in slope are penalized, which discourages
wiggles and yields a smoother, more stable mapping between anchor points. In parallel,
nonnegative first-derivative values are enforced at the spline constraint points. The final
spline is then checked for strict monotonicity on a dense validation grid so that the stored
[`RetentionMapper`](@ref RetentionMapper) supports both forward mapping and reliable
reverse mapping via the monotonic inverse.

For most users, the primary tuning parameter when applying [`fitmap`](@ref JuChrom.fitmap) 
is smoothing strength (λ). Larger values emphasize smoothness over exact fit to the anchor 
points, while smaller values track the anchors more tightly. The default `λ=3e-9` is tuned
for normalized RT -> RI calibration and is intended to avoid boundary-of-monotonicity
wiggles.

JuChrom includes visual diagnostics for a fitted mapper, letting you inspect the forward 
and inverse fits side by side. The plotting helpers load automatically once a 
[Makie](https://docs.makie.org) backend is available. 

## Example

```@example 1
# Load JuChrom
using JuChrom

# Load CairoMakie for plotting
using CairoMakie

# Known mapping points: retention times (minutes) and Kováts indices.
# JuChrom reexports Unitful, so unit literals work without `using Unitful`.
retention_times = [
    27.972, 30.043, 32.011, 33.884, 35.679, 37.398, 39.049, 40.636, 42.164,
    43.636, 45.058, 46.431, 47.759, 49.046, 50.296, 51.669, 53.266, 55.149,
]u"minute"

kovats_indices = [
    1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0,
    2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0,
]

# Fit a smooth, strictly increasing mapping function
mapper = fitmap(retention_times, kovats_indices)
```

The printed [`RetentionMapper`](@ref RetentionMapper) summary provides a compact overview. 
For a closer look, plot the forward map, its derivative, and the residuals. Pass
`reverse=true` to include the inverse map and inverse derivative.

```@example 1
# Visual diagnostics (forward map, derivative, and residuals)
fig = plot(mapper; size=(900, 600))
save("retention_mapper.svg", fig)
nothing # hide
```

![](retention_mapper.svg)

The black and red numbers at some anchor points show how much the fitted mapping departs 
from those anchors. The inferred mapping is fully satisfactory, but you may want to tune 
the smoothing strength if the anchor points are noisy or not fully trustworthy. In the 
derivative plots, for example, visible oscillation indicates that a larger `λ` may be
appropriate. Increasing `λ` reduces overfitting (less wiggle, more stability), while
decreasing `λ` tracks the anchors more tightly.

Let us fit a smoother mapper with `λ=1e-6` and examine the effect on the derivative plots.

```@example 1
# Fit mapping function with stronger smoothing
mapper_smooth = fitmap(retention_times, kovats_indices, λ=1e-6)
fig = plot(mapper_smooth; size=(900, 600))
save("retention_mapper_smooth.svg", fig)
nothing # hide
```

![](retention_mapper_smooth.svg)

In the derivative plots the mapping is noticeably smoother, but it now deviates more from 
the anchor points, as shown by the residuals. Whether the suppressed wiggles reflect real 
structure that should be modeled or are noise that should be smoothed away is a judgment 
call for the analyst. If a fit fails monotonicity validation, use a larger `λ` or inspect
the anchor points for false or noisy calibration steps.

Let’s continue with the original mapper (inferred from using the default `λ=3e-9`) and use
it to compute retention indices for a few retention times, including extrapolation beyond 
the domain.

```@example 1
ri = applymap(mapper, 41.5u"minute")  # single value
```

```@example 1
rts = [10, 29.3, 35.0]u"minute"
ri = applymap.(mapper, rts)  # dot form broadcasts over the vector
```

!!! warning "Domain limits and extrapolation"
    Mappings are defined over the anchor domain. Values outside that domain are linearly 
    extrapolated using the slope at the nearest boundary. Use `warn=true` on 
    [`applymap`](@ref)
    /
    [`invmap`](@ref)
    /
    [`derivmap`](@ref)
    /
    [`derivinvmap`](@ref)
    or the `raw*` variants to surface extrapolation during
    analysis.

To invert the mapping (e.g., from Kováts to retention time), use [`invmap`](@ref).

```@example 1
ris = [1853.2, 3137.3, 3501.0]
rts = invmap.(mapper, ris)  # dot form broadcasts over the vector
```

To transform intensities with the Jacobian, use the derivative of the mapping. If
`ri = f(t)`, then `d(ri)/dt` is the local stretch factor; to preserve area, you divide the
intensity by this slope at the same time point. However, even when intensities are
reported as unitless, they usually represent counts accumulated over a finite scan
interval, so the implied unit is typically `time^-1`. Here we assume the given intensities
are counts acquired over a 0.5‑second scan interval, and therefore divide by 0.5 seconds.
Because the intensities and `d(ri)/dt` must use the same time unit for the division to
cancel cleanly, we convert the Jacobian to the unit `s^-1`. This is done by specifying the
input‑domain unit as `u"s"`.

```@example 1
scantimes = [1802.5, 1803.0, 1803.5]u"s"
intensities = [1000, 4000, 3500] / 0.5u"s"
```

```@example 1
dridt = derivmap.(mapper, scantimes, rA_unit=u"s")
```

The derivative tells you how much the Kováts retention index changes per unit time. To
transform the intensities associated with each scan time, divide by the Jacobian (i.e.,
`d(ri)/dt`).

```@example 1
ints_transformed = intensities ./ dridt
```

The transformed intensities are now expressed per unit of retention index rather than per
unit time, reflecting the change of variables.

!!! note "Use ion-specific retention times for sequential MS data"
    The Jacobian must be evaluated at the retention coordinate of the measurement whose
    intensity is transformed. For a chromatographic trace with one intensity per scan, the
    scan retention time is the appropriate coordinate. For sequentially scanned MS data,
    however, different m/z channels in the same scan are measured at slightly different
    times. In that case, a single scan-level Jacobian applies the same correction to all
    ions in a scan and is only an approximation.

    For the most accurate RT -> RI transformation of ion traces or MS scan matrices, first
    correct the scan retention time to the acquisition time of the specific m/z channel,
    then evaluate the Jacobian at that ion-specific retention time:

    ```julia
    rt_ion = mzretention(scan_retention, mz; order=:descending, ...)
    ri_ion = applymap(mapper, rt_ion)
    jacobian = derivmap(mapper, rt_ion)
    intensity_per_ri = intensity_per_rt / jacobian
    ```

    This distinction is mainly relevant for sequential quadrupole-like acquisition. For
    simultaneous m/z acquisition, or when the m/z scan duration is negligible relative to
    chromatographic peak widths and mapper curvature, the scan-level Jacobian can be a
    reasonable approximation.

## Mapping tools at a glance

| Function | Use case |
| :--- | :--- |
| [`fitmap`](@ref) | Infer mapping function from paired points |
| [`applymap`](@ref) | Map retention domain A → retention domain B |
| [`invmap`](@ref) | Map retention domain B → retention domain A |
| [`derivmap`](@ref) | Jacobian `d(ri)/dt` for intensity scaling |
| [`derivinvmap`](@ref) | Inverse mapping derivative |
| [`rawapplymap`](@ref) | Unitless variant of [`applymap`](@ref) |
| [`rawinvmap`](@ref) | Unitless variant of [`invmap`](@ref) |
| [`rawderivmap`](@ref) | Unitless variant of [`derivmap`](@ref) |
| [`rawderivinvmap`](@ref) | Unitless variant of [`derivinvmap`](@ref) |
| [`retentions_A`](@ref), [`retentions_B`](@ref) | Anchor vectors used to fit the mapper |
| [`rawretentions_A`](@ref), [`rawretentions_B`](@ref) | Unitless anchor vectors |
| [`retentionunit_A`](@ref), [`retentionunit_B`](@ref) | Stored units for the anchor domains |
| [`extras`](@ref) | Metadata attached to the mapper |
| [`mapmin`](@ref), [`mapmax`](@ref) | Numeric minimum and maximum value of the input domain |
| [`invmapmin`](@ref), [`invmapmax`](@ref) | Numeric minimum and maximum value of the output domain |
| [`rawmapmin`](@ref), [`rawmapmax`](@ref) | Unitless variant of [`mapmin`](@ref) and [`mapmax`](@ref) |
| [`rawinvmapmin`](@ref), [`rawinvmapmax`](@ref) | Unitless variant of [`invmapmin`](@ref), [`invmapmax`](@ref) |

For full API details, see [Mapping tools](mapping_tools.md).

## JLD2 support

JuChrom also ships a `JLD2` extension so [`RetentionMapper`](@ref) objects can be stored 
and restored with [`JLD2`](https://github.com/JuliaIO/JLD2.jl). The extension loads 
automatically once `JLD2` is available.

```@example 1
using JLD2

# Save and load single mapper
save_object("retention_mapper.jld2", mapper)
mapper_loaded = load_object("retention_mapper.jld2")

# Save multiple mappers under their own names and load them back
jldsave("retention_mappers.jld2"; mapper, mapper_smooth)
mapper_reloaded = JLD2.load("retention_mappers.jld2", "mapper")
```
```@example 1
mapper_smooth_loaded = JLD2.load("retention_mappers.jld2", "mapper_smooth")
```

## References

- Kováts E (1958): Gas-Chromatographische Charakterisierung organischer Verbindungen. Teil 1: Retentionsindices aliphatischer Halogenide,Alkohole, Aldehyde und Ketone. Helvetica Chimica Acta 41: 1915-1932.
- Skoog DA, Holler FJ, Crouch SR (2007): Principles of Instrumental Analysis. 6th ed. Thomson Brooks/Cole.
