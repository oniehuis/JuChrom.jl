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
(Kováts 1958). Another approach uses a curated set of internal standards (e.g., a few stable 
compounds spiked into every run), producing matched retention pairs that anchor the mapping 
across batches (Skoog et al. 2007).

## JuChrom retention mapping

JuChrom provides [`fitmap`](@ref JuChrom.fitmap) to infer an empirical, smooth mapping 
function from paired retention values in domain A (e.g., retention times) and reference 
values in domain B (e.g., Kováts indices). The fit constructs a cubic B-spline and chooses 
the smallest smoothing penalty that still enforces a strictly increasing curve. Concretely, 
the objective minimizes squared residuals plus a curvature penalty based on the spline’s 
second derivative: large changes in slope are penalized, which discourages wiggles and 
yields a smoother, more stable mapping between anchor points. In parallel, nonnegative 
first-derivative values are enforced at a dense grid of points to guarantee monotonicity. 
This yields a continuous, differentiable, and invertible function (stored in a 
[`RetentionMapper`](@ref RetentionMapper)) that supports both forward mapping and 
reliable reverse mapping via the monotonic inverse.

For most users, the primary tuning parameter when applying [`fitmap`](@ref JuChrom.fitmap) 
is smoothing strength (λ). Larger values emphasize smoothness over exact fit to the anchor 
points, while smaller values track the anchors more tightly. The default automatically 
searches for the smallest λ that still yields a strictly monotonic fit.

JuChrom includes visual diagnostics for a fitted mapper, letting you inspect the forward 
and inverse fits side by side. The plotting helpers load automatically once a 
[Makie](https://docs.makie.org) backend is available. 

## Example

```@example 1
# Load JuChrom
using JuChrom

# Load CairoMakie for plotting
using CairoMakie
CairoMakie.activate!()

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
For a closer look, plot the forward and inverse maps (including derivatives) and the 
residuals.

```@example 1
# Visual diagnostics (forward/inverse maps and residuals)
fig = plot(mapper; reverse=true, size=(900, 600))
save("retention_mapper.svg", fig)
nothing # hide
```

![](retention_mapper.svg)

While the inferred mapping is fully satisfactory, you may want to tune the smoothing 
strength if the anchor points themselves are noisy or not fully trustworthy. Raising `λ_min` 
prevents overfitting (less wiggle, more stability), while raising `λ_max` can help recover 
monotonicity when a fit would otherwise fail. Here we increase `λ_min` from its default 
(`1e-20`) to `1e-7` and examine the effect on the derivative plots.

```@example 1
# Fit mapping function with λ_min set to 1e-7
mapper_λ_min_set = fitmap(retention_times, kovats_indices, λ_min=1e-7)
fig = plot(mapper_λ_min_set; reverse=true, size=(900, 600))
save("retention_mapper_λ_min_set.svg", fig)
nothing # hide
```

![](retention_mapper_λ_min_set.svg)

In the derivative plots the mapping is noticeably smoother, but it no longer passes exactly 
through the anchor points, as shown by the residuals. Whether the suppressed wiggles reflect 
real structure that should be modeled or are noise that should be smoothed away is a 
judgment call for the analyst.

Let's continue with the mapper inferred using the default `λ_min` and use it to compute
retention indices for a few retention times, including extrapolation beyond the domain.

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

To invert the mapping (e.g., from Kováts to retention time), use 
[`invmap`](@ref).

```@example 1
ris = [1853.2, 3137.3, 3501.0]
rts = invmap.(mapper, ris)  # dot form broadcasts over the vector
```

To transform intensities with the Jacobian, use the derivative of the mapping. If
`ri = f(t)`, then `d(ri)/dt` is the local stretch factor; to preserve area, you divide the
intensity by this slope at the same time point.

```@example 1
scantimes = [1802.5, 1803.0, 1803.5]u"s"
intensities = [1000, 4000, 3500]

dridt = derivmap.(mapper, scantimes)
```

Note that you can supply mapping inputs in any compatible time unit; values are converted
automatically.

The derivative tells you how much the Kováts retention index changes per unit time. To
transform the intensities associated with each scan time, divide by the Jacobian (i.e.,
`d(ri)/dt`).

```@example 1
ints_transformed = intensities ./ dridt
```

The transformed intensities are now expressed per unit of retention index rather than per
unit time, reflecting the change of variables.

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
jldsave("retention_mappers.jld2"; mapper, mapper_λ_min_set)
mapper_reloaded = JLD2.load("retention_mappers.jld2", "mapper")
```
```@example 1
mapper_λ_min_set_loaded = JLD2.load("retention_mappers.jld2", "mapper_λ_min_set")
```

## References

- Kováts E (1958): Gas-Chromatographische Charakterisierung organischer Verbindungen. Teil 1: Retentionsindices aliphatischer Halogenide,Alkohole, Aldehyde und Ketone. Helvetica Chimica Acta 41: 1915-1932.
- Skoog DA, Holler FJ, Crouch SR (2007): Principles of Instrumental Analysis. 6th ed. Thomson Brooks/Cole.
