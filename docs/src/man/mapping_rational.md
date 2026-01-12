# Mapping rational

Retention mapping converts an irregular axis (time or distance) into a stable, comparable 
index so chromatograms from different runs can be aligned, interpolated, and compared on a 
common grid. This makes retention behavior reproducible across batches, instruments, and 
methods while preserving ordering and monotonicity.

A continuous, differentiable mapping does more than align coordinates: it enables intensity-
aware transforms via the Jacobian. When you warp the axis, the mapping derivative provides a 
principled way to rescale intensities so areas and peak shapes remain physically meaningful. 
This yields consistent peak integrals across transformed domains, supports smooth 
interpolation, and avoids artifacts from piecewise or discontinuous mappings.

Mapping points are typically collected from standards with known reference positions. A 
common example is an n-alkane ladder, which yields paired arrays of retention times and 
Kováts retention indices (Kováts 1958). Another example is a curated set of internal 
standards (e.g., a few stable compounds spiked into every run), producing matched retention 
pairs that anchor the mapping across batches (Skoog et al. 2007).

JuChrom provides `fitmap` to infer an empirical, smooth mapping function from these paired 
points. The fit constructs a cubic B-spline and chooses the smallest smoothing penalty that 
still enforces a strictly increasing curve. Concretely, the objective minimizes squared 
residuals plus a curvature penalty based on the spline’s second derivative: large changes 
in slope are penalized, which discourages wiggles and yields a smoother, more stable mapping 
between anchor points. In parallel, nonnegative first-derivative values are enforced at a 
dense grid of points to guarantee monotonicity. This yields a continuous, differentiable, 
and invertible function (stored in a `RetentionMapper`) that supports both forward mapping 
and reliable reverse mapping via the monotonic inverse.

For most users, the primary tuning parameter when applying `fitmap` is smoothing strength 
(λ). Larger values emphasize smoothness over exact fit to the anchor points, while smaller 
values track the anchors more tightly. The default automatically searches for the smallest λ 
that still yields a strictly monotonic fit.

JuChrom also includes Makie-based diagnostics for a fitted mapper, so you can inspect the 
forward and inverse fits visually. JuChrom reexports Unitful, so unit literals like 
`u"minute"` work without an explicit `using Unitful`. Example:

```@example 1
# Load JuChrom
using JuChrom

# Load CairoMakie for plotting
using CairoMakie
CairoMakie.activate!()

# Known mapping points: retention times (minutes) and Kováts indices
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

The printed `RetentionMapper` summary provides a compact overview. For a closer look,
plot the forward and inverse maps (including derivatives) and the residuals.

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

Note that the first retention time lies below the spline's lower domain boundary
(27.972 minute). In that case the function extrapolates linearly using the slope at the
nearest boundary. Set `warn=true` to emit a warning when extrapolation occurs (default
`false`).

```@example 1
ri = applymap.(mapper, rts, warn=true)  # dot form broadcasts over the vector
```

```@example 1
ris = [1853.2, 3137.3, 3501.0]
rts = invmap.(mapper, ris, warn=true)  # dot form broadcasts over the vector
```

To transform intensities with the Jacobian, use the derivative of the mapping. If
`ri = f(t)`, then `d(ri)/dt` is the local stretch factor; to preserve area you divide the
intensity by this slope at the same time point.

```@example 1
scantimes = [1802.5, 1803.0, 1803.5]u"s"
intensities = [1000, 4000, 3500]

ris = applymap.(mapper, scantimes, warn=true)
jacobian = derivmap.(mapper, scantimes)
ints_transformed = intensities ./ jacobian
```

## References

- Kováts E (1958): Gas-Chromatographische Charakterisierung organischer Verbindungen. Teil 1: Retentionsindices aliphatischer Halogenide,Alkohole, Aldehyde und Ketone. Helvetica Chimica Acta 41: 1915-1932.
- Skoog DA, Holler FJ, Crouch SR (2007): Principles of Instrumental Analysis. 6th ed. Thomson Brooks/Cole.
