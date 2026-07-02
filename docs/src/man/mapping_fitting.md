# Fitting maps

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

For a complete example that continues into plotting, applying, and persisting the mapper,
see the [retention mapping workflow](mapping_workflow.md).

```@docs
JuChrom.fitmap
```
