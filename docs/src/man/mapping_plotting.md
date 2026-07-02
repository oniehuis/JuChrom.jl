# Diagnostic plots

The plotting methods are available after loading a Makie backend such as `CairoMakie`.
For a complete example with generated diagnostics, see the
[retention mapping workflow](mapping_workflow.md).

```julia
plot(
    mapper::RetentionMapper;
    size=(1600, 800),
    linecolor=:blue,
    linewidth=1.0,
    markersize=10.0,
    markercolor=:orange,
    rA_unit=retentionunit_A(mapper),
    rB_unit=retentionunit_B(mapper),
    reverse=false,
    digits=1,
)
```

`plot(mapper)` returns a `Makie.Figure` with diagnostic views of a fitted
`RetentionMapper`. The forward mapping from retention domain A to domain B and the
forward derivative `dB/dA` are always shown. Calibration anchors are drawn as markers,
and nonzero anchor residuals are annotated next to the corresponding anchors. Set
`reverse=true` to add the inverse mapping from domain B to domain A and the inverse
derivative `dA/dB`.

Keyword arguments control figure layout, styling, unit conversion, and residual labels:

- `size`: figure size passed to `Makie.Figure`.
- `linecolor`, `linewidth`: style for fitted mapping and derivative curves.
- `markersize`, `markercolor`: style for calibration anchor markers.
- `rA_unit`: unit used for domain-A axes and raw domain-A values. Defaults to
  `retentionunit_A(mapper)`.
- `rB_unit`: unit used for domain-B axes and raw domain-B values. Defaults to
  `retentionunit_B(mapper)`.
- `reverse`: include inverse mapping and inverse derivative panels when `true`.
- `digits`: number of decimal digits used for residual annotations.

`rA_unit` and `rB_unit` must be `nothing` for unitless mapper domains or compatible
`Unitful.Units` values for unitful domains.

```julia
using JuChrom
using CairoMakie

mapper = fitmap([1.2, 2.5, 4.1, 6.8]u"minute", [100.0, 200.0, 300.0, 400.0])
fig = plot(mapper; rA_unit=u"minute", reverse=true, size=(1200, 600))
```
