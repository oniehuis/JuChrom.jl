# Plotting

JuChrom provides Makie plotting methods through its optional Makie extension. Load a Makie
backend, for example `CairoMakie` or `GLMakie`, before plotting.

## TIC traces

`tictrace` and `tictrace!` plot a total ion chromatogram (TIC) for JuChrom scan
containers:

- `MassScanMatrix` and `VarianceMassScanMatrix`: y values are the row-wise sums of the
  intensity matrix.
- `MassScanSeries`: y values are the per-scan sums of mass-spectral intensities.
- `ChromScanSeries`: y values are the chromatographic scan intensities.

```julia
using CairoMakie
using JuChrom
using Unitful

fig = tictrace(
    msm;
    figure=(; size=(900, 350)),
    retentionunit=u"minute",
    intensityunit=u"nA",
    color=:black,
    linewidth=1.5,
)
```

The mutating form plots into an existing axis and returns Makie's `Lines` plot:

```julia
fig = Figure(; size=(900, 350))
ax = Axis(fig[1, 1])
plt = tictrace!(
    ax,
    msm;
    retentionunit=u"minute",
    intensityunit=u"nA",
    color=:black,
)
```

By default, the axis labels are set to `"Retention [$unit]"` and
`"Intensity [$unit]"` using the displayed retention and intensity units, and the title is
empty. Use `title` or `axis=(; title=...)` to set one:

```julia
fig = tictrace(
    msm;
    retentionunit=u"minute",
    intensityunit=u"nA",
    title="Total ion chromatogram",
)
```

You can also use Makie's layout-position form, which creates an axis in the selected grid
cell and returns `AxisPlot`:

```julia
fig = Figure()
axplot = tictrace(
    fig[1, 1],
    msm;
    retentionunit=u"minute",
    intensityunit=u"nA",
)
```

`retentionunit` controls the x-axis unit and `intensityunit` controls the unit used before
summing and plotting intensities. In both cases, JuChrom converts to the requested unit and
then strips the unit before passing numeric vectors to Makie. The default value `nothing`
means "strip the stored unit if there is one"; unitless data are plotted unchanged. If data
are unitless and a target unit is requested, JuChrom throws an `ArgumentError`.

All other keywords are forwarded to Makie's `lines!` plot, so standard line attributes
such as `color`, `linewidth`, `linestyle`, `linecap`, and `alpha` are available.

For manual data extraction, the same conversion rule is used by the raw getters:

```julia
x = rawretentions(msm; unit=u"minute")
y = vec(sum(rawintensities(msm; unit=u"nA"); dims=2))
lines!(ax, x, y)
```

Makie's `lines` and `lines!` also accept the same JuChrom scan containers directly. These
lower-level methods use the same TIC extraction and unit-conversion rules as `tictrace`,
but they do not set axis labels.

## Annotated alkane ladders

Use `tictrace(msm, result)` or `tictrace!(ax, msm, result)` for the higher-level TIC view
with alkane ladder annotations. These methods use the same `figure`, `axis`,
`retentionunit`, and `intensityunit` keywords:

```julia
fig = tictrace(
    msm,
    result;
    figure=(; size=(900, 350)),
    retentionunit=u"minute",
    intensityunit=u"nA",
)
```

The TIC line is converted to `intensityunit` before plotting, and ladder-step retentions
are converted to `retentionunit`. The fitted baseline is hidden by default; pass
`baseline=true` to show it.
