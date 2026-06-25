# Plotting

JuChrom provides Makie plotting methods through its optional Makie extension. Load a Makie
backend, for example `CairoMakie` or `GLMakie`, before plotting.

## TIC lines

Makie's `lines` and `lines!` accept JuChrom scan containers directly and plot a total ion
chromatogram (TIC):

- `MassScanMatrix` and `VarianceMassScanMatrix`: y values are the row-wise sums of the
  intensity matrix.
- `MassScanSeries`: y values are the per-scan sums of mass-spectral intensities.
- `ChromScanSeries`: y values are the chromatographic scan intensities.

```julia
using CairoMakie
using JuChrom
using Unitful

fig = Figure(; size=(900, 350))
ax = Axis(fig[1, 1], xlabel="Retention [minute]", ylabel="TIC [nA]")

plt = lines!(
    ax,
    msm;
    retentionunit=u"minute",
    intensityunit=u"nA",
    color=:black,
    linewidth=1.5,
)
```

The non-mutating form creates a figure and axis and returns Makie's `FigureAxisPlot`:

```julia
fap = lines(
    msm;
    figure=(; size=(900, 350)),
    axis=(; xlabel="Retention [minute]", ylabel="TIC [nA]"),
    retentionunit=u"minute",
    intensityunit=u"nA",
    color=:black,
)
```

You can also use Makie's layout-position form, which creates an axis in the selected grid
cell and returns `AxisPlot`:

```julia
fig = Figure()
axplot = lines(
    fig[1, 1],
    msm;
    axis=(; xlabel="Retention [minute]", ylabel="TIC [nA]"),
    retentionunit=u"minute",
    intensityunit=u"nA",
)
```

`retentionunit` controls the x-axis unit and `intensityunit` controls the unit used before
summing and plotting intensities. In both cases, JuChrom converts to the requested unit and
then strips the unit before passing numeric vectors to Makie. The default value `nothing`
means "strip the stored unit if there is one"; unitless data are plotted unchanged. If data
are unitless and a target unit is requested, JuChrom throws an `ArgumentError`.

All other keywords are forwarded to Makie's `lines` plot, so standard line attributes such
as `color`, `linewidth`, `linestyle`, `linecap`, and `alpha` are available.

For manual data extraction, the same conversion rule is used by the raw getters:

```julia
x = rawretentions(msm; unit=u"minute")
y = vec(sum(rawintensities(msm; unit=u"nA"); dims=2))
lines!(ax, x, y)
```

## Annotated alkane ladders

Use `tictrace(msm, result)` or `tictrace!(ax, msm, result)` for the higher-level TIC view
with baseline and alkane ladder annotations. The `lines` methods above are the lower-level
Makie primitive for plotting only the TIC line.
