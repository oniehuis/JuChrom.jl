# JuChrom.jl

JuChrom is a Julia package for loading, storing, and analyzing chromatographic data,
with an emphasis on gas chromatography-mass spectrometry (GC-MS). The core API is pure
Julia, while some vendor-specific loaders (for example, the Shimadzu MS loader) rely on
optional extensions that use PyCall and the Python `olefile` module.

## Installation

The package is not registered, so install it directly from GitHub:

```
julia> ]
pkg> add https://github.com/oniehuis/JuChrom.jl
```

After the installation finishes you can load it in the Julia REPL with:

```
julia> using JuChrom
```

## Current capabilities

- Data structures for chromatographic and mass-spectral scans (retention, intensity,
  and m/z values), scan series, and matrices, including metadata. Retention times,
  intensities, and m/z values can carry `Unitful` units.
- Loaders for proprietary vendor formats, including Agilent MassHunter GC-MS,
  Agilent ChemStation MS, Agilent FID, and Shimadzu GC-MS data (the Shimadzu loader is
  available when the PyCall extension is loaded).
- Conversions between representations, including extracted-ion and total ion
  chromatograms from MS scan series or matrices.
- Transformations such as retention trimming, retention binning, and m/z binning within
  scans.
- Quadratic variance modeling of GC-MS intensities from technical replicates.
- Baseline estimation for chromatographic traces (e.g., individual ions or TICs).
- Retention mapping to fit continuous functions (and derivatives) that convert between
  retention times and retention indices.
- Spectral similarity tools (cosine distance/similarity) and a dynamic-programming
  alignment routine for matching spectra based on similarity scores.

## Quickstart

Most users will work with GC-MS data collected on an instrument. This example loads
Agilent GC-MS data stored in ChemStation MS format.

```@example 1
# Load JuChrom and the Agilent ChemStation MS loader
using JuChrom
using JuChrom.ChemStationMSLoader

# Load CairoMakie for plotting
using CairoMakie
CairoMakie.activate!()

# Load an example Agilent ChemStation GC-MS run
file = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
mss = load(ChemStationMS(file; mode=:ms))

# Compute the total ion chromatogram (TIC)
tic = mzchrom(mss)

# Plot the TIC and save to SVG
fig_1 = Figure(; size=(1200, 400))
axis_1 = Axis(fig_1[1,1], title="Total Ion Chromatogram",
                          ylabel="Intensity [no unit]",
                          xlabel="Retention [minute]")
lines!(axis_1, rawretentions(tic, unit=u"minute"), 
               rawintensities(tic),
               color=:blue)
save("tic.svg", fig_1)
nothing # hide
```

![](tic.svg)

```@example 1
# Bin m/z values to integer bins
bmss = binmzvalues(mss)

# Extract the m/z 109 chromatogram; silence warnings for scans without m/z 109 signal
xic = mzchrom(bmss, 109, warning=false)

# Plot the m/z 109 chromatogram and save to SVG
fig_2 = Figure(; size=(1200, 400))
axis_2= Axis(fig_2[1,1], title="m/z 109",
                         ylabel="Intensity [no unit]",
                         xlabel="Retention [minute]")
lines!(axis_2, rawretentions(xic, unit=u"minute"),
               rawintensities(xic),
               color=:blue)
save("xic.svg", fig_2)
nothing # hide
```

![](xic.svg)

```@example 1
# Estimate the baseline with airPLS, assuming Poisson count variance (Var ≈ intensity)
baseline = airpls(rawretentions(xic, unit=u"minute"),
                  rawintensities(xic), 
                  variances=rawintensities(xic),
                  λ=0.01)  # λ controls baseline smoothness (higher = smoother)

# Plot the m/z 109 chromatogram with its baseline and save to SVG.
fig_3 = Figure(; size=(1200, 400))
axis_3= Axis(fig_3[1,1], title="m/z 109",
                         ylabel="Intensity [no unit]",
                         xlabel="Retention [minute]")
lines!(axis_3, rawretentions(xic, unit=u"minute"),
               rawintensities(xic),
               color=:blue)
lines!(axis_3, rawretentions(xic, unit=u"minute"),
               baseline,
               color=:red)
save("xic-baseline.svg", fig_3)
nothing # hide
```

![](xic-baseline.svg)

## Disclaimer

JuChrom is provided "as is," without warranty of any kind. Users are responsible for
independently validating all outputs and interpretations and for determining suitability
for their specific applications. The authors and contributors disclaim any liability for
errors, omissions, or any consequences arising from use of the software, including use
in regulated, clinical, or safety-critical contexts.
