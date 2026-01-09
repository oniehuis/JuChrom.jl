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
# Load JuChrom package and its optional Agilent ChemStation MS loader
using JuChrom
using JuChrom.ChemStationMSLoader

# Load CairoMakie package for plotting
using CairoMakie
CairoMakie.activate!()

# Load GC-MS run data from an example Agilent ChemStation "data.ms" file
file = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
mss = load(ChemStationMS(file; mode=:ms))

# Infer tic
tic = mzchrom(mss)

# Plot TIC into SVG file
fig_1 = Figure(; size=(1000, 400))
axis_1 = Axis(fig_1[1,1], title="Total Ion Chromatogram",
                          ylabel="Intensity [no unit]",
                          xlabel="Retention [minute]")
lines!(axis_1, rawretentions(tic, unit=u"minute"), rawintensities(tic), color=:red)
save("tic.svg", fig_1)
nothing # hide
```

![](tic.svg)

```@example 1
# Integer-bin mz values
bmss = binmzvalues(mss)

# Extract chromatogram of m/z 85, repressing warnings for scans without m/z 85 signal
t85 = mzchrom(bmss, 85, warning=false)

# Plot m/z 85 chromatogram into SVG file
fig_2 = Figure(; size=(1000, 400))
axis_2= Axis(fig_2[1,1], title="m/z 85",
                         ylabel="Intensity [no unit]",
                         xlabel="Retention [minute]")
lines!(axis_2, rawretentions(t85, unit=u"minute"), rawintensities(t85), color=:red)
save("t85.svg", fig_2)
nothing # hide
```

![](t85.svg)


## Disclaimer

JuChrom is provided "as is," without warranty of any kind. Users are responsible for
independently validating all outputs and interpretations and for determining suitability
for their specific applications. The authors and contributors disclaim any liability for
errors, omissions, or any consequences arising from use of the software, including use
in regulated, clinical, or safety-critical contexts.
