# JuChrom.jl

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][GHA-img]][GHA-url] [![][codecov-img]][codecov-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://oniehuis.github.io/JuChrom.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://oniehuis.github.io/JuChrom.jl/dev

[GHA-img]: https://github.com/oniehuis/JuChrom.jl/actions/workflows/CI.yml/badge.svg?branch=main
[GHA-url]: https://github.com/oniehuis/JuChrom.jl/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/oniehuis/JuChrom.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/oniehuis/JuChrom.jl

A Julia package for loading, processing, and analyzing chromatographic data, with
an emphasis on gas chromatography–mass spectrometry (GC-MS). Work with individual
scans, chromatograms, mass spectra, or retention-by-m/z matrices while preserving
metadata and physical units.

## Features

- Load Agilent ChemStation MS, MassHunter MS, and FID data, plus Shimadzu MS data
  through an optional Python extension.
- Extract total-ion and extracted-ion chromatograms, combine scans into mass
  spectra, and compare spectra using cosine similarity.
- Prepare data with trimming, m/z and retention binning, normalization, baseline
  estimation, and baseline subtraction.
- Detect n-alkane ladders, calibrate retention indices, and estimate and propagate
  intensity variances.
- Fit chromatographic peak shapes, perform PARAFAC2 deconvolution across samples,
  and visualize results with Makie.

## Installation

Julia 1.12 or later is required. Install the development version from GitHub in
your Julia environment:

```julia
using Pkg
Pkg.add(url="https://github.com/oniehuis/JuChrom.jl", rev="main")
```

Then load the package with `using JuChrom`.

## Quickstart

Run this in the Julia REPL. The example uses an n-alkane standard distributed with
JuChrom, so no measurement file of your own is needed.

```julia
using JuChrom
using JuChrom.ChemStationMSLoader

# Load the example run as a series of mass-spectral scans.
path = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
series = load(ChemStationMS(path; mode=:ms))

# Bin to integer m/z values and build a retention-by-m/z intensity matrix.
binned = binmzvalues(series, JuChrom.integer)
matrix = mscanmatrix(binned)

# Extract the total-ion chromatogram (TIC) and an m/z 57 chromatogram (XIC).
tic = mzchrom(matrix)
xic = mzchrom(matrix, 57)

# Extract a mass spectrum at the strongest m/z 57 signal.
peak_scan = argmax(rawintensities(xic))
spectrum = extractmassspectrum(matrix, [peak_scan])
```

`series` stores individual scans; `matrix` stores scans as rows and m/z channels
as columns. Use `retentions`, `mzvalues`, and `intensities` to access data with
their units, or the corresponding `rawretentions`, `rawmzvalues`, and
`rawintensities` getters for numeric values. JuChrom re-exports Unitful, so unit
expressions such as `u"minute"` are available after `using JuChrom`.

### Plot chromatograms and spectra

Install a Makie backend in the same environment to use JuChrom's plotting methods:

```julia
using Pkg
Pkg.add("CairoMakie")
```

Continue with the data from the quickstart:

```julia
using CairoMakie

tic_figure = tictrace(tic; retentionunit=u"minute")
save("tic.png", tic_figure)

xic_figure = tictrace(xic; retentionunit=u"minute")
save("xic.png", xic_figure)

spectrum_figure = massspectrum(spectrum)
save("spectrum.png", spectrum_figure)
```

The figures are saved in the current working directory. For interactive
inspection and comparison of runs, see
[JuChromExplorer](https://github.com/oniehuis/JuChromExplorer.jl).

### Load your own data

Replace `path` in the quickstart with the path to your ChemStation `data.ms` file
or its containing `.D` directory. Other formats use their own loader modules and
constructors; see the [loader guide](docs/src/man/loaderapi.md).
The [Shimadzu loader](docs/src/man/shimadzums.md) additionally requires PyCall and
the Python `olefile` module.

## Documentation

See the [documentation](https://oniehuis.github.io/JuChrom.jl/dev/) for the full
API and worked examples, including
[n-alkane detection](docs/src/man/alkaneseries.md),
[retention mapping](docs/src/man/mapping_workflow.md),
[variance estimation](docs/src/man/variance_model.md), and
[deconvolution](docs/src/man/deconvolution.md).
