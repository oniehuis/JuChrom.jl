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

## Disclaimer

JuChrom is research software provided "as is." You remain responsible for validating
each step and interpretation of any data analysis; the authors cannot be held liable if
the package produces misleading or incorrect results.
