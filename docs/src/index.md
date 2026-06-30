# JuChrom.jl

JuChrom is a Julia package for loading, storing, and analyzing chromatographic data,
with an emphasis on gas chromatography-mass spectrometry (GC-MS). The core API is pure
Julia; some vendor-specific loaders are provided via optional extensions.

## Installation

JuChrom is not yet registered in the General registry. Install it directly from GitHub:

```
julia> ]
pkg> add https://github.com/oniehuis/JuChrom.jl
```

After the installation finishes, you can load it in the Julia REPL with:

```
julia> using JuChrom
```

## Current capabilities

- Data structures for chromatographic scans, mass-spectral scans, scan series,
  mass-scan matrices, and variance-carrying matrices, including metadata. Retention
  times, intensities, and m/z values can carry
  [Unitful.jl](https://github.com/JuliaPhysics/Unitful.jl) units.
- Loaders for proprietary vendor formats, including Agilent MassHunter GC-MS,
  Agilent ChemStation MS, Agilent FID, and Shimadzu GC-MS data.
- Conversions between representations, including extracted-ion and total-ion
  chromatograms from MS scan series or matrices and conversion to mass-scan matrices.
- Preprocessing transformations such as retention and index trimming, m/z binning,
  retention binning and gridding, scan-timing correction, dwell normalization, TIC
  normalization, censored-value replacement, centered log-ratio transformation, and
  whitening.
- Automatic identification of n-alkane series in standard ladder runs for
  retention-index calibration, including ladder-step detection, calibration-point
  selection, and extraction of ladder mass spectra.
- Variance estimation and modeling, including count-based variance estimates and linear
  GC-MS intensity variance models fitted from n-alkane standards.
- Retention mapping for fitting continuous forward and inverse maps between retention
  coordinates and retention indices, including derivatives and Jacobian-aware intensity
  transformation.
- Baseline estimation for chromatographic traces, such as individual ion traces or TICs.
- Spectral similarity tools for cosine similarity and distance, plus a dynamic-programming
  alignment routine for matching spectra based on similarity scores.
- Deconvolution tools for chromatographic peak-shape fitting and PARAFAC2 decomposition
  across samples.
- Optional Makie plotting methods for TIC traces, mass spectra, retention mappers, and
  annotated n-alkane ladders.

## Quickstart

```@example 1
# Load JuChrom and the Agilent ChemStation MS loader
using JuChrom
using JuChrom.ChemStationMSLoader

# Load a Makie backend for plotting
using CairoMakie

# Load an example n-alkane standard run
file = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
mss = load(ChemStationMS(file; mode=:ms))

# Plot the raw total ion chromatogram (TIC); notice unitfree intensity
fig₁ = tictrace(mss; retentionunit=u"minute", figure=(; size=(1000, 350)))
save("rt_tic.svg", fig₁)
nothing # hide
```

![](rt_tic.svg)

```@example 1
# Trim to the n-alkane ladder retention window
tmss = retentiontrim(mss; start=10u"minute", stop=30u"minute")

# Bin m/z values to a common nominal-mass axis
btmss = binmzvalues(tmss)

# Convert scans to a dense retention-by-m/z matrix
msm = mscanmatrix(btmss)

# Normalize intensity by dwell time
nmsm = dwellnormalize(msm)

# Extract the trace of some ion of interest (e.g., m/z 104)
xic = mzchrom(nmsm, 104, warning=false)

# Plot m/z 104 trace; notice unitful intensity
fig₂ = tictrace(xic; retentionunit=u"minute", figure=(; size=(1000, 350)))
save("xic.svg", fig₂)
nothing # hide
```

![](xic.svg)

```@example 1
# Find the n-alkane ladder for retention-index calibration
asr = findalkanes(nmsm)
asr.success || throw(ArgumentError(
    "Failed to find a suitable ladder in standard file: $file"))

# Overlay detected ladder steps on the TIC
fig₃ = tictrace(nmsm, asr; retentionunit=u"minute", figure=(; size=(1000, 350)))
save("ladder.svg", fig₃)
nothing # hide
```

![](ladder.svg)

```@example 1
# Extract spectra for the annotated ladder steps
mss = alkaneladdermassspectra(nmsm, asr)

# Select eicosane (C20), which must have been extracted
c_count = 20
haskey(mss, c_count) || throw(ArgumentError("No spectrum extracted for C$c_count"))
ms_c20 = mss[c_count]

# Plot the selected ladder spectrum
fig₄ = massspectrum(ms_c20)
save("ms_c20.svg", fig₄)
nothing # hide
```

![](ms_c20.svg)

```@example 1
# Fit an intensity-variance model from the ladder peaks
avf = fitalkanevariancemodel(nmsm, asr)
avf.success || throw(ArgumentError(
    "Failed to infer a variance model from standard file: $file"))

# Attach predicted variances to the intensity matrix
vmsm = varpred(nmsm, avf)

# Fit a smooth map from retention time to Kováts retention index
mapper = fitmap(asr)

# Map retention times to indices and propagate intensities and variances
vmsm_ri = applymap(mapper, vmsm)

# Plot the retention-index TIC
fig₅ = tictrace(vmsm_ri; figure=(; size=(1000, 350)))
save("ri_tic.svg", fig₅)
nothing # hide
```

![](ri_tic.svg)

```@example 1
# Estimate a smooth baseline on the mapped data
baseline = arpls(vmsm_ri)

# Subtract the baseline from the mapped data
signal = subtractbaseline(vmsm_ri, baseline)

# Plot the baseline-corrected TIC
fig₆ = tictrace(signal; figure=(; size=(1000, 350)))
save("signal.svg", fig₆)
nothing # hide
```

![](signal.svg)

## Disclaimer

JuChrom is provided "as is," without warranty of any kind. Users are responsible for
independently validating all outputs and interpretations and for determining suitability
for their specific applications. The authors and contributors disclaim any liability for
errors, omissions, or any consequences arising from use of the software, including use
in regulated, clinical, or safety-critical contexts.
