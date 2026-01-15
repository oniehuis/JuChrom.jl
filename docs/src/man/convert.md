# Format conversion

Conversion tools change representation without changing the underlying data. In JuChrom
this means moving between a scan-oriented view (`MassScanSeries`) and a grid-oriented view
(`MassScanMatrix`), or extracting chromatograms from selected m/z values. These conversions
are often useful after preprocessing (for example, m/z binning to a fixed grid) because
they make the data shape explicit and unlock workflows that benefit from a matrix
representation. Conversion itself does not bin or resample; `mzchrom` can sum intensities
across selected m/z values or tolerances, but it does not merge nearby m/z into new bins.
Binning and resampling are handled by the [binning tools](@ref).

## Example

```@example 1
using JuChrom
using JuChrom.ChemStationMSLoader

file = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
mss = load(ChemStationMS(file; mode=:ms))

# Extract a chromatogram of the target m/z with tolerance on raw m/z values
xic_104_tol = mzchrom(mss, 104.0; tol=1e-2)

# Integer-bin m/z values to a fixed grid
mss_binned = binmzvalues(mss)

# Extract a chromatogram of the target m/z on the binned grid
xic_104 = mzchrom(mss_binned, 104)

# Extract a summed chromatogram across multiple target m/z
xic_43_57_71_85 = mzchrom(mss_binned, [43, 57, 71, 85])  # or 43:14:85

# Extract total ion chromatogram (across all ions)
tic = mzchrom(mss_binned)

# Convert the binned series to a matrix for matrix-based workflows
msm = mscanmatrix(mss_binned)  # dense by default; use SPARSE for sparse storage
```

```@docs
JuChrom.mscanmatrix
JuChrom.mzchrom
```
