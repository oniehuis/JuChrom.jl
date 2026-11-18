# Format conversion

Conversion tools change representation without changing the underlying data. In JuChrom this 
means moving between a scan-oriented view ([`MassScanSeries`](@ref JuChrom.MassScanSeries)) 
and a grid-oriented view ([`MassScanMatrix`](@ref JuChrom.MassScanMatrix)), or extracting 
chromatograms from selected m/z values. These conversions are often useful after 
preprocessing (for example, m/z binning to a fixed grid) because they make the data shape 
explicit and unlock workflows that benefit from a matrix representation. Conversion itself 
does not bin or resample; [`mzchrom`](@ref) can sum intensities across selected m/z values 
or tolerances, but it does not merge nearby m/z into new bins. Binning and resampling are 
handled by [m/z binning](@ref "m/z binning") and
[retention binning and gridding](@ref "Retention binning and gridding").

## Explicit matrix columns

[`mscanmatrix`](@ref) requires the output m/z axis as its second argument:

```julia
mscanmatrix(series, 29:563)
mscanmatrix(series, (57, 71, 85); unmatched=:drop)
mscanmatrix(series, [57, 71, 85]; format=SPARSE, unmatched=:drop)
```

The grid may be a vector, range, or tuple, with regular or irregular spacing. It must
be nonempty, finite, positive, and strictly increasing; duplicate values are rejected.
Every grid value becomes a column, including ions with no recorded intensities. Supplying
the same grid across runs fixes their m/z axes independently of which ions were reported.
The grid is copied into the result; the input series and grid are not modified.

Matching is exact. Floating-point grids such as `20.1:0.1:200.4` are supported, but
recorded values must equal the generated grid values. There is no implicit rounding or
tolerance, and `57.1f0` does not equal `57.1`. Use [`binmzvalues`](@ref) beforehand to
assign nearby masses to shared values. Its `validmzvalues` argument only filters bins;
it does not add missing columns.

The default `unmatched=:error` raises an error for any recorded m/z absent from the grid.
Use `unmatched=:drop` for an intentional subset. This check applies only to scans at
`target_level` and happens **before** intensity filtering, so increasing `threshold`
cannot hide a grid mismatch. Only intensities strictly greater than the threshold are
retained. No scans at the selected level is an error; an entirely zero matrix is valid.
Dense and sparse storage use the same matching and filtering rules.

For a series with a unitful m/z axis, grids may also contain quantities:

```julia
mscanmatrix(unitful_series, (57u"Th", 71u"Th"); unmatched=:drop)
mscanmatrix(unitful_series, [57000, 71000]u"mTh"; unmatched=:drop)
mscanmatrix(unitful_series, (20.1u"Th"):(0.1u"Th"):(200.4u"Th"); unmatched=:drop)
```

Compatible grid units are converted to the series' m/z unit before exact matching, and
that unit is retained in the result. Plain numeric grids use the series' m/z unit, if
present. Mixing plain numbers with quantities, using incompatible units, or supplying a
unitful grid for a unitless series is an error. Retention and intensity units and series
metadata are preserved as well.

Zero filling is a representation choice: it does not distinguish an unacquired ion from
an acquired ion omitted by the instrument. An explicit matrix grid does not by itself
establish acquisition coverage or dwell times.

## Example

```@example 1
# Load JuChrom and the Agilent ChemStation MS loader
using JuChrom
using JuChrom.ChemStationMSLoader

# Load an example Agilent ChemStation GC-MS run
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
msm = mscanmatrix(mss_binned, 29:563)  # nominal grid from the example acquisition method
```

# Conversion tools

```@docs
JuChrom.mscanmatrix
JuChrom.mzchrom
```
