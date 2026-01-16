# Trimming and filtering

Filtering tools operate on scan series to remove scans outside a region of interest while
leaving values unchanged. Support for [`MassScanMatrix`](@ref JuChrom.MassScanMatrix) is
not yet provided, because these tools target preprocessing workflows that typically
operate on [`AbstractScanSeries`](@ref JuChrom.AbstractScanSeries). Index-based trimming
([`indextrim`](@ref), [`indextrim!`](@ref)) keeps a contiguous slice by scan position 
(e.g., scans 2 through 100) regardless of retention values, while retention-based trimming 
([`retentiontrim`](@ref), [`retentiontrim!`](@ref)) keeps the scans whose retention times 
fall inside a given window. Use index trimming when you already know the scan indices to 
keep or when spacing is uniform across runs; use retention trimming when you want a fixed 
time window or when scan spacing varies across runs. Typical use cases include dropping 
solvent delay and equilibration regions, isolating a peak cluster for targeted analysis, 
and trimming tails before alignment. For MS data, [`levelscans`](@ref) filters to a target 
MS level (e.g., MS¹ only) so you can compare like with like or focus downstream processing 
on a specific acquisition layer. The default operations return new series and leave the 
original untouched; use the mutating variants when you want to edit a series in place.

## Example

```@example 1
# Load JuChrom and the Agilent ChemStation MS loader
using JuChrom
using JuChrom.ChemStationMSLoader

# Load an example Agilent ChemStation GC-MS run
file = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
series = load(ChemStationMS(file; mode=:ms))

# Index-based trimming
series_trimmed_by_index = indextrim(series; start=50, stop=150)

# Retention-based trimming
series_trimmed_by_retention = retentiontrim(series; start=10u"minute", stop=30u"minute")

# Mutating variants operate in place
series_mut = deepcopy(series)
retentiontrim!(series_mut; start=10u"minute", stop=30u"minute")
after_retention = scancount(series_mut)
indextrim!(series_mut; start=1, stop=50)
after_index = scancount(series_mut)

# Filter MS level (e.g., MS1 only)
# Note that this example series only includes MS1; use levels(series) to check levels.
series_ms1 = levelscans(series, 1)

# Quick checks
(; orig=scancount(series),
   by_index=scancount(series_trimmed_by_index),
   by_retention=scancount(series_trimmed_by_retention),
   ms1=scancount(series_ms1),
   after_retention=after_retention,
   after_index=after_index)
```

## Filtering tools

```@docs
indextrim
indextrim!
levelscans
retentiontrim
retentiontrim!
```
