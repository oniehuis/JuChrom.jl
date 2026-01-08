# ScansSeries

```@docs
AbstractScanSeries
AbstractChromScanSeries
ChromScanSeries
AbstractMassScanSeries
MassScanSeries
acquisition(::AbstractScanSeries)
extras(::AbstractScanSeries)
intensities(::AbstractChromScanSeries)
intensities(::AbstractMassScanSeries, ::Integer)
intensity(::AbstractChromScanSeries, ::Integer)
intensityunit(::AbstractScanSeries)
instrument(::AbstractScanSeries)
levels(::AbstractMassScanSeries)
mscanmatrix
mzchrom(::MassScanSeries)
mzunit(::AbstractMassScanSeries)
mzvalues(::AbstractMassScanSeries, ::Integer)
rawintensities(::AbstractChromScanSeries)
rawintensities(::AbstractMassScanSeries, ::Integer)
rawretentions(::AbstractScanSeries)
retentions(::AbstractScanSeries)
retentionunit(::AbstractScanSeries)
sample(::AbstractScanSeries)
scan(::AbstractScanSeries, ::Integer)
scancount(::AbstractScanSeries)
scans(::AbstractScanSeries)
uniquemzvalues(::AbstractMassScanSeries, ::Integer)
user(::AbstractScanSeries) 
```