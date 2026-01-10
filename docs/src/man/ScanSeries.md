# ScansSeries

```@docs
AbstractScanSeries
AbstractChromScanSeries
ChromScanSeries
ChromScanSeries(::AbstractVector{<:AbstractChromScan}; instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::Dict{<:AbstractString, <:Any}=Dict())
AbstractMassScanSeries
MassScanSeries
MassScanSeries(::AbstractVector{<:AbstractMassScan}; instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::Dict{<:AbstractString, <:Any}=Dict())
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
