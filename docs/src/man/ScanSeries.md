# ScansSeries

## AbstractScanSeries

```@docs
AbstractScanSeries
acquisition(::AbstractScanSeries)
extras(::AbstractScanSeries)
instrument(::AbstractScanSeries)
intensityunit(::AbstractScanSeries)
rawretentions(::AbstractScanSeries)
retentions(::AbstractScanSeries)
retentionunit(::AbstractScanSeries)
sample(::AbstractScanSeries)
scan(::AbstractScanSeries, ::Integer)
scancount(::AbstractScanSeries)
scans(::AbstractScanSeries)
user(::AbstractScanSeries) 
```

## AbstractChromScanSeries

```@docs
AbstractChromScanSeries
intensities(::AbstractChromScanSeries)
intensity(::AbstractChromScanSeries, ::Integer)
rawintensities(::AbstractChromScanSeries)
```

## ChromScanSeries (concrete type)

```@docs
ChromScanSeries
ChromScanSeries(::AbstractVector{<:AbstractChromScan}; instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::Dict{<:AbstractString, <:Any}=Dict())
```

## AbstractMassScanSeries

```@docs
AbstractMassScanSeries
intensities(::AbstractMassScanSeries, ::Integer)
levels(::AbstractMassScanSeries)
mzunit(::AbstractMassScanSeries)
mzvalues(::AbstractMassScanSeries, ::Integer)
rawintensities(::AbstractMassScanSeries, ::Integer)
uniquemzvalues(::AbstractMassScanSeries, ::Integer)
```

## MassScanSeries (concrete type)

```@docs
MassScanSeries
MassScanSeries(::AbstractVector{<:AbstractMassScan}; instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::Dict{<:AbstractString, <:Any}=Dict())
mscanmatrix
mzchrom(::MassScanSeries)
```
