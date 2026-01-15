# Scan series

## AbstractScanSeries

```@docs
JuChrom.AbstractScanSeries
JuChrom.acquisition(::AbstractScanSeries)
JuChrom.extras(::AbstractScanSeries)
JuChrom.instrument(::AbstractScanSeries)
JuChrom.intensityunit(::AbstractScanSeries)
JuChrom.rawretentions(::AbstractScanSeries)
JuChrom.retentions(::AbstractScanSeries)
JuChrom.retentionunit(::AbstractScanSeries)
JuChrom.sample(::AbstractScanSeries)
JuChrom.scan(::AbstractScanSeries, ::Integer)
JuChrom.scancount(::AbstractScanSeries)
JuChrom.scans(::AbstractScanSeries)
JuChrom.user(::AbstractScanSeries) 
```

## AbstractChromScanSeries

```@docs
JuChrom.AbstractChromScanSeries
JuChrom.intensities(::AbstractChromScanSeries)
JuChrom.intensity(::AbstractChromScanSeries, ::Integer)
JuChrom.rawintensity(::AbstractChromScanSeries, ::Integer)
JuChrom.rawintensities(::AbstractChromScanSeries)
```

## ChromScanSeries

```@docs
JuChrom.ChromScanSeries
JuChrom.ChromScanSeries(::AbstractVector{<:AbstractChromScan}; instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::Dict{<:AbstractString, <:Any}=Dict())
Base.:(==)(::ChromScanSeries, ::ChromScanSeries)
```

## AbstractMassScanSeries

```@docs
JuChrom.AbstractMassScanSeries
JuChrom.intensities(::AbstractMassScanSeries, ::Integer)
JuChrom.levels(::AbstractMassScanSeries)
JuChrom.mzunit(::AbstractMassScanSeries)
JuChrom.mzvalues(::AbstractMassScanSeries, ::Integer)
JuChrom.rawintensities(::AbstractMassScanSeries, ::Integer)
JuChrom.uniquemzvalues(::AbstractMassScanSeries, ::Integer)
```

## MassScanSeries

```@docs
JuChrom.MassScanSeries
JuChrom.MassScanSeries(::AbstractVector{<:AbstractMassScan}; instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::Dict{<:AbstractString, <:Any}=Dict())
Base.:(==)(::MassScanSeries, ::MassScanSeries)
```
