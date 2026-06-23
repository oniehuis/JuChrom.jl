# Retention binning and gridding

```@docs
RetentionGrid
extras(::RetentionGrid)
retentionunit(::RetentionGrid)
rawbinedges
binedges
rawbinwidth
binwidth
rawtolerance
tolerance
rawoverlapmin
rawoverlapmax
overlapmin
overlapmax
binretentions(::MassScanMatrix, ::AbstractVector{<:Number}, ::Union{AbstractVector{<:LinearObservedIntensityVarianceModel},LinearObservedIntensityVarianceModel})
binretentions(::AbstractVarianceMassScanMatrix, ::AbstractVector{<:Number}, ::Union{Real, AbstractVector{<:Real}})
densestgrid(::AbstractVector{<:AbstractMassScanMatrix})
densestgrid(::AbstractDict{<:Any,<:AbstractMassScanMatrix})
```
