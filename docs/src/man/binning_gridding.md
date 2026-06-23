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
binretentions(::AbstractVarianceMassScanMatrix, ::RetentionGrid, ::Union{AbstractVector{<:Real}, Real})
densestgrid(::AbstractVector{<:AbstractMassScanMatrix})
densestgrid(::AbstractDict{<:Any,<:AbstractMassScanMatrix})
```
