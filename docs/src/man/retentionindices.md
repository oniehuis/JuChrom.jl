# Retention Indices

```@docs
PolationMethod
Linear
PiecewiseLinear
NaturalCubicBSpline
AbstractRiMapper
RiMapper
extrapolationmethod
interpolationmethod
maxretentionindex
maxretentiontime
metadata(::AbstractRiMapper)
minretentionindex
minretentiontime
retentionindex(::RiMapper, ::Unitful.Time; ::Bool)
retentionindex(::AbstractChromatogram, ::Unitful.Time; ::Bool)
retentionindexname(::RiMapper)
retentionindices
retentiontimes
```