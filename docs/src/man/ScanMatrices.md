# ScanMatrices

```@docs
AbstractMassScanMatrix
MassScanMatrix
acquisition(::AbstractMassScanMatrix)
extras(::AbstractMassScanMatrix)
instrument(::AbstractMassScanMatrix)
intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})
intensityunit(::AbstractMassScanMatrix)
level(::AbstractMassScanMatrix)
mzcount(::AbstractMassScanMatrix)
mzunit(::AbstractMassScanMatrix)
mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})
rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})
rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})
rawretentions(::AbstractMassScanMatrix{Nothing})
retentions(::AbstractMassScanMatrix{Nothing})
retentionunit(::AbstractMassScanMatrix)
sample(::AbstractMassScanMatrix)
scancount(::AbstractMassScanMatrix)
user(::AbstractMassScanMatrix)
Base.:(-)(::MassScanMatrix, ::MassScanMatrix)
```
