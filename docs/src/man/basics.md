# Basics

```@docs
AbstractChromatogram
AbstractGC
AbstractFID
FID
RiFID
AbstractTIC
TIC
AbstractGCMS
GCMS
IonScanOrder
LinearAscending
LinearDescending
binions
cosine
integer
intensities(::AbstractGC; ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
intensities(::AbstractGCMS; ::OrdinalRange{T1, S1}, ::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
intensity(::AbstractGC, ::Integer)
intensity(::AbstractGC, ::Unitful.Time; ::Bool)
intensity(::AbstractGCMS, ::Integer, ::Integer)
intensity(::AbstractGCMS, ::Unitful.Time, ::Real; ::Bool)
ion
ioncount
ionindex
ionscantime
ionscantimeindex
ionscantimes
ionscantimeshift
ions
maxintensity
maxion
maxretentionindex
maxscantime
metadata
minintensity
minion
minretentionindex
minscantime
retentionindices
retentionindexname
runduration
scancount
scanduration
scantime
scantimeindex
scantimes
totalionchromatogram
```