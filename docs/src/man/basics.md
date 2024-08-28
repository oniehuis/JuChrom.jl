# Basics

```@docs
AbstractChromatogram
AbstractGC
AbstractFID
FID
RiFID
AbstractTIC
TIC
RiTIC
AbstractGCMS
GCMS
RiGCMS
IonScanOrder
LinearAscending
LinearDescending
binions
bsplineinterpolation
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
maxintensity(::AbstractGC; ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
maxintensity(::AbstractGCMS; ::OrdinalRange{T1, S1}, ::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
maxion
maxretentionindex
maxscantime
metadata
minintensity(::AbstractGC; ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
minintensity(::AbstractGCMS; ::OrdinalRange{T1, S1}, ::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
minion
minretentionindex
minscantime
retentionindex
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