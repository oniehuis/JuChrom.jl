# Basics

```@docs
AbstractChromatogram
AbstractGC
AbstractFID
FID
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
maxscantime
metadata
minintensity
minion
minscantime
runduration
scancount
scanduration
scantime
scantimeindex
scantimes
totalionchromatogram
```