# Basics

```@docs
AbstractChromatogram
AbstractFID
AbstractGC
AbstractGCMS
AbstractTIC
FID
GCMS
IonScanOrder
LinearAscending
LinearDescending
TIC
binions
cosine
integer
intensities(::AbstractChromatogram)
intensities(::AbstractGC, ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
intensities(::AbstractGCMS, ::OrdinalRange{<:T1, <:S1}, ::OrdinalRange{<:T2, <:S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
intensity(::AbstractGC, ::Integer)
intensity(::AbstractGC, ::Unitful.Time; ::Bool)
intensity(::AbstractGCMS, ::Integer, ::Integer)
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