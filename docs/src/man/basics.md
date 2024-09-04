# Basics

```@docs
AbstractChromatogram
AbstractChrom
AbstractChromMS
ChromMS
Chrom
IonScanOrder
LinearAscending
LinearDescending
binions
cosine
integer
intensities(::AbstractChrom; ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
intensities(::AbstractChromMS; ::OrdinalRange{T1, S1}, ::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
intensity(::AbstractChrom, ::Integer)
intensity(::AbstractChrom, ::Unitful.Time; ::Bool)
intensity(::AbstractChromMS, ::Integer, ::Integer)
intensity(::AbstractChromMS, ::Unitful.Time, ::Real; ::Bool)
ion
ioncount
ionindex
ionscantime
ionscantimeindex
ionscantimes
ionscantimeshift
ions
maxintensity(::AbstractChrom; ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
maxintensity(::AbstractChromMS; ::OrdinalRange{T1, S1}, ::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
maxion
maxscantime
metadata(::AbstractChromatogram)
minintensity(::AbstractChrom; ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
minintensity(::AbstractChromMS; ::OrdinalRange{T1, S1}, ::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
minion
minscantime
rimapper
runduration
scancount
scanduration
scantime
scantimeindex
scantimes
totalionchromatogram
```