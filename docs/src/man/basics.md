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
integer
intensities(::AbstractChrom; ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
intensities(::AbstractChromMS; ::OrdinalRange{T1, S1}, ::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
intensity(::AbstractChrom, ::Integer)
intensity(::AbstractChrom, ::Unitful.Time; ::Bool)
intensity(::AbstractChromMS, ::Integer, ::Integer)
intensity(::AbstractChromMS, ::Unitful.Time, ::Real; ::Bool)
ion(::AbstractChromMS, ::Integer)
ioncount(::AbstractChromMS)
ionindex(::AbstractChromMS, ::Real)
ionscantime
ionscantimeindex
ionscantimes
ionscantimeshift
ions(::AbstractChromMS)
maxintensity(::AbstractChrom; ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
maxintensity(::AbstractChromMS; ::OrdinalRange{T1, S1}, ::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
maxion(::AbstractChromMS)
maxscantime
metadata(::AbstractChromatogram)
minintensity(::AbstractChrom; ::OrdinalRange{T, S}) where {T<:Integer, S<:Integer}
minintensity(::AbstractChromMS; ::OrdinalRange{T1, S1}, ::OrdinalRange{T2, S2}) where {T1<:Integer, S1<:Integer, T2<:Integer, S2<:Integer}
minion(::AbstractChromMS)
minscantime
rimapper(::AbstractChromatogram)
rimapper!(::AbstractChromatogram, ::AbstractRiMapper)
runduration
scancount
scanduration
scantime
scantimeindex
scantimes
totalionchromatogram
```