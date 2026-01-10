# Scans

```@docs
AbstractScan
AbstractChromScan
ChromScan
ChromScan(::Real, ::Union{Unitful.Units, Nothing}, ::Real, ::Union{Unitful.Units, Nothing}; attrs::NamedTuple=NamedTuple())
AbstractMassScan
MassScan
MassScan(::Real, ::Union{Unitful.Units, Nothing}, ::AbstractVector{<:Real}, ::Union{Unitful.Units, Nothing}, ::AbstractVector{<:Real}, ::Union{Unitful.Units, Nothing}; level::Integer=1, attrs::NamedTuple=NamedTuple())
attrs(::AbstractScan)
intensity(::AbstractChromScan{<:Any, Nothing})
intensityunit(::AbstractScan)
intensities(::AbstractMassScan{<:Any, <:Any, Nothing})
level(::AbstractMassScan)
mzcount(::AbstractMassScan)
mzunit(::AbstractMassScan)
mzvalues(::AbstractMassScan{<:Any, Nothing, <:Any})
rawintensities(::AbstractMassScan{<:Any, <:Any, Nothing})
rawintensity(::AbstractChromScan{<:Any, Nothing})
rawmzvalues(::AbstractMassScan{<:Any, <:Any, Nothing})
rawretention(::AbstractScan{Nothing})
retention(::AbstractScan{Nothing})
retentionunit(::AbstractScan)
Base.:(==)(::ChromScan, ::ChromScan)
Base.:(==)(::MassScan, ::MassScan)
```
