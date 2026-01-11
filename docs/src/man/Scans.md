# Scans

## AbstractScan

```@docs
AbstractScan
attrs(::AbstractScan)
intensityunit(::AbstractScan)
rawretention(::AbstractScan{Nothing})
retention(::AbstractScan{Nothing})
retentionunit(::AbstractScan)
```

## AbstractChromScan

```@docs
AbstractChromScan
intensity(::AbstractChromScan{<:Any, Nothing})
rawintensity(::AbstractChromScan{<:Any, Nothing})
```

## ChromScan (concrete type)

```@docs
ChromScan
ChromScan(::Unitful.AbstractQuantity, ::Unitful.AbstractQuantity; attrs::NamedTuple=NamedTuple())
ChromScan(::Real, ::Union{Unitful.Units, Nothing}, ::Real, ::Union{Unitful.Units, Nothing}; attrs::NamedTuple=NamedTuple())
Base.:(==)(::ChromScan, ::ChromScan)
```

## AbstractMassScan

```@docs
AbstractMassScan
intensities(::AbstractMassScan{<:Any, <:Any, Nothing})
level(::AbstractMassScan)
mzcount(::AbstractMassScan)
mzunit(::AbstractMassScan)
mzvalues(::AbstractMassScan{<:Any, Nothing, <:Any})
rawintensities(::AbstractMassScan{<:Any, <:Any, Nothing})
rawmzvalues(::AbstractMassScan{<:Any, <:Any, Nothing})
```

## MassScan (concrete type)

```@docs
MassScan
MassScan(::Unitful.AbstractQuantity, ::AbstractVector{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}, ::AbstractVector{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}; level::Integer=1, attrs::NamedTuple=NamedTuple())
MassScan(::Real, ::Union{Unitful.Units, Nothing}, ::AbstractVector{<:Real}, ::Union{Unitful.Units, Nothing}, ::AbstractVector{<:Real}, ::Union{Unitful.Units, Nothing}; level::Integer=1, attrs::NamedTuple=NamedTuple())
Base.:(==)(::MassScan, ::MassScan)
```
