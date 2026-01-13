# Scans

## AbstractScan

```@docs
JuChrom.AbstractScan
JuChrom.attrs(::AbstractScan)
JuChrom.intensityunit(::AbstractScan)
JuChrom.rawretention(::AbstractScan{Nothing})
JuChrom.retention(::AbstractScan{Nothing})
JuChrom.retentionunit(::AbstractScan)
```

## AbstractChromScan

```@docs
JuChrom.AbstractChromScan
JuChrom.intensity(::AbstractChromScan{<:Any, Nothing})
JuChrom.rawintensity(::AbstractChromScan{<:Any, Nothing})
```

## ChromScan

```@docs
JuChrom.ChromScan
JuChrom.ChromScan(::Unitful.AbstractQuantity, ::Unitful.AbstractQuantity; attrs::NamedTuple=NamedTuple())
JuChrom.ChromScan(::Real, ::Union{Unitful.Units, Nothing}, ::Real, ::Union{Unitful.Units, Nothing}; attrs::NamedTuple=NamedTuple())
Base.:(==)(::JuChrom.ChromScan, ::JuChrom.ChromScan)
```

## AbstractMassScan

```@docs
JuChrom.AbstractMassScan
JuChrom.intensities(::AbstractMassScan{<:Any, <:Any, Nothing})
JuChrom.level(::AbstractMassScan)
JuChrom.mzcount(::AbstractMassScan)
JuChrom.mzunit(::AbstractMassScan)
JuChrom.mzvalues(::AbstractMassScan{<:Any, Nothing, <:Any})
JuChrom.rawintensities(::AbstractMassScan{<:Any, <:Any, Nothing})
JuChrom.rawmzvalues(::AbstractMassScan{<:Any, <:Any, Nothing})
```

## MassScan

```@docs
JuChrom.MassScan
JuChrom.MassScan(::Unitful.AbstractQuantity, ::AbstractVector{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}, ::AbstractVector{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}; level::Integer=1, attrs::NamedTuple=NamedTuple())
JuChrom.MassScan(::Real, ::Union{Unitful.Units, Nothing}, ::AbstractVector{<:Real}, ::Union{Unitful.Units, Nothing}, ::AbstractVector{<:Real}, ::Union{Unitful.Units, Nothing}; level::Integer=1, attrs::NamedTuple=NamedTuple())
Base.:(==)(::MassScan, ::MassScan)
```
