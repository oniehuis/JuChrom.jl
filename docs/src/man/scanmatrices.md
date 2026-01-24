# Scan matrices

## AbstractMassScanMatrix

```@docs
JuChrom.AbstractMassScanMatrix
JuChrom.acquisition(::AbstractMassScanMatrix)
JuChrom.extras(::AbstractMassScanMatrix)
JuChrom.instrument(::AbstractMassScanMatrix)
JuChrom.intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})
JuChrom.intensityunit(::AbstractMassScanMatrix)
JuChrom.level(::AbstractMassScanMatrix)
JuChrom.mzcount(::AbstractMassScanMatrix)
JuChrom.mzindex(::AbstractMassScanMatrix{<:Any, Nothing}, ::Real)
JuChrom.mzunit(::AbstractMassScanMatrix)
JuChrom.mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})
JuChrom.rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})
JuChrom.rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})
JuChrom.rawretentions(::AbstractMassScanMatrix{Nothing})
JuChrom.retentions(::AbstractMassScanMatrix{Nothing})
JuChrom.retentionunit(::AbstractMassScanMatrix)
JuChrom.sample(::AbstractMassScanMatrix)
JuChrom.scancount(::AbstractMassScanMatrix)
JuChrom.user(::AbstractMassScanMatrix)
```

## MassScanMatrix

```@docs
JuChrom.MassScanMatrix
JuChrom.MassScanMatrix(::AbstractVector{<:Real}, ::Union{Nothing, Unitful.Units}, ::AbstractVector{<:Real}, ::Union{Nothing, Unitful.Units}, ::AbstractMatrix{<:Real}, ::Union{Nothing, Unitful.Units}; level::Integer=1, instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::AbstractDict=Dict{String, Any}())
JuChrom.MassScanMatrix(::AbstractVector{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}, ::AbstractVector{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}, ::AbstractMatrix{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}; level::Integer=1, instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::AbstractDict=Dict{String, Any}())
Base.:(-)(::MassScanMatrix, ::MassScanMatrix)
```
