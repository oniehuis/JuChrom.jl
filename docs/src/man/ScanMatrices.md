# ScanMatrices

Mass scan matrices are commonly built from a `MassScanSeries` using [`mscanmatrix`](@ref).
This aligns the m/z grid across scans and produces a matrix representation.

```julia
# Convert a MassScanSeries to a MassScanMatrix
msm = mscanmatrix(mss)
```

```@docs
AbstractMassScanMatrix
MassScanMatrix
MassScanMatrix(::AbstractVector{<:Real}, ::Union{Nothing, Unitful.Units}, ::AbstractVector{<:Real}, ::Union{Nothing, Unitful.Units}, ::AbstractMatrix{<:Real}, ::Union{Nothing, Unitful.Units}; level::Integer=1, instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::AbstractDict=Dict{String, Any}())
MassScanMatrix(::AbstractVector{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}, ::AbstractVector{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}, ::AbstractMatrix{<:Union{Real, Unitful.AbstractQuantity{<:Real}}}; level::Integer=1, instrument::NamedTuple=NamedTuple(), acquisition::NamedTuple=NamedTuple(), user::NamedTuple=NamedTuple(), sample::NamedTuple=NamedTuple(), extras::AbstractDict=Dict{String, Any}())
acquisition(::AbstractMassScanMatrix)
extras(::AbstractMassScanMatrix)
instrument(::AbstractMassScanMatrix)
intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})
intensityunit(::AbstractMassScanMatrix)
level(::AbstractMassScanMatrix)
mzcount(::AbstractMassScanMatrix)
mzunit(::AbstractMassScanMatrix)
mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})
rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})
rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})
rawretentions(::AbstractMassScanMatrix{Nothing})
retentions(::AbstractMassScanMatrix{Nothing})
retentionunit(::AbstractMassScanMatrix)
sample(::AbstractMassScanMatrix)
scancount(::AbstractMassScanMatrix)
user(::AbstractMassScanMatrix)
Base.:(-)(::MassScanMatrix, ::MassScanMatrix)
```
