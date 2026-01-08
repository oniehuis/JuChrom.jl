# Scans

```@docs
AbstractScan
AbstractChromScan
ChromScan
AbstractMassScan
MassScan
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
```