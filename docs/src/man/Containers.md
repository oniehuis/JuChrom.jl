# Core Containers

JuChrom provides three core container types for chromatographic and mass spectrometric
data. Single scans capture one measurement, while scan series and mass-scan matrices
provide two alternative multi-scan representations.

Raw getters (e.g., [`rawretention`](@ref JuChrom.rawretention(::AbstractScan{Nothing})),
[`rawintensity`](@ref JuChrom.rawintensity(::AbstractChromScan{<:Any, Nothing})),
[`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScan{<:Any, <:Any, Nothing})))
return unitless numeric values. Unit-aware convenience functions (e.g.,
[`retention`](@ref JuChrom.retention(::AbstractScan{Nothing})),
[`intensity`](@ref JuChrom.intensity(::AbstractChromScan{<:Any, Nothing})),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScan{<:Any, Nothing, <:Any}))) preserve
or convert units when available.

```julia
rawretention(scan)               # -> 60.0
retention(scan)                  # -> 60.0u"s"
retention(scan, unit=u"minute")  # -> 1.0u"minute"
```

## Single-scan containers

Single scans store one chromatographic point (such as [`ChromScan`](@ref JuChrom.ChromScan)) or one mass
spectrometric scan (such as [`MassScan`](@ref JuChrom.MassScan)).

### AbstractScan and subtypes

The abstract scan types define the common interface for single scans. Getter functions
return the stored, unitless values; convenience functions apply unit handling or
conversions.

**AbstractScan** ([`AbstractScan`](@ref JuChrom.AbstractScan)) expects:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `retention` | `Real` | [`rawretention`](@ref JuChrom.rawretention(::AbstractScan{Nothing})) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref JuChrom.retentionunit(::AbstractScan)) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref JuChrom.intensityunit(::AbstractScan)) |
| `attrs` | `NamedTuple` | [`attrs`](@ref JuChrom.attrs(::AbstractScan)) |

**AbstractChromScan** ([`AbstractChromScan`](@ref JuChrom.AbstractChromScan)) adds:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `intensity` | `Real` | [`rawintensity`](@ref JuChrom.rawintensity(::AbstractChromScan{<:Any, Nothing})) |

**AbstractMassScan** ([`AbstractMassScan`](@ref JuChrom.AbstractMassScan)) adds:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `mzvalues` | `AbstractVector{<:Real}` | [`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScan{<:Any, <:Any, Nothing})) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref JuChrom.mzunit(::AbstractMassScan)) |
| `intensities` | `AbstractVector{<:Real}` | [`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScan{<:Any, <:Any, Nothing})) |
| `level` | `Integer` | [`level`](@ref JuChrom.level(::AbstractMassScan)) |

Convenience:

- [`intensities`](@ref JuChrom.intensities(::AbstractMassScan{<:Any, <:Any, Nothing})),
  [`intensity`](@ref JuChrom.intensity(::AbstractChromScan{<:Any, Nothing})),
  [`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScan{<:Any, Nothing, <:Any})),
  [`retention`](@ref JuChrom.retention(::AbstractScan{Nothing}))

Concrete types such as [`ChromScan`](@ref JuChrom.ChromScan) and 
[`MassScan`](@ref JuChrom.MassScan) are subtypes that implement these abstract interfaces.

If you define custom container types, subtype the corresponding abstract type. The
generic getters provide the shared API; you are free to add additional fields and
methods as needed.

### ChromScan

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `retention` | `Real` | [`rawretention`](@ref JuChrom.rawretention(::AbstractScan{Nothing})) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref JuChrom.retentionunit(::AbstractScan)) |
| `intensity` | `Real` | [`rawintensity`](@ref JuChrom.rawintensity(::AbstractChromScan{<:Any, Nothing})) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref JuChrom.intensityunit(::AbstractScan)) |
| `attrs` | `NamedTuple` | [`attrs`](@ref JuChrom.attrs(::AbstractScan)) |

Convenience:

- [`intensity`](@ref JuChrom.intensity(::AbstractChromScan{<:Any, Nothing})),
  [`retention`](@ref JuChrom.retention(::AbstractScan{Nothing}))

### MassScan

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `retention` | `Real` | [`rawretention`](@ref JuChrom.rawretention(::AbstractScan{Nothing})) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref JuChrom.retentionunit(::AbstractScan)) |
| `mzvalues` | `AbstractVector{<:Real}` | [`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScan{<:Any, <:Any, Nothing})) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref JuChrom.mzunit(::AbstractMassScan)) |
| `intensities` | `AbstractVector{<:Real}` | [`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScan{<:Any, <:Any, Nothing})) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref JuChrom.intensityunit(::AbstractScan)) |
| `level` | `Integer` | [`level`](@ref JuChrom.level(::AbstractMassScan)) |
| `attrs` | `NamedTuple` | [`attrs`](@ref JuChrom.attrs(::AbstractScan)) |

Convenience:

- [`intensities`](@ref JuChrom.intensities(::AbstractMassScan{<:Any, <:Any, Nothing})),
  [`mzcount`](@ref JuChrom.mzcount(::AbstractMassScan)),
  [`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScan{<:Any, Nothing, <:Any})),
  [`retention`](@ref JuChrom.retention(::AbstractScan{Nothing}))

## Scan series (batches of scans)

Series hold multiple scans plus shared metadata. For series-level getters, the `raw*`
functions return unitless numeric arrays, while the non-raw variants preserve or convert
units.

### AbstractScanSeries and subtypes

The abstract series types define the shared metadata and series API for scan access.

**AbstractScanSeries** ([`AbstractScanSeries`](@ref JuChrom.AbstractScanSeries)) expects:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `scans` | `AbstractVector{<:AbstractScan}` | [`scans`](@ref JuChrom.scans(::AbstractScanSeries)) |
| `instrument` | `NamedTuple` | [`instrument`](@ref JuChrom.instrument(::AbstractScanSeries)) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref JuChrom.acquisition(::AbstractScanSeries)) |
| `user` | `NamedTuple` | [`user`](@ref JuChrom.user(::AbstractScanSeries)) |
| `sample` | `NamedTuple` | [`sample`](@ref JuChrom.sample(::AbstractScanSeries)) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref JuChrom.extras(::AbstractScanSeries)) |

Convenience:

- [`scan`](@ref JuChrom.scan(::AbstractScanSeries, ::Integer)),
  [`scancount`](@ref JuChrom.scancount(::AbstractScanSeries))

**AbstractChromScanSeries** ([`AbstractChromScanSeries`](@ref JuChrom.AbstractChromScanSeries)) 
and **AbstractMassScanSeries** ([`AbstractMassScanSeries`](@ref JuChrom.AbstractMassScanSeries)) 
specialize the scan type and enable series-level convenience getters such as
[`retentions`](@ref JuChrom.retentions(::AbstractScanSeries)),
[`intensities`](@ref JuChrom.intensities(::AbstractChromScanSeries)),
[`intensities`](@ref JuChrom.intensities(::AbstractMassScanSeries, ::Integer)),
[`levels`](@ref JuChrom.levels(::AbstractMassScanSeries)), and
[`uniquemzvalues`](@ref JuChrom.uniquemzvalues(::AbstractMassScanSeries, ::Integer)).

Concrete types such as [`ChromScanSeries`](@ref JuChrom.ChromScanSeries) and 
[`MassScanSeries`](@ref JuChrom.MassScanSeries) are subtypes that implement these abstract 
interfaces.

### ChromScanSeries

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `scans` | `AbstractVector{<:ChromScan}` | [`scans`](@ref JuChrom.scans(::AbstractScanSeries)) |
| `instrument` | `NamedTuple` | [`instrument`](@ref JuChrom.instrument(::AbstractScanSeries)) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref JuChrom.acquisition(::AbstractScanSeries)) |
| `user` | `NamedTuple` | [`user`](@ref JuChrom.user(::AbstractScanSeries)) |
| `sample` | `NamedTuple` | [`sample`](@ref JuChrom.sample(::AbstractScanSeries)) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref JuChrom.extras(::AbstractScanSeries)) |

Convenience:

- [`intensities`](@ref JuChrom.intensities(::AbstractChromScanSeries)),
  [`intensity`](@ref JuChrom.intensities(::AbstractChromScanSeries, ::Integer)),
  [`rawintensities`](@ref JuChrom.rawintensities(::AbstractChromScanSeries)),
  [`rawintensity`](@ref JuChrom.rawintensities(::AbstractChromScanSeries, ::Integer)),
  [`retentions`](@ref JuChrom.retentions(::AbstractScanSeries)),
  [`scancount`](@ref JuChrom.scancount(::AbstractScanSeries))

### MassScanSeries

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `scans` | `AbstractVector{<:MassScan}` | [`scans`](@ref JuChrom.scans(::AbstractScanSeries)) |
| `instrument` | `NamedTuple` | [`instrument`](@ref JuChrom.instrument(::AbstractScanSeries)) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref JuChrom.acquisition(::AbstractScanSeries)) |
| `user` | `NamedTuple` | [`user`](@ref JuChrom.user(::AbstractScanSeries)) |
| `sample` | `NamedTuple` | [`sample`](@ref JuChrom.sample(::AbstractScanSeries)) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref JuChrom.extras(::AbstractScanSeries)) |

Convenience:

- [`intensities`](@ref JuChrom.intensities(::AbstractMassScanSeries, ::Integer)),
  [`levels`](@ref JuChrom.levels(::AbstractMassScanSeries)),
  [`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanSeries)),
  [`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScanSeries, ::Integer)),
  [`retentions`](@ref JuChrom.retentions(::AbstractScanSeries)),
  [`scancount`](@ref JuChrom.scancount(::AbstractScanSeries)),
  [`uniquemzvalues`](@ref JuChrom.uniquemzvalues(::AbstractMassScanSeries, ::Integer))

## Mass-scan matrices (aligned grid)

Mass-scan matrices store aligned mass scans on a shared m/z grid.
Use series when you need scan-level structure or heterogeneous m/z grids; use matrices
for aligned m/z operations and linear algebra workflows.

### AbstractMassScanMatrix interface

The abstract matrix type defines the aligned grid interface for mass scans.

**AbstractMassScanMatrix** ([`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix)) 
expects:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `retentions` | `AbstractVector{<:Real}` | [`rawretentions`](@ref JuChrom.rawretentions(::AbstractMassScanMatrix{Nothing})) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref JuChrom.retentionunit(::AbstractMassScanMatrix)) |
| `mzvalues` | `AbstractVector{<:Real}` | [`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanMatrix)) |
| `intensities` | `AbstractMatrix{<:Real}` | [`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassScanMatrix)) |
| `level` | `Integer` | [`level`](@ref JuChrom.level(::AbstractMassScanMatrix)) |
| `instrument` | `NamedTuple` | [`instrument`](@ref JuChrom.instrument(::AbstractMassScanMatrix)) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref JuChrom.acquisition(::AbstractMassScanMatrix)) |
| `user` | `NamedTuple` | [`user`](@ref JuChrom.user(::AbstractMassScanMatrix)) |
| `sample` | `NamedTuple` | [`sample`](@ref JuChrom.sample(::AbstractMassScanMatrix)) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref JuChrom.extras(::AbstractMassScanMatrix)) |

### MassScanMatrix

[`MassScanMatrix`](@ref JuChrom.MassScanMatrix) is a subtype of 
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix) that implements this 
interface.

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `retentions` | `AbstractVector{<:Real}` | [`rawretentions`](@ref JuChrom.rawretentions(::AbstractMassScanMatrix{Nothing})) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref JuChrom.retentionunit(::AbstractMassScanMatrix)) |
| `mzvalues` | `AbstractVector{<:Real}` | [`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanMatrix)) |
| `intensities` | `AbstractMatrix{<:Real}` | [`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassScanMatrix)) |
| `level` | `Integer` | [`level`](@ref JuChrom.level(::AbstractMassScanMatrix)) |
| `instrument` | `NamedTuple` | [`instrument`](@ref JuChrom.instrument(::AbstractMassScanMatrix)) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref JuChrom.acquisition(::AbstractMassScanMatrix)) |
| `user` | `NamedTuple` | [`user`](@ref JuChrom.user(::AbstractMassScanMatrix)) |
| `sample` | `NamedTuple` | [`sample`](@ref JuChrom.sample(::AbstractMassScanMatrix)) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref JuChrom.extras(::AbstractMassScanMatrix)) |

Convenience:

- [`intensities`](@ref JuChrom.intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})),
  [`mzcount`](@ref JuChrom.mzcount(::AbstractMassScanMatrix)),
  [`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})),
  [`retentions`](@ref JuChrom.retentions(::AbstractMassScanMatrix{Nothing})),
  [`scancount`](@ref JuChrom.scancount(::AbstractMassScanMatrix))

Use [`mscanmatrix`](@ref JuChrom.mscanmatrix(::JuChrom.MassScanSeries, ::JuChrom.AbstractMatrixFormat)) 
to convert a [`MassScanSeries`](@ref JuChrom.MassScanSeries) to a 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix):

```julia
msm = mscanmatrix(mss)
```

See the dedicated container pages for the full API reference:
[`Scans`](Scans.md), [`Scan Series`](ScanSeries.md), and [`Scan Matrices`](ScanMatrices.md).
