# Core Containers

JuChrom provides three core container types for chromatographic and mass spectrometric
data. Single scans capture one measurement, while scan series and mass-scan matrices
provide two alternative multi-scan representations.

Raw getters (e.g., [`rawretention`](@ref), [`rawintensity`](@ref), [`rawmzvalues`](@ref))
return unitless numeric values. Unit-aware convenience functions (e.g., [`retention`](@ref),
[`intensity`](@ref), [`mzvalues`](@ref)) preserve or convert units when available.

```julia
rawretention(scan)  # -> 12.3
retention(scan)     # -> 12.3u"s"
```

## Single-scan containers

Single scans store one chromatographic point (such as `ChromScan`) or one mass
spectrometric scan (such as `MassScan`).

### AbstractScan and subtypes

The abstract scan types define the common interface for single scans. Getter functions
return the stored, unitless values; convenience functions apply unit handling or
conversions.

**AbstractScan** ([`AbstractScan`](@ref)) expects:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `retention` | `Real` | [`rawretention`](@ref) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref) |
| `attrs` | `NamedTuple` | [`attrs`](@ref) |

**AbstractChromScan** ([`AbstractChromScan`](@ref)) adds:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `intensity` | `Real` | [`rawintensity`](@ref) |

**AbstractMassScan** ([`AbstractMassScan`](@ref)) adds:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `mzvalues` | `AbstractVector{<:Real}` | [`rawmzvalues`](@ref) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref) |
| `intensities` | `AbstractVector{<:Real}` | [`rawintensities`](@ref) |
| `level` | `Integer` | [`level`](@ref) |

Convenience:

- [`retention`](@ref), [`intensity`](@ref), [`mzvalues`](@ref), [`intensities`](@ref)

Concrete types such as `ChromScan` and `MassScan` are subtypes that implement these
abstract interfaces.

If you define custom container types, subtype the corresponding abstract type. The
generic getters provide the shared API; you are free to add additional fields and
methods as needed.

### ChromScan

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `retention` | `Real` | [`rawretention`](@ref) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref) |
| `intensity` | `Real` | [`rawintensity`](@ref) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref) |
| `attrs` | `NamedTuple` | [`attrs`](@ref) |

Convenience:

- [`retention`](@ref), [`intensity`](@ref), [`retentionunit`](@ref), [`intensityunit`](@ref)

### MassScan

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `retention` | `Real` | [`rawretention`](@ref) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref) |
| `mzvalues` | `AbstractVector{<:Real}` | [`rawmzvalues`](@ref) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref) |
| `intensities` | `AbstractVector{<:Real}` | [`rawintensities`](@ref) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref) |
| `level` | `Integer` | [`level`](@ref) |
| `attrs` | `NamedTuple` | [`attrs`](@ref) |

Convenience:

- [`retention`](@ref), [`mzvalues`](@ref), [`intensities`](@ref), [`mzcount`](@ref)

## Scan series (batches of scans)

Series hold multiple scans plus shared metadata. For series-level getters, the `raw*`
functions return unitless numeric arrays, while the non-raw variants preserve or convert
units.

### AbstractScanSeries and subtypes

The abstract series types define the shared metadata and scanning interface.

**AbstractScanSeries** ([`AbstractScanSeries`](@ref)) expects:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `scans` | `AbstractVector{<:AbstractScan}` | [`scans`](@ref) |
| `instrument` | `NamedTuple` | [`instrument`](@ref) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref) |
| `user` | `NamedTuple` | [`user`](@ref) |
| `sample` | `NamedTuple` | [`sample`](@ref) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref) |

**AbstractChromScanSeries** ([`AbstractChromScanSeries`](@ref)) and
**AbstractMassScanSeries** ([`AbstractMassScanSeries`](@ref)) specialize the scan type and
enable series-level convenience getters such as [`retentions`](@ref),
[`intensities`](@ref), [`levels`](@ref), and [`uniquemzvalues`](@ref).

Concrete types such as `ChromScanSeries` and `MassScanSeries` are subtypes that implement
these abstract interfaces.

### ChromScanSeries

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `scans` | `AbstractVector{<:ChromScan}` | [`scans`](@ref) |
| `instrument` | `NamedTuple` | [`instrument`](@ref) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref) |
| `user` | `NamedTuple` | [`user`](@ref) |
| `sample` | `NamedTuple` | [`sample`](@ref) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref) |

Convenience:

- [`retentions`](@ref), [`intensities`](@ref), [`scancount`](@ref)

### MassScanSeries

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `scans` | `AbstractVector{<:MassScan}` | [`scans`](@ref) |
| `instrument` | `NamedTuple` | [`instrument`](@ref) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref) |
| `user` | `NamedTuple` | [`user`](@ref) |
| `sample` | `NamedTuple` | [`sample`](@ref) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref) |

Convenience:

- [`retentions`](@ref), [`intensities`](@ref), [`mzvalues`](@ref), [`mzunit`](@ref),
  [`levels`](@ref), [`uniquemzvalues`](@ref), [`scancount`](@ref)

## Mass-scan matrices (aligned grid)

Mass-scan matrices store aligned mass scans on a shared m/z grid.
Use series when you need scan-level structure or heterogeneous m/z grids; use matrices
for aligned m/z operations and linear algebra workflows.

### AbstractMassScanMatrix

The abstract matrix type defines the aligned grid interface for mass scans.

**AbstractMassScanMatrix** ([`AbstractMassScanMatrix`](@ref)) expects:

| Expected field | Type | Getter |
| :--- | :--- | :--- |
| `retentions` | `AbstractVector{<:Real}` | [`rawretentions`](@ref) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref) |
| `mzvalues` | `AbstractVector{<:Real}` | [`rawmzvalues`](@ref) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref) |
| `intensities` | `AbstractMatrix{<:Real}` | [`rawintensities`](@ref) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref) |
| `level` | `Integer` | [`level`](@ref) |
| `instrument` | `NamedTuple` | [`instrument`](@ref) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref) |
| `user` | `NamedTuple` | [`user`](@ref) |
| `sample` | `NamedTuple` | [`sample`](@ref) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref) |

### MassScanMatrix

`MassScanMatrix` is a subtype of `AbstractMassScanMatrix` that implements this interface.

Fields and accessors:

| Field | Type | Getter |
| :--- | :--- | :--- |
| `retentions` | `AbstractVector{<:Real}` | [`rawretentions`](@ref) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref) |
| `mzvalues` | `AbstractVector{<:Real}` | [`rawmzvalues`](@ref) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref) |
| `intensities` | `AbstractMatrix{<:Real}` | [`rawintensities`](@ref) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref) |
| `level` | `Integer` | [`level`](@ref) |
| `instrument` | `NamedTuple` | [`instrument`](@ref) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref) |
| `user` | `NamedTuple` | [`user`](@ref) |
| `sample` | `NamedTuple` | [`sample`](@ref) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref) |

Convenience:

- [`retentions`](@ref), [`mzvalues`](@ref), [`intensities`](@ref), [`mzcount`](@ref),
  [`scancount`](@ref)

Use `mscanmatrix(mss)` to convert a `MassScanSeries` to a `MassScanMatrix`:

```julia
msm = mscanmatrix(mss)
```

See the dedicated container pages for the full API reference:
[`Scans`](Scans.md), [`ScanSeries`](ScanSeries.md), and [`ScanMatrices`](ScanMatrices.md).
