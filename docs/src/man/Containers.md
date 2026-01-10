# Core Containers

JuChrom centers on three container layers for chromatographic and mass-spectrometric data.
Single scans capture one measurement, scan series group multiple scans, and mass-scan
matrices store aligned mass scans on a shared m/z grid.

## Single-scan containers

Single scans store one chromatographic point (`ChromScan`) or one mass-spectrometric scan
(`MassScan`).

### ChromScan

Fields and accessors:

| Field | Type | Getter |
| --- | --- | --- |
| `retention` | `Real` | [`retention`](@ref) / [`rawretention`](@ref) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref) |
| `intensity` | `Real` | [`intensity`](@ref) / [`rawintensity`](@ref) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref) |
| `attrs` | `NamedTuple` | [`attrs`](@ref) |

Convenience:

- [`retentionunit`](@ref), [`intensityunit`](@ref)

### MassScan

Fields and accessors:

| Field | Type | Getter |
| --- | --- | --- |
| `retention` | `Real` | [`retention`](@ref) / [`rawretention`](@ref) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref) |
| `mzvalues` | `AbstractVector{<:Real}` | [`mzvalues`](@ref) / [`rawmzvalues`](@ref) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref) |
| `intensities` | `AbstractVector{<:Real}` | [`intensities`](@ref) / [`rawintensities`](@ref) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref) |
| `level` | `Integer` | [`level`](@ref) |
| `attrs` | `NamedTuple` | [`attrs`](@ref) |

Convenience:

- [`mzcount`](@ref)

## Scan series (batches of scans)

Series hold multiple scans plus shared metadata.

### ChromScanSeries

Fields and accessors:

| Field | Type | Getter |
| --- | --- | --- |
| `scans` | `AbstractVector{<:ChromScan}` | [`scans`](@ref), [`scancount`](@ref) |
| `instrument` | `NamedTuple` | [`instrument`](@ref) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref) |
| `user` | `NamedTuple` | [`user`](@ref) |
| `sample` | `NamedTuple` | [`sample`](@ref) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref) |

Convenience:

- [`retentions`](@ref), [`intensities`](@ref)

### MassScanSeries

Fields and accessors:

| Field | Type | Getter |
| --- | --- | --- |
| `scans` | `AbstractVector{<:MassScan}` | [`scans`](@ref), [`scancount`](@ref) |
| `instrument` | `NamedTuple` | [`instrument`](@ref) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref) |
| `user` | `NamedTuple` | [`user`](@ref) |
| `sample` | `NamedTuple` | [`sample`](@ref) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref) |

Convenience:

- [`retentions`](@ref), [`intensities`](@ref), [`mzvalues`](@ref), [`mzunit`](@ref),
  [`levels`](@ref), [`uniquemzvalues`](@ref)

## Mass-scan matrices (aligned grid)

Mass-scan matrices store aligned mass scans on a shared m/z grid.

Conversion from series:

```julia
msm = mscanmatrix(mss)
```

### MassScanMatrix

Fields and accessors:

| Field | Type | Getter |
| --- | --- | --- |
| `retentions` | `AbstractVector{<:Real}` | [`retentions`](@ref) / [`rawretentions`](@ref) |
| `retentionunit` | `Union{Unitful.Units, Nothing}` | [`retentionunit`](@ref) |
| `mzvalues` | `AbstractVector{<:Real}` | [`mzvalues`](@ref) / [`rawmzvalues`](@ref) |
| `mzunit` | `Union{Unitful.Units, Nothing}` | [`mzunit`](@ref) |
| `intensities` | `AbstractMatrix{<:Real}` | [`intensities`](@ref) / [`rawintensities`](@ref) |
| `intensityunit` | `Union{Unitful.Units, Nothing}` | [`intensityunit`](@ref) |
| `level` | `Integer` | [`level`](@ref) |
| `instrument` | `NamedTuple` | [`instrument`](@ref) |
| `acquisition` | `NamedTuple` | [`acquisition`](@ref) |
| `user` | `NamedTuple` | [`user`](@ref) |
| `sample` | `NamedTuple` | [`sample`](@ref) |
| `extras` | `Dict{String, Any}` | [`extras`](@ref) |

Convenience:

- [`mzcount`](@ref), [`scancount`](@ref)

See the dedicated container pages for the full API reference:
[`Scans`](Scans.md), [`ScanSeries`](ScanSeries.md), and [`ScanMatrices`](ScanMatrices.md).
