# Transformation

The transformations provided here focus on intensity normalization for downstream analysis. 
In JuChrom, the centered log-ratio (CLR) transform converts strictly positive, compositional 
data from the simplex into log-ratio coordinates in Euclidean space, which many multivariate 
projection methods assume (for example PCA, PLS, and CPPLS-DA). Whitening rescales values by 
their variances to make the data homoscedastic to first order, so features with larger 
uncertainty do not dominate distances, correlations, or model fits. Apply these steps when 
downstream analyses assume or benefit from Euclidean geometry or homoscedasticity (for 
example PCA, PLS, or CPPLS-DA).

Use [`replacecensored`](@ref) before CLR when a variance-carrying mass-scan matrix contains
non-positive or below-limit intensities. It replaces censored values by positive
below-limit estimates and updates their variances on the same final analysis grid.

[`ticnormalize`](@ref) scales a mass-scan matrix by the sum of positive
total-ion-chromatogram (TIC) values. For a `VarianceMassScanMatrix`, variances are
propagated as fixed scaling (`V ./ scale^2`). The normalization denominator is treated as
the observed scaling factor; uncertainty in that denominator is not propagated because it
would introduce covariances between all normalized cells, while `VarianceMassScanMatrix`
stores only per-cell variances.

CLR is invariant to this global TIC scaling. For strictly positive data,
`clr(ticnormalize(vmsm))` gives the same CLR intensities and propagated CLR variances as
`clr(vmsm)`. If the matrix contains zeros or negative values and you want CLR output, first
replace censored values:

```julia
clr(ticnormalize(replacecensored(vmsm)))
```

This gives the same CLR result as `clr(replacecensored(vmsm))`, with the TIC-normalized
intermediate available if you need it for inspection or other analysis.

## Example

```@example 1
using JuChrom

msm = MassScanMatrix(
    [1.0, 2.0]u"s",
    [100.0, 101.0],
    [1.0 2.0; 4.0 8.0]u"pA",
)
vmsm = VarianceMassScanMatrix(msm, [0.01 0.04; 0.16 0.64] .* u"pA"^2)

vmsm_clr = clr(vmsm)
vmsm_whitened = whiten(vmsm_clr)
```

`whiten` returns dimensionless intensities and variances. By default it infers a
standard-deviation floor from the positive propagated standard deviations
(`sigmafloor = :auto`, `floorquantile = 0.05`). Pass an explicit bare number for
dimensionless transformed data, for example `whiten(vmsm_clr; sigmafloor=0.01)`. Pass a
unitful floor for physical intensity data, for example
`whiten(vmsm; sigmafloor=0.1u"pA")`; the intensities and variances are converted to that
unit before units are stripped.

## Intensity units and dwell normalization

Use [`withintensityunit`](@ref) when the numeric intensity values already have a known
unit that is absent from the container. It attaches that unit without rescaling values.
For example, `withintensityunit(msm, u"s^-1")` labels an already normalized rate.

Use [`dwellnormalize`](@ref) when stored intensity is proportional to the duration of
measurement. Having no intensity unit does not by itself establish this relationship.
The function requires an explicit unitful duration, or explicit m/z-to-duration pairs;
acquisition-method interpretation and estimation of dwell are separate operations.

```@example 2
using JuChrom

msm = MassScanMatrix(
    [1.0, 2.0]u"s",
    [100.0, 200.0, 300.0],
    [10.0 20.0 30.0; 40.0 50.0 60.0],
)

# Known dwell duration for each ion; pair order need not follow column order.
dwell = (300 => 500u"ms", 100 => 0.2u"s", 200 => 300u"ms")
normalized = dwellnormalize(msm, dwell)
rawintensities(normalized)
```

Tuples and vectors of pairs, and dictionaries, are accepted. Every matrix m/z value
needs a duration, including columns containing only zeros. Additional ions are allowed,
so the same acquisition specification can be reused after selecting fewer columns.
Duplicate m/z keys and missing durations are errors. All entries are validated, including
ions absent from the matrix. Durations must be finite, positive quantities with time units.

Keys follow the same unit conventions as [`mscanmatrix`](@ref): plain keys use the matrix's
m/z unit, if any; unitful keys are converted to that unit before exact matching. Mixed
plain/unitful keys, incompatible units, and unitful keys for a unitless m/z axis are
rejected. Keys must be finite and positive. Compatible units may differ between keys;
duplicates are checked after conversion. There is no rounding or approximate matching.

For a known duration shared by all columns, use a scalar. This form makes no assumption
about whether acquisition was sequential or simultaneous:

```@example 2
shared = dwellnormalize(msm, 500u"ms")
rawintensities(shared)
```

All durations are converted to seconds. Results always have intensity unit `s^-1`,
independently of input time units, pair order, or additional ions. Input intensities must
have no attached unit. Both dense and sparse matrices are supported; storage is preserved.
Axes, MS level, and metadata are copied.

```@example 2
intensityunit(normalized)
```

For [`VarianceMassScanMatrix`](@ref), the same API divides variances by the squared dwell
durations, yielding variance unit `s^-2`:

```@example 2
vmsm = VarianceMassScanMatrix(msm, [1.0 4.0 9.0; 16.0 25.0 36.0])
normalized_variances = dwellnormalize(vmsm, Dict(dwell))
rawvariances(normalized_variances)
```

There is no argument-free inference, positional duration-vector form, or `acquisition`
keyword. Neither recorded ions nor adjacent retention spacing determines the normalization.
Durations are not checked against retention spacing: explicit dwell works for single
scans, selected scans, and transformed or unitless retention coordinates. Establishing
whether a duration correctly describes the acquisition is the caller's responsibility.

## Transformation tools

```@docs
CensoredReplacementInfo
clr
dwellnormalize
replacecensored
ticnormalize
withintensityunit
whiten
```
