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

Mass spectrometry intensities may be stored as signal accumulated over dwell intervals.
Use [`withintensityunit`](@ref) when the numeric intensity values are known to have a unit
that is not already recorded in the container. Raw or vendor-scaled intensity matrices are
normally left unitless until they are dwell-normalized:

```@example 2
using JuChrom

msm = MassScanMatrix(
    [1.0, 2.0]u"s",
    [100.0, 200.0, 300.0],
    [10.0 20.0 30.0; 40.0 50.0 60.0],
)

intensityunit(msm)
```

Use [`dwellnormalize`](@ref) to divide intensities by m/z-specific dwell intervals. The
function accepts inputs whose intensity unit is absent and returns a rate-like signal with
the reciprocal dwell unit. For explicit dwell intervals, their sum must not exceed the
shortest scan interval in the retention axis; this explicit-vector form is therefore a
sequential-acquisition model.

```@example 2
normalized = dwellnormalize(msm, [0.2, 0.3, 0.5], u"s")
intensityunit(normalized)
```

```@example 2
rawintensities(normalized)
```

For simultaneous acquisition with a known effective dwell time, pass a scalar dwell
interval. The scalar dwell must not exceed the shortest scan interval.

```@example 2
simultaneous = dwellnormalize(msm, 0.5, u"s")
rawintensities(simultaneous)
```

When the scan interval is uniformly divided among all m/z values, `dwellnormalize(msm)`
infers a single dwell interval from the shortest scan interval and the number of m/z
values. For simultaneous acquisition, use `dwellnormalize(msm; acquisition=:simultaneous)`
to infer the full shortest scan interval as the dwell time for every m/z value.

## Transformation tools

```@docs
CensoredReplacementInfo
clr
dwellnormalize
replacecensored
withintensityunit
whiten
```
