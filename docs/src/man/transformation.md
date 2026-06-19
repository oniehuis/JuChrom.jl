# Transformation

The transformations provided here focus on intensity normalization for downstream analysis. 
In JuChrom, the centered log-ratio (CLR) transform converts strictly positive, compositional 
data from the simplex into log-ratio coordinates in Euclidean space, which many multivariate 
projection methods assume (for example PCA, PLS, and CPPLS-DA). Whitening rescales values by 
their variances to make the data homoscedastic to first order, so features with larger 
uncertainty do not dominate distances, correlations, or model fits. Apply these steps when 
downstream analyses assume or benefit from Euclidean geometry or homoscedasticity (for 
example PCA, PLS, or CPPLS-DA).

## Example

```@example 1
using JuChrom

x = [0.1, 0.7, 0.2]
σ² = [1, 9, 3]

x_clr, σ²_clr = clr(x, σ²)
x_whitened = whiten(x_clr, σ²_clr, 0.01)
```

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
clr
dwellnormalize
withintensityunit
whiten
```
