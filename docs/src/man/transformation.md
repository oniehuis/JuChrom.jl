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

Mass spectrometry intensities may be stored as integrated ion counts over dwell intervals.
Use [`withintensityunit`](@ref) when the numeric intensity values are known to have a unit
that is not already recorded in the container. Raw ion-count matrices are normally left
unitless until they are dwell-normalized:

```@example 2
using JuChrom

msm = MassScanMatrix(
    [1.0, 2.0]u"s",
    [100.0, 200.0, 300.0],
    [10.0 20.0 30.0; 40.0 50.0 60.0],
)

intensityunit(msm)
```

Use [`dwellnormalize`](@ref) to divide count intensities by m/z-specific dwell intervals.
The function accepts inputs whose intensity unit is absent and returns a rate with the
reciprocal dwell unit.

```@example 2
normalized = dwellnormalize(msm, [1.0, 2.0, 5.0], u"s")
intensityunit(normalized)
```

```@example 2
rawintensities(normalized)
```

When the scan interval is uniformly divided among all m/z values, `dwellnormalize(msm)`
infers a single dwell interval from the mean scan interval and the number of m/z values.

## Transformation tools

```@docs
clr
dwellnormalize
withintensityunit
whiten
```
