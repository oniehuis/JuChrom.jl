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

## Transformation tools

```@docs
clr
whiten
```
