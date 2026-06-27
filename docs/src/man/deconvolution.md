# Deconvolution

## PARAFAC2 retention metadata

`parafac2` stores retention metadata per sample. If retention metadata are available,
`rawretentions(fit)` and `retentions(fit)` return one retention vector for each sample,
even when the input was a shared retention grid.

Use a three-dimensional tensor only when the data are represented as a rectangular
`samples x retentions x m/z` block. In that case, a single `retentions` vector is treated
as the common retention grid for all samples:

```julia
fit = parafac2(
    Xtensor,
    ncomponents;
    retentions=retention_grid,
    mzvalues=mz_axis,
)
```

Use a vector of sample matrices when samples have separate retention axes. Each matrix must
be `retentions x m/z`, and the retention metadata should be passed as one vector per
sample:

```julia
fit = parafac2(
    Xmatrices,
    ncomponents;
    retentions=retention_axes,
    mzvalues=mz_axis,
)
```

For `MassScanMatrix` inputs there is currently no dedicated `parafac2` overload. Extract
the matrix data and metadata explicitly:

```julia
Xmatrices = rawintensities.(msms)
retention_axes = retentions.(msms)
mz_axis = mzvalues(first(msms))

fit = parafac2(
    Xmatrices,
    ncomponents;
    retentions=retention_axes,
    mzvalues=mz_axis,
)
```

Before doing this, make sure all samples use compatible m/z axes. The retention axes do not
need to be identical when you pass `retentions` as a vector of vectors.

```@docs
JuChrom.Parafac2Fit
JuChrom.parafac2
JuChrom.parafac2areas
JuChrom.parafac2apexes
JuChrom.parafac2fitpercent
JuChrom.parafac2intensities
JuChrom.parafac2loss
JuChrom.parafac2profilediagnostics
JuChrom.parafac2profileminima
JuChrom.parafac2reconstruct
JuChrom.parafac2residuals
JuChrom.parafac2scores
JuChrom.parafac2spectra
JuChrom.unimodalfit
JuChrom.unimodalfit_apexsearch
```
