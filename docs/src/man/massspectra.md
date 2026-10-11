# Mass spectra

## AbstractMassSpectrum

```@docs
JuChrom.AbstractMassSpectrum
JuChrom.attrs(::JuChrom.AbstractMassSpectrum)
JuChrom.intensityunit(::JuChrom.AbstractMassSpectrum)
JuChrom.intensities(::JuChrom.AbstractMassSpectrum{<:Any, Nothing})
JuChrom.mzcount(::JuChrom.AbstractMassSpectrum)
JuChrom.mzunit(::JuChrom.AbstractMassSpectrum)
JuChrom.mzvalues(::JuChrom.AbstractMassSpectrum{Nothing, <:Any})
```

## MassSpectrum

```@docs
JuChrom.MassSpectrum
```

## Extraction from scans

```@docs
JuChrom.extractmassspectrum
```

## Plotting spectra

Load a Makie backend, such as CairoMakie or GLMakie, to enable spectrum plotting.

```@docs
JuChrom.massspectrum
JuChrom.massspectrum!
```

For an interactive view with a fixed intensity range, disable the recipe's
automatic Y limits. Labels can also be disabled to avoid their placement work
during frequent updates:

```julia
using GLMakie, JuChrom
using GLMakie: Makie

fig = Figure()
ax = Axis(fig[1, 1]; limits=(40, 300, 0, 100))
plot = massspectrum!(ax, [50.0, 100.0], [100.0, 20.0];
    autoylimits=false, show_mz_labels=false)

ylims!(ax, 0, 5)
Makie.update!(plot; arg1=[60.0, 90.0, 150.0], arg2=[2.0, 100.0, 1.0])
# The manual 0–5 Y range is retained.
```
