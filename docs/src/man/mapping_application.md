# Applying maps

Use [`applymap`](@ref) to transform values from retention domain A to domain B and
[`invmap`](@ref) for the inverse transformation. The mapping functions accept scalar values
and can be broadcast over vectors with Julia's dot syntax. For a complete example that uses
one mapper from fitting through plotting and application, see the
[retention mapping workflow](mapping_workflow.md).

!!! warning "Domain limits and extrapolation"
    Mappings are defined over the anchor domain. Values outside that domain are linearly 
    extrapolated using the slope at the nearest boundary. Use `warn=true` on 
    [`applymap`](@ref)
    /
    [`invmap`](@ref)
    /
    [`derivmap`](@ref)
    /
    [`derivinvmap`](@ref)
    or the `raw*` variants to surface extrapolation during
    analysis.

For `MassScanMatrix` and `VarianceMassScanMatrix`, [`applymap`](@ref) applies the
Jacobian correction and updates intensity and variance units automatically. For a vector
of intensities, use the derivative of the mapping directly. If `ri = f(t)`, then
`d(ri)/dt` is the local stretch factor; to preserve area, you divide the intensity by this
slope at the same time point. However, even when intensities are reported as unitless,
they usually represent counts accumulated over a finite scan interval, so the implied unit
is typically `time^-1`. Here we assume the given intensities are counts acquired over a
0.5‑second scan interval, and therefore divide by 0.5 seconds.
Because the intensities and `d(ri)/dt` must use the same time unit for the division to
cancel cleanly, we convert the Jacobian to the unit `s^-1`. This is done by specifying the
input‑domain unit as `u"s"`.

The transformed intensities are now expressed per unit of retention index rather than per
unit time, reflecting the change of variables.

!!! note "Use ion-specific retention times for sequential MS data"
    The Jacobian must be evaluated at the retention coordinate of the measurement whose
    intensity is transformed. For a chromatographic trace with one intensity per scan, the
    scan retention time is the appropriate coordinate. For sequentially scanned MS data,
    however, different m/z channels in the same scan are measured at slightly different
    times. In that case, a single scan-level Jacobian applies the same correction to all
    ions in a scan and is only an approximation.

    For the most accurate RT -> RI transformation of ion traces or MS scan matrices, first
    correct the scan retention time to the acquisition time of the specific m/z channel,
    then evaluate the Jacobian at that ion-specific retention time:

    ```julia
    rt_ion = mzretention(scan_retention, mz; order=:descending, ...)
    ri_ion = applymap(mapper, rt_ion)
    jacobian = derivmap(mapper, rt_ion)
    intensity_per_ri = intensity_per_rt / jacobian
    ```

    This distinction is mainly relevant for sequential quadrupole-like acquisition. For
    simultaneous m/z acquisition, or when the m/z scan duration is negligible relative to
    chromatographic peak widths and mapper curvature, the scan-level Jacobian can be a
    reasonable approximation.

## Unit-Aware Mapping

```@docs
JuChrom.applymap
JuChrom.invmap
JuChrom.derivmap
JuChrom.derivinvmap
```

## Raw Numeric Mapping

```@docs
JuChrom.rawapplymap
JuChrom.rawinvmap
JuChrom.rawderivmap
JuChrom.rawderivinvmap
```
