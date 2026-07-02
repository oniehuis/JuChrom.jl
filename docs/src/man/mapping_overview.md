# Overview

## Mapping rationale

Retention mapping converts an irregular axis (time or distance) into a stable, comparable 
index so chromatograms from different runs can be aligned, interpolated, and compared on a 
common grid. This makes retention behavior reproducible across batches, instruments, and 
methods while preserving ordering and monotonicity.

A continuous, differentiable mapping does more than align coordinates: it enables 
intensity-aware transforms via the 
[`Jacobian`](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant). When you warp 
the axis, the mapping derivative provides a principled way to rescale intensities so areas 
and peak shapes remain physically meaningful. This yields consistent peak integrals across 
transformed domains, supports smooth interpolation, and avoids artifacts from piecewise or 
discontinuous mappings.

Mapping points are typically collected from standards with known reference positions. In 
gas chromatography, a common example is an n-alkane ladder, which yields paired arrays of 
retention times and 
[Kováts retention indices](https://en.wikipedia.org/wiki/Kovats_retention_index) 
(Kováts 1958). JuChrom can infer these points from GC-MS alkane ladder data with
[`findalkanes`](@ref) and expose the accepted mapper anchors with
[`alkaneladdercalibrationpoints`](@ref). Another approach uses a curated set of internal
standards (e.g., a few stable compounds spiked into every run), producing matched
retention pairs that anchor the mapping across batches (Skoog et al. 2007).

For a complete worked example, start with the
[retention mapping workflow](mapping_workflow.md). The workflow fits one mapper, inspects
diagnostic plots, applies the mapper, and stores it with JLD2.

## Mapping tools at a glance

| Function | Use case |
| :--- | :--- |
| [`fitmap`](@ref) | Infer mapping function from paired points |
| [`applymap`](@ref) | Map retention domain A → retention domain B |
| [`invmap`](@ref) | Map retention domain B → retention domain A |
| [`derivmap`](@ref) | Jacobian `d(ri)/dt` for intensity scaling |
| [`derivinvmap`](@ref) | Inverse mapping derivative |
| [`rawapplymap`](@ref) | Unitless variant of [`applymap`](@ref) |
| [`rawinvmap`](@ref) | Unitless variant of [`invmap`](@ref) |
| [`rawderivmap`](@ref) | Unitless variant of [`derivmap`](@ref) |
| [`rawderivinvmap`](@ref) | Unitless variant of [`derivinvmap`](@ref) |
| [`retentions_A`](@ref), [`retentions_B`](@ref) | Anchor vectors used to fit the mapper |
| [`rawretentions_A`](@ref), [`rawretentions_B`](@ref) | Unitless anchor vectors |
| [`retentionunit_A`](@ref), [`retentionunit_B`](@ref) | Stored units for the anchor domains |
| [`extras`](@ref) | Metadata attached to the mapper |
| [`mapmin`](@ref), [`mapmax`](@ref) | Numeric minimum and maximum value of the input domain |
| [`invmapmin`](@ref), [`invmapmax`](@ref) | Numeric minimum and maximum value of the output domain |
| [`rawmapmin`](@ref), [`rawmapmax`](@ref) | Unitless variant of [`mapmin`](@ref) and [`mapmax`](@ref) |
| [`rawinvmapmin`](@ref), [`rawinvmapmax`](@ref) | Unitless variant of [`invmapmin`](@ref), [`invmapmax`](@ref) |

For full API details, see the dedicated reference pages for
[fitting maps](mapping_fitting.md), [diagnostic plots](mapping_plotting.md),
[applying maps](mapping_application.md), and [mapper containers](mapping_containers.md).

## References

- Kováts E (1958): Gas-Chromatographische Charakterisierung organischer Verbindungen. Teil 1: Retentionsindices aliphatischer Halogenide,Alkohole, Aldehyde und Ketone. Helvetica Chimica Acta 41: 1915-1932.
- Skoog DA, Holler FJ, Crouch SR (2007): Principles of Instrumental Analysis. 6th ed. Thomson Brooks/Cole.
