# m/z binning

## Rationale

Mass spectra are measured at high precision, with centroiding yielding continuous m/z
values rather than discrete integer bins. In many instruments, centroiding is performed
during acquisition such that the detector’s profile signal is summarized as peak centers
and intensities before storage, and scans are therefore stored as sparse peak lists with
non-uniform spacing in the m/z domain. Many downstream workflows, however, operate on a
fixed m/z grid. Binning consolidates nearby m/z values into shared bins, which reduces
spectrum sparsity, stabilizes signal across scans, and makes spectra comparable for
alignment, similarity scoring, and matrix-based processing. This discretization is
standard in low-resolution GC–MS workflows, where spectra are treated as nominal-mass
vectors rather than accurate-mass distributions.

For unit-mass GC–MS (e.g., quadrupole electron ionization (EI) instruments), spectra are 
acquired as centroided floating-point m/z values, whereas interpretation, library matching, 
and fingerprinting operate on nominal integer mass channels. As a result, centroid masses 
must be mapped onto an integer grid by assigning each measured m/z value to a nominal mass 
bin. The physically meaningful placement of such bins is governed by the systematic offset
between exact and nominal masses arising from elemental mass defects. Let $w(X)$ be the
monoisotopic atomic mass of element $X$, and define the elemental mass defect as
$\Delta(X) = w(X) - \mathrm{round}(w(X))$. For a fragment ion $Y$ with elemental
composition $Y = \sum_{i=1}^{r} k_i X_i$, where $X_i$ denotes the $i$-th element present,
$k_i$ is the number of atoms of that element, and $r$ is the number of distinct elements
in the fragment, the total mass defect is $\Delta(Y) = \sum_{i=1}^{r} k_i \Delta(X_i)$.
Because typical GC–EI fragments are composed primarily of C, H, N, and O, with only
limited incorporation of heavier heteroatoms, their exact masses are systematically
offset from nominal integers but remain well separated from adjacent nominal channels.
As a result, mapping floating-point centroid m/z values onto integer mass bins is a
natural and physically meaningful discretization step for low-resolution GC–MS.

From a theoretical perspective, there is no unique binning interval: any asymmetric
window that captures the physically plausible mass-defect range and exceeds instrumental
centroid uncertainty is acceptable. Accordingly, multiple asymmetric binning schemes are
valid in principle. In practice, however, a window spanning $[n-0.3, n+0.7)$ around each
integer $n$ has become the most widely adopted convention in nominal-mass GC–MS
workflows. This interval provides sufficient tolerance for mass defect, centroid jitter,
and modest mass-axis drift while maintaining clear separation between adjacent nominal
channels.

## JuChrom m/z binning

JuChrom provides [`binmzvalues`](@ref), which applies a binning function to the m/z values
in each scan of a `MassScanSeries` and sums intensities within each bin to produce a new
`MassScanSeries` on a consistent grid. By default it uses [`integer`](@ref), which 
implements the $[n-0.3, n+0.7)$ integer binning convention described above. However, users 
can apply [`integer`](@ref) with a different offset or supply custom binning functions to 
target a different bin width or precision (for example, rounding to a fixed decimal grid). 
To enforce a fixed nominal-mass range across runs, you can pass `validmzvalues` to restrict 
the accepted bins; any binned values outside the provided vector are discarded. This is 
useful when occasional out-of-range measurements would otherwise introduce extra bins. Once 
the m/z values of a [`MassScanSeries`](@ref JuChrom.MassScanSeries) are on a consistent 
grid, you can convert the container to a [`MassScanMatrix`](@ref JuChrom.MassScanMatrix)
with [`mscanmatrix`](@ref). This makes the fixed m/z grid explicit and enables efficient
matrix-based workflows (dense or sparse storage, linear algebra, and downstream binning or
mapping that operate on matrices).

## Example

```@example 1
# Load JuChrom and the Agilent ChemStation MS loader
using JuChrom
using JuChrom.ChemStationMSLoader

# Load an example Agilent ChemStation GC-MS run
file = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
mss = load(ChemStationMS(file; mode=:ms))

# Unique m/z values before binning
uniquemzvalues(mss)
```

```@example 1
# Bin m/z values with the [n-0.3, n+0.7) convention
mss_mz_binned = binmzvalues(mss)  # Returns new MassScanSeries object

# Unique m/z values after binning
uniquemzvalues(mss_mz_binned)

# Optional conversion to MassScanMatrix for efficient matrix-based workflows
msm = mscanmatrix(mss_mz_binned)
```

```@example 1
# Keep only nominal masses 29–500
mss_mz_binned_fixed = binmzvalues(mss; validmzvalues=29:1:500)
```

## m/z binning tools

```@docs
binmzvalues(::MassScanSeries, ::Function)
integer
```
