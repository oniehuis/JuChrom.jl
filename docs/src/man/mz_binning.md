# m/z binning

## Rationale

Mass spectra are measured at high precision, with centroiding yielding continuous m/z
values rather than discrete integer bins. In many instruments, centroiding is performed
on the fly during acquisition: the detector’s continuous profile signal is reduced to
peak centers and intensities before storage. As a result, scans are lists of peak centers
with non-uniform spacing in the m/z domain, not evenly sampled vectors. Many downstream
workflows, however, operate on a fixed m/z grid. Binning consolidates nearby m/z values
into shared bins, which reduces spectrum sparsity, stabilizes signal across scans, and
makes spectra comparable for alignment, similarity scoring, and matrix-based processing.
It also provides a simple, deterministic way to handle small instrument-to-instrument or
run-to-run shifts without overfitting to noise.

For GC-MS data, integer binning with boundaries at -0.3 (lower) and +0.7 (upper) around
each integer m/z is meaningful because fragment masses cluster around nominal integers
once you account for elemental mass defects (O'Callaghan et al. 2012). If $w(X)$ is the 
atomic weight of element $X$ and $\Delta(X) = w(X) - \mathrm{round}(w(X))$, then 
$\Delta(^{12}\mathrm{C}) = 0$, $\Delta(^{14}\mathrm{N}) \approx 0.00022$, 
$\Delta(^{16}\mathrm{O}) \approx -0.00032$, and $\Delta(^{1}\mathrm{H}) \approx 0.00783$. 
For a fragment $Y$ with composition $k_1 X_1 + \cdots + k_r X_r$, the mass defect is
$\Delta(Y) = k_1 \Delta(X_1) + \cdots + k_r \Delta(X_r)$. Typical GC-MS fragments 
(m/z ≤ 550) cannot accumulate large negative defects because they rarely contain many 
P/Si/O atoms, and those elements are the main contributors to negative mass defect. They 
also cannot accumulate large positive defects because the hydrogen count is limited, so the
feasible defect range is roughly bounded between about -0.13 and +0.63 for fragments in
the GC-MS mass range. A bin that spans $[n-0.3, n+0.7)$ around each integer $n$ thus 
comfortably contains the physically plausible defect window for a nominal mass, reducing 
the chance that a real fragment sits near a bin edge.

## JuChrom m/z binning

JuChrom provides [`binmzvalues`](@ref), which applies a binning function to the m/z values
in each scan of a `MassScanSeries` and sums intensities within each bin to produce a new
`MassScanSeries` on a consistent grid. By default it uses [`integer`](@ref), which 
implements the $[n-0.3, n+0.7)$ integer binning convention described above, but any custom 
binning function can be supplied to target a different bin width or precision (for example, 
rounding to a fixed decimal grid).

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
```

## References

- O'Callaghan S, De Souza DP, Isaac A, Wang Q, Hodkinson L, Olshansky M, Erwin T, Appelbe B, 
Tull DL, Roessner U, Bacic A, McConville MJ, Likic VA (2012): PyMS: a Python toolkit for 
processing of gas chromatography-mass spectrometry (GC-MS) data. Application and 
comparative study of selected tools. BMC Bioinformatics 13: 115.

## m/z binning tools

```@docs
binmzvalues(::MassScanSeries, ::Function)
integer
```
