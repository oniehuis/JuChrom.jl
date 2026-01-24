# Scan timing

## Background

In gas chromatography–mass spectrometry (GC–MS), a scan-level retention coordinate is not,
in general, the coordinate at which every ion in that scan is acquired. Each scan spans a
finite retention interval during which ions are sampled sequentially according to the
instrument’s acquisition scheme. Consequently, the ion-specific retention coordinate
depends on how the scan interval is referenced and how sampling is scheduled within it.

Instrument software and data formats typically associate a single retention value with each
scan. However, that value may represent the start, midpoint, or end of the scan interval,
depending on vendor conventions, acquisition mode, and export settings. Within a scan, ions
are acquired in a defined order (e.g., ascending or descending m/z), and each ion is
assigned a dwell interval whose duration may be constant across ions or vary with
instrument settings. As a result, different ions are sampled at systematically different
retention coordinates within the same scan.

For many downstream analyses—such as accurate peak-shape reconstruction and deconvolution—
treating all ions in a scan as if they were measured at the same retention coordinate
introduces avoidable systematic error. This is particularly relevant when scan durations
are non-negligible relative to chromatographic peak widths or when dwell times are
heterogeneous.

To obtain an ion-specific retention coordinate, it is therefore necessary to combine:
- the reference point associated with the scan-level retention value,
- the total scan interval,
- the dwell allocation across ions and the acquisition order within the scan, and
- the desired reference point within each ion’s dwell interval.

The JuChrom function [`mzretention`](@ref) formalizes this mapping. Given a scan-level
retention coordinate and a description of the within-scan sampling scheme, it computes the
effective retention coordinate at which a specific ion is sampled. This enables consistent,
reproducible ion-level retention calculations across different acquisition configurations
and data representations.

## Example

```@example 1
# Load JuChrom, the plotting backend, and the Agilent ChemStation MS loader
using JuChrom
using JuChrom.ChemStationMSLoader

# Load CairoMakie for plotting and Statistics for calculating the mean
using CairoMakie, Statistics

# Load an example Agilent ChemStation GC-MS run
file = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
mss = load(ChemStationMS(file; mode=:ms))

# Trim scans to the RT interval (20.25–20.40 minutes) containing the n‑nonacosane peak
retentiontrim!(mss, start=20.25u"minute", stop=20.4u"minute")

# Integer-bin m/z values and convert to a MassScanMatrix
msm = mscanmatrix(binmzvalues(mss, validmzvalues=29:562))

# Plot the chromatograms of m/z 57, 239, and 408 using scan-level retentions
fig₁ = Figure(; size=(1000, 350))
ax₁ = Axis(fig₁[1,1], title="Selected ion chromatograms using scan-level retentions",
                      ylabel="Intensity [no unit]",
                      xlabel="Retention [minute]")

for mz in [57, 239, 408]
    i = findfirst(==(mz), mzvalues(msm))        # m/z index in the binned list
    ints = vec(rawintensities(msm)[:, i])       # intensity trace for that m/z
    lines!(ax₁,
           rawretentions(msm, unit=u"minute"),
           ints / sum(ints),                    # normalize to a comparable scale
           label = "m/z = $(mz)")
end
axislegend(ax₁; position=:rt)
save("xic.svg", fig₁)
nothing
```

![](xic.svg)

From the three traces, the peaks at higher m/z appear shifted to slightly later retention
times, consistent with a **descending** m/z acquisition order. While Agilent GC–MS
instruments are known to acquire in descending order, plotting prominent ions with widely
separated m/z values allows the acquisition order to be inferred empirically when it is not
known *a priori*. We next use this information to correct for the intra-scan time shift and
align the chromatographic traces.

```@example 1
# Plot the chromatograms of m/z 57, 239, and 408 using shifted retentions
fig₂ = Figure(; size=(1000, 350))
ax₂ = Axis(fig₂[1,1], title="Selected ion chromatograms using shifted retentions",
                      ylabel="Intensity [no unit]",
                      xlabel="Retention [minute]")

rawrts = rawretentions(msm, unit=u"minute")
rawscanduration = mean(diff(rawrts))
for mz in [57, 239, 408]
    i = findfirst(==(mz), mzvalues(msm))
    rts_corrected = mzretention.(rawrts;
                                 mzindex=i,
                                 mzcount=mzcount(msm),
                                 dwell=:homogeneous,
                                 order=:descending,
                                 scan_interval=rawscanduration,
                                 retention_ref=:start,
                                 dwell_ref=:middle)

    ints = vec(rawintensities(msm)[:, i])
    lines!(ax₂,
           rts_corrected,
           ints / sum(ints),
           label = "m/z = $(mz)")
end
axislegend(ax₂; position=:rt)
save("xic_shifted.svg", fig₂)
nothing
```

As we can see, the three traces are now much better aligned and jointly provide a much
better estimate of the peak shape than any individual trace alone. Properly aligned ion
traces are therefore a prerequisite for accurate peak‑shape reconstruction and
[Deconvolution](deconvolution.md).

![](xic_shifted.svg)

## Scan timing function

```@docs
JuChrom.mzretention
```
