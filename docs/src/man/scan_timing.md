# Scan timing

## Theoretical background

In gas chromatography–mass spectrometry (GC–MS), a reported scan-level retention coordinate 
is not, in general, the coordinate at which all ions in that scan are measured. Each scan 
spans a finite retention interval during which ions are sampled sequentially according to 
the instrument’s acquisition scheme. Consequently, the effective retention coordinate 
associated with a specific ion depends on how that scan interval is interpreted and how ion 
sampling is scheduled within it.

Instrument software and data formats typically associate a single retention value with each 
scan. However, this value may represent the start, middle, or end of the scan interval, 
depending on vendor conventions, acquisition mode, and export settings. Within the scan, 
ions are measured in a defined order (e.g. ascending or descending mass-to-charge), and 
each ion is assigned a dwell interval whose duration may be constant across ions or vary 
according to instrument settings. As a result, the sampling of different ions occurs at 
systematically different retention coordinates within the same scan.

For many downstream analyses—such as precise peak shape reconstruction, retention-based 
alignment at the ion level, deconvolution, or modeling of chromatographic 
distortions—treating all ions in a scan as if they were measured at the same retention 
coordinate introduces avoidable systematic error. This is particularly relevant when scan 
durations are non-negligible relative to chromatographic peak widths or when dwell 
distributions are heterogeneous.

To obtain an ion-specific retention coordinate, it is therefore necessary to combine:
the reference point associated with the scan-level retention value,

the total scan interval,

the dwell allocation across ions,

the acquisition order of ions within the scan, and

the desired reference point within each ion’s dwell interval.

The function provided here formalizes this mapping. Given a scan-level retention coordinate 
and a description of the within-scan sampling scheme, it computes the effective retention 
coordinate at which a specific ion is sampled. This enables consistent, reproducible 
ion-level retention calculations across different acquisition configurations and data 
representations, without assuming that retention is inherently tied to physical time or to 
a particular unit system.


## Example

```@example 1
# Load JuChrom and the Agilent ChemStation MS loader
using JuChrom
using JuChrom.ChemStationMSLoader

# Load an example Agilent ChemStation GC-MS run
file = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
mss = load(ChemStationMS(file; mode=:ms))
```

## Scan timing tools

```@docs
JuChrom.mzretention
```
