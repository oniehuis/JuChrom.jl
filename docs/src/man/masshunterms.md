# Agilent MassHunter MS

Provide these files inside the `.D` folder's `AcqData` directory:

| File | Required when |
| --- | --- |
| `MSScan.bin` | Both `:ms` and `:tic` modes |
| `MSScan.xsd` | Both `:ms` and `:tic` modes |
| `MSPeak.bin` | `mode=:ms` |
| `MSTS.xml` | In either mode, when `MSScan.xsd` declares `SpectrumParamValues` within `ScanRecordType` with `maxOccurs="unbounded"` |

For that schema declaration, the loader uses the scan count in `MSTS.xml` to verify
that the index contains one spectrum block per scan. For other supported schemas,
`MSTS.xml` is optional; its scan count is checked if the file is present.

```julia
using JuChrom
using JuChrom.MassHunterMSLoader

mss = load(MassHunterMS("run.D"; mode=:ms))
tic = load(MassHunterMS("run.D"; mode=:tic))
ms1 = load(MassHunterMS("run.D"; level=1,
    scantimetype=Float32, iontype=Float32, intensitytype=Float32))
```

`mode=:ms` returns a `MassScanSeries`; `mode=:tic` returns a `ChromScanSeries`.
`level=nothing` includes all MS levels;
an integer selects that level. The type keywords convert retention times, m/z values,
and intensities to the requested floating-point types.

Available scan metadata can be inspected with `attrs(first(scans(mss)))`. Its fields
depend on the acquisition. Sample, acquisition, instrument, and user metadata are
currently empty.

Full-spectrum loading supports uncompressed spectra in `MSPeak.bin`. Profile data,
compressed spectra, multiple spectra per scan record, and empty spectra are currently
unsupported. For unsupported spectral encodings, `mode=:tic` may still be available.

```@docs
JuChrom.MassHunterMSLoader.MassHunterMSLoaderSpec
JuChrom.MassHunterMSLoader.MassHunterMS
JuChrom.MassHunterMSLoader.load(::JuChrom.MassHunterMSLoader.MassHunterMSLoaderSpec{JuChrom.MassHunterMSLoader.MassHunterMSv1})
```
