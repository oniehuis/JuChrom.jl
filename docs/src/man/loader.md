# Loader

```@docs
JuChrom.AgilentFIDLoader.AgilentFIDLoaderSpec
JuChrom.AgilentFIDLoader.AgilentFID
JuChrom.ChemStationMSLoader.ChemStationMS
JuChrom.ChemStationMSLoader.ChemStationMSLoaderSpec
JuChrom.ChemStationMSLoader.load
JuChrom.MassHunterMSLoader.MassHunterMSLoaderSpec
JuChrom.MassHunterMSLoader.MassHunterMS
```

Shimadzu loaders are available when PyCall is installed and loaded (extension).

```@eval
using JuChrom
if isdefined(JuChrom, :ShimadzuMSLoader)
    :( @docs JuChrom.ShimadzuMSLoader.ShimadzuMSLoaderSpec JuChrom.ShimadzuMSLoader.ShimadzuMS )
end
```
