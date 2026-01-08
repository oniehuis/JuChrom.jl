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
using Markdown
if isdefined(JuChrom, :ShimadzuMSLoader)
    Markdown.MD(Markdown.CodeBlock("@docs", """
JuChrom.ShimadzuMSLoader.ShimadzuMSLoaderSpec
JuChrom.ShimadzuMSLoader.ShimadzuMS
"""))
end
```
