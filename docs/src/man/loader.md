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
    bt = Char(96)
    docs_block = string(
        bt, bt, bt, "@docs\n",
        "JuChrom.ShimadzuMSLoader.ShimadzuMSLoaderSpec\n",
        "JuChrom.ShimadzuMSLoader.ShimadzuMS\n",
        bt, bt, bt
    )
    Markdown.parse(docs_block)
end
```
