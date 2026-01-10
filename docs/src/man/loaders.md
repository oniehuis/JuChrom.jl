# Loaders

```@docs
JuChrom.AgilentFIDLoader.AgilentFIDLoaderSpec
JuChrom.AgilentFIDLoader.AgilentFID
JuChrom.ChemStationMSLoader.ChemStationMS
JuChrom.ChemStationMSLoader.ChemStationMSLoaderSpec
JuChrom.ChemStationMSLoader.load
JuChrom.MassHunterMSLoader.MassHunterMSLoaderSpec
JuChrom.MassHunterMSLoader.MassHunterMS
JuChrom.ShimadzuMSLoader.ShimadzuMSLoaderSpec
JuChrom.ShimadzuMSLoader.ShimadzuMS
```

Shimadzu loaders are available when PyCall is installed and loaded (extension).
The extension uses the Python `olefile` module to read Shimadzu GC-MS files.
If the docs build does not load PyCall, the API docs below will not render,
but the loader is still available at runtime.

To enable the loader in your session:

```julia
using JuChrom
using PyCall # triggers the extension
```

To install PyCall (and configure which Python it uses):

```julia
import Pkg
Pkg.add("PyCall")
# Optional: point PyCall to a specific Python before rebuilding
# ENV["PYTHON"] = "/path/to/python"
Pkg.build("PyCall")
```

If needed, install the Python dependency:

```julia
using PyCall
pyimport_conda("olefile", "olefile")
```

You can verify the loader is available with:

```julia
isdefined(JuChrom, :ShimadzuMSLoader)
```
