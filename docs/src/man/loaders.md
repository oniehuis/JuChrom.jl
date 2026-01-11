# Loaders

JuChrom ships loader APIs for several chromatographic file formats. Some vendor readers rely
on reverse-engineered formats, so a given loader may not work with every file or file-set
version.

The extension system lets JuChrom offer multiple reader implementations side-by-side, so
different reader versions can be used for different file format versions.

## Loader extensions and specs

Loaders live in extension modules so optional dependencies stay opt-in. When an extension is 
loaded, its loader types and `load` methods become available without changing the core API.

Each module defines a loader *spec* type that encodes the reader version in its type 
parameter `F` (for example, `MassHunterMSLoaderSpec{F}` where `F` is a concrete format type 
such as `MassHunterMSv1`). The `load` method dispatches on that spec type, so multiple 
reader versions can coexist and be selected by constructor helpers in the module. This is 
how JuChrom can support different file versions with different reader implementations.

## Loading a file

Loaders are regular Julia types with a `load` method. The common pattern is to `using` the 
loader module, construct a loader spec, and call `load`.

```julia
using JuChrom
using JuChrom.ChemStationMSLoader

file = "/path/to/run.D/data.ms"
mss = load(ChemStationMS(file; mode=:ms))  # convenience constructor
```

Convenience constructors like `ChemStationMS(...)` and `MassHunterMS(...)` are thin wrappers 
that build the corresponding loader spec. To select a specific reader version, construct the 
spec type explicitly:

```julia
using JuChrom
using JuChrom.MassHunterMSLoader

spec = MassHunterMSLoaderSpec{MassHunterMSv1}(file; mode=:ms)
mss = load(spec)
```

## Agilent FID

```@docs
JuChrom.AgilentFIDLoader.AgilentFIDLoaderSpec
JuChrom.AgilentFIDLoader.AgilentFID
JuChrom.AgilentFIDLoader.load
```

## Aligent ChemStation MS

```@docs
JuChrom.ChemStationMSLoader.ChemStationMS
JuChrom.ChemStationMSLoader.ChemStationMSLoaderSpec
JuChrom.ChemStationMSLoader.load(::JuChrom.ChemStationMSLoader.ChemStationMSLoaderSpec{JuChrom.ChemStationMSLoader.ChemStationMSv2})
```

## Aligent MassHunter MS

```@docs
JuChrom.MassHunterMSLoader.MassHunterMSLoaderSpec
JuChrom.MassHunterMSLoader.MassHunterMS
JuChrom.MassHunterMSLoader.load(::JuChrom.MassHunterMSLoader.MassHunterMSLoaderSpec{JuChrom.MassHunterMSLoader.MassHunterMSv1})
```

## Shimadzu MS

The Shimadzu MS loader uses the PyCall extension and the Python `olefile` module. If you
do not see the Shimadzu docs below, ensure PyCall and the the Python `olefile` module are 
installed and the extension is loaded in your session.

To enable the loader in your session:

```julia
using JuChrom
using PyCall # triggers the extension
```

To install PyCall (and select a Python if needed):

```julia
import Pkg
Pkg.add("PyCall")
# Optional: point PyCall to a specific Python before rebuilding
# ENV["PYTHON"] = "/path/to/python"
Pkg.build("PyCall")
```

To install the Python dependency:

```julia
using PyCall
pyimport_conda("olefile", "olefile")
```

You can verify availability with `isdefined(JuChrom, :ShimadzuMSLoader)`.

```@docs
JuChrom.ShimadzuMSLoader.ShimadzuMSLoaderSpec
JuChrom.ShimadzuMSLoader.ShimadzuMS
JuChrom.ShimadzuMSLoader.load(::JuChrom.ShimadzuMSLoader.ShimadzuMSLoaderSpec{JuChrom.ShimadzuMSLoader.ShimadzuMSv1})
```
