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
reader versions can coexist and be selected by constructor helpers in the module.

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

## Loader pages

Each loader has its own page with the relevant docstrings.

- [Agilent FID](AgilentFID.md)
- [Agilent ChemStation MS](ChemStationMS.md)
- [Agilent MassHunter MS](MassHunterMS.md)
- [Shimadzu MS](ShimadzuMS.md)
