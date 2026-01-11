# Shimadzu MS

The Shimadzu MS loader uses the PyCall extension and the Python `olefile` module. If you
cannot see the Shimadzu docs below, ensure PyCall and the Python `olefile` module are
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
