# deps/build.jl

using Pkg

Pkg.add("PyCall")
ENV["PYTHON"] = ""
Pkg.build("PyCall")

Pkg.add("https://github.com/rdeits/NNLS.jl.git")
Pkg.build("https://github.com/rdeits/NNLS.jl.git")
