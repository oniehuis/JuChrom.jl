# deps/build.jl

using Pkg

Pkg.add("PyCall")
ENV["PYTHON"] = ""
Pkg.build("PyCall")
