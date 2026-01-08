module PyCallExtension

using PyCall
using JuChrom

const ROOT = joinpath(@__DIR__, "..")
@eval JuChrom include(joinpath(ROOT, "src", "IO", "Loaders", "ShimadzuMSLoader.jl"))

end
