module MakieExtension

using Makie
using JuChrom

const ROOT = joinpath(@__DIR__, "..")
include(joinpath(ROOT, "ext", "makie_extensions", "massspectrum.jl"))
include(joinpath(ROOT, "ext", "makie_extensions", "quadvarfit.jl"))
include(joinpath(ROOT, "ext", "makie_extensions", "retentionmapper.jl"))

end
