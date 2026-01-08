module PyCallExtension

using PyCall
using JuChrom

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(ROOT, "src", "IO", "Loaders", "ShimadzuMSLoader.jl"))

function __init__()
    if !isdefined(JuChrom, :ShimadzuMSLoader)
        @eval JuChrom const ShimadzuMSLoader = $(ShimadzuMSLoader)
    end
end

end
