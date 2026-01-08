module PyCallExtension

using PyCall
using JuChrom

const ROOT = joinpath(@__DIR__, "..")

function __init__()
    if !isdefined(JuChrom, :ShimadzuMSLoader)
        Base.include(JuChrom, joinpath(ROOT, "src", "IO", "Loaders", "ShimadzuMSLoader.jl"))
    end
end

end
