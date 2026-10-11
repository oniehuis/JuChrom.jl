using Pkg
using TOML

length(ARGS) == 1 || error("Usage: julia test_explorer.jl /path/to/JuChromExplorer")
const juchrom_root = realpath(joinpath(@__DIR__, "..", ".."))
const explorer_root = realpath(only(ARGS))
const test_root = joinpath(explorer_root, "test")
const test_entry = joinpath(test_root, "runtests.jl")
isfile(test_entry) || error("JuChromExplorer's test/runtests.jl is missing")

# Use the explorer's declared test dependencies in an isolated environment. Run
# the tests directly so Pkg.test cannot reintroduce its released JuChrom pin.
project = TOML.parsefile(joinpath(test_root, "Project.toml"))
sources = get!(project, "sources", Dict{String, Any}())
for source in values(sources)
    if haskey(source, "path")
        source["path"] = abspath(test_root, source["path"])
    end
end
sources["JuChrom"] = Dict("path" => juchrom_root)
sources["JuChromExplorer"] = Dict("path" => explorer_root)

Pkg.activate(; temp=true)
empty!(LOAD_PATH)
append!(LOAD_PATH, ["@", "@stdlib"])
open(Base.active_project(), "w") do io
    TOML.print(io, project)
end
Pkg.instantiate()

# Fail rather than accidentally test the published dependency instead of the PR.
using JuChrom
using JuChromExplorer
@assert realpath(pkgdir(JuChrom)) == juchrom_root "Wrong JuChrom source loaded"
@assert realpath(pkgdir(JuChromExplorer)) == explorer_root "Wrong explorer source loaded"
@info "Downstream test sources" JuChrom=pkgdir(JuChrom) JuChromExplorer=pkgdir(JuChromExplorer)
Pkg.status()
cd(explorer_root) do
    include(test_entry)
end
