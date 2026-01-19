using Pkg
using SafeTestsets

const TEST_PROJECT = normpath(joinpath(@__DIR__, "Project.toml"))

if normpath(Base.active_project()) != TEST_PROJECT
    Pkg.activate(@__DIR__)
end

# Avoid leaking packages from the user default environment (@v#.#).
if LOAD_PATH != ["@", "@stdlib"]
    empty!(LOAD_PATH)
    append!(LOAD_PATH, ["@", "@stdlib"])
end

# Ensure the local package is available when running runtests.jl directly.
if !haskey(Pkg.project().dependencies, "JuChrom")
    Pkg.develop(PackageSpec(path=normpath(joinpath(@__DIR__, ".."))))
end

Pkg.instantiate()

# Core
@safetestset "core/alignment/alignment.jl" begin
    include("./core/alignment/alignment.jl")
end
@safetestset "core/core_containers/scans.jl" begin 
    include("./core/core_containers/scans.jl")
end
@safetestset "core/core_containers/series.jl" begin
    include("./core/core_containers/series.jl")
end
@safetestset "core/core_containers/matrices.jl" begin
    include("./core/core_containers/matrices.jl")
end
@safetestset "core/core_getters/getters.jl" begin
    include("./core/core_getters/getters.jl")
end
@safetestset "core/baseline/baseline.jl" begin
    include("./core/baseline/baseline.jl")
end
@safetestset "core/convert/mscanmatrix.jl" begin
    include("./core/convert/mscanmatrix.jl")
end
@safetestset "core/convert/mzchrom.jl" begin
    include("./core/convert/mzchrom.jl")
end
@safetestset "core/utils/collections.jl" begin
    include("./core/utils/collections.jl")
end
@safetestset "core/utils/math.jl" begin
    include("./core/utils/math.jl")
end
@safetestset "core/utils/units.jl" begin
    include("./core/utils/units.jl")
end
@safetestset "core/retention_mapping/retention_mapping.jl" begin
    include("./core/retention_mapping/retention_mapping.jl")
end
@safetestset "core/transform/binning.jl" begin
    include("./core/transform/binning.jl")
end
@safetestset "core/transform/clr.jl" begin
    include("./core/transform/clr.jl")
end
@safetestset "core/transform/gridding.jl" begin
    include("./core/transform/gridding.jl")
end
@safetestset "core/transform/mscanmatrix_ops.jl" begin
    include("./core/transform/mscanmatrix_ops.jl")
end
@safetestset "core/transform/series_ops.jl" begin
    include("./core/transform/series_ops.jl")
end
@safetestset "core/transform/whitening.jl" begin
    include("./core/transform/whitening.jl")
end
@safetestset "core/quadvar_model/quadvar_model.jl" begin
    include("./core/quadvar_model/quadvar_model.jl")
end
@safetestset "core/deconvolution/unimodalfit.jl" begin
    include("./core/deconvolution/unimodalfit.jl")
end

# IO
@safetestset "IO/InputOutput" begin
    include("IO/inputoutput.jl")
end
@safetestset "IO/Loaders/AgilentFIDLoader" begin
    include("IO/Loaders/AgilentFIDLoader.jl")
end
@safetestset "IO/Loaders/ChemStationMSLoader" begin
    include("IO/Loaders/ChemStationMSLoader.jl")
end
@safetestset "IO/Loaders/MassHunterMSLoader" begin
    include("IO/Loaders/MassHunterMSLoader.jl")
end
if Base.find_package("PyCall") !== nothing
    @safetestset "IO/Loaders/ShimadzuMSLoader" begin
        include("IO/Loaders/ShimadzuMSLoader.jl")
    end
else
    @info "Skipping ShimadzuMSLoader tests because PyCall is not available."
end

# Extensions
@safetestset "ext/jld2_extensions/retentionmapper" begin
    include("ext/jld2_extensions/retentionmapper.jl")
end
@safetestset "ext/makie_extensions/retentionmapper" begin
    include("ext/makie_extensions/retentionmapper.jl")
end
@safetestset "ext/makie_extensions/quadvarfit" begin
    include("ext/makie_extensions/quadvarfit.jl")
end
@safetestset "ext/makie_extensions/massspectrum" begin
    include("ext/makie_extensions/massspectrum.jl")
end
