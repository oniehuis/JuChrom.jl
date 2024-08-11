using SafeTestsets

# Testing contents of base.jl
@safetestset "GCMS" begin include("base/GCMS.jl") end
@safetestset "FID" begin include("base/FID.jl") end
@safetestset "TIC" begin include("base/TIC.jl") end

@safetestset "AbstractGCMS" begin include("base/AbstractGCMS.jl") end
@safetestset "AbstractChromatogram" begin include("base/AbstractChromatogram.jl") end
