using SafeTestsets

# Testing contents of base.jl
@safetestset "GCMS" begin include("base/GCMS.jl") end
@safetestset "FID" begin include("base/FID.jl") end
@safetestset "RiFID" begin include("base/RiFID.jl") end
@safetestset "TIC" begin include("base/TIC.jl") end
@safetestset "RiTIC" begin include("base/RiTIC.jl") end

@safetestset "AbstractGCMS" begin include("base/AbstractGCMS.jl") end
@safetestset "AbstractChromatogram" begin include("base/AbstractChromatogram.jl") end

@safetestset "IonScanOrder" begin include("base/IonScanOrder.jl") end

@safetestset "Utilities" begin include("utilities/utilities.jl") end
