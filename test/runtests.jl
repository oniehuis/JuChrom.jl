using SafeTestsets

# Testing contents of base.jl
@safetestset "base GCMS" begin include("base/GCMS.jl") end
@safetestset "base RiGCMS" begin include("base/RiGCMS.jl") end
@safetestset "base FID" begin include("base/FID.jl") end
@safetestset "base RiFID" begin include("base/RiFID.jl") end
@safetestset "base TIC" begin include("base/TIC.jl") end
@safetestset "base RiTIC" begin include("base/RiTIC.jl") end
@safetestset "base AbstractGCMS" begin include("base/AbstractGCMS.jl") end
@safetestset "base AbstractChromatogram" begin include("base/AbstractChromatogram.jl") end
@safetestset "base IonScanOrder" begin include("base/IonScanOrder.jl") end
@safetestset "base Utilities" begin include("utilities/utilities.jl") end

@safetestset "IO Utilities" begin include("IO/utilities.jl") end
@safetestset "IO ChemStationMS" begin include("IO/ChemStationMS.jl") end
@safetestset "IO AgilentFID" begin include("IO/AgilentFID.jl") end
@safetestset "IO MassHunterMS" begin include("IO/MassHunterMS.jl") end