using SafeTestsets

# Testing contents of base.jl
@safetestset "base ChromMS" begin include("base/ChromMS.jl") end
@safetestset "base Chrom" begin include("base/Chrom.jl") end
@safetestset "base AbstractChromMS" begin include("base/AbstractChromMS.jl") end
@safetestset "base AbstractChromatogram" begin include("base/AbstractChromatogram.jl") end
@safetestset "base RiMapper" begin include("base/RiMapper.jl") end
@safetestset "base IonScanOrder" begin include("base/IonScanOrder.jl") end
@safetestset "base MassSpectrum" begin include("base/MassSpectrum.jl") end
@safetestset "base Utilities" begin include("utilities/utilities.jl") end

@safetestset "IO AgilentFID" begin include("IO/AgilentFID.jl") end
@safetestset "IO ANDI" begin include("IO/ANDI.jl") end
@safetestset "IO ChemStationMS" begin include("IO/ChemStationMS.jl") end
@safetestset "IO DelimitedTextWriter" begin include("IO/DelimitedTextWriter.jl") end
@safetestset "IO Excel" begin include("IO/Excel.jl") end
@safetestset "IO MassHunterMS" begin include("IO/MassHunterMS.jl") end
@safetestset "IO ShimadzuMS" begin include("IO/ShimadzuMS.jl") end
@safetestset "IO Utilities" begin include("IO/utilities.jl") end

# @safetestset "Deconvolution deconvolution" begin include("deconvolution/deconvolution.jl") end