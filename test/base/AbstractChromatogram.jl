using JuChrom
using Test
using Unitful: 𝐓


############################################################################################
# intensities(chrom::AbstractChromatogram)
############################################################################################
@testset "intensities FID" begin
    # Same return values as those provided as arguments to construct the object
    @test [12, 956, 23] == intensities(FID([1, 2, 3]u"s", [12, 956, 23]))

    # Same return container and element type as used to construct the object
    @test Vector{Int} == typeof(intensities(FID([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Vector{Float64} == typeof(intensities(FID([1, 2, 3]u"s", 
        Float64[12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(intensities(FID([1, 2, 3]u"s", 10:5:20)))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(FID([1, 2, 3]u"s", 10.0:5:20.0)))
end


@testset "intensities GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [0 12; 34 956; 23 1] == intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
    @test reshape([0, 956, 1], (:,1)) == intensities(GCMS([1, 2, 3]u"s", [85],
        reshape([0, 956, 1], (:,1))))

    # Same return container and element type as used to construct the object
    @test Matrix{Int64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Int64[0 12; 34 956; 23 1])))
    @test Matrix{Float64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1])))
    @test Matrix{Int64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1],
         (:,1)))))
end


@testset "intensities TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test [12, 956, 23] == intensities(TIC([1, 2, 3]u"s", [12, 956, 23]))

    # Same return container and element type as used to construct the object
    @test Vector{Int} == typeof(intensities(TIC([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Vector{Float64} == typeof(intensities(TIC([1, 2, 3]u"s", 
        Float64[12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(intensities(TIC([1, 2, 3]u"s", 10:5:20)))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(TIC([1, 2, 3]u"s", 10.0:5:20.0)))
end


############################################################################################
# maxintensity(chrom::AbstractChromatogram)
############################################################################################
@testset "maxintensity FID" begin
    # Validate the returned value
    @test 956 == maxintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 956 == maxintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]))

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(maxintensity(FID([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Float64 == typeof(maxintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23])))
end


@testset "maxintensity GCMS" begin
    # Validate the returned value
    @test 956 == maxintensity(GCMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]))
    @test 956 == maxintensity(GCMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]))

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(maxintensity(GCMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1])))
    @test Float64 == typeof(maxintensity(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1])))
end


@testset "maxintensity TIC" begin
    # Validate the returned value
    @test 956 == maxintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 956 == maxintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]))

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(maxintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Float64 == typeof(maxintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23])))
end


############################################################################################
# maxscantime(chrom::AbstractChromatogram); timeunit::Unitful.TimeUnits, 
# ustripped::Bool=false)
############################################################################################
@testset "maxscantime FID" begin
    # Same return values as those provided as arguments to construct the object
    @test 3u"s" == maxscantime(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/20)u"minute" ≈ maxscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 3 == maxscantime(FID([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/20 ≈ maxscantime(FID([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(FID(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(FID(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Int == typeof(maxscantime(FID(Int[1, 2, 3]u"s", [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))
end


@testset "maxscantime GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 3u"s" == maxscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test (1/20)u"minute" ≈ maxscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute")
    @test 3 == maxscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test 1/20 ≈ maxscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Int == typeof(maxscantime(GCMS(Int[1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true))
    @test Float64 == typeof(maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
end


@testset "maxscantime TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 3u"s" == maxscantime(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/20)u"minute" ≈ maxscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 3 == maxscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/20 ≈ maxscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(TIC(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(TIC(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Int == typeof(maxscantime(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))
end


############################################################################################
# metadata(chrom::AbstractChromatogram) -> Dict{Any, Any}
############################################################################################
@testset "metadata FID" begin
    # Same return values as those provided as arguments to construct the object
    @test Dict() == metadata(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test Dict(:id => 4, :name => "sample") == metadata(FID([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 4, :name => "sample")))
    
    # Same return container and element type as used to construct the object
    @test Dict{Any, Any} == typeof(metadata(FID([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 4, :name => "sample"))))
end


@testset "metadata GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test Dict() == metadata(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test Dict(:id => 4, :name => "sample") == metadata(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 4, :name => "sample")))

    # Same return container and element type as used to construct the object
    @test Dict{Any, Any} == typeof(metadata(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 4, :name => "sample"))))   
end


@testset "metadata TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test Dict() == metadata(TIC([1, 2, 3]u"s", [12, 956, 23]))

    # Same return container and element type as used to construct the object
    @test Dict(:id => 4, :name => "sample") == metadata(TIC([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 4, :name => "sample")))
    @test Dict{Any, Any} == typeof(metadata(TIC([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 4, :name => "sample"))))
end


############################################################################################
# minintensity(chrom::AbstractChromatogram)
############################################################################################
@testset "minintensity FID" begin
    # Validate the returned value
    @test 12 == minintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 12 == minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]))

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(minintensity(FID([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Float64 == typeof(minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23])))
end


@testset "minintensity GCMS" begin
    # Validate the returned value
    @test 0 == minintensity(GCMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]))
    @test 0 == minintensity(GCMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]))

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(minintensity(GCMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1])))
    @test Float64 == typeof(minintensity(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1])))
end


@testset "minintensity TIC" begin
    # Validate the returned value
    @test 12 == minintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 12 == minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]))

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(minintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Float64 == typeof(minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23])))
end


############################################################################################
# minscantime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
# ustripped::Bool=false)
############################################################################################
@testset "minscantime FID" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/60)u"minute" ≈ minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 1 == minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(FID(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(FID(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Int == typeof(minscantime(FID(Int[1, 2, 3]u"s", [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))
end


@testset "minscantime GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test (1/60)u"minute" ≈ minscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute")
    @test 1 == minscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test 1/60 ≈ minscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Int == typeof(minscantime(GCMS(Int[1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true))
    @test Float64 == typeof(minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
end


@testset "minscantime TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/60)u"minute" ≈ minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 1 == minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(TIC(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(TIC(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Int == typeof(minscantime(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))
end


############################################################################################
# runduration(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
# ustripped::Bool=false)
############################################################################################
@testset "runduration FID" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == runduration(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test (1//30)u"minute" == runduration(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test (1/30)u"minute" ≈ runduration(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 2 == runduration(FID([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1//30 == runduration(FID([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)
    @test 1/30 ≈ runduration(FID([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(FID(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(runduration(FID(Int64[1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Int == typeof(runduration(FID(Int[1, 2, 3]u"s", [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(runduration(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Rational{Int64} == typeof(runduration(FID(Int[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(runduration(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true))
end


@testset "runduration GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == runduration(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test (1//30)u"minute" == runduration(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test (1/30)u"minute" ≈ runduration(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test 2 == runduration(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test 1//30 == runduration(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)
    @test 1/30 ≈ runduration(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(runduration(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Int == typeof(runduration(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Float64 == typeof(runduration(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Rational{Int64} == typeof(runduration(GCMS(Int[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(runduration(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
end


@testset "runduration TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == runduration(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test (1//30)u"minute" == runduration(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test (1/30)u"minute" ≈ runduration(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 2 == runduration(TIC([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1//30 == runduration(TIC([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)
    @test 1/30 ≈ runduration(TIC([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(TIC(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(runduration(TIC(Int64[1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Int == typeof(runduration(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(runduration(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Rational{Int64} == typeof(runduration(TIC(Int[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(runduration(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true))
end


############################################################################################
# scancount(chrom::AbstractChromatogram)
############################################################################################
@testset "scancount FID" begin
    # Validate the returned value
    @test 3 == scancount(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test 1 == scancount(FID([1]u"s", [12]))

    # Return value must be an integer
    @test Int == typeof(scancount(FID([1, 2, 3]u"s", [12, 956, 23])))
end


@testset "scancount GCMS" begin
    # Validate the returned value
    @test 3 == scancount(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test 1 == scancount(GCMS([1]u"s", [85, 100], [0 12]))

    # Return value must be an integer
    @test Int == typeof(scancount(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])))
end


@testset "scancount TIC" begin
    # Validate the returned value
    @test 3 == scancount(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test 1 == scancount(TIC([1]u"s", [12]))

    # Return value must be an integer
    @test Int == typeof(scancount(TIC([1, 2, 3]u"s", [12, 956, 23])))
end


############################################################################################
# scantime(chrom::AbstractChromatogram, index::Integer; timeunit::Unitful.TimeUnits, 
# ustripped::Bool=false)
############################################################################################
@testset "scantime FID" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == scantime(FID([1, 2, 3]u"s", [12, 956, 23]), 2)
    @test (1//20)u"minute" == scantime(FID([1, 2, 3]u"s", [12, 956, 23]), 3,
        timeunit=u"minute")
    @test (1/20)u"minute" ≈ scantime(FID([1.0, 2.0, 3.0]u"s", [12, 956, 23]), 3,
        timeunit=u"minute")
    @test 2 == scantime(FID([1, 2, 3]u"s", [12, 956, 23]), 2, ustripped=true)
    @test 1//20 == scantime(FID([1, 2, 3]u"s", [12, 956, 23]), 3, timeunit=u"minute", 
        ustripped=true)
    @test 1/20 ≈ scantime(FID([1.0, 2.0, 3.0]u"s", [12, 956, 23]), 3, timeunit=u"minute", 
        ustripped=true)

    # Same return element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}} == typeof(scantime(FID(Int64[1, 2, 3]u"s", [12, 956, 23]), 2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}} == typeof(scantime(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}} == typeof(scantime(FID(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(scantime(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2, timeunit=u"minute"))
    @test Int == typeof(scantime(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 2,
        ustripped=true))
    @test Rational{Int64} == typeof(scantime(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 2,
        timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(scantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 2, 
        ustripped=true))
    @test Float64 == typeof(scantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 2, 
        timeunit=u"minute", ustripped=true))
end


@testset "scantime GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == scantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2)
    @test (1//20)u"minute" == scantime(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3, timeunit=u"minute")
    @test (1/20)u"minute" ≈ scantime(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3, timeunit=u"minute")
    @test 2 == scantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2, 
        ustripped=true)
    @test 1//20 == scantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 3, 
        timeunit=u"minute", ustripped=true)
    @test 1/20 ≈ scantime(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), 3, 
        timeunit=u"minute", ustripped=true)

    # Same return element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(scantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, timeunit=u"minute"))
    @test Int == typeof(scantime(GCMS(Int[1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        2, ustripped=true))
    @test Rational{Int64} == typeof(scantime(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(scantime(GCMS(Float64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, ustripped=true))
    @test Float64 == typeof(scantime(GCMS(Float64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, timeunit=u"minute", ustripped=true))
end


@testset "scantime TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == scantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 2)
    @test (1//20)u"minute" == scantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 3,
        timeunit=u"minute")
    @test (1/20)u"minute" ≈ scantime(TIC([1.0, 2.0, 3.0]u"s", [12, 956, 23]), 3,
        timeunit=u"minute")
    @test 2 == scantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 2, ustripped=true)
    @test 1//20 == scantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 3, timeunit=u"minute", 
        ustripped=true)
    @test 1/20 ≈ scantime(TIC([1.0, 2.0, 3.0]u"s", [12, 956, 23]), 3, timeunit=u"minute", 
        ustripped=true)

    # Same return element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}} == typeof(scantime(TIC(Int64[1, 2, 3]u"s", [12, 956, 23]), 2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}} == typeof(scantime(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}} == typeof(scantime(TIC(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(scantime(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2, timeunit=u"minute"))
    @test Int == typeof(scantime(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 2,
        ustripped=true))
    @test Rational{Int64} == typeof(scantime(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 2,
        timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(scantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 2, 
        ustripped=true))
    @test Float64 == typeof(scantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 2, 
        timeunit=u"minute", ustripped=true))
end


############################################################################################
# scantimeindex(gcms::AbstractGCMS, chrom::AbstractChromatogram, time::Unitful.Time; 
# precisetime::Bool=false) -> Int
############################################################################################
@testset "scantimeindex FID" begin
    # Validate the returned value
    @test 1 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s")
    @test 2 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 2u"s")
    @test 3 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 3u"s")
    @test 1 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=true)
    @test 2 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 2u"s", precisetime=true)
    @test 3 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 3u"s", precisetime=true)
    @test 1 == scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s", 
        precisetime=true)
    @test 2 == scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s", 
        precisetime=true)
    @test 3 == scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s", 
        precisetime=true)
    @test 1 == scantimeindex(FID(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(FID(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(FID(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(FID(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s", 
        precisetime=true)
    @test 2 == scantimeindex(FID(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s", 
        precisetime=true)
    @test 3 == scantimeindex(FID(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s", 
        precisetime=true)
    @test 1 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 0u"s")
    @test 3 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 4u"s")
    
    # Function throws ArgumentError if precisetime is set to true and the specified time 
    # does not does not exist
    @test_throws ArgumentError scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 0u"s", 
        precisetime=true)
    @test_throws ArgumentError scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 4u"s", 
        precisetime=true)
    @test_throws ArgumentError scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]),
        2.1001u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]),
        1.9999u"s", precisetime=true)

    # Return value must be an integer
    @test Int == typeof(scantimeindex(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 1u"s"))
    @test Int == typeof(scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 
        1.1u"s"))
    @test Int == typeof(scantimeindex(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=true))
    @test Int == typeof(scantimeindex(FID(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 
        1.1u"s", precisetime=true))

    # The function argument time must be a quantity with a valid time unit
    @test 1 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"ms")
    @test 3 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"minute")
    @test_throws MethodError scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"m")
    @test_throws MethodError scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"m", 
        precisetime=true)

    # The function argument precisetime must be a Boolean value
    @test 1 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=true)
    @test 1 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=false)
    @test 1 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=:true)
    @test 1 == scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=:false)
    @test_throws TypeError scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=:truly)
    @test_throws TypeError scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime="true")
    @test_throws TypeError scantimeindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=1)
end


@testset "scantimeindex GCMS" begin
    # Validate the returned value
    @test 1 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s")
    @test 2 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2u"s")
    @test 3 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 3u"s")
    @test 1 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=true)
    @test 2 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2u"s", 
        precisetime=true)
    @test 3 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 3u"s", 
        precisetime=true)
    @test 1 == scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s")
    @test 2 == scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1u"s")
    @test 3 == scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3.1u"s")
    @test 1 == scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(GCMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s")
    @test 2 == scantimeindex(GCMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1u"s")
    @test 3 == scantimeindex(GCMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3.1u"s")
    @test 1 == scantimeindex(GCMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(GCMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(GCMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 0u"s")
    @test 3 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 4u"s")
    
    # Function throws ArgumentError if precisetime is set to true and the specified time 
    # does not does not exist
    @test_throws ArgumentError scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 4u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1001u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.9999u"s", precisetime=true)

    # Return value must be an integer
    @test Int == typeof(scantimeindex(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s"))
    @test Int == typeof(scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s"))
    @test Int == typeof(scantimeindex(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=true))
    @test Int == typeof(scantimeindex(GCMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s", precisetime=true))

    # The function argument time must be a quantity with a valid time unit
    @test 1 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"ms")
    @test 3 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        1u"minute")
    @test_throws MethodError scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"m")
    @test_throws MethodError scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"m", precisetime=true)

    # The function argument precisetime must be a Boolean value
    @test 1 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=true)
    @test 1 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=false)
    @test 1 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=:true)
    @test 1 == scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=:false)
    @test_throws TypeError scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=:truly)
    @test_throws TypeError scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime="true")
    @test_throws TypeError scantimeindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=1)
end


@testset "scantimeindex TIC" begin
    # Validate the returned value
    @test 1 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s")
    @test 2 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 2u"s")
    @test 3 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 3u"s")
    @test 1 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=true)
    @test 2 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 2u"s", precisetime=true)
    @test 3 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 3u"s", precisetime=true)
    @test 1 == scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s", 
        precisetime=true)
    @test 2 == scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s", 
        precisetime=true)
    @test 3 == scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s", 
        precisetime=true)
    @test 1 == scantimeindex(TIC(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(TIC(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(TIC(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(TIC(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s", 
        precisetime=true)
    @test 2 == scantimeindex(TIC(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s", 
        precisetime=true)
    @test 3 == scantimeindex(TIC(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s", 
        precisetime=true)
    @test 1 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 0u"s")
    @test 3 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 4u"s")
    
    # Function throws ArgumentError if precisetime is set to true and the specified time 
    # does not does not exist
    @test_throws ArgumentError scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 0u"s", 
        precisetime=true)
    @test_throws ArgumentError scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 4u"s", 
        precisetime=true)
    @test_throws ArgumentError scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]),
        2.1001u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]),
        1.9999u"s", precisetime=true)

    # Return value must be an integer
    @test Int == typeof(scantimeindex(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 1u"s"))
    @test Int == typeof(scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 
        1.1u"s"))
    @test Int == typeof(scantimeindex(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=true))
    @test Int == typeof(scantimeindex(TIC(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 
        1.1u"s", precisetime=true))

    # The function argument time must be a quantity with a valid time unit
    @test 1 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"ms")
    @test 3 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"minute")
    @test_throws MethodError scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"m")
    @test_throws MethodError scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"m", 
        precisetime=true)

    # The function argument precisetime must be a Boolean value
    @test 1 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=true)
    @test 1 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=false)
    @test 1 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=:true)
    @test 1 == scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=:false)
    @test_throws TypeError scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=:truly)
    @test_throws TypeError scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime="true")
    @test_throws TypeError scantimeindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=1)
end

############################################################################################
# scantimes(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "scantimes FID" begin
    # Same return values as those provided as arguments to construct the object
    @test [1, 2, 3]u"s" == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test [1//60, 1//30, 1//20]u"minute" == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute")
    @test [1, 2, 3] == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test [1//60, 1//30, 1//20] == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}}} == typeof(scantimes(FID(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}}} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23])))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(FID(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(FID(Int[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute", ustripped=true))
end


@testset "scantimes GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [1, 2, 3]u"s" == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test [1, 2, 3] == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}}} == typeof(scantimes(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}}} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(GCMS(Int64[1, 2, 3]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Vector{Float64} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
end


@testset "scantimes TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test [1, 2, 3]u"s" == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test [1//60, 1//30, 1//20]u"minute" == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute")
    @test [1, 2, 3] == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test [1//60, 1//30, 1//20] == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}}} == typeof(scantimes(TIC(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}}} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23])))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(TIC(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(TIC(Int[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute", ustripped=true))
end
