using JuChrom
using Test
using Unitful: 𝐓


############################################################################################
# intensities(chrom::AbstractChrom; scanindexrange::OrdinalRange{Integer, Integer})
############################################################################################
@testset "intensities scanindexrange Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test [12, 956, 23] == intensities(Chrom([1, 2, 3]u"s", [12, 956, 23]))
    @test [956, 23] == intensities(Chrom([1, 2, 3]u"s", [12, 956, 23]), scanindexrange=2:3)
    @test [12] == intensities(Chrom([1, 2, 3]u"s", [12, 956, 23]), scanindexrange=1:1)

    # Same return same element type as used to construct the object
    @test Vector{Int} == typeof(intensities(Chrom([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Vector{Float64} == typeof(intensities(Chrom([1, 2, 3]u"s", 
        Float64[12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(intensities(Chrom([1, 2, 3]u"s", 10:5:20)))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(Chrom([1, 2, 3]u"s", 10.0:5:20.0)))

    @test SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(Chrom([1, 2, 3]u"s", Int[12, 956, 23]), scanindexrange=2:3))
    @test SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(Chrom([1, 2, 3]u"s", Float64[12.0, 956.0, 23.0]), scanindexrange=2:3))
    @test StepRange{Int64, Int64} == typeof(intensities(Chrom([1, 2, 3]u"s", 10:5:20), 
        scanindexrange=2:3))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(Chrom([1, 2, 3]u"s", 10.0:5:20.0), 
        scanindexrange=2:3))

    # Check for BoundsErrors
    @test_throws BoundsError intensities(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:3)
    @test_throws BoundsError intensities(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=1:4)
end


############################################################################################
# intensity(chrom::AbstractChrom, time::Unitful.Time; precisetime::Bool=false)
############################################################################################
@testset "intensity Chrom time" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s")
    @test 12 == intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=true)
    @test 956 == intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1.5u"s")
    @test 12 == intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), prevfloat(1.5)u"s")
    @test 12 == intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), (1/60)u"minute"; 
        precisetime=true)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]), 1u"s"))
    @test Float64 == typeof(intensity(Chrom([1, 2, 3]u"s", Float64[12.0, 956.0, 23.0]), 
        1u"s"))
    
    @test_throws ArgumentError intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1.5u"s", 
        precisetime=true)
end


############################################################################################
# intensity(chrom::AbstractChrom, scanindex::Integer)
############################################################################################
@testset "intensity Chrom scanindex" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1)
    @test 956 == intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]), 1))
    @test Float64 == typeof(intensity(Chrom([1, 2, 3]u"s", Float64[12.0, 956.0, 23.0]), 1))
    
    @test_throws BoundsError intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 0)
    @test_throws BoundsError intensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 4)
end


############################################################################################
# maxintensity(chrom::AbstractChrom; scanindexrange)
############################################################################################
@testset "maxintensity Chrom; scanindexrange" begin
    # Validate the returned value
    @test 956 == maxintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 12 == maxintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]), scanindexrange=1:2)
    @test 956 == maxintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]))
    @test 12 == maxintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(maxintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Int == typeof(maxintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:1))
    @test Int == typeof(maxintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:2))
    @test Float64 == typeof(maxintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23])))
    @test Float64 == typeof(maxintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 
        scanindexrange=1:1))
    @test Float64 == typeof(maxintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 
        scanindexrange=1:2))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:3)
    @test_throws BoundsError maxintensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=1:4)
    @test_throws BoundsError maxintensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:4)
end


############################################################################################
# maxretentionindex(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{T, S}]) 
# where {T<:Integer, S<:Integer}
############################################################################################
@testset "maxretentionindex ChromMS" begin
    @test_throws MethodError maxretentionindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
end


@testset "maxretentionindex Chrom" begin
    @test_throws MethodError maxretentionindex(Chrom([1, 2, 3]u"s", [12, 956, 23]))
end


############################################################################################
# maxscantime(chrom::AbstractChromatogram); timeunit::Unitful.TimeUnits, 
# ustripped::Bool=false)
############################################################################################
@testset "maxscantime ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 3u"s" == maxscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test (1/20)u"minute" ≈ maxscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute")
    @test 3 == maxscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test 1/20 ≈ maxscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Int == typeof(maxscantime(ChromMS(Int[1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true))
    @test Float64 == typeof(maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
end


@testset "maxscantime Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test 3u"s" == maxscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/20)u"minute" ≈ maxscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 3 == maxscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/20 ≈ maxscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(Chrom(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Int == typeof(maxscantime(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))
end


############################################################################################
# maxscantime(chrom::AbstractChromatogram, scanindexrange::OrdinalRange; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "maxscantime scanindexrange ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == maxscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1:2)
    @test (1/30)u"minute" ≈ maxscantime(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2, timeunit=u"minute")
    @test 2 == maxscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1:2, 
        ustripped=true)
    @test 1/30 ≈ maxscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1:2, 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), 1:2, timeunit=u"minute"))
    @test Int == typeof(maxscantime(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2, ustripped=true))
    @test Float64 == typeof(maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError maxscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "maxscantime scanindexrange Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == maxscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1:2)
    @test (1/30)u"minute" ≈ maxscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1:2, 
        timeunit=u"minute")
    @test 2 == maxscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1:2, ustripped=true)
    @test 1/30 ≈ maxscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1:2, timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23]), 1:2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(Chrom(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 1:2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:2, timeunit=u"minute"))
    @test Int == typeof(maxscantime(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), 1:2, 
        ustripped=true))
    @test Float64 == typeof(maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:2, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError maxscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0.0:4.0, timeunit=u"minute", ustripped=true)
end


############################################################################################
# scanduration(chrom::AbstractChromatogram; error::Real=0.001, 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "scanduration ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == scanduration(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test (1/60)u"minute" ≈ scanduration(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test 1 == scanduration(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test 1/60 ≈ scanduration(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Float64 == typeof(scanduration(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Float64 == typeof(scanduration(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))

    @test_throws ArgumentError scanduration(ChromMS([1]u"s", [12], reshape([0], length(1), 
        1)))
end


@testset "scanduration Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == scanduration(Chrom([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/60)u"minute" ≈ scanduration(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 1 == scanduration(Chrom([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ scanduration(Chrom([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Float64 == typeof(scanduration(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Float64 == typeof(scanduration(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), ustripped=true))

    @test_throws ArgumentError scanduration(Chrom([1]u"s", [12]))
end


############################################################################################
# metadata(chrom::AbstractChromatogram) -> Dict{Any, Any}
############################################################################################
@testset "metadata ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test Dict() == metadata(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test Dict(:id => 4, :name => "sample") == metadata(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 4, :name => "sample")))

    # Same return container and element type as used to construct the object
    @test Dict{Any, Any} == typeof(metadata(ChromMS([1, 2, 3]u"s", [85, 100], 
    [0 12; 34 956; 23 1], Dict(:id => 4, :name => "sample"))))   
end


@testset "metadata Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test Dict() == metadata(Chrom([1, 2, 3]u"s", [12, 956, 23]))

    # Same return container and element type as used to construct the object
    @test Dict(:id => 4, :name => "sample") == metadata(Chrom([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 4, :name => "sample")))
    @test Dict{Any, Any} == typeof(metadata(Chrom([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 4, :name => "sample"))))
end


############################################################################################
# minintensity(chrom::AbstractChromatogram[, greaterthan::Real; scanindexrange])
############################################################################################
@testset "minintensity Chrom [, greaterthan::Real; scanindexrange]" begin
    # Validate the returned value
    @test 12 == minintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 23 == minintensity(Chrom([1, 2, 3]u"s", Int[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(Chrom([1, 2, 3]u"s", Int[23, 12, 956]), scanindexrange=1:2)
    @test 12 == minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]))
    @test 23 == minintensity(Chrom([1, 2, 3]u"s", Float64[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(Chrom([1, 2, 3]u"s", Float64[23, 12, 956]), scanindexrange=1:2)

    @test 23 == minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 12)
    @test 956 == minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 12, 
        scanindexrange=1:2)

    @test nothing === minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 956)
    @test nothing === minintensity(Chrom([1, 2, 3]u"s", Float64[23, 12, 956]), 23, 
        scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(minintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Int == typeof(minintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:1))
    @test Int == typeof(minintensity(Chrom([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:2))
    @test Float64 == typeof(minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23])))
    @test Float64 == typeof(minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 
        scanindexrange=1:1))
    @test Float64 == typeof(minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 12, 
        scanindexrange=1:2))
    @test Nothing === typeof(minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 956))
    @test Nothing === typeof(minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 956, 
        scanindexrange=1:1))
    @test Nothing === typeof(minintensity(Chrom([1, 2, 3]u"s", Float64[12, 956, 23]), 956, 
        scanindexrange=1:1))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:3)
    @test_throws BoundsError maxintensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=1:4)
    @test_throws BoundsError maxintensity(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:4)
    # "Kovats", [100, 200, 300]
end


############################################################################################
# minretentionindex(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{T, S}]) 
# where {T<:Integer, S<:Integer}
############################################################################################
@testset "minretentionindex ChromMS" begin
    @test_throws MethodError minretentionindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
end


@testset "minretentionindex Chrom" begin
    @test_throws MethodError minretentionindex(Chrom([1, 2, 3]u"s", [12, 956, 23]))
end


############################################################################################
# minscantime(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange]; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "minscantime ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test (1/60)u"minute" ≈ minscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute")
    @test 1 == minscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test 1/60 ≈ minscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)

    @test 2u"s" == minscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3)
    @test (1/30)u"minute" ≈ minscantime(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute")
    @test 2 == minscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3, 
        ustripped=true)
    @test 1/30 ≈ minscantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3, 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Int == typeof(minscantime(ChromMS(Int[1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true))
    @test Float64 == typeof(minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))

    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Int == typeof(minscantime(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true))
    @test Float64 == typeof(minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError minscantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "minscantime Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/60)u"minute" ≈ minscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 1 == minscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ minscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    @test 2u"s" == minscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2:3)
    @test (1/30)u"minute" ≈ minscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute")
    @test 2 == minscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2:3, ustripped=true)
    @test 1/30 ≈ minscantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2:3, timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(Chrom(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Int == typeof(minscantime(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23]), 2:3))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(Chrom(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3, timeunit=u"minute"))
    @test Int == typeof(minscantime(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        ustripped=true))
    @test Float64 == typeof(minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError minscantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0.0:4.0, timeunit=u"minute", ustripped=true)
end


############################################################################################
# retentionindex(chrom::AbstractChromatogram, scanindex::Integer)
############################################################################################
@testset "retentionindex ChromMS scanindex" begin
    @test_throws MethodError retentionindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1)
end


@testset "retentionindex Chrom scanindex" begin
    @test_throws MethodError retentionindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1)
end


############################################################################################
# retentionindexname(chrom::AbstractChromatogram)
############################################################################################
@testset "retentionindexname ChromMS" begin
    @test_throws MethodError retentionindexname(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
end


@testset "retentionindexname Chrom" begin
    @test_throws MethodError retentionindexname(Chrom([1, 2, 3]u"s", [12, 956, 23]))
end


############################################################################################
# retentionindices(chrom::AbstractChromatogram)
############################################################################################
@testset "retentionindices ChromMS" begin
    @test_throws MethodError retentionindices(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
end


@testset "retentionindices Chrom" begin
    @test_throws MethodError retentionindices(Chrom([1, 2, 3]u"s", [12, 956, 23]))
end


############################################################################################
# rimapper(chrom::AbstractChromatogram)
############################################################################################
@testset "rimapper(::AbstractChromatogram)" begin
    @test nothing === rimapper(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test isa(rimapper(Chrom([1, 2, 3]u"s", [12, 956, 23], rimapper=RiMapper("Kovats", 
        (1:5)u"minute", 1000:1000:5000))), JuChrom.AbstractRiMapper)
    @test nothing === rimapper(Chrom([1, 2, 3]u"s", [12, 956, 23]))
    @test isa(rimapper(Chrom([1, 2, 3]u"s", [12, 956, 23], rimapper=RiMapper("Kovats", 
            (1:5)u"minute", 1000:1000:5000))), JuChrom.AbstractRiMapper)
end


############################################################################################
# rimapper(chrom::AbstractChromatogram, rim::AbstractRiMapper)
############################################################################################
@testset "rimapper!(::AbstractChromatogram, ::AbstractRiMapper)" begin
    chrom = Chrom([1, 2, 3]u"s", [12, 956, 23])
    ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000)
    @test nothing === chrom.rimapper
    rimapper!(chrom, ld)
    @test isa(chrom.rimapper, JuChrom.AbstractRiMapper)

    chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    @test nothing === chrom.rimapper
    rimapper!(chrom, ld)
    @test isa(chrom.rimapper, JuChrom.AbstractRiMapper)
end


############################################################################################
# runduration(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
# ustripped::Bool=false)
############################################################################################
@testset "runduration ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == runduration(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test (1//30)u"minute" == runduration(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test (1/30)u"minute" ≈ runduration(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test 2 == runduration(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test 1//30 == runduration(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)
    @test 1/30 ≈ runduration(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(runduration(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Int == typeof(runduration(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Float64 == typeof(runduration(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Rational{Int64} == typeof(runduration(ChromMS(Int[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(runduration(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
end


@testset "runduration Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == runduration(Chrom([1, 2, 3]u"s", [12, 956, 23]))
    @test (1//30)u"minute" == runduration(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test (1/30)u"minute" ≈ runduration(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 2 == runduration(Chrom([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1//30 == runduration(Chrom([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)
    @test 1/30 ≈ runduration(Chrom([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(runduration(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Int == typeof(runduration(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(runduration(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Rational{Int64} == typeof(runduration(Chrom(Int[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(runduration(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true))
end


############################################################################################
# scancount(chrom::AbstractChromatogram)
############################################################################################
@testset "scancount ChromMS" begin
    # Validate the returned value
    @test 3 == scancount(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test 1 == scancount(ChromMS([1]u"s", [85, 100], [0 12]))

    # Return value must be an integer
    @test Int == typeof(scancount(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])))
end


@testset "scancount Chrom" begin
    # Validate the returned value
    @test 3 == scancount(Chrom([1, 2, 3]u"s", [12, 956, 23]))
    @test 1 == scancount(Chrom([1]u"s", [12]))

    # Return value must be an integer
    @test Int == typeof(scancount(Chrom([1, 2, 3]u"s", [12, 956, 23])))
end


############################################################################################
# scantime(chrom::AbstractChromatogram, index::Integer; timeunit::Unitful.TimeUnits, 
# ustripped::Bool=false)
############################################################################################
@testset "scantime ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == scantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2)
    @test (1//20)u"minute" == scantime(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3, timeunit=u"minute")
    @test (1/20)u"minute" ≈ scantime(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3, timeunit=u"minute")
    @test 2 == scantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2, 
        ustripped=true)
    @test 1//20 == scantime(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 3, 
        timeunit=u"minute", ustripped=true)
    @test 1/20 ≈ scantime(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), 3, 
        timeunit=u"minute", ustripped=true)

    # Same return element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(scantime(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scantime(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, timeunit=u"minute"))
    @test Int == typeof(scantime(ChromMS(Int[1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        2, ustripped=true))
    @test Rational{Int64} == typeof(scantime(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(scantime(ChromMS(Float64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, ustripped=true))
    @test Float64 == typeof(scantime(ChromMS(Float64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, timeunit=u"minute", ustripped=true))
end


@testset "scantime Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == scantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2)
    @test (1//20)u"minute" == scantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 3,
        timeunit=u"minute")
    @test (1/20)u"minute" ≈ scantime(Chrom([1.0, 2.0, 3.0]u"s", [12, 956, 23]), 3,
        timeunit=u"minute")
    @test 2 == scantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2, ustripped=true)
    @test 1//20 == scantime(Chrom([1, 2, 3]u"s", [12, 956, 23]), 3, timeunit=u"minute", 
        ustripped=true)
    @test 1/20 ≈ scantime(Chrom([1.0, 2.0, 3.0]u"s", [12, 956, 23]), 3, timeunit=u"minute", 
        ustripped=true)

    # Same return element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}} == typeof(scantime(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23]), 2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}} == typeof(scantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}} == typeof(scantime(Chrom(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(scantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2, timeunit=u"minute"))
    @test Int == typeof(scantime(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), 2,
        ustripped=true))
    @test Rational{Int64} == typeof(scantime(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), 2,
        timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(scantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 2, 
        ustripped=true))
    @test Float64 == typeof(scantime(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 2, 
        timeunit=u"minute", ustripped=true))
end


############################################################################################
# scantimeindex(chrom::AbstractChromMS, chrom::AbstractChromatogram, time::Unitful.Time; 
# precisetime::Bool=false) -> Int
############################################################################################
@testset "scantimeindex ChromMS" begin
    # Validate the returned value
    @test 1 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s")
    @test 2 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2u"s")
    @test 3 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 3u"s")
    @test 1 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=true)
    @test 2 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2u"s", 
        precisetime=true)
    @test 3 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 3u"s", 
        precisetime=true)
    @test 1 == scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s")
    @test 2 == scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1u"s")
    @test 3 == scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3.1u"s")
    @test 1 == scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(ChromMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s")
    @test 2 == scantimeindex(ChromMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1u"s")
    @test 3 == scantimeindex(ChromMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3.1u"s")
    @test 1 == scantimeindex(ChromMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(ChromMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(ChromMS(Float32[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 0u"s")
    @test 3 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 4u"s")
    
    # Function throws ArgumentError if precisetime is set to true and the specified time 
    # does not does not exist
    @test_throws ArgumentError scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 4u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1001u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.9999u"s", precisetime=true)

    # Return value must be an integer
    @test Int == typeof(scantimeindex(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s"))
    @test Int == typeof(scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s"))
    @test Int == typeof(scantimeindex(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=true))
    @test Int == typeof(scantimeindex(ChromMS(Float64[1.1, 2.1, 3.1]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.1u"s", precisetime=true))

    # The function argument time must be a quantity with a valid time unit
    @test 1 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"ms")
    @test 3 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        1u"minute")
    @test_throws MethodError scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"m")
    @test_throws MethodError scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"m", precisetime=true)

    # The function argument precisetime must be a Boolean value
    @test 1 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=true)
    @test 1 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=false)
    @test 1 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=:true)
    @test 1 == scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1u"s", 
        precisetime=:false)
    @test_throws TypeError scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=:truly)
    @test_throws TypeError scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime="true")
    @test_throws TypeError scantimeindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=1)
end


@testset "scantimeindex Chrom" begin
    # Validate the returned value
    @test 1 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s")
    @test 2 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2u"s")
    @test 3 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 3u"s")
    @test 1 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=true)
    @test 2 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2u"s", precisetime=true)
    @test 3 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 3u"s", precisetime=true)
    @test 1 == scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s", 
        precisetime=true)
    @test 2 == scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s", 
        precisetime=true)
    @test 3 == scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s", 
        precisetime=true)
    @test 1 == scantimeindex(Chrom(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(Chrom(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(Chrom(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(Chrom(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 1.1u"s", 
        precisetime=true)
    @test 2 == scantimeindex(Chrom(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 2.1u"s", 
        precisetime=true)
    @test 3 == scantimeindex(Chrom(Float32[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 3.1u"s", 
        precisetime=true)
    @test 1 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 0u"s")
    @test 3 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 4u"s")
    
    # Function throws ArgumentError if precisetime is set to true and the specified time 
    # does not does not exist
    @test_throws ArgumentError scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 0u"s", 
        precisetime=true)
    @test_throws ArgumentError scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 4u"s", 
        precisetime=true)
    @test_throws ArgumentError scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]),
        2.1001u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]),
        1.9999u"s", precisetime=true)

    # Return value must be an integer
    @test Int == typeof(scantimeindex(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), 1u"s"))
    @test Int == typeof(scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 
        1.1u"s"))
    @test Int == typeof(scantimeindex(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=true))
    @test Int == typeof(scantimeindex(Chrom(Float64[1.1, 2.1, 3.1]u"s", [12, 956, 23]), 
        1.1u"s", precisetime=true))

    # The function argument time must be a quantity with a valid time unit
    @test 1 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"ms")
    @test 3 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"minute")
    @test_throws MethodError scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"m")
    @test_throws MethodError scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"m", 
        precisetime=true)

    # The function argument precisetime must be a Boolean value
    @test 1 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=true)
    @test 1 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=false)
    @test 1 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=:true)
    @test 1 == scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s", precisetime=:false)
    @test_throws TypeError scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=:truly)
    @test_throws TypeError scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime="true")
    @test_throws TypeError scantimeindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=1)
end


############################################################################################
# scantimes(chrom::AbstractChromatogram; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "scantimes ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [1, 2, 3]u"s" == scantimes(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test [1, 2, 3] == scantimes(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true)

    # Check if reference to data structure is returned
    chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    @test scantimes(chrom) === scantimes(chrom)
    
    # Same return container and element type as used to construct the object
    @test Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}}} == typeof(scantimes(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}}} == typeof(scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(ChromMS(Int64[1, 2, 3]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Vector{Float64} == typeof(scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
end


@testset "scantimes Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test [1, 2, 3]u"s" == scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]))
    @test [1//60, 1//30, 1//20]u"minute" == scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute")
    @test [1, 2, 3] == scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test [1//60, 1//30, 1//20] == scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true)

    # Check if reference to data structure is returned
    tic = Chrom([1, 2, 3]u"s", [12, 956, 23])
    @test scantimes(tic) === scantimes(tic)

    # Same return container and element type as used to construct the object
    @test Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}}} == typeof(scantimes(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}}} == typeof(scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23])))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(Chrom(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(Chrom(Int[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute", ustripped=true))
end


############################################################################################
# scantimes(chrom::AbstractChromatogram, scanindexrange::OrdinalRange; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "scantimes with scanindexrange ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [2, 3]u"s" == scantimes(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3)
    @test [1/30, 1/20]u"minute" ≈ scantimes(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute")
    @test [2, 3] == scantimes(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3, 
        ustripped=true)
    @test [1/30, 1/20] ≈ scantimes(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        2:3, timeunit=u"minute", ustripped=true)

    # Check if reference to data structure is returned
    chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    @test scantimes(chrom, 2:3) === scantimes(chrom, 2:3)
    
    # Same return container and element type as used to construct the object
    @test SubArray{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3))
    @test SubArray{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(ChromMS(Int64[1, 2, 3]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true))
    @test Vector{Float64} == typeof(scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(ChromMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute", ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError scantimes(ChromMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "scantimes with scanindexrange Chrom" begin
    # Same return values as those provided as arguments to construct the object
    @test [2, 3]u"s" == scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2:3)
    @test [1//30, 1//20]u"minute" == scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute")
    @test [1/30, 1/20]u"minute" ≈ scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3, timeunit=u"minute")
    @test [2, 3] == scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2:3, ustripped=true)
    @test [1//30, 1//20] == scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute", ustripped=true)
    @test [1/30, 1/20] ≈ scantimes(Chrom([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute", ustripped=true)
    
    # Check if reference to data structure is returned
    tic = Chrom([1, 2, 3]u"s", [12, 956, 23])
    @test scantimes(tic, 2:3) === scantimes(tic, 2:3)
    
    # Same return container and element type as used to construct the object
    @test SubArray{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true}  == typeof(scantimes(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23]), 2:3))
    @test SubArray{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(Chrom(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(Chrom(Int[1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(Chrom(Int[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 0:3, 
        timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 1:4, 
        timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 0:4, 
        timeunit=u"minute", ustripped=true)
    @test_throws MethodError scantimes(Chrom(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0.0:4.0, timeunit=u"minute", ustripped=true)
end
