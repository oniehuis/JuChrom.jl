using JuChrom
using Test
using Unitful: 𝐓


############################################################################################
# intensities(chrom::AbstractGC; scanindexrange::OrdinalRange{Integer, Integer})
############################################################################################
@testset "intensities scanindexrange FID" begin
    # Same return values as those provided as arguments to construct the object
    @test [12, 956, 23] == intensities(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test [956, 23] == intensities(FID([1, 2, 3]u"s", [12, 956, 23]), scanindexrange=2:3)
    @test [12] == intensities(FID([1, 2, 3]u"s", [12, 956, 23]), scanindexrange=1:1)

    # Same return same element type as used to construct the object
    @test Vector{Int} == typeof(intensities(FID([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Vector{Float64} == typeof(intensities(FID([1, 2, 3]u"s", 
        Float64[12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(intensities(FID([1, 2, 3]u"s", 10:5:20)))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(FID([1, 2, 3]u"s", 10.0:5:20.0)))
    @test SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(FID([1, 2, 3]u"s", Int[12, 956, 23]), scanindexrange=2:3))
    @test SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(FID([1, 2, 3]u"s", Float64[12.0, 956.0, 23.0]), scanindexrange=2:3))
    @test StepRange{Int64, Int64} == typeof(intensities(FID([1, 2, 3]u"s", 10:5:20), 
        scanindexrange=2:3))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(FID([1, 2, 3]u"s", 10.0:5:20.0), 
        scanindexrange=2:3))

    # Check for BoundsErrors
    @test_throws BoundsError intensities(FID([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:3)
    @test_throws BoundsError intensities(FID([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=1:4)
end


@testset "intensities scanindexrange RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test [12, 956, 23] == intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test [956, 23] == intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=2:3)
    @test [12] == intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=1:1)

    # Same return same element type as used to construct the object
    @test Vector{Int} == typeof(intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23])))
    @test Vector{Float64} == typeof(intensities(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], Float64[12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(intensities(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], 10:5:20)))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
            10.0:5:20.0)))
    @test SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], Int[12, 956, 23]), 
        scanindexrange=2:3))
    @test SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12.0, 956.0, 23.0]), scanindexrange=2:3))
    @test StepRange{Int64, Int64} == typeof(intensities(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], 10:5:20), scanindexrange=2:3))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        10.0:5:20.0), scanindexrange=2:3))

    # Check for BoundsErrors
    @test_throws BoundsError intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:3)
    @test_throws BoundsError intensities(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=1:4)
end


@testset "intensities scanindexrange TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test [12, 956, 23] == intensities(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test [956, 23] == intensities(TIC([1, 2, 3]u"s", [12, 956, 23]), scanindexrange=2:3)
    @test [12] == intensities(TIC([1, 2, 3]u"s", [12, 956, 23]), scanindexrange=1:1)

    # Same return same element type as used to construct the object
    @test Vector{Int} == typeof(intensities(TIC([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Vector{Float64} == typeof(intensities(TIC([1, 2, 3]u"s", 
        Float64[12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(intensities(TIC([1, 2, 3]u"s", 10:5:20)))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(TIC([1, 2, 3]u"s", 10.0:5:20.0)))

    @test SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(TIC([1, 2, 3]u"s", Int[12, 956, 23]), scanindexrange=2:3))
    @test SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(TIC([1, 2, 3]u"s", Float64[12.0, 956.0, 23.0]), scanindexrange=2:3))
    @test StepRange{Int64, Int64} == typeof(intensities(TIC([1, 2, 3]u"s", 10:5:20), 
        scanindexrange=2:3))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(TIC([1, 2, 3]u"s", 10.0:5:20.0), 
        scanindexrange=2:3))

    # Check for BoundsErrors
    @test_throws BoundsError intensities(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:3)
    @test_throws BoundsError intensities(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=1:4)
end


@testset "intensities scanindexrange RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test [12, 956, 23] == intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test [956, 23] == intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=2:3)
    @test [12] == intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=1:1)

    # Same return same element type as used to construct the object
    @test Vector{Int} == typeof(intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23])))
    @test Vector{Float64} == typeof(intensities(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], Float64[12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(intensities(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], 10:5:20)))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
            10.0:5:20.0)))
    @test SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], Int[12, 956, 23]), 
        scanindexrange=2:3))
    @test SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} == typeof(
        intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12.0, 956.0, 23.0]), scanindexrange=2:3))
    @test StepRange{Int64, Int64} == typeof(intensities(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], 10:5:20), scanindexrange=2:3))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        10.0:5:20.0), scanindexrange=2:3))

    # Check for BoundsErrors
    @test_throws BoundsError intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:3)
    @test_throws BoundsError intensities(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=1:4)
end


############################################################################################
# intensity(chrom::AbstractGC, time::Unitful.Time; precisetime::Bool=false)
############################################################################################
@testset "intensity FID time" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s")
    @test 12 == intensity(FID([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=true)
    @test 956 == intensity(FID([1, 2, 3]u"s", [12, 956, 23]), 1.5u"s")
    @test 12 == intensity(FID([1, 2, 3]u"s", [12, 956, 23]), prevfloat(1.5)u"s")
    @test 12 == intensity(FID([1, 2, 3]u"s", [12, 956, 23]), (1/60)u"minute"; 
        precisetime=true)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(FID([1, 2, 3]u"s", Int[12, 956, 23]), 1u"s"))
    @test Float64 == typeof(intensity(FID([1, 2, 3]u"s", Float64[12.0, 956.0, 23.0]), 
        1u"s"))
    
    @test_throws ArgumentError intensity(FID([1, 2, 3]u"s", [12, 956, 23]), 1.5u"s", 
        precisetime=true)
end


@testset "intensity RiFID time" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1u"s")
    @test 12 == intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1u"s", precisetime=true)
    @test 956 == intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1.5u"s")
    @test 12 == intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        prevfloat(1.5)u"s")
    @test 12 == intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        (1/60)u"minute"; precisetime=true)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), 1u"s"))
    @test Float64 == typeof(intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12.0, 956.0, 23.0]), 1u"s"))
    
    @test_throws ArgumentError intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.5u"s", precisetime=true)
end


@testset "intensity TIC time" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s")
    @test 12 == intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 1u"s", 
        precisetime=true)
    @test 956 == intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 1.5u"s")
    @test 12 == intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), prevfloat(1.5)u"s")
    @test 12 == intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), (1/60)u"minute"; 
        precisetime=true)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]), 1u"s"))
    @test Float64 == typeof(intensity(TIC([1, 2, 3]u"s", Float64[12.0, 956.0, 23.0]), 
        1u"s"))
    
    @test_throws ArgumentError intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 1.5u"s", 
        precisetime=true)
end


@testset "intensity RiTIC time" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1u"s")
    @test 12 == intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1u"s", precisetime=true)
    @test 956 == intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1.5u"s")
    @test 12 == intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        prevfloat(1.5)u"s")
    @test 12 == intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        (1/60)u"minute"; precisetime=true)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), 1u"s"))
    @test Float64 == typeof(intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12.0, 956.0, 23.0]), 1u"s"))
    
    @test_throws ArgumentError intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.5u"s", precisetime=true)
end


############################################################################################
# intensity(chrom::AbstractGC, scanindex::Integer)
############################################################################################
@testset "intensity FID scanindex" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(FID([1, 2, 3]u"s", [12, 956, 23]), 1)
    @test 956 == intensity(FID([1, 2, 3]u"s", [12, 956, 23]), 2)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(FID([1, 2, 3]u"s", Int[12, 956, 23]), 1))
    @test Float64 == typeof(intensity(FID([1, 2, 3]u"s", Float64[12.0, 956.0, 23.0]), 1))
    
    @test_throws BoundsError intensity(FID([1, 2, 3]u"s", [12, 956, 23]), 0)
    @test_throws BoundsError intensity(FID([1, 2, 3]u"s", [12, 956, 23]), 4)
end


@testset "intensity RiFID scanindex" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1)
    @test 956 == intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        2)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), 1))
    @test Float64 == typeof(intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12.0, 956.0, 23.0]), 1))
    
    @test_throws BoundsError intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 0)
    @test_throws BoundsError intensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 4)
end


@testset "intensity TIC scanindex" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 1)
    @test 956 == intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 2)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]), 1))
    @test Float64 == typeof(intensity(TIC([1, 2, 3]u"s", Float64[12.0, 956.0, 23.0]), 1))
    
    @test_throws BoundsError intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 0)
    @test_throws BoundsError intensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 4)
end


@testset "intensity RiTIC scanindex" begin
    # Same return values as those provided as arguments to construct the object
    @test 12 == intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1)
    @test 956 == intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        2)

    # Same return container and element type as used to construct the object
    @test Int == typeof(intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), 1))
    @test Float64 == typeof(intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12.0, 956.0, 23.0]), 1))
    
    @test_throws BoundsError intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 0)
    @test_throws BoundsError intensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 4)
end


############################################################################################
# maxintensity(chrom::AbstractGC; scanindexrange)
############################################################################################
@testset "maxintensity FID; scanindexrange" begin
    # Validate the returned value
    @test 956 == maxintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 12 == maxintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]), scanindexrange=1:2)
    @test 956 == maxintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]))
    @test 12 == maxintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(maxintensity(FID([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Int == typeof(maxintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:1))
    @test Int == typeof(maxintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:2))
    @test Float64 == typeof(maxintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23])))
    @test Float64 == typeof(maxintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 
        scanindexrange=1:1))
    @test Float64 == typeof(maxintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 
        scanindexrange=1:2))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(FID([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:3)
    @test_throws BoundsError maxintensity(FID([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=1:4)
    @test_throws BoundsError maxintensity(FID([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:4)
end


@testset "maxintensity RiFID; scanindexrange" begin
    # Validate the returned value
    @test 956 == maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]))
    @test 12 == maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:2)
    @test 956 == maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]))
    @test 12 == maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23])))
    @test Int == typeof(maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:1))
    @test Int == typeof(maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:2))
    @test Float64 == typeof(maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23])))
    @test Float64 == typeof(maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:1))
    @test Float64 == typeof(maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:2))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:3)
    @test_throws BoundsError maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=1:4)
    @test_throws BoundsError maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:4)
end


@testset "maxintensity TIC; scanindexrange" begin
    # Validate the returned value
    @test 956 == maxintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 12 == maxintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]), scanindexrange=1:2)
    @test 956 == maxintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]))
    @test 12 == maxintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(maxintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Int == typeof(maxintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:1))
    @test Int == typeof(maxintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:2))
    @test Float64 == typeof(maxintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23])))
    @test Float64 == typeof(maxintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 
        scanindexrange=1:1))
    @test Float64 == typeof(maxintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 
        scanindexrange=1:2))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:3)
    @test_throws BoundsError maxintensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=1:4)
    @test_throws BoundsError maxintensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:4)
end


@testset "maxintensity RiTIC; scanindexrange" begin
    # Validate the returned value
    @test 956 == maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]))
    @test 12 == maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:2)
    @test 956 == maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]))
    @test 12 == maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:1)
    @test 956 == maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23])))
    @test Int == typeof(maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:1))
    @test Int == typeof(maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:2))
    @test Float64 == typeof(maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23])))
    @test Float64 == typeof(maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:1))
    @test Float64 == typeof(maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:2))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:3)
    @test_throws BoundsError maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=1:4)
    @test_throws BoundsError maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:4)
end


############################################################################################
# maxretentionindex(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{T, S}]) 
# where {T<:Integer, S<:Integer}
############################################################################################
@testset "maxretentionindex FID" begin
    @test_throws MethodError maxretentionindex(FID([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "maxretentionindex RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test 300 == maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test 200 == maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2)
    @test 100 == maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:1)

    # Same return same element type as used to construct the object
    @test Int == typeof(maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23])))
    @test Float64 == typeof(maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0])))
    @test Float64 == typeof(maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [missing, 200.0, 300.0], [12.0, 956.0, 23.0])))
    @test Int64 == typeof(maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 10:5:20, 
        [12.0, 956.0, 23.0])))
    @test Float64 == typeof(maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0])))
    @test Int64 == typeof(maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23]), 2:3))
    @test Float64 == typeof(maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0]), 2:3))
    @test Int64 == typeof(maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 10:5:20, 
        [12.0, 956.0, 23.0]), 2:3))
    @test Float64 == typeof(maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0]), 2:3))

    # Check for BoundsErrors
    @test_throws BoundsError maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3)
    @test_throws BoundsError maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4)

    # Check for ArgumentError when specified range contains no numercial values
    @test_throws ArgumentError maxretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, missing, missing], [12, 956, 23]), 2:3)
end


@testset "maxretentionindex GCMS" begin
    @test_throws MethodError maxretentionindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
end


@testset "maxretentionindex TIC" begin
    @test_throws MethodError maxretentionindex(TIC([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "maxretentionindex RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 300 == maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test 200 == maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2)
    @test 100 == maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:1)

    # Same return same element type as used to construct the object
    @test Int == typeof(maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23])))
    @test Float64 == typeof(maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0])))
    @test Float64 == typeof(maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [missing, 200.0, 300.0], [12.0, 956.0, 23.0])))
    @test Int64 == typeof(maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 10:5:20, 
        [12.0, 956.0, 23.0])))
    @test Float64 == typeof(maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0])))
    @test Int64 == typeof(maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23]), 2:3))
    @test Float64 == typeof(maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0]), 2:3))
    @test Int64 == typeof(maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 10:5:20, 
        [12.0, 956.0, 23.0]), 2:3))
    @test Float64 == typeof(maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0]), 2:3))

    # Check for BoundsErrors
    @test_throws BoundsError maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3)
    @test_throws BoundsError maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4)

    # Check for ArgumentError when specified range contains no numercial values
    @test_throws ArgumentError maxretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, missing, missing], [12, 956, 23]), 2:3)
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


@testset "maxscantime RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test 3u"s" == maxscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test (1/20)u"minute" ≈ maxscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute")
    @test 3 == maxscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true)
    @test 1/20 ≈ maxscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Int == typeof(maxscantime(RiFID(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
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


@testset "maxscantime RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 3u"s" == maxscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]))
    @test (1/20)u"minute" ≈ maxscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test 3 == maxscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true)
    @test 1/20 ≈ maxscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Int == typeof(maxscantime(RiGCMS(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), ustripped=true))
    @test Float64 == typeof(maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), ustripped=true))
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


@testset "maxscantime RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 3u"s" == maxscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test (1/20)u"minute" ≈ maxscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute")
    @test 3 == maxscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true)
    @test 1/20 ≈ maxscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Int == typeof(maxscantime(RiTIC(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
end


############################################################################################
# maxscantime(chrom::AbstractChromatogram, scanindexrange::OrdinalRange; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "maxscantime scanindexrange FID" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == maxscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 1:2)
    @test (1/30)u"minute" ≈ maxscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 1:2, 
        timeunit=u"minute")
    @test 2 == maxscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 1:2, ustripped=true)
    @test 1/30 ≈ maxscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 1:2, timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(FID(Int64[1, 2, 3]u"s", [12, 956, 23]), 1:2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(FID(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 1:2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:2, timeunit=u"minute"))
    @test Int == typeof(maxscantime(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 1:2, 
        ustripped=true))
    @test Float64 == typeof(maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:2, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError maxscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "maxscantime scanindexrange RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == maxscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2)
    @test (1/30)u"minute" ≈ maxscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2, timeunit=u"minute")
    @test 2 == maxscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1:2, ustripped=true)
    @test 1/30 ≈ maxscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]),
         1:2, timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiFID(Int64[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:2, timeunit=u"minute"))
    @test Int == typeof(maxscantime(RiFID(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2, ustripped=true))
    @test Float64 == typeof(maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:2, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError maxscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "maxscantime scanindexrange GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == maxscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1:2)
    @test (1/30)u"minute" ≈ maxscantime(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2, timeunit=u"minute")
    @test 2 == maxscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1:2, 
        ustripped=true)
    @test 1/30 ≈ maxscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1:2, 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), 1:2, timeunit=u"minute"))
    @test Int == typeof(maxscantime(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2, ustripped=true))
    @test Float64 == typeof(maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:2, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError maxscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "maxscantime scanindexrange RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == maxscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1:2)
    @test (1/30)u"minute" ≈ maxscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1:2, timeunit=u"minute")
    @test 2 == maxscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1:2, ustripped=true)
    @test 1/30 ≈ maxscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1:2, timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1:2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1:2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1:2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1:2, 
        timeunit=u"minute"))
    @test Int == typeof(maxscantime(RiGCMS(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1:2, ustripped=true))
    @test Float64 == typeof(maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1:2, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0:3, timeunit=u"minute", 
        ustripped=true)
    @test_throws BoundsError maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1:4, timeunit=u"minute", 
        ustripped=true)
    @test_throws BoundsError maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0:4, timeunit=u"minute", 
        ustripped=true)
    @test_throws MethodError maxscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0.0:4.0, timeunit=u"minute", 
        ustripped=true)
end


@testset "maxscantime scanindexrange TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == maxscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 1:2)
    @test (1/30)u"minute" ≈ maxscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 1:2, 
        timeunit=u"minute")
    @test 2 == maxscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 1:2, ustripped=true)
    @test 1/30 ≈ maxscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 1:2, timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(TIC(Int64[1, 2, 3]u"s", [12, 956, 23]), 1:2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(TIC(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 1:2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:2, timeunit=u"minute"))
    @test Int == typeof(maxscantime(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 1:2, 
        ustripped=true))
    @test Float64 == typeof(maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:2, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError maxscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "maxscantime scanindexrange RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == maxscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2)
    @test (1/30)u"minute" ≈ maxscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2, timeunit=u"minute")
    @test 2 == maxscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        1:2, ustripped=true)
    @test 1/30 ≈ maxscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]),
         1:2, timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiTIC(Int64[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(maxscantime(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:2, timeunit=u"minute"))
    @test Int == typeof(maxscantime(RiTIC(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:2, ustripped=true))
    @test Float64 == typeof(maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:2, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError maxscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end

############################################################################################
# scanduration(chrom::AbstractChromatogram; error::Real=0.001, 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "scanduration FID" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == scanduration(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/60)u"minute" ≈ scanduration(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 1 == scanduration(FID([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ scanduration(FID([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)
    @test 1u"s" ≈ scanduration(FID([1.0, 1.99, 3.0]u"s", [12, 956, 1]), error=0.02)

    # Same return container and element type as used to construct the object
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(FID(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(FID(Int64[1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Float64 == typeof(scanduration(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Float64 == typeof(scanduration(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        ustripped=true))

    @test_throws ArgumentError scanduration(FID([1]u"s", [12]))
    @test_throws ArgumentError scanduration(scanduration(FID([1.0, 1.99, 3.0]u"s", 
        [12, 956, 1])))
end


@testset "scanduration RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == scanduration(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test (1/60)u"minute" ≈ scanduration(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300],
        [12, 956, 23]), timeunit=u"minute")
    @test 1 == scanduration(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ scanduration(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300],
        [12, 956, 23]), timeunit=u"minute", ustripped=true)
    @test 1u"s" ≈ scanduration(RiFID([1.0, 1.99, 3.0]u"s", "Kovats", [100, 200, 300],
        [12, 956, 1]), error=0.02)

    # Same return container and element type as used to construct the object
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Float64 == typeof(scanduration(RiFID(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(scanduration(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))

    @test_throws ArgumentError scanduration(RiFID([1]u"s", "Kovats", [100], 
        [12]))
    @test_throws ArgumentError scanduration(scanduration(RiFID([1.0, 1.99, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 1])))
end


@testset "scanduration GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == scanduration(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test (1/60)u"minute" ≈ scanduration(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test 1 == scanduration(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test 1/60 ≈ scanduration(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Float64 == typeof(scanduration(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Float64 == typeof(scanduration(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))

    @test_throws ArgumentError scanduration(GCMS([1]u"s", [12], reshape([0], length(1), 
        1)))
end


@testset "scanduration RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == scanduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]))
    @test (1/60)u"minute" ≈ scanduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test 1 == scanduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true)
    @test 1/60 ≈ scanduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Float64 == typeof(scanduration(RiGCMS(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), ustripped=true))
    @test Float64 == typeof(scanduration(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), ustripped=true))

    @test_throws ArgumentError scanduration(RiGCMS([1]u"s", "Kovats", [100], [12], 
        reshape([0], length(1), 1)))
end


@testset "scanduration TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == scanduration(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/60)u"minute" ≈ scanduration(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 1 == scanduration(TIC([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ scanduration(TIC([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(TIC(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(TIC(Int64[1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Float64 == typeof(scanduration(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Float64 == typeof(scanduration(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), ustripped=true))

    @test_throws ArgumentError scanduration(TIC([1]u"s", [12]))
end


@testset "scanduration RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == scanduration(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test (1/60)u"minute" ≈ scanduration(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300],
        [12, 956, 23]), timeunit=u"minute")
    @test 1 == scanduration(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ scanduration(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300],
        [12, 956, 23]), timeunit=u"minute", ustripped=true)
    @test 1u"s" ≈ scanduration(RiTIC([1.0, 1.99, 3.0]u"s", "Kovats", [100, 200, 300],
        [12, 956, 1]), error=0.02)

    # Same return container and element type as used to construct the object
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scanduration(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Float64 == typeof(scanduration(RiTIC(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(scanduration(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))

    @test_throws ArgumentError scanduration(RiTIC([1]u"s", "Kovats", [100], 
        [12]))
    @test_throws ArgumentError scanduration(scanduration(RiTIC([1.0, 1.99, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 1])))
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


@testset "metadata RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test Dict() == metadata(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test Dict(:id => 4, :name => "sample") == metadata(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23], Dict(:id => 4, :name => "sample")))
    
    # Same return container and element type as used to construct the object
    @test Dict{Any, Any} == typeof(metadata(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23], Dict(:id => 4, :name => "sample"))))
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


@testset "metadata RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test Dict() == metadata(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]))
    @test Dict(:id => 4, :name => "sample") == metadata(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1], Dict(:id => 4, 
        :name => "sample")))

    # Same return container and element type as used to construct the object
    @test Dict{Any, Any} == typeof(metadata(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1], Dict(:id => 4, 
        :name => "sample"))))   
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


@testset "metadata RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test Dict() == metadata(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test Dict(:id => 4, :name => "sample") == metadata(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23], Dict(:id => 4, :name => "sample")))
    
    # Same return container and element type as used to construct the object
    @test Dict{Any, Any} == typeof(metadata(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23], Dict(:id => 4, :name => "sample"))))
end

############################################################################################
# minintensity(chrom::AbstractChromatogram[, greaterthan::Real; scanindexrange])
############################################################################################
@testset "minintensity FID [, greaterthan::Real; scanindexrange]" begin
    # Validate the returned value
    @test 12 == minintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 23 == minintensity(FID([1, 2, 3]u"s", Int[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(FID([1, 2, 3]u"s", Int[23, 12, 956]), scanindexrange=1:2)
    @test 12 == minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]))
    @test 23 == minintensity(FID([1, 2, 3]u"s", Float64[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(FID([1, 2, 3]u"s", Float64[23, 12, 956]), scanindexrange=1:2)

    @test 23 == minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 12)
    @test 956 == minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 12, 
        scanindexrange=1:2)

    @test nothing === minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 956)
    @test nothing === minintensity(FID([1, 2, 3]u"s", Float64[23, 12, 956]), 23, 
        scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(minintensity(FID([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Int == typeof(minintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:1))
    @test Int == typeof(minintensity(FID([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:2))
    @test Float64 == typeof(minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23])))
    @test Float64 == typeof(minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 
        scanindexrange=1:1))
    @test Float64 == typeof(minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 12, 
        scanindexrange=1:2))
    @test Nothing === typeof(minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 956))
    @test Nothing === typeof(minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 956, 
        scanindexrange=1:1))
    @test Nothing === typeof(minintensity(FID([1, 2, 3]u"s", Float64[12, 956, 23]), 956, 
        scanindexrange=1:1))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(FID([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:3)
    @test_throws BoundsError maxintensity(FID([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=1:4)
    @test_throws BoundsError maxintensity(FID([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:4)
    # "Kovats", [100, 200, 300]
end


@testset "minintensity RiFID [, greaterthan::Real; scanindexrange]" begin
    # Validate the returned value
    @test 12 == minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]))
    @test 23 == minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[23, 12, 956]), scanindexrange=1:2)
    @test 12 == minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]))
    @test 23 == minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[23, 12, 956]), scanindexrange=1:2)

    @test 23 == minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 12)
    @test 956 == minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 12, scanindexrange=1:2)

    @test nothing === minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 956)
    @test nothing === minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[23, 12, 956]), 23, scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23])))
    @test Int == typeof(minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:1))
    @test Int == typeof(minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:2))
    @test Float64 == typeof(minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23])))
    @test Float64 == typeof(minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:1))
    @test Float64 == typeof(minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 12, scanindexrange=1:2))
    @test Nothing === typeof(minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 956))
    @test Nothing === typeof(minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 956, scanindexrange=1:1))
    @test Nothing === typeof(minintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 956, scanindexrange=1:1))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:3)
    @test_throws BoundsError maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=1:4)
    @test_throws BoundsError maxintensity(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:4)
end


@testset "minintensity TIC [, greaterthan::Real; scanindexrange]" begin
    # Validate the returned value
    @test 12 == minintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]))
    @test 23 == minintensity(TIC([1, 2, 3]u"s", Int[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(TIC([1, 2, 3]u"s", Int[23, 12, 956]), scanindexrange=1:2)
    @test 12 == minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]))
    @test 23 == minintensity(TIC([1, 2, 3]u"s", Float64[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(TIC([1, 2, 3]u"s", Float64[23, 12, 956]), scanindexrange=1:2)

    @test 23 == minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 12)
    @test 956 == minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 12, 
        scanindexrange=1:2)

    @test nothing === minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 956)
    @test nothing === minintensity(TIC([1, 2, 3]u"s", Float64[23, 12, 956]), 23, 
        scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(minintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Int == typeof(minintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:1))
    @test Int == typeof(minintensity(TIC([1, 2, 3]u"s", Int[12, 956, 23]), 
        scanindexrange=1:2))
    @test Float64 == typeof(minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23])))
    @test Float64 == typeof(minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 
        scanindexrange=1:1))
    @test Float64 == typeof(minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 12, 
        scanindexrange=1:2))
    @test Nothing === typeof(minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 956))
    @test Nothing === typeof(minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 956, 
        scanindexrange=1:1))
    @test Nothing === typeof(minintensity(TIC([1, 2, 3]u"s", Float64[12, 956, 23]), 956, 
        scanindexrange=1:1))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:3)
    @test_throws BoundsError maxintensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=1:4)
    @test_throws BoundsError maxintensity(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        scanindexrange=0:4)
    # "Kovats", [100, 200, 300]
end


@testset "minintensity RiTIC [, greaterthan::Real; scanindexrange]" begin
    # Validate the returned value
    @test 12 == minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]))
    @test 23 == minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[23, 12, 956]), scanindexrange=1:2)
    @test 12 == minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]))
    @test 23 == minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[23, 12, 956]), scanindexrange=1:1)
    @test 12 == minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[23, 12, 956]), scanindexrange=1:2)

    @test 23 == minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 12)
    @test 956 == minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 12, scanindexrange=1:2)

    @test nothing === minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 956)
    @test nothing === minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[23, 12, 956]), 23, scanindexrange=1:2)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23])))
    @test Int == typeof(minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:1))
    @test Int == typeof(minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[12, 956, 23]), scanindexrange=1:2))
    @test Float64 == typeof(minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23])))
    @test Float64 == typeof(minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), scanindexrange=1:1))
    @test Float64 == typeof(minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 12, scanindexrange=1:2))
    @test Nothing === typeof(minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 956))
    @test Nothing === typeof(minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 956, scanindexrange=1:1))
    @test Nothing === typeof(minintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]), 956, scanindexrange=1:1))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:3)
    @test_throws BoundsError maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=1:4)
    @test_throws BoundsError maxintensity(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), scanindexrange=0:4)
end


############################################################################################
# minretentionindex(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange{T, S}]) 
# where {T<:Integer, S<:Integer}
############################################################################################
@testset "minretentionindex FID" begin
    @test_throws MethodError minretentionindex(FID([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "minretentionindex RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test 100 == minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test 200 == minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3)
    @test 100 == minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:1)

    # Same return same element type as used to construct the object
    @test Int == typeof(minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23])))
    @test Float64 == typeof(minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0])))
    @test Float64 == typeof(minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [missing, 200.0, 300.0], [12.0, 956.0, 23.0])))
    @test Int64 == typeof(minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 10:5:20, 
        [12.0, 956.0, 23.0])))
    @test Float64 == typeof(minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0])))
    @test Int64 == typeof(minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23]), 2:3))
    @test Float64 == typeof(minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0]), 2:3))
    @test Int64 == typeof(minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 10:5:20, 
        [12.0, 956.0, 23.0]), 2:3))
    @test Float64 == typeof(minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0]), 2:3))

    # Check for BoundsErrors
    @test_throws BoundsError minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3)
    @test_throws BoundsError minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4)

    # Check for ArgumentError when specified range contains no numercial values
    @test_throws ArgumentError minretentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, missing, missing], [12, 956, 23]), 2:3)
end


@testset "minretentionindex GCMS" begin
    @test_throws MethodError minretentionindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
end


@testset "minretentionindex TIC" begin
    @test_throws MethodError minretentionindex(TIC([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "minretentionindex RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 100 == minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test 200 == minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3)
    @test 100 == minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:1)

    # Same return same element type as used to construct the object
    @test Int == typeof(minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23])))
    @test Float64 == typeof(minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0])))
    @test Float64 == typeof(minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [missing, 200.0, 300.0], [12.0, 956.0, 23.0])))
    @test Int64 == typeof(minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 10:5:20, 
        [12.0, 956.0, 23.0])))
    @test Float64 == typeof(minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0])))
    @test Int64 == typeof(minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23]), 2:3))
    @test Float64 == typeof(minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0]), 2:3))
    @test Int64 == typeof(minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 10:5:20, 
        [12.0, 956.0, 23.0]), 2:3))
    @test Float64 == typeof(minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0]), 2:3))

    # Check for BoundsErrors
    @test_throws BoundsError minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3)
    @test_throws BoundsError minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4)

    # Check for ArgumentError when specified range contains no numercial values
    @test_throws ArgumentError minretentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, missing, missing], [12, 956, 23]), 2:3)
end


############################################################################################
# minscantime(chrom::AbstractChromatogram[, scanindexrange::OrdinalRange]; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "minscantime FID" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/60)u"minute" ≈ minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 1 == minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    @test 2u"s" == minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 2:3)
    @test (1/30)u"minute" ≈ minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute")
    @test 2 == minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 2:3, ustripped=true)
    @test 1/30 ≈ minscantime(FID([1, 2, 3]u"s", [12, 956, 23]), 2:3, timeunit=u"minute", 
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

    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(FID(Int64[1, 2, 3]u"s", [12, 956, 23]), 2:3))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(FID(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3, timeunit=u"minute"))
    @test Int == typeof(minscantime(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        ustripped=true))
    @test Float64 == typeof(minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError minscantime(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "minscantime RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test (1/60)u"minute" ≈ minscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute")
    @test 1 == minscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        ustripped=true)
    @test 1/60 ≈ minscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)

    @test 2u"s" == minscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3)
    @test (1/30)u"minute" ≈ minscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, timeunit=u"minute")
    @test 2 == minscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        2:3, ustripped=true)
    @test 1/30 ≈ minscantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Int == typeof(minscantime(RiFID(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))

    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(RiFID(Int64[1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Int == typeof(minscantime(RiFID(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, ustripped=true))
    @test Float64 == typeof(minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300],  [12, 956, 23]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError minscantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300],  [12, 956, 23]), 0.0:4.0, timeunit=u"minute", ustripped=true)
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

    @test 2u"s" == minscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3)
    @test (1/30)u"minute" ≈ minscantime(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute")
    @test 2 == minscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3, 
        ustripped=true)
    @test 1/30 ≈ minscantime(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3, 
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

    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Int == typeof(minscantime(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true))
    @test Float64 == typeof(minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError minscantime(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "minscantime RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]))
    @test (1/60)u"minute" ≈ minscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test 1 == minscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true)
    @test 1/60 ≈ minscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true)

    @test 2u"s" == minscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2:3)
    @test (1/30)u"minute" ≈ minscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute")
    @test 2 == minscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true)
    @test 1/30 ≈ minscantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(RiGCMS(Int64[1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Int == typeof(minscantime(RiGCMS(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), ustripped=true))
    @test Float64 == typeof(minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), ustripped=true))

    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Int == typeof(minscantime(RiGCMS(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2:3, ustripped=true))
    @test Float64 == typeof(minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0:3, timeunit=u"minute", 
        ustripped=true)
    @test_throws BoundsError minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1:4, timeunit=u"minute", 
        ustripped=true)
    @test_throws BoundsError minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0:4, timeunit=u"minute", 
        ustripped=true)
    @test_throws MethodError minscantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0.0:4.0, timeunit=u"minute", 
        ustripped=true)
end


@testset "minscantime TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test (1/60)u"minute" ≈ minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test 1 == minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test 1/60 ≈ minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), timeunit=u"minute", 
        ustripped=true)

    @test 2u"s" == minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 2:3)
    @test (1/30)u"minute" ≈ minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute")
    @test 2 == minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 2:3, ustripped=true)
    @test 1/30 ≈ minscantime(TIC([1, 2, 3]u"s", [12, 956, 23]), 2:3, timeunit=u"minute", 
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

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(TIC(Int64[1, 2, 3]u"s", [12, 956, 23]), 2:3))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(TIC(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3, timeunit=u"minute"))
    @test Int == typeof(minscantime(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        ustripped=true))
    @test Float64 == typeof(minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError minscantime(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "minscantime RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 1u"s" == minscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test (1/60)u"minute" ≈ minscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute")
    @test 1 == minscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        ustripped=true)
    @test 1/60 ≈ minscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)
    @test 2u"s" == minscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3)
    @test (1/30)u"minute" ≈ minscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, timeunit=u"minute")
    @test 2 == minscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        2:3, ustripped=true)
    @test 1/30 ≈ minscantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))

    @test Int == typeof(minscantime(RiTIC(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))

    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(minscantime(RiTIC(Int64[1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Int == typeof(minscantime(RiTIC(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, ustripped=true))
    @test Float64 == typeof(minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300],  [12, 956, 23]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError minscantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300],  [12, 956, 23]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end


############################################################################################
# retentionindex(chrom::AbstractChromatogram, scanindex::Integer)
############################################################################################
@testset "retentionindex FID scanindex" begin
    @test_throws MethodError retentionindex(FID([1, 2, 3]u"s", [12, 956, 23]), 1)
end


@testset "retentionindex RiFID scanindex" begin
    @test 100 == retentionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, missing], 
        [12, 956, 23]), 1)
    @test missing === retentionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, missing], 
        [12, 956, 23]), 3)

    @test Int == typeof(retentionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1))
    @test Float64 == typeof(retentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100.0, 200.0, missing], [12, 956, 23]), 1))
    @test Missing == typeof(retentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100.0, 200.0, missing], [12, 956, 23]), 3))

    @test_throws BoundsError retentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, missing], [12, 956, 23]), 0)
    @test_throws BoundsError retentionindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, missing], [12, 956, 23]), 4)
end


@testset "retentionindex GCMS scanindex" begin
    @test_throws MethodError retentionindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1)
end


@testset "retentionindex RiGCMS scanindex" begin
    @test 100 == retentionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, missing], 
        [85, 100], [0 12; 34 956; 23 1]), 1)
    @test missing === retentionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, missing], 
        [85, 100], [0 12; 34 956; 23 1]), 3)

    @test Int == typeof(retentionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1))
    @test Float64 == typeof(retentionindex(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100.0, 200.0, missing], [85, 100], [0 12; 34 956; 23 1]), 1))
    @test Missing == typeof(retentionindex(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100.0, 200.0, missing], [85, 100], [0 12; 34 956; 23 1]), 3))

    @test_throws BoundsError retentionindex(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, missing], [85, 100], [0 12; 34 956; 23 1]), 0)
    @test_throws BoundsError retentionindex(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, missing], [85, 100], [0 12; 34 956; 23 1]), 4)
end


@testset "retentionindex TIC scanindex" begin
    @test_throws MethodError retentionindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 1)
end


@testset "retentionindex RiTIC scanindex" begin
    @test 100 == retentionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, missing], 
        [12, 956, 23]), 1)
    @test missing === retentionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, missing], 
        [12, 956, 23]), 3)

    @test Int == typeof(retentionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1))
    @test Float64 == typeof(retentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100.0, 200.0, missing], [12, 956, 23]), 1))
    @test Missing == typeof(retentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100.0, 200.0, missing], [12, 956, 23]), 3))

    @test_throws BoundsError retentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, missing], [12, 956, 23]), 0)
    @test_throws BoundsError retentionindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, missing], [12, 956, 23]), 4)
end


############################################################################################
# retentionindexname(chrom::AbstractChromatogram)
############################################################################################
@testset "retentionindexname FID" begin
    @test_throws MethodError retentionindexname(FID([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "retentionindexname RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test "Kovats" == retentionindexname(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))

    # Same return same element type as used to construct the object
    @test String == typeof(retentionindexname(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test SubString{String} == typeof(retentionindexname(RiFID([1, 2, 3]u"s", 
        SubString("Kovats", 1:6), [100, 200, 300], [12.0, 956.0, 23.0])))
end


@testset "retentionindexname GCMS" begin
    @test_throws MethodError retentionindexname(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
end


@testset "retentionindexname TIC" begin
    @test_throws MethodError retentionindexname(TIC([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "retentionindexname RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test "Kovats" == retentionindexname(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))

    # Same return same element type as used to construct the object
    @test String == typeof(retentionindexname(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test SubString{String} == typeof(retentionindexname(RiTIC([1, 2, 3]u"s", 
        SubString("Kovats", 1:6), [100, 200, 300], [12.0, 956.0, 23.0])))
end


############################################################################################
# retentionindices(chrom::AbstractChromatogram)
############################################################################################
@testset "retentionindices FID" begin
    @test_throws MethodError retentionindices(FID([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "retentionindices RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test [100, 200, 300] == retentionindices(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]))
    @test [200, 300] == retentionindices(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3)
    @test [100] == retentionindices(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:1)

    # Same return same element type as used to construct the object
    @test Vector{Int} == typeof(retentionindices(RiFID([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23])))
    @test Vector{Float64} == typeof(retentionindices(RiFID([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0])))
    @test Vector{Union{Float64, Missing}} == typeof(retentionindices(RiFID([1, 2, 3]u"s", "Kovats", 
        [missing, 200.0, 300.0], [12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(retentionindices(RiFID([1, 2, 3]u"s", 
        "Kovats", 10:5:20, [12.0, 956.0, 23.0])))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(retentionindices(RiFID([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0])))
    @test SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true} == typeof(
        retentionindices(RiFID([1, 2, 3]u"s", "Kovats", Int[100, 200, 300], [12, 956, 23]), 
        2:3))
    @test SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} == typeof(
        retentionindices(RiFID([1, 2, 3]u"s", "Kovats", Float64[100, 200, 300], 
        [12.0, 956.0, 23.0]), 2:3))
    @test StepRange{Int64, Int64} == typeof(retentionindices(RiFID([1, 2, 3]u"s", "Kovats", 
        10:5:20, [12.0, 956.0, 23.0]), 2:3))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(retentionindices(RiFID([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0]), 2:3))

    # Check for BoundsErrors
    @test_throws BoundsError retentionindices(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3)
    @test_throws BoundsError retentionindices(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4)
end


@testset "retentionindices GCMS" begin
    @test_throws MethodError retentionindices(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
end


@testset "retentionindices TIC" begin
    @test_throws MethodError retentionindices(TIC([1, 2, 3]u"s", [12, 956, 23]))
end

@testset "retentionindices RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test [100, 200, 300] == retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]))
    @test [200, 300] == retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3)
    @test [100] == retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1:1)

    # Same return same element type as used to construct the object
    @test Vector{Int} == typeof(retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", 
        Int[100, 200, 300], [12, 956, 23])))
    @test Vector{Float64} == typeof(retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", 
        Float64[100, 200, 300], [12.0, 956.0, 23.0])))
    @test Vector{Union{Float64, Missing}} == typeof(retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", 
        [missing, 200.0, 300.0], [12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(retentionindices(RiTIC([1, 2, 3]u"s", 
        "Kovats", 10:5:20, [12.0, 956.0, 23.0])))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0])))
    @test SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true} == typeof(
        retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", Int[100, 200, 300], [12, 956, 23]), 
        2:3))
    @test SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} == typeof(
        retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", Float64[100, 200, 300], 
        [12.0, 956.0, 23.0]), 2:3))
    @test StepRange{Int64, Int64} == typeof(retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", 
        10:5:20, [12.0, 956.0, 23.0]), 2:3))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", 10.0:5:20.0, 
        [12.0, 956.0, 23.0]), 2:3))

    # Check for BoundsErrors
    @test_throws BoundsError retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3)
    @test_throws BoundsError retentionindices(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4)
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


@testset "runduration RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == runduration(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test (1//30)u"minute" == runduration(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute")
    @test (1/30)u"minute" ≈ runduration(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute")
    @test 2 == runduration(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        ustripped=true)
    @test 1//30 == runduration(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)
    @test 1/30 ≈ runduration(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(runduration(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Int == typeof(runduration(RiFID(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(runduration(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
    @test Rational{Int64} == typeof(runduration(RiFID(Int[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(runduration(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true))
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


@testset "runduration RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == runduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]))
    @test (1//30)u"minute" == runduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test (1/30)u"minute" ≈ runduration(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test 2 == runduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true)
    @test 1//30 == runduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true)
    @test 1/30 ≈ runduration(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(runduration(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Int == typeof(runduration(RiGCMS(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), ustripped=true))
    @test Float64 == typeof(runduration(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), ustripped=true))
    @test Rational{Int64} == typeof(runduration(RiGCMS(Int[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute", 
        ustripped=true))
    @test Float64 == typeof(runduration(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute", 
        ustripped=true))
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


@testset "runduration RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == runduration(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test (1//30)u"minute" == runduration(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute")
    @test (1/30)u"minute" ≈ runduration(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute")
    @test 2 == runduration(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        ustripped=true)
    @test 1//30 == runduration(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)
    @test 1/30 ≈ runduration(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)

    # Same return container and element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(runduration(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(runduration(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Int == typeof(runduration(RiTIC(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true))
    @test Float64 == typeof(runduration(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
    @test Rational{Int64} == typeof(runduration(RiTIC(Int[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(runduration(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true))
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


@testset "scancount RiFID" begin
    # Validate the returned value
    @test 3 == scancount(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]))
    @test 1 == scancount(RiFID([1]u"s", "Kovats", [100], [12]))

    # Return value must be an integer
    @test Int == typeof(scancount(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
    [12, 956, 23])))
end


@testset "scancount GCMS" begin
    # Validate the returned value
    @test 3 == scancount(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test 1 == scancount(GCMS([1]u"s", [85, 100], [0 12]))

    # Return value must be an integer
    @test Int == typeof(scancount(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])))
end


@testset "scancount RiGCMS" begin
    # Validate the returned value
    @test 3 == scancount(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]))
    @test 1 == scancount(RiGCMS([1]u"s", "Kovats", [100], [85, 100], [0 12]))

    # Return value must be an integer
    @test Int == typeof(scancount(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1])))
end


@testset "scancount TIC" begin
    # Validate the returned value
    @test 3 == scancount(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test 1 == scancount(TIC([1]u"s", [12]))

    # Return value must be an integer
    @test Int == typeof(scancount(TIC([1, 2, 3]u"s", [12, 956, 23])))
end


@testset "scancount RiTIC" begin
    # Validate the returned value
    @test 3 == scancount(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]))
    @test 1 == scancount(RiTIC([1]u"s", "Kovats", [100], [12]))

    # Return value must be an integer
    @test Int == typeof(scancount(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
    [12, 956, 23])))
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


@testset "scantime RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == scantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2)
    @test (1//20)u"minute" == scantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3, timeunit=u"minute")
    @test (1/20)u"minute" ≈ scantime(RiFID([1.0, 2.0, 3.0]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3, timeunit=u"minute")
    @test 2 == scantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 2, 
        ustripped=true)
    @test 1//20 == scantime(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        3, timeunit=u"minute", ustripped=true)
    @test 1/20 ≈ scantime(RiFID([1.0, 2.0, 3.0]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3, timeunit=u"minute", ustripped=true)

    # Same return element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}} == typeof(scantime(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}} == typeof(scantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}} == typeof(scantime(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(scantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2, timeunit=u"minute"))
    @test Int == typeof(scantime(RiFID(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2, ustripped=true))
    @test Rational{Int64} == typeof(scantime(RiFID(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2, timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(scantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300],[12, 956, 23]), 2, ustripped=true))
    @test Float64 == typeof(scantime(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300],[12, 956, 23]), 2, timeunit=u"minute", ustripped=true))
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


@testset "scantime RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == scantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2)
    @test (1//20)u"minute" == scantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 3, timeunit=u"minute")
    @test (1/20)u"minute" ≈ scantime(RiGCMS([1.0, 2.0, 3.0]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 3, timeunit=u"minute")
    @test 2 == scantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2, ustripped=true)
    @test 1//20 == scantime(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 3, timeunit=u"minute", ustripped=true)
    @test 1/20 ≈ scantime(RiGCMS([1.0, 2.0, 3.0]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 3, timeunit=u"minute", ustripped=true)

    # Same return element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scantime(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(scantime(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(scantime(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2, timeunit=u"minute"))
    @test Int == typeof(scantime(RiGCMS(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2, ustripped=true))
    @test Rational{Int64} == typeof(scantime(RiGCMS(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2, timeunit=u"minute", 
        ustripped=true))
    @test Float64 == typeof(scantime(RiGCMS(Float64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2, ustripped=true))
    @test Float64 == typeof(scantime(RiGCMS(Float64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2, timeunit=u"minute", 
        ustripped=true))
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


@testset "scantime RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test 2u"s" == scantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2)
    @test (1//20)u"minute" == scantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3, timeunit=u"minute")
    @test (1/20)u"minute" ≈ scantime(RiTIC([1.0, 2.0, 3.0]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3, timeunit=u"minute")
    @test 2 == scantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 2, 
        ustripped=true)
    @test 1//20 == scantime(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        3, timeunit=u"minute", ustripped=true)
    @test 1/20 ≈ scantime(RiTIC([1.0, 2.0, 3.0]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3, timeunit=u"minute", ustripped=true)

    # Same return element type as used to construct the object
    @test Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}} == typeof(scantime(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}} == typeof(scantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2))
    @test Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}} == typeof(scantime(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2, timeunit=u"minute"))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(scantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2, timeunit=u"minute"))
    @test Int == typeof(scantime(RiTIC(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2, ustripped=true))
    @test Rational{Int64} == typeof(scantime(RiTIC(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2, timeunit=u"minute", ustripped=true))
    @test Float64 == typeof(scantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300],[12, 956, 23]), 2, ustripped=true))
    @test Float64 == typeof(scantime(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300],[12, 956, 23]), 2, timeunit=u"minute", ustripped=true))
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


@testset "scantimeindex RiFID" begin
    # Validate the returned value
    @test 1 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s")
    @test 2 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2u"s")
    @test 3 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3u"s")
    @test 1 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=true)
    @test 2 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2u"s", precisetime=true)
    @test 3 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3u"s", precisetime=true)
    @test 1 == scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(RiFID(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(RiFID(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(RiFID(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(RiFID(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(RiFID(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(RiFID(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 0u"s")
    @test 3 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 4u"s")
    
    # Function throws ArgumentError if precisetime is set to true and the specified time 
    # does not does not exist
    @test_throws ArgumentError scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 4u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2.1001u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 1.9999u"s", precisetime=true)

    # Return value must be an integer
    @test Int == typeof(scantimeindex(RiFID(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s"))
    @test Int == typeof(scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1.1u"s"))
    @test Int == typeof(scantimeindex(RiFID(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=true))
    @test Int == typeof(scantimeindex(RiFID(Float64[1.1, 2.1, 3.1]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1.1u"s", precisetime=true))

    # The function argument time must be a quantity with a valid time unit
    @test 1 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"ms")
    @test 3 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"minute")
    @test_throws MethodError scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"m")
    @test_throws MethodError scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"m", precisetime=true)

    # The function argument precisetime must be a Boolean value
    @test 1 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=true)
    @test 1 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=false)
    @test 1 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=:true)
    @test 1 == scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=:false)
    @test_throws TypeError scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=:truly)
    @test_throws TypeError scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime="true")
    @test_throws TypeError scantimeindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=1)
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


@testset "scantimeindex RiGCMS" begin
    # Validate the returned value
    @test 1 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s")
    @test 2 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s")
    @test 3 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 3u"s")
    @test 1 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=true)
    @test 2 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s", precisetime=true)
    @test 3 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 3u"s", precisetime=true)
    @test 1 == scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1.1u"s")
    @test 2 == scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2.1u"s")
    @test 3 == scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 3.1u"s")
    @test 1 == scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(RiGCMS(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1.1u"s")
    @test 2 == scantimeindex(RiGCMS(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2.1u"s")
    @test 3 == scantimeindex(RiGCMS(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 3.1u"s")
    @test 1 == scantimeindex(RiGCMS(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(RiGCMS(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(RiGCMS(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 0u"s")
    @test 3 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 4u"s")
    
    # Function throws ArgumentError if precisetime is set to true and the specified time 
    # does not does not exist
    @test_throws ArgumentError scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 4u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2.1001u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1.9999u"s", precisetime=true)

    # Return value must be an integer
    @test Int == typeof(scantimeindex(RiGCMS(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1u"s"))
    @test Int == typeof(scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1.1u"s"))
    @test Int == typeof(scantimeindex(RiGCMS(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1u"s", precisetime=true))
    @test Int == typeof(scantimeindex(RiGCMS(Float64[1.1, 2.1, 3.1]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1.1u"s", precisetime=true))

    # The function argument time must be a quantity with a valid time unit
    @test 1 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1u"ms")
    @test 3 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1u"minute")
    @test_throws MethodError scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1u"m")
    @test_throws MethodError scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1u"m", precisetime=true)

    # The function argument precisetime must be a Boolean value
    @test 1 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=true)
    @test 1 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=false)
    @test 1 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=:true)
    @test 1 == scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1u"s", precisetime=:false)
    @test_throws TypeError scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1u"s", precisetime=:truly)
    @test_throws TypeError scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1u"s", precisetime="true")
    @test_throws TypeError scantimeindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1u"s", precisetime=1)
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


@testset "scantimeindex RiTIC" begin
    # Validate the returned value
    @test 1 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s")
    @test 2 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2u"s")
    @test 3 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3u"s")
    @test 1 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=true)
    @test 2 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2u"s", precisetime=true)
    @test 3 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3u"s", precisetime=true)
    @test 1 == scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(RiTIC(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.1u"s")
    @test 2 == scantimeindex(RiTIC(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2.1u"s")
    @test 3 == scantimeindex(RiTIC(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3.1u"s")
    @test 1 == scantimeindex(RiTIC(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1.1u"s", precisetime=true)
    @test 2 == scantimeindex(RiTIC(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2.1u"s", precisetime=true)
    @test 3 == scantimeindex(RiTIC(Float32[1.1, 2.1, 3.1]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 3.1u"s", precisetime=true)
    @test 1 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 0u"s")
    @test 3 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 4u"s")
    
    # Function throws ArgumentError if precisetime is set to true and the specified time 
    # does not does not exist
    @test_throws ArgumentError scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 4u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2.1001u"s", precisetime=true)
    @test_throws ArgumentError scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 1.9999u"s", precisetime=true)

    # Return value must be an integer
    @test Int == typeof(scantimeindex(RiTIC(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s"))
    @test Int == typeof(scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1.1u"s"))
    @test Int == typeof(scantimeindex(RiTIC(Int[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=true))
    @test Int == typeof(scantimeindex(RiTIC(Float64[1.1, 2.1, 3.1]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1.1u"s", precisetime=true))

    # The function argument time must be a quantity with a valid time unit
    @test 1 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"ms")
    @test 3 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"minute")
    @test_throws MethodError scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"m")
    @test_throws MethodError scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"m", precisetime=true)

    # The function argument precisetime must be a Boolean value
    @test 1 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=true)
    @test 1 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=false)
    @test 1 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=:true)
    @test 1 == scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=:false)
    @test_throws TypeError scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=:truly)
    @test_throws TypeError scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime="true")
    @test_throws TypeError scantimeindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1u"s", precisetime=1)
end

############################################################################################
# scantimes(chrom::AbstractChromatogram; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
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
    
    # Check if reference to data structure is returned
    fid = FID([1, 2, 3]u"s", [12, 956, 23])
    @test scantimes(fid) === scantimes(fid)
    
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


@testset "scantimes RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test [1, 2, 3]u"s" == scantimes(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test [1//60, 1//30, 1//20]u"minute" == scantimes(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute")
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute")
    @test [1, 2, 3] == scantimes(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true)
    @test [1//60, 1//30, 1//20] == scantimes(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)
    
    # Check if reference to data structure is returned
    fid = RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23])
    @test scantimes(fid) === scantimes(fid)
    
    # Same return container and element type as used to construct the object
    @test Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}}} == typeof(scantimes(RiFID(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}}} == typeof(scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(RiFID(Int64[1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(RiFID(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(RiFID(Int[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true))
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

    # Check if reference to data structure is returned
    gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    @test scantimes(gcms) === scantimes(gcms)
    
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


@testset "scantimes RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [1, 2, 3]u"s" == scantimes(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]))
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test [1, 2, 3] == scantimes(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true)

    # Check if reference to data structure is returned
    gcms = RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])
    @test scantimes(gcms) === scantimes(gcms)
    
    # Same return container and element type as used to construct the object
    @test Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}}} == typeof(scantimes(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}}} == typeof(scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(RiGCMS(Int64[1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(RiGCMS(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(RiGCMS(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true))
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

    # Check if reference to data structure is returned
    tic = TIC([1, 2, 3]u"s", [12, 956, 23])
    @test scantimes(tic) === scantimes(tic)

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


@testset "scantimes RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test [1, 2, 3]u"s" == scantimes(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test [1//60, 1//30, 1//20]u"minute" == scantimes(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute")
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute")
    @test [1, 2, 3] == scantimes(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), ustripped=true)
    @test [1//60, 1//30, 1//20] == scantimes(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), timeunit=u"minute", ustripped=true)
    
    # Check if reference to data structure is returned
    fid = RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23])
    @test scantimes(fid) === scantimes(fid)
    
    # Same return container and element type as used to construct the object
    @test Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓, nothing}}} == typeof(scantimes(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓, nothing}}} == typeof(scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23])))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(RiTIC(Int64[1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(RiTIC(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(RiTIC(Int[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), timeunit=u"minute", ustripped=true))
end


############################################################################################
# scantimes(chrom::AbstractChromatogram, scanindexrange::OrdinalRange; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "scantimes with scanindexrange FID" begin
    # Same return values as those provided as arguments to construct the object
    @test [2, 3]u"s" == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 2:3)
    @test [1//30, 1//20]u"minute" == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute")
    @test [1/30, 1/20]u"minute" ≈ scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3, timeunit=u"minute")
    @test [2, 3] == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 2:3, ustripped=true)
    @test [1//30, 1//20] == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute", ustripped=true)
    @test [1/30, 1/20] ≈ scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute", ustripped=true)
    
    # Check if reference to data structure is returned
    fid = FID([1, 2, 3]u"s", [12, 956, 23])
    @test scantimes(fid, 2:3) === scantimes(fid, 2:3)
    
    # Same return container and element type as used to construct the object
    @test SubArray{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(FID(Int64[1, 2, 3]u"s", [12, 956, 23]), 2:3))
    @test SubArray{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(FID(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(FID(Int[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 0:3, 
        timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 1:4, 
        timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 0:4, 
        timeunit=u"minute", ustripped=true)
    @test_throws MethodError scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "scantimes with scanindexrange RiFID" begin
    # Same return values as those provided as arguments to construct the object
    @test [2, 3]u"s" == scantimes(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3)
    @test [1//30, 1//20]u"minute" == scantimes(RiFID([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute")
    @test [1/30, 1/20]u"minute" ≈ scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute")
    @test [2, 3] == scantimes(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, ustripped=true)
    @test [1//30, 1//20] == scantimes(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true)
    @test [1/30, 1/20] ≈ scantimes(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true)
    
    # Check if reference to data structure is returned
    fid = RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23])
    @test scantimes(fid, 2:3) === scantimes(fid, 2:3)
    
    # Same return container and element type as used to construct the object
    @test SubArray{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(RiFID(Int64[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3))
    @test SubArray{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(RiFID(Int64[1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(RiFID(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(RiFID(Int[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute", 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute", 
        ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError scantimes(RiFID(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "scantimes with scanindexrange GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [2, 3]u"s" == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3)
    @test [1/30, 1/20]u"minute" ≈ scantimes(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute")
    @test [2, 3] == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2:3, 
        ustripped=true)
    @test [1/30, 1/20] ≈ scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        2:3, timeunit=u"minute", ustripped=true)

    # Check if reference to data structure is returned
    gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    @test scantimes(gcms, 2:3) === scantimes(gcms, 2:3)
    
    # Same return container and element type as used to construct the object
    @test SubArray{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3))
    @test SubArray{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(GCMS(Int64[1, 2, 3]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true))
    @test Vector{Float64} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute", ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end


@testset "scantimes with scanindexrange RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [2, 3]u"s" == scantimes(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2:3)
    @test [1/30, 1/20]u"minute" ≈ scantimes(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute")
    @test [2, 3] == scantimes(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2:3, ustripped=true)
    @test [1/30, 1/20] ≈ scantimes(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute", ustripped=true)

    # Check if reference to data structure is returned
    gcms = RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])
    @test scantimes(gcms, 2:3) === scantimes(gcms, 2:3)
    
    # Same return container and element type as used to construct the object
    @test SubArray{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(RiGCMS(Int64[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2:3))
    @test SubArray{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(RiGCMS(Int64[1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, 
        timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(RiGCMS(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(RiGCMS(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute", 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 2:3, timeunit=u"minute", 
        ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0:3, timeunit=u"minute", 
        ustripped=true)
    @test_throws BoundsError scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 1:4, timeunit=u"minute", 
        ustripped=true)
    @test_throws BoundsError scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0:4, timeunit=u"minute", 
        ustripped=true)
    @test_throws MethodError scantimes(RiGCMS(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), 0.0:4.0, timeunit=u"minute", 
        ustripped=true)
end


@testset "scantimes with scanindexrange TIC" begin
    # Same return values as those provided as arguments to construct the object
    @test [2, 3]u"s" == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 2:3)
    @test [1//30, 1//20]u"minute" == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute")
    @test [1/30, 1/20]u"minute" ≈ scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        2:3, timeunit=u"minute")
    @test [2, 3] == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 2:3, ustripped=true)
    @test [1//30, 1//20] == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute", ustripped=true)
    @test [1/30, 1/20] ≈ scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        timeunit=u"minute", ustripped=true)
    
    # Check if reference to data structure is returned
    tic = TIC([1, 2, 3]u"s", [12, 956, 23])
    @test scantimes(tic, 2:3) === scantimes(tic, 2:3)
    
    # Same return container and element type as used to construct the object
    @test SubArray{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true}  == typeof(scantimes(TIC(Int64[1, 2, 3]u"s", [12, 956, 23]), 2:3))
    @test SubArray{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(TIC(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 2:3, 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(TIC(Int[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true))
    @test Vector{Float64} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 0:3, 
        timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 1:4, 
        timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 0:4, 
        timeunit=u"minute", ustripped=true)
    @test_throws MethodError scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", [12, 956, 23]), 
        0.0:4.0, timeunit=u"minute", ustripped=true)
end

@testset "scantimes with scanindexrange RiTIC" begin
    # Same return values as those provided as arguments to construct the object
    @test [2, 3]u"s" == scantimes(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3)
    @test [1//30, 1//20]u"minute" == scantimes(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute")
    @test [1/30, 1/20]u"minute" ≈ scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute")
    @test [2, 3] == scantimes(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, ustripped=true)
    @test [1//30, 1//20] == scantimes(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true)
    @test [1/30, 1/20] ≈ scantimes(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3, timeunit=u"minute", ustripped=true)
    
    # Check if reference to data structure is returned
    fid = RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23])
    @test scantimes(fid, 2:3) === scantimes(fid, 2:3)
    
    # Same return container and element type as used to construct the object
    @test SubArray{Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(RiTIC(Int64[1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 2:3))
    @test SubArray{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}, 1, Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{
        (Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, nothing}}}, Tuple{UnitRange{Int64}}, 
        true} == typeof(scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3))
    @test Vector{Quantity{Rational{Int64}, 𝐓, Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓, nothing}}} == typeof(scantimes(RiTIC(Int64[1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(RiTIC(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 2:3, ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, ustripped=true))
    @test Vector{Rational{Int64}} == typeof(scantimes(RiTIC(Int[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute", 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", 
        "Kovats", [100, 200, 300], [12, 956, 23]), 2:3, timeunit=u"minute", 
        ustripped=true))

    # Check for BoundErrors and MethodErrors
    @test_throws BoundsError scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:3, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 1:4, timeunit=u"minute", ustripped=true)
    @test_throws BoundsError scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0:4, timeunit=u"minute", ustripped=true)
    @test_throws MethodError scantimes(RiTIC(Float64[1.0, 2.0, 3.0]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), 0.0:4.0, timeunit=u"minute", ustripped=true)
end