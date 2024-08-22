using JuChrom
using Test
using Unitful: 𝐓


############################################################################################
# binions(gcms::AbstractGCMS; ionbin::Function=integerion)
############################################################################################
@testset "binions GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test (1:3)u"s" == binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])).scantimes
    @test [85, 101] == binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])).ions
    @test [24 12; 0 956; 23 1] == binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])).intensities
    @test [24.1 12.2; 1.0 956.7; 23.9 1.5] ≈ binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5])).intensities
    @test Dict(:id => 4) == binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1], Dict(:id => 4))).metadata
    @test (1:3)u"s" == binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1]), ionbin=ion->integer(ion, start=0.9)).scantimes
    @test [84, 85, 101] == binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1]), ionbin=ion->integer(ion, start=0.9)).ions
    @test [0 24 12; 0 0 956; 23 0 1] == binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1]), ionbin=ion->integer(ion, start=0.9)).intensities
    @test [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5] ≈ binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5]), 
        ionbin=ion->integer(ion, start=0.9)).intensities
    @test Dict(:id => 4) == binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1], Dict(:id => 4)), 
        ionbin=ion->integer(ion, start=0.9)).metadata

    # Check the returned type and supertypes
    @test isa(binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])), 
        GCMS)
    @test isa(binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])), 
        AbstractGCMS)
    @test isa(binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])), 
        AbstractChromatogram)

    # Binned intensities have the same type as used to construct the object
    @test Matrix{Int64} == typeof(binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])).intensities)
    @test Matrix{Float64} == typeof(binions(GCMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5])).intensities)

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError binions(FID([1, 2, 3]u"s", [12, 956, 23]), 1)
    @test_throws MethodError binions(TIC([1, 2, 3]u"s", [12, 956, 23]), 1)
end


@testset "binions RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test (1:3)u"s" == binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], 
        [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])).scantimes
    @test [85, 101] == binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], 
        [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])).ions
    @test [24 12; 0 956; 23 1] == binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], 
        [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])).intensities
    @test [24.1 12.2; 1.0 956.7; 23.9 1.5] ≈ binions(RiGCMS((1:3)u"s", "Kovats", 
        [100, 200, 300], [84.8, 85.2, 100.9], 
        [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5])).intensities
    @test Dict(:id => 4) == binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], 
        [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1], Dict(:id => 4))).metadata
    @test (1:3)u"s" == binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], 
        [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1]), 
        ionbin=ion->integer(ion, start=0.9)).scantimes
    @test [84, 85, 101] == binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], 
        [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1]), 
        ionbin=ion->integer(ion, start=0.9)).ions
    @test [0 24 12; 0 0 956; 23 0 1] == binions(RiGCMS((1:3)u"s", "Kovats", 
        [100, 200, 300], [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1]), 
        ionbin=ion->integer(ion, start=0.9)).intensities
    @test [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5] ≈ binions(RiGCMS((1:3)u"s", 
        "Kovats", [100, 200, 300], [84.8, 85.2, 100.9], 
        [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5]), 
        ionbin=ion->integer(ion, start=0.9)).intensities
    @test Dict(:id => 4) == binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], 
        [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1], Dict(:id => 4)), 
        ionbin=ion->integer(ion, start=0.9)).metadata

    # Check the returned type and supertypes
    @test isa(binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])), RiGCMS)
    @test isa(binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])), AbstractGCMS)
    @test isa(binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])), AbstractChromatogram)

    # Binned intensities have the same type as used to construct the object
    @test Matrix{Int64} == typeof(binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], 
        [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])).intensities)
    @test Matrix{Float64} == typeof(binions(RiGCMS((1:3)u"s", "Kovats", [100, 200, 300], 
        [84.8, 85.2, 100.9], [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5])).intensities)

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError binions(RiFID([1, 2, 3]u"s", [12, 956, 23]), 1)
    @test_throws MethodError binions(RiTIC([1, 2, 3]u"s", [12, 956, 23]), 1)
end


############################################################################################
# ion(gcms::AbstractGCMS, index::Integer)
############################################################################################
@testset "ion GCMS index" begin
    # Same return values as those provided as arguments to construct the object
    @test 85 == ion(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1)
    @test 100 == ion(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2)
    
    # Provoke a BoundsError by specifying an index that does not exist
    @test_throws BoundsError ion(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 3)

    # Same return container and element type as used to construct the object
    @test Int64 == typeof(ion(GCMS([1, 2, 3]u"s", Int64[85, 100], [0 12; 34 956; 23 1]), 1))
    @test Float64 == typeof(ion(GCMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1]), 1))
    @test Int64 == typeof(ion(GCMS([1, 2, 3]u"s", 85:86, [0 12; 34 956; 23 1]), 1))
    @test Int64 == typeof(ion(GCMS([1, 2, 3]u"s", 85:85, reshape([0, 956, 1], (:,1))), 1))
    @test Int64 == typeof(ion(GCMS([1, 2, 3]u"s", 85:15:100, [0 12; 34 956; 23 1]), 1))
    @test Float64 == typeof(ion(GCMS([1, 2, 3]u"s", 85.0:15:100.0, [0 12; 34 956; 23 1]), 
        1))
    @test Float64 == typeof(ion(GCMS([1, 2, 3]u"s", 85.0:1.0:85.0, reshape([0, 956, 1], 
        (:,1))), 1))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError ion(FID([1, 2, 3]u"s", [12, 956, 23]), 1)
    @test_throws MethodError ion(TIC([1, 2, 3]u"s", [12, 956, 23]), 1)
end


@testset "ion RiGCMS index" begin
    # Same return values as those provided as arguments to construct the object
    @test 85 == ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1)
    @test 100 == ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2)
    
    # Provoke a BoundsError by specifying an index that does not exist
    @test_throws BoundsError ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 3)

    # Same return container and element type as used to construct the object
    @test Int64 == typeof(ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int64[85, 100], [0 12; 34 956; 23 1]), 1))
    @test Float64 == typeof(ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[85, 100], [0 12; 34 956; 23 1]), 1))
    @test Int64 == typeof(ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 85:86, 
        [0 12; 34 956; 23 1]), 1))
    @test Int64 == typeof(ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 85:85, 
        reshape([0, 956, 1], (:,1))), 1))
    @test Int64 == typeof(ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 85:15:100, 
        [0 12; 34 956; 23 1]), 1))
    @test Float64 == typeof(ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        85.0:15:100.0, [0 12; 34 956; 23 1]), 1))
    @test Float64 == typeof(ion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        85.0:1.0:85.0, reshape([0, 956, 1], (:,1))), 1))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError ion(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1)
    @test_throws MethodError ion(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 1)
end


############################################################################################
# ioncount(gcms::AbstractGCMS) -> Int
############################################################################################
@testset "ioncount GCMS" begin
    # Validate the returned value
    @test 2 == ioncount(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test 1 == ioncount(GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1], (:,1))))

    # Return value must be an integer
    @test Int == typeof(ioncount(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError ioncount(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test_throws MethodError ioncount(TIC([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "ioncount RiGCMS" begin
    # Validate the returned value
    @test 2 == ioncount(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]))
    @test 1 == ioncount(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85], 
        reshape([0, 956, 1], (:,1))))

    # Return value must be an integer
    @test Int == typeof(ioncount(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1])))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError ioncount(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test_throws MethodError ioncount(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
end


############################################################################################
# ionindex(gcms::AbstractGCMS, ion::Real) -> Int
############################################################################################
@testset "ionindex GCMS" begin
    # Validate the returned value
    @test 1 == ionindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 85)
    @test 2 == ionindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 100)
    @test_throws ArgumentError ionindex(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 200)
    @test 1 == ionindex(GCMS([1, 2, 3]u"s", Float32[85, 100], [0 12; 34 956; 23 1]), 85)
    @test 2 == ionindex(GCMS([1, 2, 3]u"s", Float64[85, 100], [0 12; 34 956; 23 1]), 100)
    @test 1 == ionindex(GCMS([1, 2, 3]u"s", Float32[85.5, 100.1], [0 12; 34 956; 23 1]), 
        85.5)
    @test 2 == ionindex(GCMS([1, 2, 3]u"s", Float32[85.5, 100.1], [0 12; 34 956; 23 1]), 
        100.1)

    # Return value must be an integer
    @test Int == typeof(ionindex(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 85))
    @test Int == typeof(ionindex(GCMS([1, 2, 3]u"s", Float32[85.5, 100.1], 
        [0 12; 34 956; 23 1]), 100.1))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError ionindex(FID([1, 2, 3]u"s", [12, 956, 23]), 85)
    @test_throws MethodError ionindex(TIC([1, 2, 3]u"s", [12, 956, 23]), 85)
end


@testset "ionindex RiGCMS" begin
    # Validate the returned value
    @test 1 == ionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 85)
    @test 2 == ionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 100)
    @test_throws ArgumentError ionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 200)
    @test 1 == ionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], Float32[85, 100], 
        [0 12; 34 956; 23 1]), 85)
    @test 2 == ionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], Float64[85, 100], 
        [0 12; 34 956; 23 1]), 100)
    @test 1 == ionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float32[85.5, 100.1], [0 12; 34 956; 23 1]), 85.5)
    @test 2 == ionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float32[85.5, 100.1], [0 12; 34 956; 23 1]), 100.1)

    # Return value must be an integer
    @test Int == typeof(ionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 85))
    @test Int == typeof(ionindex(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float32[85.5, 100.1], [0 12; 34 956; 23 1]), 100.1))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError ionindex(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 85)
    @test_throws MethodError ionindex(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]), 85)
end


############################################################################################
# ions(gcms::AbstractGCMS)
############################################################################################
@testset "ions GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [85, 100] == ions(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test [85] == ions(GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1], (:,1))))
    
    # Same return container and element type as used to construct the object
    @test Vector{Int64} == typeof(ions(GCMS([1, 2, 3]u"s", Int64[85, 100], 
        [0 12; 34 956; 23 1])))
    @test Vector{Float64} == typeof(ions(GCMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1])))
    @test UnitRange{Int64} == typeof(ions(GCMS([1, 2, 3]u"s", 85:86, [0 12; 34 956; 23 1])))
    @test UnitRange{Int64} == typeof(ions(GCMS([1, 2, 3]u"s", 85:85, 
        reshape([0, 956, 1], (:,1)))))
    @test StepRange{Int64, Int64} == typeof(ions(GCMS([1, 2, 3]u"s", 85:15:100, 
        [0 12; 34 956; 23 1])))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(ions(GCMS([1, 2, 3]u"s", 85.0:15:100.0, [0 12; 34 956; 23 1])))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(ions(GCMS([1, 2, 3]u"s", 85.0:1.0:85.0, 
        reshape([0, 956, 1], (:,1)))))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError ions(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test_throws MethodError ions(TIC([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "ions RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [85, 100] == ions(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]))
    @test [85] == ions(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85], 
        reshape([0, 956, 1], (:,1))))
    
    # Same return container and element type as used to construct the object
    @test Vector{Int64} == typeof(ions(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int64[85, 100], [0 12; 34 956; 23 1])))
    @test Vector{Float64} == typeof(ions(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[85, 100], [0 12; 34 956; 23 1])))
    @test UnitRange{Int64} == typeof(ions(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        85:86, [0 12; 34 956; 23 1])))
    @test UnitRange{Int64} == typeof(ions(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        85:85, reshape([0, 956, 1], (:,1)))))
    @test StepRange{Int64, Int64} == typeof(ions(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], 85:15:100, [0 12; 34 956; 23 1])))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(ions(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        85.0:15:100.0, [0 12; 34 956; 23 1])))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(ions(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        85.0:1.0:85.0, reshape([0, 956, 1], (:,1)))))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError ions(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test_throws MethodError ions(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
end


############################################################################################
# maxion(gcms::AbstractGCMS)
############################################################################################
@testset "maxion GCMS" begin
    @test 100 == maxion(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test Int == typeof(maxion(GCMS([1, 2, 3]u"s", Int[85, 100], [0 12; 34 956; 23 1])))
    @test Float64 == typeof(maxion(GCMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1])))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError maxion(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test_throws MethodError maxion(TIC([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "maxion RiGCMS" begin
    @test 100 == maxion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]))
    @test Int == typeof(maxion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[85, 100], [0 12; 34 956; 23 1])))
    @test Float64 == typeof(maxion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[85, 100], [0 12; 34 956; 23 1])))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError maxion(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test_throws MethodError maxion(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
end


############################################################################################
# minion(gcms::AbstractGCMS)
############################################################################################
@testset "minion GCMS" begin
    @test 85 == minion(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test Int == typeof(minion(GCMS([1, 2, 3]u"s", Int[85, 100], [0 12; 34 956; 23 1])))
    @test Float64 == typeof(minion(GCMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1])))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError minion(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test_throws MethodError minion(TIC([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "minion RiGCMS" begin
    @test 85 == minion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]))
    @test Int == typeof(minion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int[85, 100], [0 12; 34 956; 23 1])))
    @test Float64 == typeof(minion(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[85, 100], [0 12; 34 956; 23 1])))

    # Function does not work with objects that are not subtypes of AbstractGCMS
    @test_throws MethodError minion(RiFID([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
    @test_throws MethodError minion(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]))
end


############################################################################################
# intensities(gcms::AbstractGCMS; scanindexrange::OrdinalRange{Integer, Integer}, 
# ionindexrange::OrdinalRange{Integer, Integer})
############################################################################################
@testset "intensities GCMS scanindexrange ionindexrange" begin
    # Same return values as those provided as arguments to construct the object
    @test [0 12; 34 956; 23 1] == intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
    @test reshape([0, 956, 1], (:,1)) == intensities(GCMS([1, 2, 3]u"s", [85],
        reshape([0, 956, 1], (:,1))))
    @test [34 956; 23 1] == intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2)
    @test reshape([12], (1, 1)) == intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:1, ionindexrange=2:2)

    # Same return predictable container and element types
    @test Matrix{Int64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Int64[0 12; 34 956; 23 1])))
    @test Matrix{Float64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1])))
    @test Matrix{Int64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1],
         (:,1)))))
    @test SubArray{Int64, 2, Matrix{Int64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, 
        false} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2))
    @test SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, 
        false} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2))
    @test Vector{Int64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2)[:])
    @test Vector{Float64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2)[:])
    
    # Check for BoundsErrors
    @test_throws BoundsError intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:3, ionindexrange=1:2)
    @test_throws BoundsError intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:4, ionindexrange=1:2)
    @test_throws BoundsError intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:3, ionindexrange=0:2)
    @test_throws BoundsError intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:3, ionindexrange=1:3)
end


@testset "intensities RiGCMS scanindexrange ionindexrange" begin
    # Same return values as those provided as arguments to construct the object
    @test [0 12; 34 956; 23 1] == intensities(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]))
    @test reshape([0, 956, 1], (:,1)) == intensities(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85], reshape([0, 956, 1], (:,1))))
    @test [34 956; 23 1] == intensities(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2)
    @test reshape([12], (1, 1)) == intensities(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), scanindexrange=1:1, 
        ionindexrange=2:2)

    # Same return predictable container and element types
    @test Matrix{Int64} == typeof(intensities(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], Int64[0 12; 34 956; 23 1])))
    @test Matrix{Float64} == typeof(intensities(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], Float64[0 12; 34 956; 23 1])))
    @test Matrix{Int64} == typeof(intensities(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85], reshape([0, 956, 1], (:,1)))))
    @test SubArray{Int64, 2, Matrix{Int64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, 
        false} == typeof(intensities(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], Int[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2))
    @test SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, 
        false} == typeof(intensities(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], Float64[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2))
    @test Vector{Int64} == typeof(intensities(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], Int[0 12; 34 956; 23 1]), scanindexrange=2:3, 
        ionindexrange=1:2)[:])
    @test Vector{Float64} == typeof(intensities(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], Float64[0 12; 34 956; 23 1]), scanindexrange=2:3, 
        ionindexrange=1:2)[:])
    
    # Check for BoundsErrors
    @test_throws BoundsError intensities(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), scanindexrange=0:3, ionindexrange=1:2)
    @test_throws BoundsError intensities(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), scanindexrange=1:4, ionindexrange=1:2)
    @test_throws BoundsError intensities(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), scanindexrange=1:3, ionindexrange=0:2)
    @test_throws BoundsError intensities(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), scanindexrange=1:3, ionindexrange=1:3)
end


############################################################################################
# intensity(gcms::AbstractGCMS, scanindex::Integer, ionindex::Integer)
############################################################################################
@testset "intensity GCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 34 == intensity(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2, 1)
    @test 12 == intensity(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1, 2)

    # Returned value has same type as element type was used to construct the object
    @test Int == typeof(intensity(GCMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        1, 2))
    @test Float64 == typeof(intensity(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), 1, 2))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError intensity(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0, 1)
    @test_throws BoundsError intensity(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 4, 1)
    @test_throws BoundsError intensity(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, 0)
    @test_throws BoundsError intensity(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, 3)
end


@testset "intensity RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 34 == intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2, 1)
    @test 12 == intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1, 2)

    # Returned value has same type as element type was used to construct the object
    @test Int == typeof(intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], Int[0 12; 34 956; 23 1]), 1, 2))
    @test Float64 == typeof(intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], Float64[0 12; 34 956; 23 1]), 1, 2))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 0, 1)
    @test_throws BoundsError intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 4, 1)
    @test_throws BoundsError intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2, 0)
    @test_throws BoundsError intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 2, 3)
end


############################################################################################
# intensity(gcms::AbstractGCMS, time::Unitful.Time, ion::Real; precisetime::Bool=false)
############################################################################################
@testset "intensity GCMS time ion" begin
    # Same return values as those provided as arguments to construct the object
    @test 34 == intensity(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1.9u"s", 
        85)
    @test 34 == intensity(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2u"s", 
        85)
    @test 34 == intensity(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2.0u"s", 
        85, precisetime=true)

    # Same return predictable container and element types
    @test Int == typeof(intensity(GCMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        1.9u"s", 85))
    @test Float64 == typeof(intensity(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), 1.9u"s", 85))
    
    @test_throws ArgumentError intensity(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.9u"s", 101)
    @test_throws ArgumentError intensity(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.9u"s", 101, precisetime=true)
end


@testset "intensity RiGCMS time ion" begin
    # Same return values as those provided as arguments to construct the object
    @test 34 == intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 1.9u"s", 85)
    @test 34 == intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s", 85)
    @test 34 == intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), 2.0u"s", 85, precisetime=true)

    # Same return predictable container and element types
    @test Int == typeof(intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], Int[0 12; 34 956; 23 1]), 1.9u"s", 85))
    @test Float64 == typeof(intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], Float64[0 12; 34 956; 23 1]), 1.9u"s", 85))
    
    @test_throws ArgumentError intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1.9u"s", 101)
    @test_throws ArgumentError intensity(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]), 1.9u"s", 101, precisetime=true)
end


############################################################################################
# ionscantime(δtᵢ::Function, gcms::AbstractGCMS, scanindex::Integer, ionindex::Integer;
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "ionscantime GCMS" begin
    # Same return values as those provided as arguments to construct the object
    gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    δtᵢ = ionscantimeshift(gcms, LinearDescending())
    @test 2u"s" ≈ ionscantime(δtᵢ, gcms, 2, 1)
    @test 1.5u"s" ≈ ionscantime(δtᵢ, gcms, 2, 2)
    @test 0.03333333333333333u"minute" ≈ ionscantime(δtᵢ, gcms, 2, 1, timeunit=u"minute")
    @test 0.025u"minute" ≈ ionscantime(δtᵢ, gcms, 2, 2, timeunit=u"minute")
    @test 2 ≈ ionscantime(δtᵢ, gcms, 2, 1, ustripped=true)
    @test 1.5 ≈ ionscantime(δtᵢ, gcms, 2, 2, ustripped=true)
    @test 0.03333333333333333 ≈ ionscantime(δtᵢ, gcms, 2, 1, timeunit=u"minute", 
        ustripped=true)
    @test 0.025 ≈ ionscantime(δtᵢ, gcms, 2, 2, timeunit=u"minute", ustripped=true)

    # Returned element type is always Float64
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(ionscantime(δtᵢ, gcms, 2, 1))
    @test Float64 == typeof(ionscantime(δtᵢ, gcms, 2, 1, ustripped=true))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(ionscantime(δtᵢ, gcms, 2, 1, timeunit=u"minute"))
    @test Float64 == typeof(ionscantime(δtᵢ, gcms, 2, 1, timeunit=u"minute", ustripped=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError ionscantime(δtᵢ, gcms, 2, 0)
    @test_throws BoundsError ionscantime(δtᵢ, gcms, 2, 3)
    @test_throws BoundsError ionscantime(δtᵢ, gcms, 0, 2)
    @test_throws BoundsError ionscantime(δtᵢ, gcms, 4, 2)
end


@testset "ionscantime RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    gcms = RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])
    δtᵢ = ionscantimeshift(gcms, LinearDescending())
    @test 2u"s" ≈ ionscantime(δtᵢ, gcms, 2, 1)
    @test 1.5u"s" ≈ ionscantime(δtᵢ, gcms, 2, 2)
    @test 0.03333333333333333u"minute" ≈ ionscantime(δtᵢ, gcms, 2, 1, timeunit=u"minute")
    @test 0.025u"minute" ≈ ionscantime(δtᵢ, gcms, 2, 2, timeunit=u"minute")
    @test 2 ≈ ionscantime(δtᵢ, gcms, 2, 1, ustripped=true)
    @test 1.5 ≈ ionscantime(δtᵢ, gcms, 2, 2, ustripped=true)
    @test 0.03333333333333333 ≈ ionscantime(δtᵢ, gcms, 2, 1, timeunit=u"minute", 
        ustripped=true)
    @test 0.025 ≈ ionscantime(δtᵢ, gcms, 2, 2, timeunit=u"minute", ustripped=true)

    # Returned element type is always Float64
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(ionscantime(δtᵢ, gcms, 2, 1))
    @test Float64 == typeof(ionscantime(δtᵢ, gcms, 2, 1, ustripped=true))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(ionscantime(δtᵢ, gcms, 2, 1, timeunit=u"minute"))
    @test Float64 == typeof(ionscantime(δtᵢ, gcms, 2, 1, timeunit=u"minute", ustripped=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError ionscantime(δtᵢ, gcms, 2, 0)
    @test_throws BoundsError ionscantime(δtᵢ, gcms, 2, 3)
    @test_throws BoundsError ionscantime(δtᵢ, gcms, 0, 2)
    @test_throws BoundsError ionscantime(δtᵢ, gcms, 4, 2)
end


############################################################################################
# ionscantimes(δtᵢ::Function, gcms::AbstractGCMS, ionindex::Integer; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "ionscantimes GCMS" begin
    # Same return values as those provided as arguments to construct the object
    gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    δtᵢ = ionscantimeshift(gcms, LinearDescending())
    @test [1, 2, 3]u"s" ≈ ionscantimes(δtᵢ, gcms, 1)
    @test [0.5, 1.5, 2.5]u"s" ≈ ionscantimes(δtᵢ, gcms, 2)
    @test [0.016666666666666666, 0.03333333333333333, 0.05]u"minute" ≈ ionscantimes(δtᵢ, 
        gcms, 1, timeunit=u"minute")
    @test [0.008333333333333333, 0.025, 0.041666666666666664]u"minute" ≈ ionscantimes(δtᵢ, 
        gcms, 2, timeunit=u"minute")
    @test [1, 2, 3] ≈ ionscantimes(δtᵢ, gcms, 1, ustripped=true)
    @test [0.5, 1.5, 2.5] ≈ ionscantimes(δtᵢ, gcms, 2, ustripped=true)
    @test [0.016666666666666666, 0.03333333333333333, 0.05] ≈ ionscantimes(δtᵢ, gcms, 1, 
        timeunit=u"minute", ustripped=true)
    @test [0.008333333333333333, 0.025, 0.041666666666666664] ≈ ionscantimes(δtᵢ, gcms, 2, 
        timeunit=u"minute", ustripped=true)

    # Same return element type includes always Float64
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(ionscantimes(δtᵢ, gcms, 1))
    @test Vector{Float64} == typeof(ionscantimes(δtᵢ, gcms, 1, ustripped=true))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(ionscantimes(δtᵢ, gcms, 1, timeunit=u"minute"))
    @test Vector{Float64} == typeof(ionscantimes(δtᵢ, gcms, 1, timeunit=u"minute", 
        ustripped=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError ionscantimes(δtᵢ, gcms, 0)
    @test_throws BoundsError ionscantimes(δtᵢ, gcms, 3)
end


@testset "ionscantimes RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    gcms = RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])
    δtᵢ = ionscantimeshift(gcms, LinearDescending())
    @test [1, 2, 3]u"s" ≈ ionscantimes(δtᵢ, gcms, 1)
    @test [0.5, 1.5, 2.5]u"s" ≈ ionscantimes(δtᵢ, gcms, 2)
    @test [0.016666666666666666, 0.03333333333333333, 0.05]u"minute" ≈ ionscantimes(δtᵢ, 
        gcms, 1, timeunit=u"minute")
    @test [0.008333333333333333, 0.025, 0.041666666666666664]u"minute" ≈ ionscantimes(δtᵢ, 
        gcms, 2, timeunit=u"minute")
    @test [1, 2, 3] ≈ ionscantimes(δtᵢ, gcms, 1, ustripped=true)
    @test [0.5, 1.5, 2.5] ≈ ionscantimes(δtᵢ, gcms, 2, ustripped=true)
    @test [0.016666666666666666, 0.03333333333333333, 0.05] ≈ ionscantimes(δtᵢ, gcms, 1, 
        timeunit=u"minute", ustripped=true)
    @test [0.008333333333333333, 0.025, 0.041666666666666664] ≈ ionscantimes(δtᵢ, gcms, 2, 
        timeunit=u"minute", ustripped=true)

    # Same return element type includes always Float64
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(ionscantimes(δtᵢ, gcms, 1))
    @test Vector{Float64} == typeof(ionscantimes(δtᵢ, gcms, 1, ustripped=true))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(ionscantimes(δtᵢ, gcms, 1, timeunit=u"minute"))
    @test Vector{Float64} == typeof(ionscantimes(δtᵢ, gcms, 1, timeunit=u"minute", 
        ustripped=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError ionscantimes(δtᵢ, gcms, 0)
    @test_throws BoundsError ionscantimes(δtᵢ, gcms, 3)
end


############################################################################################
# ionscantimeindex(δtᵢ::Function, gcms::AbstractGCMS, ionindex::Integer, 
# time::Unitful.Time; precisetime::Bool=false) -> Int
############################################################################################
@testset "ionscantimeindex GCMS" begin
    # Same return values as those provided as arguments to construct the object
    gcms = GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    δtᵢ = ionscantimeshift(gcms, LinearDescending())
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s")
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 0.025u"minute")
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 1.6u"s")
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s", precisetime=true)
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 0.025u"minute", precisetime=true)
    @test 1 == ionscantimeindex(δtᵢ, gcms, 1, 1.4u"s")
    @test 1 == ionscantimeindex(δtᵢ, gcms, 1, 1u"s")
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 1.0u"s")
    @test 1 == ionscantimeindex(δtᵢ, gcms, 2, prevfloat(1.0)*u"s")

    # Same return element type includes always Int
    @test Int == typeof(ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s"))
    @test Int == typeof(ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s", precisetime=true))
    @test Int == typeof(ionscantimeindex(δtᵢ, gcms, 2, 0.025u"minute", precisetime=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws ArgumentError ionscantimeindex(δtᵢ, gcms, 2, 1.6u"s", precisetime=true)
    @test_throws BoundsError ionscantimeindex(δtᵢ, gcms, 3, 1.5u"s")
    @test_throws BoundsError ionscantimeindex(δtᵢ, gcms, 0, 1.5u"s")
end


@testset "ionscantimeindex RiGCMS" begin
    # Same return values as those provided as arguments to construct the object
    gcms = RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])
    δtᵢ = ionscantimeshift(gcms, LinearDescending())
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s")
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 0.025u"minute")
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 1.6u"s")
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s", precisetime=true)
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 0.025u"minute", precisetime=true)
    @test 1 == ionscantimeindex(δtᵢ, gcms, 1, 1.4u"s")
    @test 1 == ionscantimeindex(δtᵢ, gcms, 1, 1u"s")
    @test 2 == ionscantimeindex(δtᵢ, gcms, 2, 1.0u"s")
    @test 1 == ionscantimeindex(δtᵢ, gcms, 2, prevfloat(1.0)*u"s")

    # Same return element type includes always Int
    @test Int == typeof(ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s"))
    @test Int == typeof(ionscantimeindex(δtᵢ, gcms, 2, 1.5u"s", precisetime=true))
    @test Int == typeof(ionscantimeindex(δtᵢ, gcms, 2, 0.025u"minute", precisetime=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws ArgumentError ionscantimeindex(δtᵢ, gcms, 2, 1.6u"s", precisetime=true)
    @test_throws BoundsError ionscantimeindex(δtᵢ, gcms, 3, 1.5u"s")
    @test_throws BoundsError ionscantimeindex(δtᵢ, gcms, 0, 1.5u"s")
end


############################################################################################
# totalionchromatogram(gcms::AbstractGCMS)
############################################################################################
@testset "totalionchromatogram GCMS" begin
    # Check the returned type and associated supertypes
    @test isa(totalionchromatogram(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])), 
        TIC)
    @test isa(totalionchromatogram(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])), 
        AbstractTIC)
    @test isa(totalionchromatogram(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])), 
        AbstractChromatogram)

    # Check the contents of the TIC for plausibility
    @test [1, 2, 3]u"s" == scantimes(totalionchromatogram(GCMS([1, 2, 3]u"s", [85, 100],
        [0 12; 34 956; 23 1])))
    @test [12, 990, 24] == intensities(totalionchromatogram(GCMS([1, 2, 3]u"s", [85, 100],
        [0 12; 34 956; 23 1])))
    @test Dict(:id => 4, "name" => "sample") == metadata(totalionchromatogram(GCMS(
        [1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1], Dict(:id => 4, "name" => "sample")
        )))
    
    # Verify that the element type of the TIC intensities is the same as the GCMS
    @test Vector{Int} == typeof(intensities(totalionchromatogram(GCMS([1, 2, 3]u"s", 
        [85, 100], Int[0 12; 34 956; 23 1]))))
    @test Vector{Float64} == typeof(intensities(totalionchromatogram(GCMS([1, 2, 3]u"s", 
        [85, 100], Float64[0 12; 34 956; 23 1]))))
end


@testset "totalionchromatogram RiGCMS" begin
    # Check the returned type and associated supertypes
    @test isa(totalionchromatogram(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1])), RiTIC)
    @test isa(totalionchromatogram(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1])), AbstractTIC)
    @test isa(totalionchromatogram(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1])), AbstractChromatogram)

    # Check the contents of the TIC for plausibility
    @test [1, 2, 3]u"s" == scantimes(totalionchromatogram(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test [12, 990, 24] == intensities(totalionchromatogram(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test Dict(:id => 4, "name" => "sample") == metadata(totalionchromatogram(RiGCMS(
        [1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1], 
        Dict(:id => 4, "name" => "sample"))))
    
    # Verify that the element type of the TIC intensities is the same as the GCMS
    @test Vector{Int} == typeof(intensities(totalionchromatogram(RiGCMS([1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [85, 100], Int[0 12; 34 956; 23 1]))))
    @test Vector{Float64} == typeof(intensities(totalionchromatogram(RiGCMS([1, 2, 3]u"s", 
        "Kovats", [100, 200, 300], [85, 100], Float64[0 12; 34 956; 23 1]))))
end