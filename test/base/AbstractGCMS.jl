using JuChrom
using Test
using Unitful: 𝐓


############################################################################################
# binions(gcms::AbstractGCMS; ionbin::Function=integerion)
############################################################################################
@testset "ion GCMS" begin
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


############################################################################################
# ion(gcms::AbstractGCMS, index::Integer)
############################################################################################
@testset "ion GCMS" begin
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


############################################################################################
# ionscantime(timeshift::Function, gcms::AbstractGCMS, scanindex::Integer, 
# ionindex::Integer; timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "ionscantime" begin
    # Same return values as those provided as arguments to construct the object
    gcms = GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])
    δt = timeshift(gcms, LinearDescending())
    @test 2u"s" ≈ ionscantime(δt, gcms, 1, 2)
    @test 1.5u"s" ≈ ionscantime(δt, gcms, 2, 2)
    @test 0.03333333333333333u"minute" ≈ ionscantime(δt, gcms, 1, 2, timeunit=u"minute")
    @test 0.025u"minute" ≈ ionscantime(δt, gcms, 2, 2, timeunit=u"minute")
    @test 2 ≈ ionscantime(δt, gcms, 1, 2, ustripped=true)
    @test 1.5 ≈ ionscantime(δt, gcms, 2, 2, ustripped=true)
    @test 0.03333333333333333 ≈ ionscantime(δt, gcms, 1, 2, timeunit=u"minute", ustripped=true)
    @test 0.025 ≈ ionscantime(δt, gcms, 2, 2, timeunit=u"minute", ustripped=true)

    # Same return element type as used to construct the object
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(ionscantime(δt, gcms, 1, 2))
    @test Float64 == typeof(ionscantime(δt, gcms, 1, 2, ustripped=true))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(ionscantime(δt, gcms, 1, 2, timeunit=u"minute"))
    @test Float64 == typeof(ionscantime(δt, gcms, 1, 2, timeunit=u"minute", ustripped=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError ionscantime(δt, gcms, 0, 2)
    @test_throws BoundsError ionscantime(δt, gcms, 3, 2)
    @test_throws BoundsError ionscantime(δt, gcms, 2, 0)
    @test_throws BoundsError ionscantime(δt, gcms, 2, 4)
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