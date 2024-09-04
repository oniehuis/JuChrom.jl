using JuChrom
using Test
using Unitful: 𝐓


############################################################################################
# binions(chrom::AbstractChromMS; ionbin::Function=integerion)
############################################################################################
@testset "binions ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test (1:3)u"s" == binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])).scantimes
    @test [85, 101] == binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])).ions
    @test [24 12; 0 956; 23 1] == binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])).intensities
    @test [24.1 12.2; 1.0 956.7; 23.9 1.5] ≈ binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5])).intensities
    @test Dict(:id => 4) == binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1], Dict(:id => 4))).metadata
    @test (1:3)u"s" == binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1]), ionbin=ion->integer(ion, start=0.9)).scantimes
    @test [84, 85, 101] == binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1]), ionbin=ion->integer(ion, start=0.9)).ions
    @test [0 24 12; 0 0 956; 23 0 1] == binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1]), ionbin=ion->integer(ion, start=0.9)).intensities
    @test [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5] ≈ binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5]), 
        ionbin=ion->integer(ion, start=0.9)).intensities
    @test Dict(:id => 4) == binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1], Dict(:id => 4)), 
        ionbin=ion->integer(ion, start=0.9)).metadata

    # Check the returned type and supertypes
    @test isa(binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])), 
        ChromMS)
    @test isa(binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])), 
        AbstractChromMS)
    @test isa(binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])), 
        AbstractChromatogram)

    # Binned intensities have the same type as used to construct the object
    @test Matrix{Int64} == typeof(binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0 24 12; 0 0 956; 23 0 1])).intensities)
    @test Matrix{Float64} == typeof(binions(ChromMS((1:3)u"s", [84.8, 85.2, 100.9], 
        [0.0 24.1 12.2; 0.1 0.9 956.7; 23.1 0.8 1.5])).intensities)

    # Function does not work with objects that are not subtypes of AbstractChromMS
    @test_throws MethodError binions(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1)
end


############################################################################################
# ion(chrom::AbstractChromMS, index::Integer)
############################################################################################
@testset "ion ChromMS index" begin
    # Same return values as those provided as arguments to construct the object
    @test 85 == ion(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1)
    @test 100 == ion(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2)
    
    # Provoke a BoundsError by specifying an index that does not exist
    @test_throws BoundsError ion(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 3)

    # Same return container and element type as used to construct the object
    @test Int64 == typeof(ion(ChromMS([1, 2, 3]u"s", Int64[85, 100], [0 12; 34 956; 23 1]), 1))
    @test Float64 == typeof(ion(ChromMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1]), 1))
    @test Int64 == typeof(ion(ChromMS([1, 2, 3]u"s", 85:86, [0 12; 34 956; 23 1]), 1))
    @test Int64 == typeof(ion(ChromMS([1, 2, 3]u"s", 85:85, reshape([0, 956, 1], (:,1))), 1))
    @test Int64 == typeof(ion(ChromMS([1, 2, 3]u"s", 85:15:100, [0 12; 34 956; 23 1]), 1))
    @test Float64 == typeof(ion(ChromMS([1, 2, 3]u"s", 85.0:15:100.0, [0 12; 34 956; 23 1]), 
        1))
    @test Float64 == typeof(ion(ChromMS([1, 2, 3]u"s", 85.0:1.0:85.0, reshape([0, 956, 1], 
        (:,1))), 1))

    # Function does not work with objects that are not subtypes of AbstractChromMS
    @test_throws MethodError ion(Chrom([1, 2, 3]u"s", [12, 956, 23]), 1)
end


############################################################################################
# ioncount(chrom::AbstractChromMS) -> Int
############################################################################################
@testset "ioncount ChromMS" begin
    # Validate the returned value
    @test 2 == ioncount(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test 1 == ioncount(ChromMS([1, 2, 3]u"s", [85], reshape([0, 956, 1], (:,1))))

    # Return value must be an integer
    @test Int == typeof(ioncount(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])))

    # Function does not work with objects that are not subtypes of AbstractChromMS
    @test_throws MethodError ioncount(Chrom([1, 2, 3]u"s", [12, 956, 23]))
end


############################################################################################
# ionindex(chrom::AbstractChromMS, ion::Real) -> Int
############################################################################################
@testset "ionindex ChromMS" begin
    # Validate the returned value
    @test 1 == ionindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 85)
    @test 2 == ionindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 100)
    @test_throws ArgumentError ionindex(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 200)
    @test 1 == ionindex(ChromMS([1, 2, 3]u"s", Float32[85, 100], [0 12; 34 956; 23 1]), 85)
    @test 2 == ionindex(ChromMS([1, 2, 3]u"s", Float64[85, 100], [0 12; 34 956; 23 1]), 100)
    @test 1 == ionindex(ChromMS([1, 2, 3]u"s", Float32[85.5, 100.1], [0 12; 34 956; 23 1]), 
        85.5)
    @test 2 == ionindex(ChromMS([1, 2, 3]u"s", Float32[85.5, 100.1], [0 12; 34 956; 23 1]), 
        100.1)

    # Return value must be an integer
    @test Int == typeof(ionindex(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 85))
    @test Int == typeof(ionindex(ChromMS([1, 2, 3]u"s", Float32[85.5, 100.1], 
        [0 12; 34 956; 23 1]), 100.1))

    # Function does not work with objects that are not subtypes of AbstractChromMS
    @test_throws MethodError ionindex(Chrom([1, 2, 3]u"s", [12, 956, 23]), 85)
end


############################################################################################
# ions(chrom::AbstractChromMS)
############################################################################################
@testset "ions ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test [85, 100] == ions(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test [85] == ions(ChromMS([1, 2, 3]u"s", [85], reshape([0, 956, 1], (:,1))))
    
    # Same return container and element type as used to construct the object
    @test Vector{Int64} == typeof(ions(ChromMS([1, 2, 3]u"s", Int64[85, 100], 
        [0 12; 34 956; 23 1])))
    @test Vector{Float64} == typeof(ions(ChromMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1])))
    @test UnitRange{Int64} == typeof(ions(ChromMS([1, 2, 3]u"s", 85:86, [0 12; 34 956; 23 1])))
    @test UnitRange{Int64} == typeof(ions(ChromMS([1, 2, 3]u"s", 85:85, 
        reshape([0, 956, 1], (:,1)))))
    @test StepRange{Int64, Int64} == typeof(ions(ChromMS([1, 2, 3]u"s", 85:15:100, 
        [0 12; 34 956; 23 1])))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(ions(ChromMS([1, 2, 3]u"s", 85.0:15:100.0, [0 12; 34 956; 23 1])))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(ions(ChromMS([1, 2, 3]u"s", 85.0:1.0:85.0, 
        reshape([0, 956, 1], (:,1)))))

    # Function does not work with objects that are not subtypes of AbstractChromMS
    @test_throws MethodError ions(Chrom([1, 2, 3]u"s", [12, 956, 23]))
end


############################################################################################
# maxion(chrom::AbstractChromMS)
############################################################################################
@testset "maxion ChromMS" begin
    @test 100 == maxion(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test Int == typeof(maxion(ChromMS([1, 2, 3]u"s", Int[85, 100], [0 12; 34 956; 23 1])))
    @test Float64 == typeof(maxion(ChromMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1])))

    # Function does not work with objects that are not subtypes of AbstractChromMS
    @test_throws MethodError maxion(Chrom([1, 2, 3]u"s", [12, 956, 23]))
end


############################################################################################
# minion(chrom::AbstractChromMS)
############################################################################################
@testset "minion ChromMS" begin
    @test 85 == minion(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test Int == typeof(minion(ChromMS([1, 2, 3]u"s", Int[85, 100], [0 12; 34 956; 23 1])))
    @test Float64 == typeof(minion(ChromMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1])))

    # Function does not work with objects that are not subtypes of AbstractChromMS
    @test_throws MethodError minion(Chrom([1, 2, 3]u"s", [12, 956, 23]))
end


############################################################################################
# intensities(chrom::AbstractChromMS; scanindexrange::OrdinalRange{Integer, Integer}, 
# ionindexrange::OrdinalRange{Integer, Integer})
############################################################################################
@testset "intensities ChromMS scanindexrange ionindexrange" begin
    # Same return values as those provided as arguments to construct the object
    @test [0 12; 34 956; 23 1] == intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
    @test reshape([0, 956, 1], (:,1)) == intensities(ChromMS([1, 2, 3]u"s", [85],
        reshape([0, 956, 1], (:,1))))
    @test [34 956; 23 1] == intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2)
    @test reshape([12], (1, 1)) == intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:1, ionindexrange=2:2)

    # Same return predictable container and element types
    @test Matrix{Int64} == typeof(intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int64[0 12; 34 956; 23 1])))
    @test Matrix{Float64} == typeof(intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1])))
    @test Matrix{Int64} == typeof(intensities(ChromMS([1, 2, 3]u"s", [85], reshape([0, 956, 1],
         (:,1)))))
    @test SubArray{Int64, 2, Matrix{Int64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, 
        false} == typeof(intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2))
    @test SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, 
        false} == typeof(intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2))
    @test Vector{Int64} == typeof(intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2)[:])
    @test Vector{Float64} == typeof(intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=2:3, ionindexrange=1:2)[:])
    
    # Check for BoundsErrors
    @test_throws BoundsError intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:3, ionindexrange=1:2)
    @test_throws BoundsError intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:4, ionindexrange=1:2)
    @test_throws BoundsError intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:3, ionindexrange=0:2)
    @test_throws BoundsError intensities(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:3, ionindexrange=1:3)
end


############################################################################################
# intensity(chrom::AbstractChromMS, scanindex::Integer, ionindex::Integer)
############################################################################################
@testset "intensity ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    @test 34 == intensity(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2, 1)
    @test 12 == intensity(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1, 2)

    # Returned value has same type as element type was used to construct the object
    @test Int == typeof(intensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        1, 2))
    @test Float64 == typeof(intensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), 1, 2))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError intensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0, 1)
    @test_throws BoundsError intensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 4, 1)
    @test_throws BoundsError intensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, 0)
    @test_throws BoundsError intensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, 3)
end


############################################################################################
# intensity(chrom::AbstractChromMS, time::Unitful.Time, ion::Real; precisetime::Bool=false)
############################################################################################
@testset "intensity ChromMS time ion" begin
    # Same return values as those provided as arguments to construct the object
    @test 34 == intensity(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 1.9u"s", 
        85)
    @test 34 == intensity(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2u"s", 
        85)
    @test 34 == intensity(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 2.0u"s", 
        85, precisetime=true)

    # Same return predictable container and element types
    @test Int == typeof(intensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        1.9u"s", 85))
    @test Float64 == typeof(intensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), 1.9u"s", 85))
    
    @test_throws ArgumentError intensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.9u"s", 101)
    @test_throws ArgumentError intensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 1.9u"s", 101, precisetime=true)
end


############################################################################################
# ionscantime(δtᵢ::Function, chrom::AbstractChromMS, scanindex::Integer, ionindex::Integer;
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "ionscantime ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    δtᵢ = ionscantimeshift(chrom, LinearDescending())
    @test 2u"s" ≈ ionscantime(δtᵢ, chrom, 2, 1)
    @test 1.5u"s" ≈ ionscantime(δtᵢ, chrom, 2, 2)
    @test 0.03333333333333333u"minute" ≈ ionscantime(δtᵢ, chrom, 2, 1, timeunit=u"minute")
    @test 0.025u"minute" ≈ ionscantime(δtᵢ, chrom, 2, 2, timeunit=u"minute")
    @test 2 ≈ ionscantime(δtᵢ, chrom, 2, 1, ustripped=true)
    @test 1.5 ≈ ionscantime(δtᵢ, chrom, 2, 2, ustripped=true)
    @test 0.03333333333333333 ≈ ionscantime(δtᵢ, chrom, 2, 1, timeunit=u"minute", 
        ustripped=true)
    @test 0.025 ≈ ionscantime(δtᵢ, chrom, 2, 2, timeunit=u"minute", ustripped=true)

    # Returned element type is always Float64
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(ionscantime(δtᵢ, chrom, 2, 1))
    @test Float64 == typeof(ionscantime(δtᵢ, chrom, 2, 1, ustripped=true))
    @test Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓, 
        nothing}} == typeof(ionscantime(δtᵢ, chrom, 2, 1, timeunit=u"minute"))
    @test Float64 == typeof(ionscantime(δtᵢ, chrom, 2, 1, timeunit=u"minute", ustripped=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError ionscantime(δtᵢ, chrom, 2, 0)
    @test_throws BoundsError ionscantime(δtᵢ, chrom, 2, 3)
    @test_throws BoundsError ionscantime(δtᵢ, chrom, 0, 2)
    @test_throws BoundsError ionscantime(δtᵢ, chrom, 4, 2)
end


############################################################################################
# ionscantimes(δtᵢ::Function, chrom::AbstractChromMS, ionindex::Integer; 
# timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "ionscantimes ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    δtᵢ = ionscantimeshift(chrom, LinearDescending())
    @test [1, 2, 3]u"s" ≈ ionscantimes(δtᵢ, chrom, 1)
    @test [0.5, 1.5, 2.5]u"s" ≈ ionscantimes(δtᵢ, chrom, 2)
    @test [0.016666666666666666, 0.03333333333333333, 0.05]u"minute" ≈ ionscantimes(δtᵢ, 
        chrom, 1, timeunit=u"minute")
    @test [0.008333333333333333, 0.025, 0.041666666666666664]u"minute" ≈ ionscantimes(δtᵢ, 
        chrom, 2, timeunit=u"minute")
    @test [1, 2, 3] ≈ ionscantimes(δtᵢ, chrom, 1, ustripped=true)
    @test [0.5, 1.5, 2.5] ≈ ionscantimes(δtᵢ, chrom, 2, ustripped=true)
    @test [0.016666666666666666, 0.03333333333333333, 0.05] ≈ ionscantimes(δtᵢ, chrom, 1, 
        timeunit=u"minute", ustripped=true)
    @test [0.008333333333333333, 0.025, 0.041666666666666664] ≈ ionscantimes(δtᵢ, chrom, 2, 
        timeunit=u"minute", ustripped=true)

    # Same return element type includes always Float64
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(ionscantimes(δtᵢ, chrom, 1))
    @test Vector{Float64} == typeof(ionscantimes(δtᵢ, chrom, 1, ustripped=true))
    @test Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}}} == typeof(ionscantimes(δtᵢ, chrom, 1, timeunit=u"minute"))
    @test Vector{Float64} == typeof(ionscantimes(δtᵢ, chrom, 1, timeunit=u"minute", 
        ustripped=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws BoundsError ionscantimes(δtᵢ, chrom, 0)
    @test_throws BoundsError ionscantimes(δtᵢ, chrom, 3)
end


############################################################################################
# ionscantimeindex(δtᵢ::Function, chrom::AbstractChromMS, ionindex::Integer, 
# time::Unitful.Time; precisetime::Bool=false) -> Int
############################################################################################
@testset "ionscantimeindex ChromMS" begin
    # Same return values as those provided as arguments to construct the object
    chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    δtᵢ = ionscantimeshift(chrom, LinearDescending())
    @test 2 == ionscantimeindex(δtᵢ, chrom, 2, 1.5u"s")
    @test 2 == ionscantimeindex(δtᵢ, chrom, 2, 0.025u"minute")
    @test 2 == ionscantimeindex(δtᵢ, chrom, 2, 1.6u"s")
    @test 2 == ionscantimeindex(δtᵢ, chrom, 2, 1.5u"s", precisetime=true)
    @test 2 == ionscantimeindex(δtᵢ, chrom, 2, 0.025u"minute", precisetime=true)
    @test 1 == ionscantimeindex(δtᵢ, chrom, 1, 1.4u"s")
    @test 1 == ionscantimeindex(δtᵢ, chrom, 1, 1u"s")
    @test 2 == ionscantimeindex(δtᵢ, chrom, 2, 1.0u"s")
    @test 1 == ionscantimeindex(δtᵢ, chrom, 2, prevfloat(1.0)*u"s")

    # Same return element type includes always Int
    @test Int == typeof(ionscantimeindex(δtᵢ, chrom, 2, 1.5u"s"))
    @test Int == typeof(ionscantimeindex(δtᵢ, chrom, 2, 1.5u"s", precisetime=true))
    @test Int == typeof(ionscantimeindex(δtᵢ, chrom, 2, 0.025u"minute", precisetime=true))

    # Check BoundsErrors for scanindex and ionindex
    @test_throws ArgumentError ionscantimeindex(δtᵢ, chrom, 2, 1.6u"s", precisetime=true)
    @test_throws BoundsError ionscantimeindex(δtᵢ, chrom, 3, 1.5u"s")
    @test_throws BoundsError ionscantimeindex(δtᵢ, chrom, 0, 1.5u"s")
end


############################################################################################
# maxintensity(chrom::AbstractChromMS[; scanindexrange, ionindexrange])
############################################################################################
@testset "maxintensity ChromMS" begin
    # Validate the returned value
    @test 956 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]))
    @test 12 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        scanindexrange=1:1)
    @test 956 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        scanindexrange=1:2)
    @test 956 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]))
    @test 12 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        scanindexrange=1:1)
    @test 956 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        scanindexrange=1:2)

    @test 34 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        ionindexrange=1:1)
    @test 956 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        ionindexrange=1:2)
    @test 34 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        ionindexrange=1:1)
    @test 956 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        ionindexrange=1:2)

    @test 0 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        scanindexrange=1:1, ionindexrange=1:1)
    @test 12 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        scanindexrange=1:1, ionindexrange=1:2)
    @test 956 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        scanindexrange=1:2, ionindexrange=2:2)
    @test 34 == maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        scanindexrange=2:2, ionindexrange=1:1)

    # Same element return type as the one that was used to construct the object
    @test Int == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1])))
    @test Int == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=1:1))
    @test Int == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=1:2))
    @test Float64 == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1])))
    @test Float64 == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=1:1))
    @test Float64 == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=1:2))
    @test Int == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), ionindexrange=1:1))
    @test Int == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), ionindexrange=1:2))
    @test Float64 == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), ionindexrange=1:1))
    @test Float64 == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), ionindexrange=1:2))
    @test Int == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=1:1, ionindexrange=1:1))
    @test Int == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=1:2, ionindexrange=1:2))
    @test Float64 == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=1:1, ionindexrange=1:1))
    @test Float64 == typeof(maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=1:2, ionindexrange=1:2))

    # Check for BoundsErrors
    @test_throws BoundsError maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:3)
    @test_throws BoundsError maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:4)
    @test_throws BoundsError maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:4)
    @test_throws BoundsError maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ionindexrange=0:2)
    @test_throws BoundsError maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ionindexrange=1:3)
    @test_throws BoundsError maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ionindexrange=0:3)
    @test_throws BoundsError maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:4, ionindexrange=0:3)
    @test_throws BoundsError maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:4, ionindexrange=1:3)
    @test_throws BoundsError maxintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:3, ionindexrange=0:2)  
end


############################################################################################
# minintensity(chrom::AbstractChromMS[; scanindexrange, ionindexrange])
############################################################################################
@testset "minintensity ChromMS [, greaterthan::Real; scanindexrange, ionindexrange]" begin
    # Validate the returned value
    @test 0 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]))
    @test 34 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        scanindexrange=2:2)
    @test 0 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        scanindexrange=1:2)
    @test 0 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]))
    @test 34 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        scanindexrange=2:2)
    @test 0 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        scanindexrange=1:2)

    @test 0 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        ionindexrange=1:1)
    @test 1 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        ionindexrange=2:2)
    @test 0 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        ionindexrange=1:1)
    @test 1 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        ionindexrange=2:2)

    @test 0 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        scanindexrange=1:1, ionindexrange=1:1)
    @test 34 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        scanindexrange=2:2, ionindexrange=1:2)
    @test 12 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Int[0 12; 34 956; 23 1]), 
        scanindexrange=1:2, ionindexrange=2:2)
    @test 34 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 
        scanindexrange=2:2, ionindexrange=1:1)

    @test 1 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 0)
    @test 12 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 0, 
        scanindexrange=1:2)
    @test 23 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 0, 
        ionindexrange=1:1)
    @test 34 == minintensity(ChromMS([1, 2, 3]u"s", [85, 100], Float64[0 12; 34 956; 23 1]), 0, 
        scanindexrange=1:2, ionindexrange=1:1)
    
    @test nothing === minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 23 1; 34 956]), 956)
    @test nothing === minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 23 1; 34 956]), 23, scanindexrange=1:2)
    @test nothing === minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 23 1; 34 956]), 34, ionindexrange=1:1)
    @test nothing === minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 23 1; 34 956]), 23, scanindexrange=1:2, ionindexrange=1:1)
    
    # Same element return type as the one that was used to construct the object
    @test Int == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1])))
    @test Int == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=1:1))
    @test Int == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=1:2))
    @test Float64 == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1])))
    @test Float64 == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=1:1))
    @test Float64 == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=1:2))
    @test Int == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), ionindexrange=1:1))
    @test Int == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), ionindexrange=1:2))
    @test Float64 == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), ionindexrange=1:1))
    @test Float64 == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), ionindexrange=1:2))
    @test Int == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=1:1, ionindexrange=1:1))
    @test Int == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int[0 12; 34 956; 23 1]), scanindexrange=1:2, ionindexrange=1:2))
    @test Float64 == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=1:1, ionindexrange=1:1))
    @test Float64 == typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]), scanindexrange=1:2, ionindexrange=1:2))
    @test Nothing === typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 23 1; 34 956]), 956))
    @test Nothing === typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 23 1; 34 956]), 23, scanindexrange=1:2))
    @test Nothing === typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 23 1; 34 956]), 34, ionindexrange=1:1))
    @test Nothing === typeof(minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 23 1; 34 956]), 23, scanindexrange=1:2, ionindexrange=1:1))

    # Check for BoundsErrors
    @test_throws BoundsError minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:3)
    @test_throws BoundsError minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:4)
    @test_throws BoundsError minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:4)
    @test_throws BoundsError minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ionindexrange=0:2)
    @test_throws BoundsError minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ionindexrange=1:3)
    @test_throws BoundsError minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ionindexrange=0:3)
    @test_throws BoundsError minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:4, ionindexrange=0:3)
    @test_throws BoundsError minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=1:4, ionindexrange=1:3)
    @test_throws BoundsError minintensity(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), scanindexrange=0:3, ionindexrange=0:2)
end


############################################################################################
# totalionchromatogram(chrom::AbstractChromMS)
############################################################################################
@testset "totalionchromatogram ChromMS" begin
    # Check the returned type and associated supertypes
    @test isa(totalionchromatogram(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])), 
        Chrom)
    @test isa(totalionchromatogram(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])), 
        AbstractChrom)
    @test isa(totalionchromatogram(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])), 
        AbstractChromatogram)

    # Check the contents of the Chrom for plausibility
    @test [1, 2, 3]u"s" == scantimes(totalionchromatogram(ChromMS([1, 2, 3]u"s", [85, 100],
        [0 12; 34 956; 23 1])))
    @test [12, 990, 24] == intensities(totalionchromatogram(ChromMS([1, 2, 3]u"s", [85, 100],
        [0 12; 34 956; 23 1])))
    @test Dict(:id => 4, "name" => "sample") == metadata(totalionchromatogram(ChromMS(
        [1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1], Dict(:id => 4, "name" => "sample")
        )))
    
    # Verify that the element type of the Chrom intensities is the same as the ChromMS
    @test Vector{Int} == typeof(intensities(totalionchromatogram(ChromMS([1, 2, 3]u"s", 
        [85, 100], Int[0 12; 34 956; 23 1]))))
    @test Vector{Float64} == typeof(intensities(totalionchromatogram(ChromMS([1, 2, 3]u"s", 
        [85, 100], Float64[0 12; 34 956; 23 1]))))
end
