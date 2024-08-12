using JuChrom
using Test
using Unitful: 𝐓


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
