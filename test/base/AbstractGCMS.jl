using JuChrom
using Test
using Unitful: 𝐓

@testset "ions GCMS" begin
    @test [85, 100] == ions(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test [85] == ions(GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1], (:,1))))
    
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

    @test_throws MethodError ions(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test_throws MethodError ions(TIC([1, 2, 3]u"s", [12, 956, 23]))
end


@testset "intensities GCMS" begin
    @test [0 12; 34 956; 23 1] == intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
    @test reshape([0, 956, 1], (:,1)) == intensities(GCMS([1, 2, 3]u"s", [85],
        reshape([0, 956, 1], (:,1))))
    @test Matrix{Int64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Int64[0 12; 34 956; 23 1])))
    @test Matrix{Float64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1])))
    @test Matrix{Int64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1],
         (:,1)))))
end
