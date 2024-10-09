using JuChrom
using Test
using Unitful: 𝐓


############################################################################################
# Base.iterate(LM::LocalMaxima, currentindex=LM.startindex)
############################################################################################
@testset "LocalMaxima" begin
    maxima = []
    for lm in JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7])
        push!(maxima, lm)
    end
    @test [3:3, 5:7] == maxima

    empty!(maxima)
    for lm in JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, stopindex=4)
        push!(maxima, lm)
    end
    @test [3:3] == maxima
end


############################################################################################
# JuChrom.LocalMaxima(values::AbstractVector{<:Real}; 
# startindex::Integer=firstindex(values), stopindex::Integer=lastindex(values))
############################################################################################
@testset "LocalMaxima" begin
    @test [4, 3, 5, 3, 6, 6, 6, 4, 7] == JuChrom.LocalMaxima(
        [4, 3, 5, 3, 6, 6, 6, 4, 7]).values
    @test 1 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7]).startindex
    @test 9 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7]).stopindex
    @test 9 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7]).length
    @test [4, 3, 5, 3, 6, 6, 6, 4, 7] == JuChrom.LocalMaxima(
        [4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, stopindex=8).values
    @test 2 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, 
        stopindex=8).startindex
    @test 8 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, 
        stopindex=8).stopindex
    @test 9 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, 
        stopindex=8).length

    @test Vector{Int} == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7]).values)
    @test Int == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7]).startindex)
    @test UInt8 == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0x1).startindex)
    @test UInt8 == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=0x9).stopindex)
    @test Int == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7]).length)

    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=10)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=0)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=10)
    @test_throws ArgumentError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=1, stopindex=2)
    @test_throws ArgumentError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=8)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0, stopindex=1)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=9, stopindex=10)
 end


############################################################################################
# JuChrom.nextlocalmaximum(values::AbstractVector{<:Real}; 
# startindex::Integer=firstindex(values), stopindex::Integer=lastindex(values))
############################################################################################
@testset "nextlocalmaximum" begin
    @test 3:3 == JuChrom.nextlocalmaximum(Int[4, 3, 5, 3, 6, 6, 6, 4, 7])
    @test 3:3 == JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7])
    @test 5:7 == JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=4)
    @test nothing === JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=4, stopindex=7)
    @test 3:3 === JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=2, stopindex=7)

    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=10)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=0)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=10)
    @test_throws ArgumentError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=1, stopindex=2)
    @test_throws ArgumentError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=8)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0, stopindex=1)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=9, stopindex=10)
 end
 

############################################################################################
# JuChrom.lsfit(xs, ys)
############################################################################################
@testset "lsfit" begin
    # Basic linear relationship
    xs1 = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys1 = [2.0, 3.0, 4.0, 5.0, 6.0]
    intercept1, slope1 = JuChrom.Deconvolution.lsfit(xs1, ys1)
    @test intercept1 ≈ 1.0 atol=1e-10
    @test slope1 ≈ 1.0 atol=1e-10

    # Negative slope
    xs2 = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys2 = [5.0, 4.0, 3.0, 2.0, 1.0]
    intercept2, slope2 = JuChrom.Deconvolution.lsfit(xs2, ys2)
    @test intercept2 ≈ 6.0 atol=1e-10
    @test slope2 ≈ -1.0 atol=1e-10

    # Zero slope (horizontal line)
    xs3 = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys3 = [3.0, 3.0, 3.0, 3.0, 3.0]
    intercept3, slope3 = JuChrom.Deconvolution.lsfit(xs3, ys3)
    @test intercept3 ≈ 3.0 atol=1e-10
    @test slope3 ≈ 0.0 atol=1e-10

    # Insufficient elements in vectors (should trigger an error)
    xs4 = [1.0]
    ys4 = [2.0]
    @test_throws ArgumentError JuChrom.Deconvolution.lsfit(xs4, ys4)

    # Vectors of different lengths (should throw an error)
    xs5 = [1.0, 2.0, 3.0]
    ys5 = [1.0, 2.0]
    @test_throws ArgumentError JuChrom.Deconvolution.lsfit(xs5, ys5)
end


############################################################################################
# JuChrom.nnls(A, b)
############################################################################################
@testset "nnls" begin
    # Test with a basic system (simple solution)
    x1 = JuChrom.nnls([1.0 2.0; 3.0 4.0], [1.0, 1.0])
    @test x1 ≈ [0.0, 0.30] atol=1e-10

    # Basic non-negative solution
    A2 = [1.0 2.0; 3.0 4.0]
    b2 = [1.0, 2.0]
    x2 = JuChrom.nnls(A2, b2)
    @test all(x2 .>= 0)  # check non-negativity of solution
    @test A2 * x2 ≈ b2 atol=1e-10

    # Zero vector as b
    A3 = [1.0 2.0; 3.0 4.0]
    x3 = JuChrom.nnls(A3, [0.0, 0.0])
    @test x3 ≈ zeros(size(A3, 2)) atol=1e-10 # solution should be 0 since b is 0

    # Infeasible solution (forcing zero in the solution)
    A4 = [1.0 0.0; 0.0 1.0]
    x4 = JuChrom.nnls(A4, [-1.0, -1.0])
    @test x4 ≈ zeros(size(A4, 2)) atol=1e-10

    # Random matrix and vector test
    x5 = JuChrom.nnls(randn(100, 200), randn(100))
    @test all(x5 .>= 0)  # check non-negativity of solution
end
