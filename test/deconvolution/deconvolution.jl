using JuChrom
using Test
using Unitful: 𝐓

# ############################################################################################
# # Base.iterate(LM::LocalMaxima, currentindex=LM.startindex)
# ############################################################################################
# @testset "LocalMaxima" begin
#     maxima = []
#     for lm in JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7])
#         push!(maxima, lm)
#     end
#     @test [3:3, 5:7] == maxima

#     empty!(maxima)
#     for lm in JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, stopindex=4)
#         push!(maxima, lm)
#     end
#     @test [3:3] == maxima
# end


# ############################################################################################
# # JuChrom.LocalMaxima(values::AbstractVector{<:Real}; 
# # startindex::Integer=firstindex(values), stopindex::Integer=lastindex(values))
# ############################################################################################
# @testset "LocalMaxima" begin
#     @test [4, 3, 5, 3, 6, 6, 6, 4, 7] == JuChrom.LocalMaxima(
#         [4, 3, 5, 3, 6, 6, 6, 4, 7]).values
#     @test 1 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7]).startindex
#     @test 9 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7]).stopindex
#     @test 9 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7]).length
#     @test [4, 3, 5, 3, 6, 6, 6, 4, 7] == JuChrom.LocalMaxima(
#         [4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, stopindex=8).values
#     @test 2 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, 
#         stopindex=8).startindex
#     @test 8 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, 
#         stopindex=8).stopindex
#     @test 9 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, 
#         stopindex=8).length

#     @test Vector{Int} == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7]).values)
#     @test Int == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7]).startindex)
#     @test UInt8 == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=0x1).startindex)
#     @test UInt8 == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         stopindex=0x9).stopindex)
#     @test Int == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7]).length)

#     @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=0)
#     @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=10)
#     @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         stopindex=0)
#     @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         stopindex=10)
#     @test_throws ArgumentError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=1, stopindex=2)
#     @test_throws ArgumentError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=8)
#     @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=0, stopindex=1)
#     @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=9, stopindex=10)
#  end


# ############################################################################################
# # JuChrom.nextlocalmaximum(values::AbstractVector{<:Real}; 
# # startindex::Integer=firstindex(values), stopindex::Integer=lastindex(values))
# ############################################################################################
# @testset "nextlocalmaximum" begin
#     @test 3:3 == JuChrom.nextlocalmaximum(Int[4, 3, 5, 3, 6, 6, 6, 4, 7])
#     @test 3:3 == JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7])
#     @test 5:7 == JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=4)
#     @test nothing === JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=4, stopindex=7)
#     @test 3:3 === JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=2, stopindex=7)

#     @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=0)
#     @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=10)
#     @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         stopindex=0)
#     @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         stopindex=10)
#     @test_throws ArgumentError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=1, stopindex=2)
#     @test_throws ArgumentError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=8)
#     @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=0, stopindex=1)
#     @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
#         startindex=9, stopindex=10)
#  end
 

@testset "test_ls_fit" begin
    # Test 1: Basic linear relationship
    xs1 = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys1 = [2.0, 3.0, 4.0, 5.0, 6.0]
    intercept1, slope1 = JuChrom.ls_fit(xs1, ys1)
    @test intercept1 ≈ 1.0 atol=1e-10
    @test slope1 ≈ 1.0 atol=1e-10

    # Test 2: Negative slope
    xs2 = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys2 = [5.0, 4.0, 3.0, 2.0, 1.0]
    intercept2, slope2 = JuChrom.ls_fit(xs2, ys2)
    @test intercept2 ≈ 6.0 atol=1e-10
    @test slope2 ≈ -1.0 atol=1e-10

    # Test 3: Zero slope (horizontal line)
    xs3 = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys3 = [3.0, 3.0, 3.0, 3.0, 3.0]
    intercept3, slope3 = JuChrom.ls_fit(xs3, ys3)
    @test intercept3 ≈ 3.0 atol=1e-10
    @test slope3 ≈ 0.0 atol=1e-10

    # Test 4: Insufficient elements in vectors (should trigger an error)
    xs4 = [1.0]
    ys4 = [2.0]
    @test_throws ArgumentError JuChrom.ls_fit(xs4, ys4)

    # Test 5: Vectors of different lengths (should throw an error)
    xs5 = [1.0, 2.0, 3.0]
    ys5 = [1.0, 2.0]
    @test_throws ArgumentError JuChrom.ls_fit(xs5, ys5)
end

@testset "nnls" begin
    # 1. Test with a basic system (simple solution)
    A = [1.0 2.0; 3.0 4.0]
    b = [1.0, 1.0]

    x = JuChrom.nnls(A, b)
    expected_x = [0.0, 0.30]  # known solution
    @test isapprox(x, expected_x, atol=1e-10)
end

# Test 1: Basic non-negative solution
@testset "Basic NILS Test" begin
    A = [1.0 2.0; 3.0 4.0]
    b = [1.0, 2.0]
    
    # Run the NILS function
    x = JuChrom.nnls(A, b)
    
    # Expected solution should be non-negative
    @test all(x .>= 0)  # Check non-negativity of solution

    # Validate if Ax is approximately equal to b
    @test isapprox(A * x, b, atol=1e-10)
end

# Test 2: Zero vector as b
@testset "Zero Vector NILS Test" begin
    A = [1.0 2.0; 3.0 4.0]
    b = [0.0, 0.0]
    
    # Run the NILS function
    x = JuChrom.nnls(A, b)
    
    # Expected solution should be zero since b is zero
    @test x ≈ zeros(size(A, 2))
end

# Test 3: Infeasible solution (forcing zero in the solution)
@testset "Infeasible NILS Test" begin
    A = [1.0 0.0; 0.0 1.0]
    b = [-1.0, -1.0]
    
    # Run the NILS function
    x = JuChrom.nnls(A, b)
    
    # Expected solution should be zero as no negative solutions are allowed
    @test x ≈ zeros(size(A, 2))
end

# Test 4: Random matrix and vector test
@testset "Random Matrix nnls Test" begin    
    A = randn(100, 200)
    b = randn(100)

    # Run the NILS function
    x = JuChrom.nnls(A, b)
    
    # Solution must satisfy non-negativity constraint
    @test all(x .>= 0)
end
