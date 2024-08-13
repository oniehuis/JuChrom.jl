using JuChrom
using Test



############################################################################################
# JuChrom.copy_with_eltype(array::AbstractArray, elementtype::Type)
############################################################################################
@testset "copy_with_eltype" begin
   @test [1, 2, 3, 4, 5, 6] == JuChrom.copy_with_eltype(Int[1, 2, 3, 4, 5, 6], Float64)
   @test [1, 2, 3, 4, 5, 6] == JuChrom.copy_with_eltype(Float64[1, 2, 3, 4, 5, 6], Int)
   @test Vector{Float64} == typeof(JuChrom.copy_with_eltype(Int[1, 2, 3, 4, 5, 6], Float64))
   @test Vector{Float64} == typeof(JuChrom.copy_with_eltype(1:6, Float64))
   @test Vector{Int} == typeof(JuChrom.copy_with_eltype(1.0:6.0, Int))
   @test Vector{Int} == typeof(JuChrom.copy_with_eltype(Float64[1, 2, 3, 4, 5, 6], Int))
   @test Matrix{Int} == typeof(JuChrom.copy_with_eltype(Float64[1 2; 3 4; 5 6], Int))
   @test Matrix{Float64} == typeof(JuChrom.copy_with_eltype(Int[1 2; 3 4; 5 6], Float64))
   @test_throws InexactError JuChrom.copy_with_eltype(Float64[1.1, 2, 3, 4, 5, 6], Int)
end


############################################################################################
# cosine(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
############################################################################################
@testset "cosine" begin
    @test 0.8978872704229618 ≈ cosine([100, 500, 250], [200, 1000, 0])
    @test 1 ≈ cosine([1, 1, 1], [1, 1, 1])
    @test 0 ≈ cosine([10.5, 10.5, 10.5], [-10.5, -10.5, -10.5])
    @test Float64 == typeof(cosine([100, 500, 250], [200, 1000, 0]))
    @test Float64 == typeof(cosine([1, 1, 1], [1, 1, 1]))
    @test Float64 == typeof(cosine([10.5, 10.5, 10.5], [-10.5, -10.5, -10.5]))
    @test_throws DimensionMismatch cosine([100, 500, 250], [200, 1000, 0, 10])
    @test_throws ArgumentError cosine(Int[], Int[])
    @test_throws ArgumentError cosine([0, 0, 0], [200, 1000, 0])
    @test_throws ArgumentError cosine([100, 500, 250], [0, 0, 0])
 end


############################################################################################
# JuChrom.findclosest(A::AbstractArray{<:Number}, x::Number) -> Int
############################################################################################
@testset "copy_with_eltype" begin
    @test 3 == JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 0)
    @test 5 == JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 1.5)
    @test 2 == JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], -1.5)
    @test 1 == JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], -2.5)
    @test 4 == JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 10//10)
    @test 4 == JuChrom.findclosest(Float64[-2.1, -1.1, 0.1, 1.1, 2.1, 3.1, 4.1, 5.1], 
        10//10)
    @test 8 == JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 6)
    @test 8 == JuChrom.findclosest(-2:5, 6)
    @test 8 == JuChrom.findclosest(-2.1:1.0:5.1, 6)
    @test 8 == JuChrom.findclosest(-2.1:1.0:5.1, 12//2)

    @test Int == typeof(JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 0))
    @test Int == typeof(JuChrom.findclosest(Float32[-2, -1, 0, 1, 2, 3, 4, 5], 1.5))

    @test_throws MethodError JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], :a)
    @test 4 == JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], true)
 end


############################################################################################
# integer(value:::Real; start::Real=0.7) -> Int
############################################################################################
@testset "integer" begin
    @test 29 == integer(29)
    @test 30 == integer(29.7)
    @test 29 == integer(prevfloat(29.7))
    @test 30 == integer(29.7f0)
    @test 29 == integer(prevfloat(29.7f0))
    @test 29 == integer(29, start=0.8)
    @test 30 == integer(29.8, start=0.8)
    @test 29 == integer(prevfloat(29.8), start=0.8)
    @test 29 == integer(29, start=0)
    @test 29 == integer(29.0, start=0)
    @test 28 == integer(prevfloat(29.0), start=0)

    @test Int == typeof(integer(29))
    @test Int == typeof(integer(29.7))
    @test Int == typeof(integer(29.7f0))
    @test Int == typeof(integer(29, start=0.8))
    @test Int == typeof(integer(29.7, start=0.8))
    @test Int == typeof(integer(29.7f0, start=0.8))
    @test Int == typeof(integer(29, start=0))
    @test Int == typeof(integer(29.7, start=0))
    @test Int == typeof(integer(29.7f0, start=0))

    @test_throws ArgumentError integer(30.79, start=-0.1)
    @test_throws ArgumentError integer(30.0, start=1)
    @test_throws ArgumentError integer(30.0, start=1.1)
 end