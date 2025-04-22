using JuChrom
using Test


############################################################################################
# JuChrom.copy_with_eltype(array::AbstractArray, elementtype::Type)
############################################################################################
@testset "copy_with_eltype" begin
   @test [1, 2, 3, 4, 5, 6] == JuChrom.copy_with_eltype(Int[1, 2, 3, 4, 5, 6], Float64)
   @test [1, 2, 3, 4, 5, 6] == JuChrom.copy_with_eltype(Float64[1, 2, 3, 4, 5, 6], Int)
   @test Vector{Float64} == typeof(JuChrom.copy_with_eltype(Int[1, 2, 3, 4, 5, 6], 
    Float64))
   @test Vector{Float64} == typeof(JuChrom.copy_with_eltype(1:6, Float64))
   @test Vector{Int} == typeof(JuChrom.copy_with_eltype(1.0:6.0, Int))
   @test Vector{Int} == typeof(JuChrom.copy_with_eltype(Float64[1, 2, 3, 4, 5, 6], Int))
   @test Matrix{Int} == typeof(JuChrom.copy_with_eltype(Float64[1 2; 3 4; 5 6], Int))
   @test Matrix{Float64} == typeof(JuChrom.copy_with_eltype(Int[1 2; 3 4; 5 6], Float64))
   @test_throws InexactError JuChrom.copy_with_eltype(Float64[1.1, 2, 3, 4, 5, 6], Int)
end


############################################################################################
# JuChrom.findclosest(A::AbstractArray{<:Number}, x::Number) -> Int
############################################################################################
@testset "findclosest" begin
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

 
############################################################################################
# JuChrom.invert(dictionary::Dict)
############################################################################################
@testset "invert" begin
    @test Dict(1.0 => [:a, :d], 2.0 => [:b, :c]) == JuChrom.invert(
        Dict(:a => 1.0, :b => 2.0, :c => 2.0, :d => 1.0))
    @test Dict(1.0 => [4, 1], 2.0 => [2, 3]) == JuChrom.invert(
        Dict(1 => 1.0, 2 => 2.0, 3 => 2.0, 4 => 1.0))
    @test Dict(:a => [4, 1], :b => [2, 3]) == JuChrom.invert(
        Dict(1 => :a, 2 => :b, 3 => :b, 4 => :a))

    @test Dict{Float64, Vector{Symbol}} == typeof(JuChrom.invert(
        Dict(:a => 1.0, :b => 2.0, :c => 2.0, :d => 1.0)))
    @test Dict{Float64, Vector{Int}} == typeof(JuChrom.invert(
        Dict{Int, Float64}(1 => 1.0, 2 => 2.0, 3 => 2.0, 4 => 1.0)))
    @test Dict{Symbol, Vector{Int}} == typeof(JuChrom.invert(
        Dict{Int, Symbol}(1 => :a, 2 => :b, 3 => :b, 4 => :a)))
end


############################################################################################
# JuChrom.name(::Type)
############################################################################################
@testset "name" begin
    @test Chrom == JuChrom.name(typeof(Chrom([1, 2, 3]u"s", [12, 956, 23])))
    @test ChromMS == JuChrom.name(typeof(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Array == JuChrom.name(typeof(Int[1, 2, 3]))
end
