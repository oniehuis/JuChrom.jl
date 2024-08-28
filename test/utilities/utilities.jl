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
# JuChrom.metricprefix(number::Real) -> Tuple
############################################################################################
@testset "metricprefix" begin
    @test ("Q", 30) == JuChrom.metricprefix(10.0^30)
    @test ("R", 27) == JuChrom.metricprefix(prevfloat(10.0^30))
    @test ("Q", 30) == JuChrom.metricprefix(-10.0^30)
    @test ("", 0) == JuChrom.metricprefix(1.0)
    @test ("m", -3) == JuChrom.metricprefix(prevfloat(1.0))
    @test ("", 0) == JuChrom.metricprefix(-1.0)
    @test ("m", -3) == JuChrom.metricprefix(prevfloat(-1.0))
    @test ("R", 27) == JuChrom.metricprefix(-prevfloat(10.0^30))
    @test ("r", -27) == JuChrom.metricprefix(10.0^-27)
    @test ("q", -30) == JuChrom.metricprefix(prevfloat(10.0^-27))
    @test ("Q", 30) == JuChrom.metricprefix(-10.0^30)
    @test ("R", 27) == JuChrom.metricprefix(-prevfloat(10.0^30))
    @test ("Q", 30) == JuChrom.metricprefix(10.0^234)
    @test ("q", -30) == JuChrom.metricprefix(10.0^-234)
    @test ("", 0) == JuChrom.metricprefix(-Inf)
    @test ("", 0) == JuChrom.metricprefix(+Inf)
    @test ("", 0) == JuChrom.metricprefix(0)
    @test ("", 0) == JuChrom.metricprefix(NaN)

    @test String == typeof(first(JuChrom.metricprefix(10^0)))
    @test String == typeof(first(JuChrom.metricprefix(10^3)))
    @test String == typeof(first(JuChrom.metricprefix(10^234)))
    @test String == typeof(first(JuChrom.metricprefix(10^-234)))
    @test String == typeof(first(JuChrom.metricprefix(-Inf)))
    @test String == typeof(first(JuChrom.metricprefix(+Inf)))
    @test String == typeof(first(JuChrom.metricprefix(NaN)))

    @test typeof(last(JuChrom.metricprefix(10^0))) <: Int
    @test typeof(last(JuChrom.metricprefix(10^3))) <: Int
    @test typeof(last(JuChrom.metricprefix(10^234))) <: Int
    @test typeof(last(JuChrom.metricprefix(10^-234))) <: Int
    @test typeof(last(JuChrom.metricprefix(-Inf))) <: Int
    @test typeof(last(JuChrom.metricprefix(+Inf))) <: Int
    @test typeof(last(JuChrom.metricprefix(NaN))) <: Int

end


############################################################################################
# JuChrom.name(::Type)
############################################################################################
@testset "name" begin
    @test FID == JuChrom.name(typeof(FID([1, 2, 3]u"s", [12, 956, 23])))
    @test GCMS == JuChrom.name(typeof(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Array == JuChrom.name(typeof(Int[1, 2, 3]))
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
 