using JuChrom
using Test
using Unitful: 𝐓

# @testset "copy_with_eltype" begin
#    @test [1, 2, 3, 4, 5, 6] == JuChrom.copy_with_eltype(Int[1, 2, 3, 4, 5, 6], Float64)
#    @test [1, 2, 3, 4, 5, 6] == JuChrom.copy_with_eltype(Float64[1, 2, 3, 4, 5, 6], Int)
#    @test Vector{Float64} == typeof(JuChrom.copy_with_eltype(Int[1, 2, 3, 4, 5, 6], Float64))
#    @test Vector{Float64} == typeof(JuChrom.copy_with_eltype(1:6, Float64))
#    @test Vector{Int} == typeof(JuChrom.copy_with_eltype(1.0:6.0, Int))
#    @test Vector{Int} == typeof(JuChrom.copy_with_eltype(Float64[1, 2, 3, 4, 5, 6], Int))
#    @test Matrix{Int} == typeof(JuChrom.copy_with_eltype(Float64[1 2; 3 4; 5 6], Int))
#    @test Matrix{Float64} == typeof(JuChrom.copy_with_eltype(Int[1 2; 3 4; 5 6], Float64))
#    @test_throws InexactError JuChrom.copy_with_eltype(Float64[1.1, 2, 3, 4, 5, 6], Int)
# end