using JuChrom
using Test
using Unitful: 𝐓

@testset "buildxic" begin
    # Verify that the function returns teh expected contents
    @test ([84.9, 85.1, 85.2, 99.9, 100.2, 100.6, 112.1], 
        [0 12 0 0 234 0 0; 0 0 23 324 0 0 45422; 21 0 0 0 0 523 0]
        ) == JuChrom.InputOutput.buildxic([2, 3, 2], 
        [85.1, 100.2, 85.2, 99.9, 112.1, 84.9, 100.6], [12, 234, 23, 324, 45422, 21, 523])
    @test ([84, 85, 99, 100, 112], [0 12 0 234 0; 0 23 324 0 45422; 21 0 0 523 0]
        ) == JuChrom.InputOutput.buildxic([2, 3, 2], [85, 100, 85, 99, 112, 84, 100], 
        [12, 234, 23, 324, 45422, 21, 523])

    # Verify that the function returns teh expected types
    @test Tuple{Vector{Float64}, Matrix{Int}} == typeof(JuChrom.InputOutput.buildxic([2, 3, 
        2], [85.1, 100.2, 85.2, 99.9, 112.1, 84.9, 100.6], Int[12, 234, 23, 324, 45422, 21, 
        523]))
    @test Tuple{Vector{Float64}, Matrix{Float64}} == typeof(JuChrom.InputOutput.buildxic(
        [2, 3, 2], [85.1, 100.2, 85.2, 99.9, 112.1, 84.9, 100.6], Float64[12, 234, 23, 324, 
        45422, 21, 523]))
    @test Tuple{Vector{Int}, Matrix{Float64}} == typeof(JuChrom.InputOutput.buildxic(
        [2, 3, 2], [85, 100, 85, 99, 112, 84, 100], Float64[12, 234, 23, 324, 
        45422, 21, 523]))
    
    # Verify that the function errors of vector counts are incompatible
    # ArgumentError: ion and intensity values differ in their counts
    @test_throws ArgumentError JuChrom.InputOutput.buildxic([2, 3, 2], [85, 100, 85, 99, 
        112, 84, 100], [12, 234, 23, 324, 45422, 21])
    @test_throws ArgumentError JuChrom.InputOutput.buildxic([2, 3, 2], [85, 100, 85, 99, 
        112, 84], [12, 234, 23, 324, 45422, 21, 523])

    # ArgumentError: point counts incompatible with ion and intensity value counts
    @test_throws ArgumentError JuChrom.InputOutput.buildxic([2, 3], [85, 100, 85, 99, 112, 
        84, 100], [12, 234, 23, 324, 45422, 21, 523])
    @test_throws ArgumentError JuChrom.InputOutput.buildxic([2, 3, 2], [85, 100, 85, 99, 
        112, 84], [12, 234, 23, 324, 45422, 21])
    @test_throws ArgumentError JuChrom.InputOutput.buildxic([2, 3, 3], [85, 100, 85, 99, 
        112, 84, 100], [12, 234, 23, 324, 45422, 21, 523])
    @test_throws ArgumentError JuChrom.InputOutput.buildxic([1, 3, 2], [85, 100, 85, 99, 
        112, 84, 100], [12, 234, 23, 324, 45422, 21, 523])

    # Returned list of ions has ions in ascending order
    mzs, xic = JuChrom.InputOutput.buildxic([2, 3, 2], [85, 100, 85, 99, 112, 84, 100], 
        [12, 234, 23, 324, 45422, 21, 523])
    @test issorted(mzs)
    mzs, xic = JuChrom.InputOutput.buildxic([2, 3, 2], [100, 85, 112, 99, 85, 100, 84], 
        [12, 234, 23, 324, 45422, 21, 523])
    @test issorted(mzs)
end
