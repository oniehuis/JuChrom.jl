using JuChrom
using Test


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


@testset "File" begin
    # Test object construction
    @test JuChrom.InputOutput.File() == JuChrom.InputOutput.File()

    # Test Supertype
    @test isa(JuChrom.InputOutput.File(), JuChrom.InputOutput.Source)
end


@testset "Path" begin
    # Test object construction
    @test JuChrom.InputOutput.Path() == JuChrom.InputOutput.Path()

    # Test Supertype
    @test isa(JuChrom.InputOutput.Path(), JuChrom.InputOutput.Source)
end


@testset "IOError" begin
    # Test object construction
    @test JuChrom.InputOutput.IOError() == JuChrom.InputOutput.IOError()

    # Check if content is properly returned
    @test "IOError: " == JuChrom.InputOutput.IOError().msg
    @test "IOError: message" == JuChrom.InputOutput.IOError("message").msg

    # Check supertype
    @test isa(JuChrom.InputOutput.IOError(), JuChrom.InputOutput.IOError)

    # Test if proper Error is thrown
    @test_throws JuChrom.InputOutput.IOError throw(JuChrom.InputOutput.IOError())
    @test_throws JuChrom.InputOutput.IOError throw(JuChrom.InputOutput.IOError("message"))

    io = IOBuffer()
    Base.showerror(io, JuChrom.InputOutput.IOError("message"))
    @test String(take!(io)) == "IOError: message"
end


@testset "FileExistsError" begin
    # Test object construction
    @test JuChrom.InputOutput.FileExistsError() == JuChrom.InputOutput.FileExistsError()

    # Check if content is properly returned
    @test "FileExistsError: " == JuChrom.InputOutput.FileExistsError().msg
    @test "FileExistsError: message" == JuChrom.InputOutput.FileExistsError("message").msg

    # Check supertype
    @test isa(JuChrom.InputOutput.FileExistsError(), JuChrom.InputOutput.FileExistsError)

    # Test if proper Error is thrown
    @test_throws JuChrom.InputOutput.FileExistsError throw(
        JuChrom.InputOutput.FileExistsError())
    @test_throws JuChrom.InputOutput.FileExistsError throw(
        JuChrom.InputOutput.FileExistsError("message"))

    io = IOBuffer()
    Base.showerror(io, JuChrom.InputOutput.FileExistsError("message"))
    @test String(take!(io)) == "FileExistsError: message"
end
