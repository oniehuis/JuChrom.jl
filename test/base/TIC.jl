using JuChrom
using Test
using Unitful: 𝐓

@testset "TIC constructor" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == TIC([1, 2, 3]u"s", [12, 956, 23]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(TIC(Int64[1, 2, 3]u"s", [12, 956, 23]).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(TIC(Float64[1, 2, 3]u"s", 
            [12, 956, 23]).scantimes)

    # Intensities
    @test [12, 956, 23] == TIC([1, 2, 3]u"s", [12, 956, 23]).intensities
    @test Vector{Int64} == typeof(TIC([1, 2, 3]u"s", Int64[12, 956, 23]).intensities)
    @test Vector{Float64} == typeof(TIC([1, 2, 3]u"s",  Float64[12, 956, 23]).intensities)

     # Metadata
     @test Dict{Any, Any} == typeof(TIC([1, 2, 3]u"s", [12, 956, 23]).metadata)
     @test Dict() == TIC([1, 2, 3]u"s", [12, 956, 23]).metadata
     @test Dict(:id => 1, "name" => "sample") == TIC([1, 2, 3]u"s", [12, 956, 23], 
         Dict(:id => 1, "name" => "sample")).metadata

    # Check the associated supertypes
    @test isa(TIC([1, 2, 3]u"s", [12, 956, 23]), AbstractChromatogram)
    @test isa(TIC([1, 2, 3]u"s", [12, 956, 23]), AbstractGC)
    @test isa(TIC([1, 2, 3]u"s", [12, 956, 23]), AbstractTIC)

    # Scan time vector accepts only time values
    @test_throws MethodError TIC([1, 2, 3], [12, 956, 23])
    @test_throws MethodError TIC([1, 2, 3]u"m", [12, 956, 23])
    
    # Intensity vector must not accept a unit
    @test_throws MethodError TIC([1, 2, 3], [12, 956, 23]u"s")

    # Metadata cannot be anything other than a dictionary.
    @test_throws MethodError TIC([1, 2, 3]u"s", [12, 956, 23], :sample_name)
    @test_throws MethodError TIC([1, 2, 3]u"s", [12, 956, 23], 1)
    @test_throws MethodError TIC([1, 2, 3]u"s", [12, 956, 23], 'S')

    # Number of scantimes and intensities must be identical
    @test_throws DimensionMismatch TIC([1, 2, 3, 4]u"s", [12, 956, 23])
    @test_throws DimensionMismatch TIC([1, 2]u"s", [12, 956, 23])

    # Scan times must be in ascending order
    @test_throws ArgumentError TIC([2, 1, 3]u"s", [12, 956, 23])

    # Intensity values cannot be less than zero
    @test_throws ArgumentError TIC([1, 2, 3]u"s", [-12, 956, 23])

    # TIC is broadcastable
    @test [1, 2]u"s" == ((tic, i) -> tic.scantimes[i]).(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        [1, 2])
end


@testset "show TIC" begin
    io = IOBuffer()
    show(io, TIC(Int64[1, 2, 3]u"s", Int64[12, 956, 23], Dict(:id => 1, "name" => "sample")))
    @test String(take!(io)) == string(
        "TIC {scantimes: Int64, intensities: Int64}\n",
        "3 scans; time range: 1 s - 3 s\n",
        "intensity range: 12 - 956\n",
        "metadata: 2")
end
