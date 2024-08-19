using JuChrom
using Test
using Unitful: 𝐓

@testset "FID constructor" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == FID([1, 2, 3]u"s", [12, 956, 23]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(FID(Int64[1, 2, 3]u"s", [12, 956, 23]).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(FID(Float64[1, 2, 3]u"s", 
            [12, 956, 23]).scantimes)

    # Intensities
    @test [12, 956, 23] == FID([1, 2, 3]u"s", [12, 956, 23]).intensities
    @test Vector{Int64} == typeof(FID([1, 2, 3]u"s", Int64[12, 956, 23]).intensities)
    @test Vector{Float64} == typeof(FID([1, 2, 3]u"s",  Float64[12, 956, 23]).intensities)

    # Metadata
    @test Dict{Any, Any} == typeof(FID([1, 2, 3]u"s", [12, 956, 23]).metadata)
    @test Dict() == FID([1, 2, 3]u"s", [12, 956, 23]).metadata
    @test Dict(:id => 4, "name" => "sample") == FID([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 4, "name" => "sample")).metadata

    # Check the associated supertypes
    @test isa(FID([1, 2, 3]u"s", [12, 956, 23]), AbstractChromatogram)
    @test isa(FID([1, 2, 3]u"s", [12, 956, 23]), AbstractGC)
    @test isa(FID([1, 2, 3]u"s", [12, 956, 23]), AbstractFID)

    # Scan time vector accepts only time values
    @test_throws MethodError FID([1, 2, 3], [12, 956, 23])
    @test_throws MethodError FID([1, 2, 3]u"m", [12, 956, 23])
    
    # Intensity vector must not accept a unit
    @test_throws MethodError FID([1, 2, 3], [12, 956, 23]u"s")

    # Metadata cannot be anything other than a dictionary
    @test_throws MethodError FID([1, 2, 3]u"s", [12, 956, 23], :sample_name)
    @test_throws MethodError FID([1, 2, 3]u"s", [12, 956, 23], 1)
    @test_throws MethodError FID([1, 2, 3]u"s", [12, 956, 23], 'S')

    # Must provide at least on scantime
    @test_throws ArgumentError FID(Int[]u"s", Int[])
    @test_throws ArgumentError FID(Int[]u"s", [12])
    @test [1]u"s" == FID([1]u"s", [12]).scantimes

    # Must provide at least on intensity
    @test_throws ArgumentError FID([1]u"s", Int[])

    # Number of scantimes and intensities must be identical
    @test_throws DimensionMismatch FID([1, 2, 3, 4]u"s", [12, 956, 23])
    @test_throws DimensionMismatch FID([1, 2]u"s", [12, 956, 23])

    # Scan times must be in ascending order
    @test_throws ArgumentError FID([2, 1, 3]u"s", [12, 956, 23])

    # Intensity values cannot be less than zero
    @test_throws ArgumentError FID([1, 2, 3]u"s", [-1, 956, 23])

    # FID is broadcastable
    @test [1, 2]u"s" == ((fid, i) -> fid.scantimes[i]).(FID([1, 2, 3]u"s", [12, 956, 23]), 
        [1, 2])

end


@testset "show FID" begin
    io = IOBuffer()
    show(io, FID(Int64[1]u"s", Int64[12]))
    @test String(take!(io)) == string(
        "FID {scan times: Int64, intensities: Int64}\n",
        "1 scan; scan time: 1 s\n",
        "intensity: 12\n",
        "metadata: 0 entries")
    show(io, FID(Int64[1, 2]u"s", Int64[12, 13]))
    @test String(take!(io)) == string(
        "FID {scan times: Int64, intensities: Int64}\n",
        "2 scans; scan times: 1 s, 2 s\n",
        "intensity range: 12 - 13\n",
        "metadata: 0 entries")
    show(io, FID(convert(Vector{Int64}, (1:10))u"s", convert(Vector{Int64}, (12:21))))
    @test String(take!(io)) == string(
        "FID {scan times: Int64, intensities: Int64}\n",
        "10 scans; scan times: 1 s, 2 s, 3 s, 4 s, 5 s, 6 s, 7 s, 8 s, 9 s, 10 s\n",
        "intensity range: 12 - 21\n",
        "metadata: 0 entries")
    show(io, FID(convert(Vector{Int64}, (1:11))u"s", convert(Vector{Int64}, (12:22))))
    @test String(take!(io)) == string(
        "FID {scan times: Int64, intensities: Int64}\n",
        "11 scans; scan time range: 1 s - 11 s\n",
        "intensity range: 12 - 22\n",
        "metadata: 0 entries")
    show(io, FID(Int64[1, 2, 3]u"s", Int64[12, 956, 23], Dict(:id => 4)))
    @test String(take!(io)) == string(
        "FID {scan times: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "intensity range: 12 - 956\n",
        "metadata: 1 entry")
    show(io, FID(Int64[1, 2, 3]u"s", Int64[12, 956, 23], Dict(:id => 4, 
        "name" => "sample")))
    @test String(take!(io)) == string(
        "FID {scan times: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "intensity range: 12 - 956\n",
        "metadata: 2 entries")
end
