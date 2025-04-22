using JuChrom
using Test
using Unitful: 𝐓

@testset "Chrom constructor" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == Chrom([1, 2, 3]u"s", [12, 956, 23]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(Chrom(Int64[1, 2, 3]u"s", [12, 956, 23]
        ).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(Chrom(Float64[1, 2, 3]u"s", 
            [12, 956, 23]).scantimes)

    # Intensities
    @test [12, 956, 23] == Chrom([1, 2, 3]u"s", [12, 956, 23]).intensities
    @test Vector{Int64} == typeof(Chrom([1, 2, 3]u"s", Int64[12, 956, 23]).intensities)
    @test Vector{Float64} == typeof(Chrom([1, 2, 3]u"s",  Float64[12, 956, 23]
        ).intensities)

    # Metadata
    @test Dict{Any, Any} == typeof(Chrom([1, 2, 3]u"s", [12, 956, 23]).metadata)
    @test Dict() == Chrom([1, 2, 3]u"s", [12, 956, 23]).metadata
    @test Dict(:id => 4, "name" => "sample") == Chrom([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 4, "name" => "sample")).metadata

    # RiMapper, incl. its mutability
    @test isa(Chrom([1, 2, 3]u"s", [12, 956, 23]).rimapper, Nothing)
    @test isa(Chrom([1, 2, 3]u"s", [12, 956, 23], rimapper=RiMapper("Kovats", 
        (1:5)u"minute", 1000:1000:5000)).rimapper, JuChrom.AbstractRiMapper)

    chrom = Chrom([1, 2, 3]u"s", [12, 956, 23])
    @test nothing === chrom.rimapper
    chrom.rimapper = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000)
    @test isa(chrom.rimapper, JuChrom.AbstractRiMapper)

    # Check the associated supertypes
    @test isa(Chrom([1, 2, 3]u"s", [12, 956, 23]), AbstractChromatogram)
    @test isa(Chrom([1, 2, 3]u"s", [12, 956, 23]), AbstractChrom)

    # Scan time vector accepts only time values
    @test_throws MethodError Chrom([1, 2, 3], [12, 956, 23])
    @test_throws MethodError Chrom([1, 2, 3]u"m", [12, 956, 23])
    
    # Intensity vector must not accept a unit
    @test_throws MethodError Chrom([1, 2, 3], [12, 956, 23]u"s")

    # Must provide at least one scantime
    @test_throws ArgumentError Chrom(Int[]u"s", Int[])
    @test_throws ArgumentError Chrom(Int[]u"s", [12])
    @test [1]u"s" == Chrom([1]u"s", [12]).scantimes

    # Must provide at least on intensity
    @test_throws ArgumentError Chrom([1]u"s", Int[])

    # Metadata cannot be anything other than a dictionary.
    @test_throws MethodError Chrom([1, 2, 3]u"s", [12, 956, 23], :sample_name)
    @test_throws MethodError Chrom([1, 2, 3]u"s", [12, 956, 23], 1)
    @test_throws MethodError Chrom([1, 2, 3]u"s", [12, 956, 23], 'S')

    # Number of scantimes and intensities must be identical
    @test_throws DimensionMismatch Chrom([1, 2, 3, 4]u"s", [12, 956, 23])
    @test_throws DimensionMismatch Chrom([1, 2]u"s", [12, 956, 23])

    # Scan times must be in ascending order
    @test_throws ArgumentError Chrom([2, 1, 3]u"s", [12, 956, 23])
    @test_throws ArgumentError Chrom([2, 2, 3]u"s", [12, 956, 23])

    # Chrom is broadcastable
    @test [1, 2]u"s" == ((chrom, i) -> chrom.scantimes[i]).(Chrom([1, 2, 3]u"s", 
        [12, 956, 23]), [1, 2])
end


@testset "show Chrom" begin
    io = IOBuffer()
    show(io, Chrom(Int64[1]u"s", Int64[12]))
    @test String(take!(io)) == string(
        "JuChrom.Chrom {scan times: Int64, intensities: Int64}\n",
        "1 scan; scan time: 1 s\n",
        "intensity: 12\n",
        "metadata: 0 entries")
    show(io, Chrom(Int64[1, 2]u"s", Int64[12, 13]))
    @test String(take!(io)) == string(
        "JuChrom.Chrom {scan times: Int64, intensities: Int64}\n",
        "2 scans; scan times: 1 s, 2 s\n",
        "intensity range: 12 - 13\n",
        "metadata: 0 entries")
    show(io, Chrom(convert(Vector{Int64}, (1:10))u"s", convert(Vector{Int64}, (12:21))))
    @test String(take!(io)) == string(
        "JuChrom.Chrom {scan times: Int64, intensities: Int64}\n",
        "10 scans; scan times: 1 s, 2 s, 3 s, 4 s, 5 s, 6 s, 7 s, 8 s, 9 s, 10 s\n",
        "intensity range: 12 - 21\n",
        "metadata: 0 entries")
    show(io, Chrom(convert(Vector{Int64}, (1:11))u"s", convert(Vector{Int64}, (12:22))))
    @test String(take!(io)) == string(
        "JuChrom.Chrom {scan times: Int64, intensities: Int64}\n",
        "11 scans; scan time range: 1 s - 11 s\n",
        "intensity range: 12 - 22\n",
        "metadata: 0 entries")
    show(io, Chrom(Int64[1, 2, 3]u"s", Int64[12, 956, 23], Dict(:id => 4)))
    @test String(take!(io)) == string(
        "JuChrom.Chrom {scan times: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "intensity range: 12 - 956\n",
        "metadata: 1 entry")
    show(io, Chrom(Int64[1, 2, 3]u"s", Int64[12, 956, 23], Dict(:id => 4, 
        "name" => "sample")))
    @test String(take!(io)) == string(
        "JuChrom.Chrom {scan times: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "intensity range: 12 - 956\n",
        "metadata: 2 entries")
    
    # RiMapper
    show(io, Chrom(Int64[1, 2, 3]u"s", Int64[12, 956, 23], rimapper=RiMapper("Kovats", 
        (1:5)u"minute", 1000:1000:5000)))
    @test String(take!(io)) == string(
        "JuChrom.Chrom {scan times: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "intensity range: 12 - 956\n",
        "metadata: 0 entries\n",
        "retention index mapper: Kovats")
end
