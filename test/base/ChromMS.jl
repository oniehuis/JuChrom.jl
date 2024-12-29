using JuChrom
using Test
using Unitful: 𝐓

@testset "ChromMS constructor" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(ChromMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(ChromMS(Float64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).scantimes)

    # Ions
    @test [85, 100] == ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]).ions
    @test [85] == ChromMS([1, 2, 3]u"s", [85], reshape([0, 956, 1], (:,1))).ions
    
    @test Vector{Int64} == typeof(ChromMS([1, 2, 3]u"s", Int64[85, 100], 
        [0 12; 34 956; 23 1]).ions)
    @test Vector{Float64} == typeof(ChromMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1]).ions)
    
    @test UnitRange{Int64} == typeof(ChromMS([1, 2, 3]u"s", 85:86, 
        [0 12; 34 956; 23 1]).ions)
    @test UnitRange{Int64} == typeof(ChromMS([1, 2, 3]u"s", 85:85, 
        reshape([0, 956, 1], (:,1))).ions)
    
    @test StepRange{Int64, Int64} == typeof(ChromMS([1, 2, 3]u"s", 85:15:100, 
        [0 12; 34 956; 23 1]).ions)
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(ChromMS([1, 2, 3]u"s", 85.0:15:100.0, [0 12; 34 956; 23 1]).ions)
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(ChromMS([1, 2, 3]u"s", 85.0:1.0:85.0, 
        reshape([0, 956, 1], (:,1))).ions)

    # Intensities
    @test [0 12; 34 956; 23 1] == ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).intensities
    @test reshape([0, 956, 1], (:,1)) == ChromMS([1, 2, 3]u"s", [85],
        reshape([0, 956, 1], (:,1))).intensities
    @test Matrix{Int64} == typeof(ChromMS([1, 2, 3]u"s", [85, 100], 
        Int64[0 12; 34 956; 23 1]).intensities)
    @test Matrix{Float64} == typeof(ChromMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]).intensities)
    @test Matrix{Int64} == typeof(ChromMS([1, 2, 3]u"s", [85], reshape([0, 956, 1],
         (:,1))).intensities)

    # Metadata
    @test Dict{Any, Any} == typeof(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).metadata)
    @test Dict() == ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]).metadata
    @test Dict(:id => 4, "name" => "sample") == ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 4, "name" => "sample")).metadata

    # RiMapper, incl. its mutability
    @test isa(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]).rimapper, Nothing)
    @test isa(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1], rimapper=RiMapper(
        "Kovats", (1:5)u"minute", 1000:1000:5000)).rimapper, JuChrom.AbstractRiMapper)

    chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    @test nothing === chrom.rimapper
    chrom.rimapper = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000)
    @test isa(chrom.rimapper, JuChrom.AbstractRiMapper)

    # Check the associated supertypes
    @test isa(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        AbstractChromatogram)
    @test isa(ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), AbstractChromMS)

    # Scan time vector accepts only time values
    @test_throws MethodError ChromMS([1, 2, 3], [85, 100], [0 12; 34 956; 23 1])
    @test_throws MethodError ChromMS([1, 2, 3]u"m", [85, 100], [0 12; 34 956; 23 1])

    # Ion vector must not accept a unit
    @test_throws MethodError ChromMS([1, 2, 3], [85, 100]u"s", [0 12; 34 956; 23 1])
    
    # Intensity matrix must not accept a unit
    @test_throws MethodError ChromMS([1, 2, 3], [85, 100], [0 12; 34 956; 23 1]u"s")

    # Must provide at least on scantime
    @test_throws ArgumentError ChromMS(Int[]u"s", [85, 100], Matrix{Int}(undef, (0,0)))
    @test_throws ArgumentError ChromMS(Int[]u"s", [85, 100], [0 12])
    @test [1]u"s" == ChromMS([1]u"s", [85, 100], [0 12]).scantimes

    # Must provide at least on intensity
    @test_throws ArgumentError ChromMS([1]u"s", [85, 100], Matrix{Int}(undef, (0,0)))

    # Metadata cannot be anything other than dictionary
    @test_throws MethodError ChromMS([1, 2, 3]u"s", [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], :sample_name)
    @test_throws MethodError ChromMS([1, 2, 3]u"s", [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], 1)
    @test_throws MethodError ChromMS([1, 2, 3]u"s", [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], 'S')

    # Check size compatibility of intensities, scantimes, and ions
    # Too many scan times in comparison to intensity matrix
    @test_throws DimensionMismatch ChromMS([1, 2, 3, 4]u"s", [85, 100], 
        [0 12; 34 956; 23 1])
    # Too few scan times in comparison to intensity matrix
    @test_throws DimensionMismatch ChromMS([1, 2]u"s", [85, 100], [0 12; 34 956; 23 1])
    # Too many ions in comparison to intensity matrix
    @test_throws DimensionMismatch ChromMS([1, 2, 3]u"s", [85, 100, 200], 
        [0 12; 34 956; 23 1])
    # Too few ions in comparison to intensity matrix
    @test_throws DimensionMismatch ChromMS([1, 2, 3]u"s", [85], [0 12; 34 956; 23 1])

    # Scan times must be in ascending order
    @test_throws ArgumentError ChromMS([2, 1, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    @test_throws ArgumentError ChromMS([2, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])

    # Ions must be in ascending order
    @test_throws ArgumentError ChromMS([1, 2, 3]u"s", [100, 85], [0 12; 34 956; 23 1])
    @test_throws ArgumentError ChromMS([1, 2, 3]u"s", [85, 85], [0 12; 34 956; 23 1])

    # ChromMS is broadcastable
    @test [1, 2]u"s" == ((chrom, i) -> chrom.scantimes[i]).(ChromMS([1, 2, 3]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), [1, 2])
end


@testset "show ChromMS" begin
    io = IOBuffer()
    show(io, ChromMS(Int64[1]u"s", Int64[85], reshape(Int64[12], (:,1))))
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Int64, ions: Int64, intensities: Int64}\n",
        "1 scan; scan time: 1 s\n",
        "1 ion: m/z 85\n",
        "intensity: 12\n",
        "metadata: 0 entries")
    show(io, ChromMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956]))
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Int64, ions: Int64, intensities: Int64}\n",
        "2 scans; scan times: 1 s, 2 s\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 956\n",
        "metadata: 0 entries")
    show(io, ChromMS(convert(Vector{Int64}, (1:10))u"s", Int64[85, 100], 
        reshape(convert(Vector{Int64}, (0:19)), (10,2))))
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Int64, ions: Int64, intensities: Int64}\n",
        "10 scans; scan times: 1 s, 2 s, 3 s, 4 s, 5 s, 6 s, 7 s, 8 s, 9 s, 10 s\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 19\n",
        "metadata: 0 entries")
    show(io, ChromMS(convert(Vector{Int64}, (1:11))u"s", Int64[85, 100], 
        reshape(convert(Vector{Int64}, (0:21)), (11,2))))
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Int64, ions: Int64, intensities: Int64}\n",
        "11 scans; scan time range: 1 s - 11 s\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 21\n",
        "metadata: 0 entries")
    show(io, ChromMS(Int64[1, 2, 3]u"s", Int64[85, 100], Int64[0 12; 34 956; 23 1], 
        Dict(:id => 4)))
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Int64, ions: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 956\n",
        "metadata: 1 entry")
    show(io, ChromMS(Int64[1, 2, 3]u"s", convert(Vector{Int64}, 85:94), 
        reshape(convert(Vector{Int64}, (0:29)), (3,10))))
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Int64, ions: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "10 ions: m/z 85, 86, 87, 88, 89, 90, 91, 92, 93, 94\n",
        "intensity range: 0 - 29\n",
        "metadata: 0 entries")
    show(io, ChromMS(Int64[1, 2, 3]u"s", convert(Vector{Int64}, 85:95), 
        reshape(convert(Vector{Int64}, (0:32)), (3,11))))
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Int64, ions: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "11 ions; range: m/z 85 - 95\n",
        "intensity range: 0 - 32\n",
        "metadata: 0 entries")

    # RiMapper
    show(io, ChromMS((1:4)u"minute", [85, 100], [0 12; 34 956; 23 1; 0 0], 
        rimapper=RiMapper("Kovats", (2:5)u"minute", 2000:1000:5000)))
    @test String(take!(io)) == string(
    "JuChrom.ChromMS {scan times: Int64, ions: Int64, intensities: Int64}\n",
    "4 scans; scan times: 1 minute, 2 minute, 3 minute, 4 minute\n",
    "2 ions: m/z 85, 100\n",
    "intensity range: 0 - 956\n",
    "metadata: 0 entries\n",
    "retention index mapper: Kovats")
end