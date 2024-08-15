using JuChrom
using Test
using Unitful: 𝐓

@testset "GCMS constructor" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(GCMS(Float64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).scantimes)

    # Ions
    @test [85, 100] == GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]).ions
    @test [85] == GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1], (:,1))).ions
    
    @test Vector{Int64} == typeof(GCMS([1, 2, 3]u"s", Int64[85, 100], 
        [0 12; 34 956; 23 1]).ions)
    @test Vector{Float64} == typeof(GCMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1]).ions)
    
    @test UnitRange{Int64} == typeof(GCMS([1, 2, 3]u"s", 85:86, [0 12; 34 956; 23 1]).ions)
    @test UnitRange{Int64} == typeof(GCMS([1, 2, 3]u"s", 85:85, 
        reshape([0, 956, 1], (:,1))).ions)
    
    @test StepRange{Int64, Int64} == typeof(GCMS([1, 2, 3]u"s", 85:15:100, 
        [0 12; 34 956; 23 1]).ions)
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(GCMS([1, 2, 3]u"s", 85.0:15:100.0, [0 12; 34 956; 23 1]).ions)
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(GCMS([1, 2, 3]u"s", 85.0:1.0:85.0, 
        reshape([0, 956, 1], (:,1))).ions)

    # Intensities
    @test [0 12; 34 956; 23 1] == GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).intensities
    @test reshape([0, 956, 1], (:,1)) == GCMS([1, 2, 3]u"s", [85],
        reshape([0, 956, 1], (:,1))).intensities
    @test Matrix{Int64} == typeof(GCMS([1, 2, 3]u"s", [85, 100], 
        Int64[0 12; 34 956; 23 1]).intensities)
    @test Matrix{Float64} == typeof(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]).intensities)
    @test Matrix{Int64} == typeof(GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1],
         (:,1))).intensities)

    # Metadata
    @test Dict{Any, Any} == typeof(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).metadata)
    @test Dict() == GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]).metadata
    @test Dict(:id => 4, "name" => "sample") == GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 4, "name" => "sample")).metadata

    # Check the associated supertypes
    @test isa(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), AbstractChromatogram)
    @test isa(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), AbstractGCMS)

    # Scan time vector accepts only time values
    @test_throws MethodError GCMS([1, 2, 3], [85, 100], [0 12; 34 956; 23 1])
    @test_throws MethodError GCMS([1, 2, 3]u"m", [85, 100], [0 12; 34 956; 23 1])

    # Ion vector must not accept a unit
    @test_throws MethodError GCMS([1, 2, 3], [85, 100]u"s", [0 12; 34 956; 23 1])
    
    # Intensity matrix must not accept a unit
    @test_throws MethodError GCMS([1, 2, 3], [85, 100], [0 12; 34 956; 23 1]u"s")

    # Must provide at least on scantime
    @test_throws ArgumentError GCMS(Int[]u"s", [85, 100], Matrix{Int}(undef, (0,0)))
    @test_throws ArgumentError GCMS(Int[]u"s", [85, 100], [0 12])
    @test [1]u"s" == GCMS([1]u"s", [85, 100], [0 12]).scantimes

    # Must provide at least on intensity
    @test_throws ArgumentError GCMS([1]u"s", [85, 100], Matrix{Int}(undef, (0,0)))

    # Metadata cannot be anything other than dictionary
    @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], :sample_name)
    @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], 1)
    @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], 'S')

    # Check size compatibility of intensities, scantimes, and ions
    # Too many scan times in comparison to intensity matrix
    @test_throws DimensionMismatch GCMS([1, 2, 3, 4]u"s", [85, 100], [0 12; 34 956; 23 1])
    # Too few scan times in comparison to intensity matrix
    @test_throws DimensionMismatch GCMS([1, 2]u"s", [85, 100], [0 12; 34 956; 23 1])
    # Too many ions in comparison to intensity matrix
    @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85, 100, 200], 
        [0 12; 34 956; 23 1])
    # Too few ions in comparison to intensity matrix
    @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85], [0 12; 34 956; 23 1])

    # Scan times must be in ascending order
    @test_throws ArgumentError GCMS([2, 1, 3]u"s", [85, 100], [0 12; 34 956; 23 1])

    # Ions must be in ascending order
    @test_throws ArgumentError GCMS([1, 2, 3]u"s", [100, 85], [0 12; 34 956; 23 1])

    # Intensity values cannot be less than zero
    @test_throws ArgumentError GCMS([1, 2, 3]u"s", [85, 100], [-1 12; 34 956; 23 1])

    # GCMS is broadcastable
    @test [1, 2]u"s" == ((gcms, i) -> gcms.scantimes[i]).(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), [1, 2])
end


@testset "show GCMS" begin
    io = IOBuffer()
    show(io, GCMS(Int64[1]u"s", Int64[85], reshape(Int64[12], (:,1))))
    @test String(take!(io)) == string(
        "GCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n",
        "1 scan; scantime: 1 s\n",
        "1 ion: m/z 85\n",
        "intensity: 12\n",
        "metadata: 0 entries")
    show(io, GCMS(Int64[1, 2]u"s", Int64[85, 100], Int64[0 12; 34 956]))
    @test String(take!(io)) == string(
        "GCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n",
        "2 scans; scantimes: 1 s, 2 s\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 956\n",
        "metadata: 0 entries")
    show(io, GCMS(convert(Vector{Int64}, (1:10))u"s", Int64[85, 100], 
        reshape(convert(Vector{Int64}, (0:19)), (10,2))))
    @test String(take!(io)) == string(
        "GCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n",
        "10 scans; scantimes: 1 s, 2 s, 3 s, 4 s, 5 s, 6 s, 7 s, 8 s, 9 s, 10 s\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 19\n",
        "metadata: 0 entries")
    show(io, GCMS(convert(Vector{Int64}, (1:11))u"s", Int64[85, 100], 
        reshape(convert(Vector{Int64}, (0:21)), (11,2))))
    @test String(take!(io)) == string(
        "GCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n",
        "11 scans; scantime range: 1 s - 11 s\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 21\n",
        "metadata: 0 entries")
    show(io, GCMS(Int64[1, 2, 3]u"s", Int64[85, 100], Int64[0 12; 34 956; 23 1], 
        Dict(:id => 4)))
    @test String(take!(io)) == string(
        "GCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n",
        "3 scans; scantimes: 1 s, 2 s, 3 s\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 956\n",
        "metadata: 1 entry")
    show(io, GCMS(Int64[1, 2, 3]u"s", Int64[85, 100], Int64[0 12; 34 956; 23 1], 
        Dict(:id => 4, "name" => "sample")))
    @test String(take!(io)) == string(
        "GCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n",
        "3 scans; scantimes: 1 s, 2 s, 3 s\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 956\n",
        "metadata: 2 entries")
end