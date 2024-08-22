using JuChrom
using Test
using Unitful: 𝐓

@testset "RiGCMS constructor" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(RiGCMS(Int[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(RiGCMS(Float64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]).scantimes)

    # Retention index name
    @test "Kovats" == RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]).retentionindexname
    @test String == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]).retentionindexname)
    @test SubString{String} == typeof(RiGCMS([1, 2, 3]u"s", SubString("Kovats", 1:6), 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]).retentionindexname)

    # Retention indices
    @test [100, 200, 300] == RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]).retentionindices
    @test Vector{Int64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", Int64[100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]).retentionindices)
    @test Vector{Float64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", Float64[100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]).retentionindices)
    @test Vector{Union{Missing, Float64}} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [missing, 200.0, 300.0], [85, 100], [0 12; 34 956; 23 1]).retentionindices)

    # Ions
    @test [85, 100] == RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]).ions
    @test [85] == RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85], 
        reshape([0, 956, 1], (:,1))).ions
    
    @test Vector{Int64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int64[85, 100], [0 12; 34 956; 23 1]).ions)
    @test Vector{Float64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Float64[85, 100], [0 12; 34 956; 23 1]).ions)
    
    @test UnitRange{Int64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        85:86, [0 12; 34 956; 23 1]).ions)
    @test UnitRange{Int64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        85:85, reshape([0, 956, 1], (:,1))).ions)
    
    @test StepRange{Int64, Int64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], 85:15:100, [0 12; 34 956; 23 1]).ions)
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 85.0:15:100.0, 
        [0 12; 34 956; 23 1]).ions)
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 85.0:1.0:85.0, 
        reshape([0, 956, 1], (:,1))).ions)

    # Intensities
    @test [0 12; 34 956; 23 1] == RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]).intensities
    @test reshape([0, 956, 1], (:,1)) == RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300],
        [85], reshape([0, 956, 1], (:,1))).intensities
    @test Matrix{Int64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300],
        [85, 100], Int64[0 12; 34 956; 23 1]).intensities)
    @test Matrix{Float64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300],
        [85, 100], Float64[0 12; 34 956; 23 1]).intensities)
    @test Matrix{Int64} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85], 
        reshape([0, 956, 1], (:,1))).intensities)

    # Metadata
    @test Dict{Any, Any} == typeof(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1]).metadata)
    @test Dict() == RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]).metadata
    @test Dict(:id => 4, "name" => "sample") == RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1], Dict(:id => 4, 
        "name" => "sample")).metadata

    # Check the associated supertypes
    @test isa(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), AbstractChromatogram)
    @test isa(RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]), AbstractGCMS)

    # Check the associated trait
    @test JuChrom.HasRetentionIndexData() == JuChrom.RetentionIndexStyle(typeof(RiGCMS(
        [1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))
    @test JuChrom.HasNoRetentionIndexData() ≠ JuChrom.RetentionIndexStyle(typeof(RiGCMS(
        [1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], [0 12; 34 956; 23 1])))

    # Scan time vector accepts only time values
    @test_throws MethodError RiGCMS([1, 2, 3], "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])
    @test_throws MethodError RiGCMS([1, 2, 3]u"m", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])

    # Retention index vector must not accept a unit
    @test_throws MethodError RiFID([1, 2, 3], "Kovats", [100, 200, 300]u"s", 
        [85, 100], [0 12; 34 956; 23 1])

    # Ion vector must not accept a unit
    @test_throws MethodError RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100]u"s", [0 12; 34 956; 23 1])
    
    # Intensity matrix must not accept a unit
    @test_throws MethodError RiGCMS([1, 2, 3]u"m", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1]u"s")

    # # Must provide at least on scantime
    @test_throws ArgumentError RiGCMS(Int[]u"s", "Kovats", Int[], [85, 100], 
        Matrix{Int}(undef, (0,0)))
    @test_throws ArgumentError RiGCMS(Int[]u"s", "Kovats", [100], [85, 100], 
        [0 12])
    @test [1]u"s" == RiGCMS([1]u"s", "Kovats", [100], [85, 100], [0 12]).scantimes

    # Must provide a retention index name
    @test_throws ArgumentError RiGCMS([1, 2, 3]u"s", "", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])

    # Must provide at least on retention index
    @test_throws ArgumentError RiGCMS([1]u"s", "Kovats", Int[], Int[85], 
        reshape(Int64[12], (:,1)))
    @test_throws ArgumentError RiGCMS([1]u"s", "Kovats", [missing], Int[85], 
        reshape(Int64[12], (:,1)))

    # Must provide at least on intensity
    @test_throws ArgumentError RiGCMS([1]u"s", "Kovats", [100], [85, 100], 
        Matrix{Int}(undef, (0,0)))

    # # Metadata cannot be anything other than dictionary
    @test_throws MethodError RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], :sample_name)
    @test_throws MethodError RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], 1)
    @test_throws MethodError RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], 'S')

    # Check size compatibility of intensities, scantimes, and ions
    # Too many scan times in comparison to retention indices
    @test_throws DimensionMismatch RiGCMS([1, 2, 3, 4]u"s", "Kovats", [100, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1])
    # Too few scan times in comparison to retention indices
    @test_throws DimensionMismatch RiGCMS([1, 2]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])
    # Too many ions in comparison to intensity matrix
    @test_throws DimensionMismatch RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [85, 100, 200], [0 12; 34 956; 23 1])
    # Too few ions in comparison to intensity matrix
    @test_throws DimensionMismatch RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85], 
        [0 12; 34 956; 23 1])

    # Scan times must be in ascending order
    @test_throws ArgumentError RiGCMS([2, 1, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])
    @test_throws ArgumentError RiGCMS([2, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [0 12; 34 956; 23 1])

    # Retention Indices must be in ascending order
    @test_throws ArgumentError RiGCMS([1, 2, 3]u"s", "Kovats", [200, 100, 300], [85, 100], 
        [0 12; 34 956; 23 1])
    @test_throws ArgumentError RiGCMS([1, 2, 3]u"s", "Kovats", [300, missing, 100], 
        [85, 100], [0 12; 34 956; 23 1])
    @test_throws ArgumentError RiGCMS([1, 2, 3]u"s", "Kovats", [200, 200, 300], 
        [85, 100], [0 12; 34 956; 23 1])

    # Ions must be in ascending order
    @test_throws ArgumentError RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [100, 85], 
        [0 12; 34 956; 23 1])
    @test_throws ArgumentError RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 85], 
        [0 12; 34 956; 23 1])

    # Numerical retention indices cannot be interrupted by missing values
    @test_throws ArgumentError RiGCMS([1, 2, 3]u"s", "Kovats", [100, missing, 200], 
        [85, 100], [-1 12; 34 956; 23 1])

    # Intensity values cannot be less than zero
    @test_throws ArgumentError RiGCMS([1, 2, 3]u"s", "Kovats", [100, 200, 300], [85, 100], 
        [-1 12; 34 956; 23 1])

    # GCMS is broadcastable
    @test [1, 2]u"s" == ((gcms, i) -> gcms.scantimes[i]).(RiGCMS([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [85, 100], [0 12; 34 956; 23 1]), [1, 2])
end


@testset "show RiGCMS" begin
    io = IOBuffer()
    show(io, RiGCMS(Int64[1]u"s", "Kovats", Float64[1000], Int64[85], 
        reshape(Int64[12], (:,1))))
    @test String(take!(io)) == string(
        "JuChrom.RiGCMS {scan times: Int64, retention indices: Float64, ions: Int64, ", 
            "intensities: Int64}\n",
        "1 scan; scan time: 1 s\n",
        "1 retention index: 1000.0 (≘ 1 s)\n",
        "1 ion: m/z 85\n",
        "intensity: 12\n",
        "metadata: 0 entries")
    show(io, RiGCMS(Int64[1, 2]u"s", "Kovats", Float64[1000, 2000], Int64[85, 100], 
        Int64[0 12; 34 956]))
    @test String(take!(io)) == string(
        "JuChrom.RiGCMS {scan times: Int64, retention indices: Float64, ions: Int64, ",
            "intensities: Int64}\n",
        "2 scans; scan times: 1 s, 2 s\n",
        "2 retention indices: 1000.0 (≘ 1 s), 2000.0 (≘ 2 s)\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 956\n",
        "metadata: 0 entries")
    show(io, RiGCMS(convert(Vector{Int64}, (1:10))u"s", "Kovats", 1000:100.0:1900, 
        Int64[85, 100], reshape(convert(Vector{Int64}, (0:19)), (10,2))))
    @test String(take!(io)) == string(
        "JuChrom.RiGCMS {scan times: Int64, retention indices: Float64, ions: Int64, ",
            "intensities: Int64}\n",
        "10 scans; scan times: 1 s, 2 s, 3 s, 4 s, 5 s, 6 s, 7 s, 8 s, 9 s, 10 s\n",
        "10 retention indices: 1000.0 (≘ 1 s), 1100.0 (≘ 2 s), 1200.0 (≘ 3 s), ",
            "1300.0 (≘ 4 s), 1400.0 (≘ 5 s), 1500.0 (≘ 6 s), 1600.0 (≘ 7 s), ",
            "1700.0 (≘ 8 s), 1800.0 (≘ 9 s), 1900.0 (≘ 10 s)\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 19\n",
        "metadata: 0 entries")
    show(io, RiGCMS(convert(Vector{Int64}, (1:11))u"s", "Kovats", 1000:100.0:2000, 
        Int64[85, 100], reshape(convert(Vector{Int64}, (0:21)), (11,2))))
    @test String(take!(io)) == string(
        "JuChrom.RiGCMS {scan times: Int64, retention indices: Float64, ions: Int64, ",
            "intensities: Int64}\n",
        "11 scans; scan time range: 1 s - 11 s\n",
        "11 retention indices; retention index range: 1000.0 (≘ 1 s) – 2000.0 (≘ 11 s)\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 21\n",
        "metadata: 0 entries")
    show(io, RiGCMS(Int64[1, 2, 3]u"s", "Kovats", [missing, 2000, 3000], Int64[85, 100], 
        Int64[0 12; 34 956; 23 1], Dict(:id => 4)))
    @test String(take!(io)) == string(
        "JuChrom.RiGCMS {scan times: Int64, retention indices: Union{Missing, Int64}, ", 
            "ions: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "2 retention indices: 2000 (≘ 2 s), 3000 (≘ 3 s)\n",
        "2 ions: m/z 85, 100\n",
        "intensity range: 0 - 956\n",
        "metadata: 1 entry")
    show(io, RiGCMS(Int64[1, 2, 3]u"s", "Kovats", [missing, 2000, missing], 
        convert(Vector{Int64}, 85:94), reshape(convert(Vector{Int64}, (0:29)), (3,10))))
    @test String(take!(io)) == string(
        "JuChrom.RiGCMS {scan times: Int64, retention indices: Union{Missing, Int64}, ",
            "ions: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "1 retention index: 2000 (≘ 2 s)\n",
        "10 ions: m/z 85, 86, 87, 88, 89, 90, 91, 92, 93, 94\n",
        "intensity range: 0 - 29\n",
        "metadata: 0 entries")
    show(io, RiGCMS(Int64[1, 2, 3]u"s", "Kovats", [missing, 2000, missing], 
        convert(Vector{Int64}, 85:95), reshape(convert(Vector{Int64}, (0:32)), (3,11))))
    @test String(take!(io)) == string(
        "JuChrom.RiGCMS {scan times: Int64, retention indices: Union{Missing, Int64}, ",
        "ions: Int64, intensities: Int64}\n",
        "3 scans; scan times: 1 s, 2 s, 3 s\n",
        "1 retention index: 2000 (≘ 2 s)\n",
        "11 ions; range: m/z 85 - 95\n",
        "intensity range: 0 - 32\n",
        "metadata: 0 entries")
end