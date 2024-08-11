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
    @test Dict(:id => 1, :name => "sample") == GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 1, :name => "sample")).metadata

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
    show(io, GCMS(Int64[1, 2, 3]u"s", Int64[85, 100], Int64[0 12; 34 956; 23 1], Dict(:id => 1, 
        :name => "sample")))
    @test String(take!(io)) == string(
        "GCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n", 
        "3 scans; time range: 1 s - 3 s\n",
        "2 ions; range: m/z 85 - 100\n",
        "intensity range: 0 - 956\n", 
        "metadata: :id, :name")
end
