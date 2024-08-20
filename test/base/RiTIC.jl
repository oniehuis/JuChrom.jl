using JuChrom
using Test
using Unitful: 𝐓

@testset "RiTIC constructor" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(RiTIC(Int64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(RiTIC(Float64[1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]).scantimes)

    # Intensities
    @test [12, 956, 23] == RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]).intensities
    @test Vector{Int64} == typeof(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        Int64[12, 956, 23]).intensities)
    @test Vector{Float64} == typeof(RiTIC([1, 2, 3]u"s",  "Kovats", [100, 200, 300], 
        Float64[12, 956, 23]).intensities)

    # Retention index name
    @test "Kovats" == RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]
        ).retentionindexname
    @test String == typeof(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]
        ).retentionindexname)
    @test SubString{String} == typeof(RiTIC([1, 2, 3]u"s", SubString("Kovats", 1:6), 
        [100, 200, 300], Float64[12, 956, 23]).retentionindexname)

    # # Retention indices
    @test [100, 200, 300] == RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]).retentionindices
    @test Vector{Int64} == typeof(RiTIC([1, 2, 3]u"s", "Kovats", Int64[100, 200, 300], 
        [12, 956, 23]).retentionindices)
    @test Vector{Float64} == typeof(RiTIC([1, 2, 3]u"s", "Kovats", Float64[100, 200, 300], 
        [12, 956, 23]).retentionindices)
    @test Vector{Union{Missing, Float64}} == typeof(RiTIC([1, 2, 3]u"s", "Kovats", 
        [missing, 200.0, 300.0], [12, 956, 23]).retentionindices)

    # Metadata
    @test Dict{Any, Any} == typeof(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]).metadata)
    @test Dict() == RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23]).metadata
    @test Dict(:id => 4, "name" => "sample") == RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23], Dict(:id => 4, "name" => "sample")).metadata

    # Check the associated supertypes
    @test isa(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), 
        AbstractChromatogram)
    @test isa(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), AbstractGC)
    @test isa(RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23]), AbstractFID)

    # Check the associated trait
    @test JuChrom.HasRetentionIndexData() == JuChrom.RetentionIndexStyle(typeof(RiTIC(
        [1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23])))
    @test JuChrom.HasNoRetentionIndexData() ≠ JuChrom.RetentionIndexStyle(typeof(RiTIC(
        [1, 2, 3]u"s", "Kovats", [100, 200, 300], [12, 956, 23])))
    
    # Scan time vector accepts only time values
    @test_throws MethodError RiTIC([1, 2, 3], "Kovats", [100, 200, 300], [12, 956, 23])
    @test_throws MethodError RiTIC([1, 2, 3]u"m", "Kovats", [100, 200, 300], 
        [12, 956, 23])
    
    # Intensity vector must not accept a unit
    @test_throws MethodError RiTIC([1, 2, 3], "Kovats", [100, 200, 300], [12, 956, 23]u"s")

    # Retention index vector must not accept a unit
    @test_throws MethodError RiTIC([1, 2, 3], "Kovats", [100, 200, 300]u"s", 
        [12, 956, 23]u"s")

    # Metadata cannot be anything other than a dictionary
    @test_throws MethodError RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23], :sample_name)
    @test_throws MethodError RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23], 1)
    @test_throws MethodError RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [12, 956, 23], 'S')

    # Must provide at least on scantime
    @test_throws ArgumentError RiTIC(Int[]u"s", "Kovats", Int[], Int[])
    @test_throws ArgumentError RiTIC(Int[]u"s", "Kovats", Int[100], [12])
    @test [1]u"s" == RiTIC([1]u"s", "Kovats", Int[100], [12]).scantimes

    # Must provide a reteion index name
    @test_throws ArgumentError RiTIC([1, 2, 3]u"s", "", [100, 200, 300], [12, 956, 23])

    # Must provide at least on retention index
    @test_throws ArgumentError RiTIC([1]u"s", "Kovats", Int[], Int[12])
    @test_throws ArgumentError RiTIC([1]u"s", "Kovats", [missing], Int[12])

    # Must provide at least on intensity
    @test_throws ArgumentError RiTIC([1]u"s", "Kovats", Int[100], Int[])

    # Number of scan times and retention indices must be identical
    @test_throws DimensionMismatch RiTIC([1, 2, 3]u"s", "Kovats", [100, 200], 
        [12, 956, 23])
    @test_throws DimensionMismatch RiTIC([1, 2]u"s", "Kovats", [100, 200, 300], [12, 956])

    # Number of scan times and intensities must be identical
    @test_throws DimensionMismatch RiTIC([1, 2, 3, 4]u"s", "Kovats", [100, 200, 300, 400], 
        [12, 956, 23])
    @test_throws DimensionMismatch RiTIC([1, 2]u"s", "Kovats", [100, 200], [12, 956, 23])

    # Scan times must be in ascending order
    @test_throws ArgumentError RiTIC([2, 1, 3]u"s", "Kovats", [100, 300, 200], 
        [12, 956, 23])

    # Retention indices must be in ascending order
    @test_throws ArgumentError RiTIC([1, 2, 3]u"s", "Kovats", [100, 300, 200], 
        [12, 956, 23])
    @test_throws ArgumentError RiTIC([1, 2, 3]u"s", "Kovats", [missing, 300, 200], 
        [12, 956, 23])

    # Intensity values cannot be less than zero
    @test_throws ArgumentError RiTIC([1, 2, 3]u"s", "Kovats", [100, 200, 300], 
        [-1, 956, 23])

    # RiTIC is broadcastable
    @test [1, 2]u"s" == ((fid, i) -> fid.scantimes[i]).(RiTIC([1, 2, 3]u"s", "Kovats", 
        [100, 200, 300], [12, 956, 23]), [1, 2])
end


# @testset "show RiTIC" begin
#     io = IOBuffer()
#     show(io, RiTIC(Int64[1]u"s", Int64[12]))
#     @test String(take!(io)) == string(
#         "RiTIC {scan times: Int64, intensities: Int64}\n",
#         "1 scan; scan time: 1 s\n",
#         "intensity: 12\n",
#         "metadata: 0 entries")
#     show(io, RiTIC(Int64[1, 2]u"s", Int64[12, 13]))
#     @test String(take!(io)) == string(
#         "RiTIC {scan times: Int64, intensities: Int64}\n",
#         "2 scans; scan times: 1 s, 2 s\n",
#         "intensity range: 12 - 13\n",
#         "metadata: 0 entries")
#     show(io, RiTIC(convert(Vector{Int64}, (1:10))u"s", convert(Vector{Int64}, (12:21))))
#     @test String(take!(io)) == string(
#         "RiTIC {scan times: Int64, intensities: Int64}\n",
#         "10 scans; scan times: 1 s, 2 s, 3 s, 4 s, 5 s, 6 s, 7 s, 8 s, 9 s, 10 s\n",
#         "intensity range: 12 - 21\n",
#         "metadata: 0 entries")
#     show(io, RiTIC(convert(Vector{Int64}, (1:11))u"s", convert(Vector{Int64}, (12:22))))
#     @test String(take!(io)) == string(
#         "RiTIC {scan times: Int64, intensities: Int64}\n",
#         "11 scans; scan time range: 1 s - 11 s\n",
#         "intensity range: 12 - 22\n",
#         "metadata: 0 entries")
#     show(io, RiTIC(Int64[1, 2, 3]u"s", Int64[12, 956, 23], Dict(:id => 4)))
#     @test String(take!(io)) == string(
#         "RiTIC {scan times: Int64, intensities: Int64}\n",
#         "3 scans; scan times: 1 s, 2 s, 3 s\n",
#         "intensity range: 12 - 956\n",
#         "metadata: 1 entry")
#     show(io, RiTIC(Int64[1, 2, 3]u"s", Int64[12, 956, 23], Dict(:id => 4, 
#         "name" => "sample")))
#     @test String(take!(io)) == string(
#         "RiTIC {scan times: Int64, intensities: Int64}\n",
#         "3 scans; scan times: 1 s, 2 s, 3 s\n",
#         "intensity range: 12 - 956\n",
#         "metadata: 2 entries")
# end
