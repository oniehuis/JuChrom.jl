using JuChrom
using Test
using Unitful: 𝐓

############################################################################################
# MassSpectrum(MassSpectrum(ions::AbstractVector{<:Real}, 
# intensities::AbstractVector{<:Real};
# retentiontime::Union{Unitful.Time, Nothing}=nothing,
# retentionindexname::Union{AbstractString, Nothing}=nothing,
# retentionindex::Union{<:Real, Nothing}=nothing, metadata::Dict=Dict())
############################################################################################
@testset "MassSpectrum" begin

    # Verify object construction and field types depending on constructor arguments
    # Ions
    @test [85.1, 112.2, 124.1] == MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]).ions
    @test Vector{Float64} == typeof(MassSpectrum(Float64[85.1, 112.2, 124.1], 
        [13, 0, 67]).ions)
    @test Vector{Int} == typeof(MassSpectrum(Int[85, 112, 124], 
        [13, 0, 67]).ions)
            
    # Intensities     
    @test [13, 0, 67] == MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]).intensities
    @test Vector{Int} == typeof(MassSpectrum([85.1, 112.2, 124.1], Int[13, 0, 67]
        ).intensities)
    @test Vector{Float64} == typeof(MassSpectrum([85, 112, 124], Float64[13, 0, 67]
        ).intensities)

    # Retention time
    @test 3.2u"minute" == MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentiontime=3.2u"minute").retentiontime
    @test Unitful.Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentiontime=3.2u"minute").retentiontime)
    @test Unitful.Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentiontime=Int64(3)*u"minute").retentiontime)

    # Retention index name
    @test "Kovats" == MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats", retentionindex=131).retentionindexname
    @test String == typeof(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats", retentionindex=131).retentionindexname)
    @test SubString{String} == typeof(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname=SubString("Kovats", 1, 6), retentionindex=131
        ).retentionindexname)       

    # Retention index name
    @test 131 == MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats", retentionindex=131).retentionindex
    @test Int == typeof(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats", retentionindex=Int(131)).retentionindex)
    @test Float64 == typeof(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats", retentionindex=131.0).retentionindex)

    # Metadata
    @test Dict(:test => 1) == MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        metadata=Dict(:test => 1)).metadata
    @test Dict{Any, Any} == typeof(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        metadata=Dict(:test => 1)).metadata)
    
    # Check the associated supertypes
    @test isa(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), AbstractMassSpectrum)

    # Must provide at least one ion
    @test_throws ArgumentError MassSpectrum(Int[], Int[])
    @test [1] == MassSpectrum([1], [1]).ions

    # Number of scantimes and intensities must be identical
    @test_throws DimensionMismatch MassSpectrum([1], [1, 2])
    @test_throws DimensionMismatch MassSpectrum([1, 2], [1])

    # Ions must be in ascending order and cannot have identical values
    @test_throws ArgumentError MassSpectrum([2, 1], [1, 2])
    @test_throws ArgumentError MassSpectrum([1, 1], [1, 2])

    # Intensity values cannot be less than zero
    @test_throws ArgumentError MassSpectrum([1, 2, 3], [-12, 956, 23])

    # If the retention index is given, the retention index name must also be given
    @test_throws ArgumentError MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindex=131)

    # If the retention index name is given, the retention index  must also be given
    @test_throws ArgumentError MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats")

    # If the retention index name must be some type of AstractString
    @test_throws MethodError MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname=:Kovats)

    # Metadata cannot be anything other than a dictionary.
    @test_throws TypeError MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        metadata=:sample_name)
    @test_throws TypeError MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        metadata=1)
    @test_throws TypeError MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        metadata='S')

    # MassSpectrum is broadcastable
    @test [85.1, 112.2] == ((ms, i) -> ms.ions[i]).(MassSpectrum([85.1, 112.2, 124.1], 
        [13, 0, 67]), [1, 2])
end

@testset "show AbstractMassSpectrum" begin
    io = IOBuffer()
    show(io, MassSpectrum([85.1, 112.2, 124.1], Int64[13, 0, 67]))
    @test String(take!(io)) == string(
        "JuChrom.MassSpectrum {ions: Float64, intensities: Int64}\n",
        "3 ions: m/z 85.1, 112.2, 124.1\n",
        "intensity range: 0 - 67\n",
        "metadata: 0 entries")
    show(io, MassSpectrum([85.1], Int64[13]))
    @test String(take!(io)) == string(
        "JuChrom.MassSpectrum {ions: Float64, intensities: Int64}\n",
        "1 ion: m/z 85.1\n",
        "intensity: 13\n",
        "metadata: 0 entries")
    show(io, MassSpectrum(convert(Vector{Int64}, collect(1:11)), convert(Vector{Int64}, 
        collect(1:11))))
    @test String(take!(io)) == string(
        "JuChrom.MassSpectrum {ions: Int64, intensities: Int64}\n",
        "11 ions; range: m/z 1 - 11\n",
        "intensity range: 1 - 11\n",
        "metadata: 0 entries")
    show(io, MassSpectrum([85.1], Int64[13], retentiontime=1.0u"s"))
    @test String(take!(io)) == string(
        "JuChrom.MassSpectrum {ions: Float64, intensities: Int64}\n",
        "1 ion: m/z 85.1\n",
        "intensity: 13\n",
        "retention time: 1.0 s\n",
        "metadata: 0 entries")
    show(io, MassSpectrum([85.1], Int64[13], retentiontime=1.0u"s", 
        retentionindexname="Kovats", retentionindex=123.2))
    @test String(take!(io)) == string(
        "JuChrom.MassSpectrum {ions: Float64, intensities: Int64}\n",
        "1 ion: m/z 85.1\n",
        "intensity: 13\n",
        "retention time: 1.0 s\n",
        "retention index: 123.2 (Kovats)\n",
        "metadata: 0 entries")
    show(io, MassSpectrum([85.1], Int64[13], retentiontime=1.0u"s", 
        retentionindexname="Kovats", retentionindex=123.2, metadata=Dict(:test => 1)))
    @test String(take!(io)) == string(
        "JuChrom.MassSpectrum {ions: Float64, intensities: Int64}\n",
        "1 ion: m/z 85.1\n",
        "intensity: 13\n",
        "retention time: 1.0 s\n",
        "retention index: 123.2 (Kovats)\n",
        "metadata: 1 entry")
    show(io, MassSpectrum([85.1], Int64[13], retentiontime=1.0u"s", 
        retentionindexname="Kovats", retentionindex=123.2, metadata=Dict(:test => 1, 
        2 => 3)))
    @test String(take!(io)) == string(
        "JuChrom.MassSpectrum {ions: Float64, intensities: Int64}\n",
        "1 ion: m/z 85.1\n",
        "intensity: 13\n",
        "retention time: 1.0 s\n",
        "retention index: 123.2 (Kovats)\n",
        "metadata: 2 entries")
end