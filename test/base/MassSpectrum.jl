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
@testset "MassSpectrum()" begin

    # Verify object construction and field types depending on constructor arguments
    # Ions
    @test [85.1, 112.2, 124.1] ≈ MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]).ions
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


############################################################################################
# cosine(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
############################################################################################
@testset "cosine" begin
    @test 0.8978872704229618 ≈ cosine([100, 500, 250], [200, 1000, 0])
    @test 1 ≈ cosine([1, 1, 1], [1, 1, 1])
    @test 0 ≈ cosine([10.5, 10.5, 10.5], [-10.5, -10.5, -10.5])
    @test Float64 == typeof(cosine([100, 500, 250], [200, 1000, 0]))
    @test Float64 == typeof(cosine([1, 1, 1], [1, 1, 1]))
    @test Float64 == typeof(cosine([10.5, 10.5, 10.5], [-10.5, -10.5, -10.5]))
    @test_throws DimensionMismatch cosine([100, 500, 250], [200, 1000, 0, 10])
    @test_throws ArgumentError cosine(Int[], Int[])
    @test_throws ArgumentError cosine([0, 0, 0], [200, 1000, 0])
    @test_throws ArgumentError cosine([100, 500, 250], [0, 0, 0])
 end


############################################################################################
# intensities(ms::AbstractMassSpectrum)
############################################################################################
@testset "intensities(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test [13, 0, 67] == intensities(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))

    # Same return same element type as used to construct the object
    @test Vector{Int} == typeof(intensities(MassSpectrum([85.1, 112.2, 124.1], 
        Int[13, 0, 67])))
    @test Vector{Float64} == typeof(intensities(MassSpectrum([85.1, 112.2, 124.1], 
        Float64[13, 0, 67])))
end


############################################################################################
# intensity(ms::AbstractMassSpectrum, ionindex::Integer)
############################################################################################
@testset "intensity(ms, ionindex)" begin
    # Same return values as those provided as arguments to construct the object
    @test 13 == intensity(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 1)
    @test 0 == intensity(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 2)

    # Same element type as used to construct the object
    @test Int == typeof(intensity(MassSpectrum([85.1, 112.2, 124.1], Int[13, 0, 67]), 1))
    @test Float64 == typeof(intensity(MassSpectrum([85.1, 112.2, 124.1], 
        Float64[13, 0, 67]), 1))
    
    # Provoke a BoundsError by specifying an index that does not exist
    @test_throws BoundsError intensity(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 0)
    @test_throws BoundsError intensity(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 4)
end


############################################################################################
# intensitydifferences(::AbstractMassSpectrum, ::AbstractMassSpectrum)
############################################################################################
@testset "intensitydifferences(::AbstractMassSpectrum, ::AbstractMassSpectrum)" begin
    # Same return values as those provided as arguments to construct the object
    @test [10, 70, 5] == intensitydifferences(MassSpectrum([80, 85, 90], [10, 80, 5]), 
        MassSpectrum([80, 85, 90], [0, 10, 0]))
    @test [0, 0, 0] == intensitydifferences(MassSpectrum([80, 85, 90], [0, 10, 0]), 
        MassSpectrum([80, 85, 90], [10, 80, 5]))
    @test [4.5, 59.9, 0.0] == intensitydifferences(MassSpectrum([80, 85, 90], [10, 80, 5]), 
        MassSpectrum([80, 85, 90], [5.5, 20.1, 5]))

    # Same element type as used to construct the object
    @test Vector{Int} == typeof(intensitydifferences(MassSpectrum([80, 85, 90], 
        [10, 80, 5]), MassSpectrum(Int[80, 85, 90], Int[0, 10, 0])))
    @test Vector{Int} == typeof(intensitydifferences(MassSpectrum([80, 85, 90], 
        [0, 10, 0]), MassSpectrum(Int[80, 85, 90], Int[10, 80, 5])))
    @test Vector{Float64}== typeof(intensitydifferences(MassSpectrum([80, 85, 90], 
        [10, 80, 5]), MassSpectrum(Int[80, 85, 90], Float64[5.5, 20.1, 5])))

    # Provoke a BoundsError by specifying an index that does not exist
    @test_throws ArgumentError intensitydifferences(MassSpectrum([80], [5]), 
        MassSpectrum([81], [5]))
end


############################################################################################
# intensitysums(::AbstractMassSpectrum, ::AbstractMassSpectrum)
############################################################################################
@testset "intensitysums(::AbstractMassSpectrum, ::AbstractMassSpectrum)" begin
    # Same return values as those provided as arguments to construct the object
    @test [10, 90, 5] == intensitysums(MassSpectrum([80, 85, 90], [10, 80, 5]), 
        MassSpectrum([80, 85, 90], [0, 10, 0]))
    @test [10, 90, 5] == intensitysums(MassSpectrum([80, 85, 90], [0, 10, 0]), 
        MassSpectrum([80, 85, 90], [10, 80, 5]))
    @test [15.5, 100.1, 10.0] == intensitysums(MassSpectrum([80, 85, 90], [10, 80, 5]), 
        MassSpectrum([80, 85, 90], [5.5, 20.1, 5]))

    # Same element type as used to construct the object
    @test Vector{Int} == typeof(intensitysums(MassSpectrum([80, 85, 90], 
        [10, 80, 5]), MassSpectrum(Int[80, 85, 90], Int[0, 10, 0])))
    @test Vector{Int} == typeof(intensitysums(MassSpectrum([80, 85, 90], 
        [0, 10, 0]), MassSpectrum(Int[80, 85, 90], Int[10, 80, 5])))
    @test Vector{Float64}== typeof(intensitysums(MassSpectrum([80, 85, 90], 
        [10, 80, 5]), MassSpectrum(Int[80, 85, 90], Float64[5.5, 20.1, 5])))

    # Provoke a BoundsError by specifying an index that does not exist
    @test_throws ArgumentError intensitysums(MassSpectrum([80], [5]), 
        MassSpectrum([81], [5]))
end

############################################################################################
# ion(gcms::AbstractMassSpectrum, ionindex::Integer)
############################################################################################
@testset "ion(gcms, ionindex)" begin
    # Same return values as those provided as arguments to construct the object
    @test 85.1 ≈ ion(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 1)
    @test 112.2 ≈ ion(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 2)
    
    # Same return element type as used to construct the object
    @test Float64 == typeof(ion(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 1))
    @test Int == typeof(ion(MassSpectrum(Int[85, 112, 124], [13, 0, 67]), 1))

    # Provoke a BoundsError by specifying an index that does not exist
    @test_throws BoundsError ion(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 0)
    @test_throws BoundsError ion(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 4)
end


############################################################################################
# ioncount(ms::AbstractMassSpectrum)
############################################################################################
@testset "ioncount(ms)" begin
    # Validate the returned value
    @test 3 == ioncount(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))
    @test 1 == ioncount(MassSpectrum([85.1], [13]))

    # Return value must be an integer
    @test Int == typeof(ioncount(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67])))
end


############################################################################################
# ionindex(ms::AbstractMassSpectrum, ion::Real)
############################################################################################
@testset "ionindex(ms, ion)" begin
    # Validate the returned value
    @test 1 == ionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 85.1)
    @test 2 == ionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 112.2)
    @test 1 == ionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), Float32(85.1))
    @test 2 == ionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), Float32(112.2))

    # Ion must exist
    @test_throws ArgumentError ionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 20)

    # Return value must be an integer
    @test Int == typeof(ionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]), 85.1))
end


############################################################################################
# ions(ms::AbstractMassSpectrum)
############################################################################################
@testset "ions(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test [85.1, 112.2, 124.1] == ions(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))
    
    # Same return container and element type as used to construct the object
    @test Vector{Float64} == typeof(ions(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67])))
    @test Vector{Int} == typeof(ions(MassSpectrum(Int[85, 112, 124], [13, 0, 67])))
end


############################################################################################
# massspectrum(gcms::AbstractChromMS, scanindex::Integer)
############################################################################################
@testset "massspectrum(gcms, scanindex)" begin
    @test [85, 100] == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2).ions
    @test [34, 956] == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2).intensities
    @test 2u"s" == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2).retentiontime
    @test nothing === massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2).retentionindexname
    @test nothing === massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2).retentionindex
    @test "Kovats" == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, retentionindexname="Kovats", retentionindex=123.21
        ).retentionindexname
    @test 123.21 == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, retentionindexname="Kovats", retentionindex=123.21
        ).retentionindex
    @test Dict() == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2).metadata
    @test Dict(:test => 1) == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2, metadata=Dict(:test => 1)).metadata

    @test_throws BoundsError massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0)
    @test_throws BoundsError massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 4)
end


############################################################################################
# massspectrum(chrom::AbstractChromMS, time::Unitful.Time; precisetime::Bool=false)
############################################################################################
@testset "massspectrum(chrom, time; precisetime)" begin
    @test [85, 100] == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s").ions
    @test [34, 956] == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s").intensities
    @test [23, 1] == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.5u"s").intensities
    @test [0, 12] == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 0.4u"s").intensities
    @test [23, 1] == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 4u"s").intensities
    @test 2u"s" == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
         [0 12; 34 956; 23 1]), 2u"s").retentiontime
    @test 3u"s" == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
         [0 12; 34 956; 23 1]), 2.5u"s").retentiontime
    @test nothing === massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s").retentionindexname
    @test nothing === massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s").retentionindex
    @test "Kovats" == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s", retentionindexname="Kovats", retentionindex=123.21
        ).retentionindexname
    @test 123.21 == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s", retentionindexname="Kovats", retentionindex=123.21
        ).retentionindex
    @test Dict() == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s").metadata
    @test Dict(:test => 1) == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2u"s", metadata=Dict(:test => 1)).metadata

    @test [34, 956] == massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.0u"s", precisetime=true).intensities
    @test_throws ArgumentError massspectrum(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), 2.1u"s", precisetime=true)
end

############################################################################################
# maxintensity(ms::AbstractMassSpectrum)
############################################################################################
@testset "maxintensity(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test 67 == maxintensity(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))

    # Same element type as used to construct the object
    @test Int == typeof(maxintensity(MassSpectrum([85.1, 112.2, 124.1], Int[13, 0, 67])))
    @test Float64 == typeof(maxintensity(MassSpectrum([85.1, 112.2, 124.1], 
        Float64[13, 0, 67])))
end


############################################################################################
# maxion(ms::AbstractMassSpectrum)
############################################################################################
@testset "maxion(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test 124.1 == maxion(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))

    # Same element type as used to construct the object
    @test Float64 == typeof(maxion(MassSpectrum(Float64[85.1, 112.2, 124.1], [13, 0, 67])))
    @test Int == typeof(maxion(MassSpectrum([85, 112, 124], [13, 0, 67])))
end


############################################################################################
# meanintensities(mss::AbstractVector{<:AbstractMassSpectrum})
############################################################################################
@testset "meanintensities(::AbstractVector{<:AbstractMassSpectrum})" begin
    # Same return values as those provided as arguments to construct the object
    @test [5.0, 45.0, 2.5] ≈ meanintensities([MassSpectrum([80, 85, 90], [0, 10, 0]), 
        MassSpectrum([80, 85, 90], [10, 80, 5])])
    @test [5.0, 45.0, 2.5] ≈ meanintensities([MassSpectrum([80, 85, 90], [0, 10, 0]), 
        MassSpectrum([80, 85, 90], Float64[10, 80, 5])])
    @test [5.0] ≈ meanintensities([MassSpectrum([80], [0]), MassSpectrum([80], [10])])
    @test [5.0] ≈ meanintensities([MassSpectrum([80], [0]), MassSpectrum([80], 
        Float64[10])])

    # Same element type as used to construct the object
    @test Vector{Float64} == typeof(meanintensities([MassSpectrum([80, 85, 90], 
        [0, 10, 0]), MassSpectrum([80, 85, 90], [10, 80, 5])]))
    @test Vector{Float64} == typeof(meanintensities([MassSpectrum([80, 85, 90], 
        [0, 10, 0]), MassSpectrum([80, 85, 90], Float64[10, 80, 5])]))
    @test Vector{Float64} == typeof(meanintensities([MassSpectrum([80], [0]), 
        MassSpectrum([80], [10])]))
    @test Vector{Float64} == typeof(meanintensities([MassSpectrum([80], [0]), 
        MassSpectrum([80], Float64[10])]))

    @test_throws ArgumentError meanintensities([MassSpectrum([80, 85, 90], [0, 10, 0]), 
        MassSpectrum([80, 84, 90], [10, 80, 5])])
    @test_throws ArgumentError meanintensities([MassSpectrum([80], [0]), 
        MassSpectrum([81], [10])])
end


############################################################################################
# metadata(ms::AbstractMassSpectrum)
############################################################################################
@testset "metadata(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test Dict() == metadata(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))
    @test Dict(:id => 4, :name => "sample") == metadata(MassSpectrum([85.1, 112.2, 124.1], 
        [13, 0, 67], metadata=Dict(:id => 4, :name => "sample")))

    # Same return container and element type as used to construct the object
    @test Dict{Any, Any} == typeof(metadata(ChromMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 4, :name => "sample"))))   
end


############################################################################################
# minintensity(ms::AbstractMassSpectrum)
############################################################################################
@testset "minintensity(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test 0 == minintensity(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))

    # Same element type as used to construct the object
    @test Int == typeof(minintensity(MassSpectrum([85.1, 112.2, 124.1], Int[13, 0, 67])))
    @test Float64 == typeof(minintensity(MassSpectrum([85.1, 112.2, 124.1], 
        Float64[13, 0, 67])))
end


############################################################################################
# minion(ms::AbstractMassSpectrum)
############################################################################################
@testset "minion(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test 85.1 == minion(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))

    # Same element type as used to construct the object
    @test Float64 == typeof(minion(MassSpectrum(Float64[85.1, 112.2, 124.1], [13, 0, 67])))
    @test Int == typeof(minion(MassSpectrum([85, 112, 124], [13, 0, 67])))
end


############################################################################################
# retentionindex(ms::AbstractMassSpectrum)
############################################################################################
@testset "retentionindex(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test nothing === retentionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))
    @test 123.2 == retentionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats", retentionindex=123.2))

    # Same element type as used to construct the object
    @test Nothing == typeof(retentionindex(MassSpectrum([85.1, 112.2, 124.1], 
        [13, 0, 67])))
    @test Float64 == typeof(retentionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats", retentionindex=123.2)))
    @test Int == typeof(retentionindex(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats", retentionindex=Int(123))))
end


############################################################################################
# retentionindexname(ms::AbstractMassSpectrum)
############################################################################################
@testset "retentionindexname(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test nothing === retentionindexname(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))
    @test "Kovats" == retentionindexname(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentionindexname="Kovats", retentionindex=123.2))

    # Same element type as used to construct the object
    @test Nothing == typeof(retentionindexname(MassSpectrum([85.1, 112.2, 124.1], 
        [13, 0, 67])))
    @test String == typeof(retentionindexname(MassSpectrum([85.1, 112.2, 124.1], 
        [13, 0, 67], retentionindexname="Kovats", retentionindex=123.2)))
    @test SubString{String} == typeof(retentionindexname(MassSpectrum([85.1, 112.2, 124.1], 
        [13, 0, 67], retentionindexname=SubString("Kovats", 1, 6), 
        retentionindex=Int(123))))
end


############################################################################################
# retentiontime(ms::AbstractMassSpectrum; timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "retentiontime(ms)" begin
    # Same return values as those provided as arguments to construct the object
    @test nothing === retentiontime(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67]))
    @test 123.2u"s" ≈ retentiontime(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentiontime=123.2u"s"))
    @test 2.0533333333333332u"minute" ≈ retentiontime(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentiontime=123.2u"s"), timeunit=u"minute")
    @test 2.0533333333333332 ≈ retentiontime(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentiontime=123.2u"s"), timeunit=u"minute", ustripped=true)
    @test 123.2 ≈ retentiontime(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentiontime=123.2u"s"), ustripped=true)

    # Same element type as used to construct the object
    @test Nothing == typeof(retentiontime(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67])))
    @test Unitful.Quantity{Float64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(retentiontime(MassSpectrum([85.1, 112.2, 124.1], 
        [13, 0, 67], retentiontime=123.2u"minute")))
    @test Unitful.Quantity{Int64, 𝐓, Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓, nothing}} == typeof(retentiontime(MassSpectrum([85.1, 112.2, 124.1], 
        [13, 0, 67], retentiontime=Int(123)*u"minute")))
    @test Int == typeof(retentiontime(MassSpectrum([85.1, 112.2, 124.1], [13, 0, 67], 
        retentiontime=Int(123)*u"minute"), ustripped=true))
end


############################################################################################
# sharedions(ms₁::AbstractMassSpectrum, ms₂::AbstractMassSpectrum)
############################################################################################
@testset "sharedions(ms₁, ms₂)" begin
    # Same return the correct values
    @test [1, 2] == sharedions(MassSpectrum([1, 2, 3], [1, 2, 3]), 
        MassSpectrum([1, 2, 4], [1, 2, 3]))
    @test [] == sharedions(MassSpectrum([1, 2, 3], [1, 2, 3]), 
        MassSpectrum([4, 5, 6], [1, 2, 3]))

    # Return the same element type as used to construct the object or their promoted type
    @test Vector{Int} == typeof(sharedions(MassSpectrum(Int[1, 2, 3], [1, 2, 3]), 
        MassSpectrum(Int[1, 2, 4], [1, 2, 3])))
    @test Vector{Float64} == typeof(sharedions(MassSpectrum(Float64[1, 2, 3], [1, 2, 3]), 
        MassSpectrum(Int[1, 2, 4], [1, 2, 3])))
    @test Vector{Float64} == typeof(sharedions(MassSpectrum(Int[1, 2, 3], [1, 2, 3]), 
        MassSpectrum(Float64[1, 2, 4], [1, 2, 3])))
end


############################################################################################
# similarity(ms₁::AbstractMassSpectrum, ms₂::AbstractMassSpectrum, f::Function)
############################################################################################
@testset "sharedions(ms₁, ms₂)" begin
    @test 1.0 ≈ similarity(MassSpectrum([80, 85, 90], [100, 500, 250]), 
        MassSpectrum([80, 85], [200, 1000]), cosine)
    @test 0.8978872704229618 ≈ similarity(MassSpectrum([80, 85, 90], [100, 500, 250]), 
        MassSpectrum([80, 85, 90], [200, 1000, 0]), cosine)

    @test Float64 == typeof(similarity(MassSpectrum([80, 85, 90], [100, 500, 250]), 
        MassSpectrum([80, 85, 90], [200, 1000, 0]), cosine))
    
    @test_throws ArgumentError typeof(similarity(MassSpectrum(Int[], Int[]), 
        MassSpectrum([80, 85, 90], [200, 1000, 0]), cosine))
end
