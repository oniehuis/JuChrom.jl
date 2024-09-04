using JuChrom
using Test
using Unitful: 𝐓

############################################################################################
# RiMapper(retentionindexname::AbstractString, 
# retentiontimes::AbstractVector{<:Unitful.Time},
# retentionindices::AbstractVector{<:Real};
# interpolationmethod::PolationMethod=NaturalCubicBSpline(), 
# extrapolationmethod::Union{Nothing, <:PolationMethod}=Linear(),
# metadata::Dict=Dict())
############################################################################################
@testset "RiMapper outer constructor" begin
    # Verify object construction and field types depending on constructor arguments
    # Retention index name
    @test "Kovats" == RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000).retentionindexname
    @test isa(RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000).retentionindexname, 
        AbstractString)

    # Retention times
    @test (1:5)u"minute" == RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000
        ).retentiontimes
    @test isa(RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000).retentiontimes, 
        AbstractVector)

    # Retention indices
    @test 1000:1000:5000 == RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000
        ).retentionindices
    @test isa(RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000).retentionindices, 
        AbstractVector)

    # Interpolation method
    @test isa(RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000).interpolationmethod, 
        PolationMethod)

    # Extrapolation method
    @test isa(RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000).extrapolationmethod, 
        PolationMethod)

    # rt2ri
    @test isa(RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000).rt2ri, Function)

    # Metadata
    @test Dict{Any, Any} == typeof(RiMapper("Kovats", (1:5)u"minute", 
        1000:1000:5000).metadata)
    @test Dict() == RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000).metadata
    @test Dict(:id => 4, "name" => "sample") == RiMapper("Kovats", (1:5)u"minute", 
        1000:1000:5000, metadata=Dict(:id => 4, "name" => "sample")).metadata

    # Check the associated supertypes
    @test isa(RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000), AbstractRiMapper)

    # Retention time vector accepts only time values
    @test_throws MethodError RiMapper("Kovats", 1:5, 1000:1000:5000)
    @test_throws MethodError RiMapper("Kovats", (1:5)u"m", 1000:1000:5000)

    # Retention index vector accepts no unit
    @test_throws MethodError RiMapper("Kovats", (1:5)u"minute", (1000:1000:5000)u"m")
    
    # Must provide at least two retention times and retention indices
    @test_throws ArgumentError RiMapper("Kovats", (1:1)u"minute", (1000:1000:1000))

    # Numer of retention times and of retention indices must be identical
    @test_throws DimensionMismatch RiMapper("Kovats", (1:2)u"minute", (1000:1000:3000))
    @test_throws DimensionMismatch RiMapper("Kovats", (1:3)u"minute", (1000:1000:2000))

    # Retention times must be in ascending order
    @test_throws ArgumentError RiMapper("Kovats", [2, 1]u"minute", (1000:1000:2000))

    # No identical retention times allowed
    @test_throws ArgumentError RiMapper("Kovats", [1, 1]u"minute", (1000:1000:2000))

    # Retention indices must be in ascending order
    @test_throws ArgumentError RiMapper("Kovats", (1:2)u"minute", [2000, 1000])
    
    # No identical retention indices allowed
    @test_throws ArgumentError RiMapper("Kovats", (1:2)u"minute", [1000, 1000])

    # Must provide retentzion index name
    @test_throws ArgumentError RiMapper("", (1:2)u"minute", (1000:1000:2000))
    @test_throws MethodError RiMapper(nothing, (1:2)u"minute", (1000:1000:2000))

    # Metadata cannot be anything other than a dictionary.
    @test_throws TypeError RiMapper("Kovats", (1:2)u"minute", (1000:1000:2000), 
        metadata=:sample_name)
    @test_throws TypeError RiMapper("Kovats", (1:2)u"minute", (1000:1000:2000), 
        metadata=1)
    @test_throws TypeError RiMapper("Kovats", (1:2)u"minute", (1000:1000:2000), 
        metadata='S')

    # RiMapper is broadcastable
    @test [1, 2] == ((ld, i) -> i).(RiMapper("Kovats", (1:2)u"minute", 1:2), [1, 2])
end


@testset "show RiMapper" begin
    io = IOBuffer()
    show(io, RiMapper("Kovats", (1:2)u"minute", 1000:1000:2000))
    @test String(take!(io)) == string(
        "JuChrom.RiMapper {index name: Kovats, calibration points: 2}\n",
        "retention times: 1 minute, 2 minute\n",
        "retention indices: 1000, 2000\n",
        "interpolation method: JuChrom.NaturalCubicBSpline()\n",
        "extrapolation method: JuChrom.Linear()\n",
        "metadata: 0 entries")
    show(io, RiMapper("Kovats", (1:11)u"minute", 1000:1000:11000))
    @test String(take!(io)) == string(
        "JuChrom.RiMapper {index name: Kovats, calibration points: 11}\n",
        "retention time range: 1 minute - 11 minute\n",
        "retention index range: 1000 - 11000\n",
        "interpolation method: JuChrom.NaturalCubicBSpline()\n",
        "extrapolation method: JuChrom.Linear()\n",
        "metadata: 0 entries")
    show(io, RiMapper("Kovats", (1:11)u"minute", 1000:1000:11000, metadata=Dict(:a => 1)))
    @test String(take!(io)) == string(
        "JuChrom.RiMapper {index name: Kovats, calibration points: 11}\n",
        "retention time range: 1 minute - 11 minute\n",
        "retention index range: 1000 - 11000\n",
        "interpolation method: JuChrom.NaturalCubicBSpline()\n",
        "extrapolation method: JuChrom.Linear()\n",
        "metadata: 1 entry")
    show(io, RiMapper("Kovats", (1:11)u"minute", 1000:1000:11000, metadata=Dict(:a => 1, 
        :b => 2)))
    @test String(take!(io)) == string(
        "JuChrom.RiMapper {index name: Kovats, calibration points: 11}\n",
        "retention time range: 1 minute - 11 minute\n",
        "retention index range: 1000 - 11000\n",
        "interpolation method: JuChrom.NaturalCubicBSpline()\n",
        "extrapolation method: JuChrom.Linear()\n",
        "metadata: 2 entries")
end


############################################################################################
# interpolationmethod(mapper::RiMapper)
############################################################################################
@testset "interpolationmethod" begin
    @test NaturalCubicBSpline() == interpolationmethod(RiMapper("Kovats", 
        (1:2)u"minute", 1:2))
end


############################################################################################
# extrapolationmethod(mapper::RiMapper)
############################################################################################
@testset "extrapolationmethod" begin
    @test Linear() == extrapolationmethod(RiMapper("Kovats", (1:2)u"minute", 1:2))
    @test nothing === extrapolationmethod(RiMapper("Kovats", (1:2)u"minute", 1:2, 
        extrapolationmethod=nothing))
end


############################################################################################
# maxretentionindex(mapper::RiMapper)
############################################################################################
@testset "maxretentionindex" begin
    @test 3 == maxretentionindex(RiMapper("Kovats", (1:3)u"minute", 1:3))
end


############################################################################################
# maxretentiontime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
#   ustripped::Bool=false)
############################################################################################
@testset "maxretentiontime" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);
    @test 3.8u"minute" ≈ maxretentiontime(ld)
    @test 228.0u"s" ≈ maxretentiontime(ld, timeunit=u"s")
    @test 228.0 ≈ maxretentiontime(ld, timeunit=u"s", ustripped=true)
end


############################################################################################
# metadata(mapper::AbstractRiMapper)
############################################################################################
@testset "metadata" begin
    @test Dict() == metadata(RiMapper("Kovats", (1:5)u"minute", 10:10:50))
    @test Dict(:id => 7) == metadata(RiMapper("Kovats", (1:5)u"minute", 10:10:50, 
        metadata=Dict(:id => 7)))
    @test Dict{Any, Any} == typeof(metadata(RiMapper("Kovats", (1:5)u"minute", 10:10:50)))
    @test Dict{Any, Any} == typeof(metadata(RiMapper("Kovats", (1:5)u"minute", 10:10:50, 
        metadata=Dict(:id => 7))))
end


############################################################################################
# minretentionindex(mapper::RiMapper)
############################################################################################
@testset "minretentionindex" begin
    @test 1 == minretentionindex(RiMapper("Kovats", (1:3)u"minute", 1:3))
end


############################################################################################
# minretentiontime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
#   ustripped::Bool=false)
############################################################################################
@testset "minretentiontime" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);
    @test 1.2u"minute" ≈ minretentiontime(ld)
    @test 72.0u"s" ≈ minretentiontime(ld, timeunit=u"s")
    @test 72.0 ≈ minretentiontime(ld, timeunit=u"s", ustripped=true)
end


############################################################################################
# retentionindex(mapper::RiMapper, retentiontime::Unitful.Time; info::Bool=false)
############################################################################################
@testset "retentionindex(::RiMapper, ...)" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000])
    @test 1512.3626373626375 ≈ retentionindex(ld, 1.8u"minute") 
    @test 913.9194139194141 ≈ retentionindex(ld, 1.1u"minute")
    @test 913.9194139194141 ≈ retentionindex(ld, 1.1u"minute", info=false)
    @test [827.8388278388279, 1678.876678876679] ≈ retentionindex.(ld, [1, 2]u"minute")
end


############################################################################################
# retentionindex(chrom::AbstractChromatogram, retentiontime::Unitful.Time; info::Bool=false)
############################################################################################
@testset "retentionindex(::AbstractChromatogram, ...)" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000])
    chrom = Chrom(Int64[1, 2, 3]u"s", Int32[12, 956, 1], rimapper=ld)
    @test 1512.3626373626375 ≈ retentionindex(chrom, 1.8u"minute") 
    @test 913.9194139194141 ≈ retentionindex(chrom, 1.1u"minute")
    @test 913.9194139194141 ≈ retentionindex(chrom, 1.1u"minute", info=false)
    @test [827.8388278388279, 1678.876678876679] ≈ retentionindex.(chrom, [1, 2]u"minute")
end


############################################################################################
# retentionindexname(mapper::RiMapper)
############################################################################################
@testset "retentionindexname" begin
    "Kovats" == retentionindexname(RiMapper("Kovats", (1:2)u"s", 1:2))
end


############################################################################################
# retentionindices(mapper::RiMapper)
############################################################################################
@testset "retentionindices" begin
    @test 1:2 == retentionindices(RiMapper("Kovats", (1:2)u"s", 1:2))
end


############################################################################################
# retentiontimes(mapper::RiMapper; timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "retentiontimes" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000])
    @test [1.2, 2.4, 3.8]u"minute" ≈ retentiontimes(ld)
    @test [72, 144, 228]u"s" ≈ retentiontimes(ld, timeunit=u"s")
    @test [72, 144, 228] ≈ retentiontimes(ld, timeunit=u"s", ustripped=true)
end


############################################################################################
# JuChrom.rt2ri(mapper::RiMapper)
############################################################################################
@testset "JuChrom.rt2ri" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 4.3]u"minute", [1000, 2000, 3000])
    @test 1266.712648556876 ≈ JuChrom.rt2ri(ld)(1.5u"minute") 
    @test 821.4487832484438 ≈ JuChrom.rt2ri(ld)(1u"minute")
end



############################################################################################
# bsplineinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
# retentionindices::AbstractVector{<:Real}; extrapolation::Bool=true)
############################################################################################
@testset "bsplineinterpolation" begin
    rt2ri = bsplineinterpolation([1, 2, 3, 4, 5, 6, 7, 8]*u"s", 
        [1000, 1800, 3050, 3800, 5500, 6600, 6900, 7400])
    @test 1000.0 ≈ rt2ri(1u"s")  
    @test 1333.7469941600825 ≈ rt2ri(1.5u"s")
    @test 1800.0 ≈ rt2ri((1//30)u"minute")
    @test 7950.0 ≈ rt2ri(9.1u"s")
    
    rt2ri = bsplineinterpolation([1, 2, 3, 4, 5, 6, 7, 8]*u"s", 
        [1000, 1800, 3050, 3800, 5500, 6600, 6900, 7400]; extrapolation=false)
    @test rt2ri(9.1u"s") === nothing

    @test_throws ArgumentError bsplineinterpolation([1]u"s", [1000])
    @test_throws DimensionMismatch bsplineinterpolation([1, 2, 3, 4, 5]u"s", 
        [1000, 1900, 3100, 3900])
    @test_throws DimensionMismatch bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1900, 3100, 3900, 4000])
    @test_throws ArgumentError bsplineinterpolation([1, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900])
    @test_throws ArgumentError bsplineinterpolation([2, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900])
    @test_throws ArgumentError bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1000, 3100, 3900])
    @test_throws ArgumentError bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1900, 1000, 3100, 3900])
    @test_throws ArgumentError bsplineinterpolation([1]u"s", [1000], 
        extrapolation=true)
    @test_throws DimensionMismatch bsplineinterpolation([1, 2, 3, 4, 5]u"s", 
        [1000, 1900, 3100, 3900], extrapolation=true)
    @test_throws DimensionMismatch bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1900, 3100, 3900, 4000], extrapolation=true)
    @test_throws ArgumentError bsplineinterpolation([1, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900], extrapolation=true)
    @test_throws ArgumentError bsplineinterpolation([2, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900], extrapolation=true)
    @test_throws ArgumentError bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1000, 3100, 3900], extrapolation=true)
    @test_throws ArgumentError bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1900, 1000, 3100, 3900], extrapolation=true)
end
