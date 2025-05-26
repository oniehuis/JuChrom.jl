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
        Nothing)

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
    
    # Must provide at least four retention times and retention indices
    @test_throws ArgumentError RiMapper("Kovats", (1:1)u"minute", (1000:1000:1000))
    @test_throws ArgumentError RiMapper("Kovats", (1:2)u"minute", (1000:1000:2000))
    @test_throws ArgumentError RiMapper("Kovats", (1:3)u"minute", (1000:1000:3000))

    # Numer of retention times and of retention indices must be identical
    @test_throws DimensionMismatch RiMapper("Kovats", (1:4)u"minute", (1000:1000:5000))
    @test_throws DimensionMismatch RiMapper("Kovats", (1:5)u"minute", (1000:1000:4000))

    # Retention times must be in ascending order
    @test_throws ArgumentError RiMapper("Kovats", [2, 1, 3, 4]u"minute", (1000:1000:4000))

    # No identical retention times allowed
    @test_throws ArgumentError RiMapper("Kovats", [1, 1, 2, 3]u"minute", (1000:1000:4000))

    # Retention indices must be in ascending order
    @test_throws ArgumentError RiMapper("Kovats", (1:4)u"minute", [20, 10, 30, 40])
    
    # No identical retention indices allowed
    @test_throws ArgumentError RiMapper("Kovats", (1:4)u"minute", [10, 10, 20, 30])

    # Must provide retention index name
    @test_throws ArgumentError RiMapper("", (1:4)u"minute", (1000:1000:4000))
    @test_throws MethodError RiMapper(nothing, (1:4)u"minute", (1000:1000:4000))

    # Metadata cannot be anything other than a dictionary.
    @test_throws TypeError RiMapper("Kovats", (1:4)u"minute", (1000:1000:4000), 
        metadata=:sample_name)
    @test_throws TypeError RiMapper("Kovats", (1:4)u"minute", (1000:1000:4000), 
        metadata=1)
    @test_throws TypeError RiMapper("Kovats", (1:4)u"minute", (1000:1000:4000), 
        metadata='S')

    # RiMapper is broadcastable
    @test [1, 2] == ((ld, i) -> i).(RiMapper("Kovats", (1:4)u"minute", 1:4), [1, 2])
end


@testset "show RiMapper" begin
    io = IOBuffer()
    show(io, RiMapper("Kovats", (1:4)u"minute", 1000:1000:4000))
    @test String(take!(io)) == string(
        "JuChrom.RiMapper {index name: Kovats, calibration points: 4}\n",
        "retention times: 1 minute, 2 minute, 3 minute, 4 minute\n",
        "retention indices: 1000, 2000, 3000, 4000\n",
        "interpolation method: JuChrom.NaturalCubicBSpline(false)\n",
        "extrapolation method: nothing\n",
        "metadata: 0 entries")
    show(io, RiMapper("Kovats", (1:11)u"minute", 1000:1000:11000))
    @test String(take!(io)) == string(
        "JuChrom.RiMapper {index name: Kovats, calibration points: 11}\n",
        "retention time range: 1 minute - 11 minute\n",
        "retention index range: 1000 - 11000\n",
        "interpolation method: JuChrom.NaturalCubicBSpline(false)\n",
        "extrapolation method: nothing\n",
        "metadata: 0 entries")
    show(io, RiMapper("Kovats", (1:11)u"minute", 1000:1000:11000, metadata=Dict(:a => 1)))
    @test String(take!(io)) == string(
        "JuChrom.RiMapper {index name: Kovats, calibration points: 11}\n",
        "retention time range: 1 minute - 11 minute\n",
        "retention index range: 1000 - 11000\n",
        "interpolation method: JuChrom.NaturalCubicBSpline(false)\n",
        "extrapolation method: nothing\n",
        "metadata: 1 entry")
    show(io, RiMapper("Kovats", (1:11)u"minute", 1000:1000:11000, metadata=Dict(:a => 1, 
        :b => 2)))
    @test String(take!(io)) == string(
        "JuChrom.RiMapper {index name: Kovats, calibration points: 11}\n",
        "retention time range: 1 minute - 11 minute\n",
        "retention index range: 1000 - 11000\n",
        "interpolation method: JuChrom.NaturalCubicBSpline(false)\n",
        "extrapolation method: nothing\n",
        "metadata: 2 entries")
end


############################################################################################
# interpolationmethod(mapper::RiMapper)
############################################################################################
@testset "interpolationmethod" begin
    @test NaturalCubicBSpline(force=false) == interpolationmethod(RiMapper("Kovats", 
        (1:4)u"minute", 1:4))
    @test NaturalCubicBSpline(force=false) == interpolationmethod(RiMapper("Kovats", 
        (1:4)u"minute", 1:4, interpolationmethod=NaturalCubicBSpline()))
    @test NaturalCubicBSpline(force=false) == interpolationmethod(RiMapper("Kovats", 
        (1:4)u"minute", 1:4, interpolationmethod=NaturalCubicBSpline(), 
        extrapolationmethod=Linear()))
    @test NaturalCubicBSpline(force=false) == interpolationmethod(RiMapper("Kovats", 
        (1:4)u"minute", 1:4, interpolationmethod=NaturalCubicBSpline(), 
        extrapolationmethod=nothing))
    
    @test_throws ArgumentError interpolationmethod(RiMapper("Kovats", 
       (1:6)u"s", [1, 2, 2.1, 5, 5.1, 9], interpolationmethod=NaturalCubicBSpline(), 
       extrapolationmethod=nothing))
    @test NaturalCubicBSpline(force=true) == interpolationmethod(RiMapper("Kovats", 
        (1:6)u"s", [1, 2, 2.1, 5, 5.1, 9], 
        interpolationmethod=NaturalCubicBSpline(force=true), extrapolationmethod=nothing))
end


############################################################################################
# extrapolationmethod(mapper::RiMapper)
############################################################################################
@testset "extrapolationmethod" begin
    @test nothing === extrapolationmethod(RiMapper("Kovats", (1:4)u"minute", 1:4))
    @test nothing === extrapolationmethod(RiMapper("Kovats", (1:4)u"minute", 1:4, 
        interpolationmethod=NaturalCubicBSpline()))
    @test nothing === extrapolationmethod(RiMapper("Kovats", (1:4)u"minute", 1:4, 
        interpolationmethod=PiecewiseLinear()))
    @test Linear() == extrapolationmethod(RiMapper("Kovats", (1:4)u"minute", 1:4, 
        extrapolationmethod=Linear()))
    @test Linear() == extrapolationmethod(RiMapper("Kovats", (1:4)u"minute", 1:4, 
        interpolationmethod=NaturalCubicBSpline(), extrapolationmethod=Linear()))
    @test Linear() == extrapolationmethod(RiMapper("Kovats", (1:4)u"minute", 1:4, 
        interpolationmethod=PiecewiseLinear(), extrapolationmethod=Linear()))
end


############################################################################################
# maxretentionindex(mapper::RiMapper)
############################################################################################
@testset "maxretentionindex" begin
    @test 4 == maxretentionindex(RiMapper("Kovats", (1:4)u"minute", 1:4))
end


############################################################################################
# maxretentiontime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
#   ustripped::Bool=false)
############################################################################################
@testset "maxretentiontime" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8, 4.6]u"minute", [1000, 2000, 3000, 4000]);
    @test 4.6u"minute" ≈ maxretentiontime(ld)
    @test 276.0u"s" ≈ maxretentiontime(ld, timeunit=u"s")
    @test 276.0 ≈ maxretentiontime(ld, timeunit=u"s", ustripped=true)
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
    @test 1 == minretentionindex(RiMapper("Kovats", (1:4)u"minute", 1:4))
end


############################################################################################
# minretentiontime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
#   ustripped::Bool=false)
############################################################################################
@testset "minretentiontime" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8, 4.6]u"minute", [1000, 2000, 3000, 4000]);
    @test 1.2u"minute" ≈ minretentiontime(ld)
    @test 72.0u"s" ≈ minretentiontime(ld, timeunit=u"s")
    @test 72.0 ≈ minretentiontime(ld, timeunit=u"s", ustripped=true)
end


############################################################################################
# retentionindex(mapper::RiMapper, retentiontime::Unitful.Time; info::Bool=false)
############################################################################################
@testset "retentionindex(::RiMapper, ...)" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8, 4.5]u"minute", [1000, 2000, 3000, 4000])
    @test 1540.7444668008052 ≈ retentionindex(ld, 1.8u"minute")
    @test_throws ArgumentError retentionindex(ld, 1.1u"minute")
    @test 1000.0 ≈ retentionindex(ld, 1.2u"minute")
    @test_throws ArgumentError retentionindex(ld, 4.6u"minute")
    @test [1706.9081153588197, 2351.619923623373] ≈ retentionindex.(ld, [2, 3]u"minute")
    
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8, 4.5]u"minute", [1000, 2000, 3000, 4000], 
        extrapolationmethod=Linear())
    @test 1540.7444668008052 ≈ retentionindex(ld, 1.8u"minute") 
    @test 907.6123407109322 ≈ retentionindex(ld, 1.1u"minute")
    @test 1000.0 ≈ retentionindex(ld, 1.2u"minute")
    @test 3252.2481829754033 ≈ retentionindex(ld, 4.0u"minute")
    @test [1275.4652917505032, 2000] ≈ retentionindex.(ld, [1.5, 2.4]u"minute")
    @test [815.2246814218645, 1706.9081153588197] ≈ retentionindex.(ld, [1, 2]u"minute")
end


############################################################################################
# retentionindex(chrom::AbstractChromatogram, retentiontime::Unitful.Time; info::Bool=false)
############################################################################################
@testset "retentionindex(::AbstractChromatogram, ...)" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8, 4.5]u"minute", [1000, 2000, 3000, 4000])
    chrom = Chrom(Int64[1, 2, 3]u"s", Int32[12, 956, 1], rimapper=ld)
    @test 1540.7444668008052 ≈ retentionindex(chrom, 1.8u"minute")
    @test_throws ArgumentError retentionindex(chrom, 1.1u"minute")
    @test 1000.0 ≈ retentionindex(chrom, 1.2u"minute")
    @test_throws ArgumentError retentionindex(chrom, 5.0u"minute")
    @test [1706.9081153588197, 2351.619923623373] ≈ retentionindex.(chrom, [2, 3]u"minute")

    ld = RiMapper("Kovats", [1.2, 2.4, 3.8, 4.5]u"minute", [1000, 2000, 3000, 4000], 
        extrapolationmethod=Linear())
    chrom = Chrom(Int64[1, 2, 3]u"s", Int32[12, 956, 1], rimapper=ld)
    @test 1540.7444668008052 ≈ retentionindex(chrom, 1.8u"minute")
    @test 907.6123407109322 ≈ retentionindex(chrom, 1.1u"minute")
    @test 1000.0 ≈ retentionindex(chrom, 1.2u"minute")
    @test 4782.6123407109335 ≈ retentionindex(chrom, 5.0u"minute")
    @test [815.2246814218645, 1706.9081153588197] ≈ retentionindex.(chrom, [1, 2]u"minute")
    @test [1275.4652917505032, 2000.0] ≈ retentionindex.(chrom, [1.5, 2.4]u"minute")
end


############################################################################################
# retentionindexname(mapper::RiMapper)
############################################################################################
@testset "retentionindexname" begin
    @test "Kovats" == retentionindexname(RiMapper("Kovats", (1:4)u"s", 1:4))
end


############################################################################################
# retentionindices(mapper::RiMapper)
############################################################################################
@testset "retentionindices" begin
    @test 1:4 == retentionindices(RiMapper("Kovats", (1:4)u"s", 1:4))
end


############################################################################################
# retentiontimes(mapper::RiMapper; timeunit::Unitful.TimeUnits, ustripped::Bool=false)
############################################################################################
@testset "retentiontimes" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 3.8, 4.5]u"minute", [1000, 2000, 3000, 4000])
    @test [1.2, 2.4, 3.8, 4.5]u"minute" ≈ retentiontimes(ld)
    @test [72, 144, 228, 270]u"s" ≈ retentiontimes(ld, timeunit=u"s")
    @test [72, 144, 228, 270] ≈ retentiontimes(ld, timeunit=u"s", ustripped=true)
    @test [1.2, 2.4, 3.8, 4.5] ≈ retentiontimes(ld, ustripped=true)
end


############################################################################################
# JuChrom.rt2ri(mapper::RiMapper)
############################################################################################
@testset "JuChrom.rt2ri" begin
    ld = RiMapper("Kovats", [1.2, 2.4, 4.3, 5.9]u"minute", [1000, 2000, 3000, 4000])
    @test 1531.7108900675921 ≈ retentionindex(ld, 1.8u"minute")
    @test_throws ArgumentError retentionindex(ld, 1.1u"minute")
    @test 1000.0 ≈ retentionindex(ld, 1.2u"minute")
    @test 3053.7704734395047 ≈ retentionindex(ld, 4.4u"minute")
    @test [1697.9860642642884, 2362.1065650064497] ≈ retentionindex.(ld, [2, 3]u"minute")
    
    ld = RiMapper("Kovats", [1.2, 2.4, 4.3, 5.9]u"minute", [1000, 2000, 3000, 4000], 
        extrapolationmethod=Linear())
    @test 1531.7108900675921 ≈ retentionindex(ld, 1.8u"minute")
    @test 909.619802207202 ≈ retentionindex(ld, 1.1u"minute")
    @test 1000.0 ≈ retentionindex(ld, 1.2u"minute")
    @test 3053.7704734395047 ≈ retentionindex(ld, 4.4u"minute")
    @test [819.2396044144039, 1697.9860642642884] ≈ retentionindex.(ld, [1, 2]u"minute")
end


############################################################################################
# JuChrom.bsplineinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
# retentionindices::AbstractVector{<:Real}; extrapolation::Bool=true)
############################################################################################
@testset "bsplineinterpolation" begin
    rt2ri, jacobian_rt2ri = JuChrom.bsplineinterpolation([1, 2, 3, 4, 5, 6, 7, 8]*u"s", 
        [1000, 1800, 3050, 3800, 5500, 6600, 6900, 7400])
    @test 1000.0 ≈ rt2ri(1u"s")
    @test 1333.7469941600825 ≈ rt2ri(1.5u"s")
    @test 1800.0 ≈ rt2ri((1//30)u"minute")
    @test_throws ArgumentError rt2ri(9.1u"s")
    @test_throws ArgumentError rt2ri(0.5u"s")
    @test [1000.0, 1800.0] ≈ rt2ri.([1, 2]u"s")

    rt2ri, jacobian_rt2ri = JuChrom.bsplineinterpolation([1, 2, 3, 4, 5, 6, 7, 8]*u"s", 
        [1000, 1800, 3050, 3800, 5500, 6600, 6900, 7400]; extrapolation=true)
    @test 1000.0 ≈ rt2ri(1u"s")
    @test 1333.7469941600825 ≈ rt2ri(1.5u"s")
    @test 1800.0 ≈ rt2ri((1//30)u"minute")
    @test rt2ri(9.1u"s") ≈ 8053.1226382686355
    @test rt2ri(0.5u"s") ≈ 688.3373411198902
    @test rt2ri.([0.5, 1]u"s") ≈ [688.3373411198902, 1000.0000000000001]
    @test [1000.0, 1800.0] ≈ rt2ri.([1, 2]u"s")

    @test_throws ArgumentError JuChrom.bsplineinterpolation([1, 2, 3]u"s", 
        [1000, 2000, 3000])
    @test_throws DimensionMismatch JuChrom.bsplineinterpolation([1, 2, 3, 4, 5]u"s", 
        [1000, 1900, 3100, 3900])
    @test_throws DimensionMismatch JuChrom.bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1900, 3100, 3900, 4000])
    @test_throws ArgumentError JuChrom.bsplineinterpolation([1, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900])
    @test_throws ArgumentError JuChrom.bsplineinterpolation([2, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900])
    @test_throws ArgumentError JuChrom.bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1000, 3100, 3900])
    @test_throws ArgumentError JuChrom.bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1900, 1000, 3100, 3900])
    @test_throws ArgumentError JuChrom.bsplineinterpolation([1]u"s", [1000], 
        extrapolation=true)
    @test_throws DimensionMismatch JuChrom.bsplineinterpolation([1, 2, 3, 4, 5]u"s", 
        [1000, 1900, 3100, 3900], extrapolation=true)
    @test_throws DimensionMismatch JuChrom.bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1900, 3100, 3900, 4000], extrapolation=true)
    @test_throws ArgumentError JuChrom.bsplineinterpolation([1, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900], extrapolation=true)
    @test_throws ArgumentError JuChrom.bsplineinterpolation([2, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900], extrapolation=true)
    @test_throws ArgumentError JuChrom.bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1000, 3100, 3900], extrapolation=true)
    @test_throws ArgumentError JuChrom.bsplineinterpolation([1, 2, 3, 4]u"s", 
        [1900, 1000, 3100, 3900], extrapolation=true)
end


############################################################################################
# JuChrom.piecewiselinearinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
# retentionindices::AbstractVector{<:Real}; extrapolation::Bool=true)
############################################################################################
@testset "piecewiselinearinterpolation" begin
    rt2ri = JuChrom.piecewiselinearinterpolation([1, 2, 3, 4, 5, 6, 7, 8]*u"s", 
        [1000, 1800, 3050, 3800, 5500, 6600, 6900, 7400])
    @test 1000 ≈ rt2ri(1u"s")
    @test 1400 ≈ rt2ri(1.5u"s")
    @test 1800 ≈ rt2ri((1//30)u"minute")
    @test_throws ArgumentError rt2ri(9.1u"s")
    @test_throws ArgumentError rt2ri(0.5u"s")
    @test [1000, 2425] ≈ rt2ri.([1, 2.5]u"s")

    rt2ri = JuChrom.piecewiselinearinterpolation([1, 2, 3, 4, 5, 6, 7, 8]*u"s", 
        [1000, 1800, 3050, 3800, 5500, 6600, 6900, 7400]; extrapolation=true)
    @test 1000 ≈ rt2ri(1u"s")
    @test 1400 ≈ rt2ri(1.5u"s")
    @test 1800 ≈ rt2ri((1//30)u"minute")
    @test 7950.0 ≈ rt2ri(9.1u"s")
    @test 600.0 ≈ rt2ri(0.5u"s")
    @test [1000, 2425] ≈ rt2ri.([1, 2.5]u"s")

    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([1]u"s", [1000])
    @test_throws DimensionMismatch JuChrom.piecewiselinearinterpolation([1, 2, 3, 4, 5]u"s", 
        [1000, 1900, 3100, 3900])
    @test_throws DimensionMismatch JuChrom.piecewiselinearinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1900, 3100, 3900, 4000])
    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([1, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900])
    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([2, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900])
    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1000, 3100, 3900])
    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([1, 2, 3, 4]u"s", 
        [1900, 1000, 3100, 3900])
    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([1]u"s", [1000], 
        extrapolation=true)
    @test_throws DimensionMismatch JuChrom.piecewiselinearinterpolation([1, 2, 3, 4, 5]u"s", 
        [1000, 1900, 3100, 3900], extrapolation=true)
    @test_throws DimensionMismatch JuChrom.piecewiselinearinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1900, 3100, 3900, 4000], extrapolation=true)
    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([1, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900], extrapolation=true)
    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([2, 1, 3, 4]u"s", 
        [1000, 1900, 3100, 3900], extrapolation=true)
    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([1, 2, 3, 4]u"s", 
        [1000, 1000, 3100, 3900], extrapolation=true)
    @test_throws ArgumentError JuChrom.piecewiselinearinterpolation([1, 2, 3, 4]u"s", 
        [1900, 1000, 3100, 3900], extrapolation=true)
end
