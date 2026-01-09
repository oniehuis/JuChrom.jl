module TestScans

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core/containers/scans.jl  (containers only; no getters)
# ─────────────────────────────────────────────────────────────────────────────

using Test
using Unitful
using Unitful: s, ms, μs, minute, pA

using JuChrom  # assumes ChromScan, MassScan are exported

const ATTRS_TUPLE = (column="1", flowrate=0.1)

# ─────────────────────────────────────────────────────────────────────────────
# ChromScan
# ─────────────────────────────────────────────────────────────────────────────

@testset "ChromScan – inner constructor" begin
    # Valid finite values
    cs = ChromScan{Float64, typeof(u"s"), Float64, typeof(pA), NamedTuple}(
        1.0, u"s", 100.0, pA, NamedTuple())
    @test cs.retention === 1.0
    @test cs.retention_unit == u"s"
    @test cs.intensity === 100.0
    @test cs.intensity_unit == pA
    @test cs.attrs == NamedTuple()

    # Attrs preserved
    cs2 = ChromScan{Float64, typeof(u"minute"), Float64, Nothing, typeof(ATTRS_TUPLE)}(
        5.0, u"minute", 42.0, nothing, ATTRS_TUPLE)
    @test cs2.attrs == ATTRS_TUPLE

    # Rejections: non-finite retention / intensity
    @test_throws ArgumentError ChromScan{Float64, typeof(u"s"), Float64, Nothing, NamedTuple}(
        NaN, u"s", 1.0, nothing, NamedTuple())
    @test_throws ArgumentError ChromScan{Float64, typeof(u"s"), Float64, Nothing, NamedTuple}(
        Inf, u"s", 1.0, nothing, NamedTuple())
    @test_throws ArgumentError ChromScan{Float64, typeof(u"s"), Float64, Nothing, NamedTuple}(
        1.0, u"s", NaN, nothing, NamedTuple())
    @test_throws ArgumentError ChromScan{Float64, typeof(u"s"), Float64, Nothing, NamedTuple}(
        1.0, u"s", -Inf, nothing, NamedTuple())
end

@testset "ChromScan – outer constructors" begin
    # Unitful retention & intensity → units stripped into *_unit fields
    c1 = ChromScan(5u"minute", 100u"pA")
    @test c1.retention == 5
    @test c1.retention_unit == u"minute"
    @test c1.intensity == 100
    @test c1.intensity_unit == pA
    @test c1.attrs == NamedTuple()

    # Unitful retention, unitless intensity
    c2 = ChromScan(2u"s", 10.0; attrs=ATTRS_TUPLE)
    @test c2.retention == 2
    @test c2.retention_unit == u"s"
    @test c2.intensity == 10.0
    @test c2.intensity_unit === nothing
    @test c2.attrs == ATTRS_TUPLE

    # Unitless retention, unitful intensity
    c3 = ChromScan(12.5, 3u"pA")
    @test c3.retention === 12.5
    @test c3.retention_unit === nothing
    @test c3.intensity === 3
    @test c3.intensity_unit == pA

    # Unitless retention & intensity
    c4 = ChromScan(100.0, 200.0)
    @test c4.retention === 100.0
    @test c4.retention_unit === nothing
    @test c4.intensity === 200.0
    @test c4.intensity_unit === nothing

    # Explicit (unit-stripped) constructor
    c5 = ChromScan(1.5, u"s", 7.0, pA; attrs=(id=1,))
    @test c5.retention === 1.5
    @test c5.retention_unit == u"s"
    @test c5.intensity === 7.0
    @test c5.intensity_unit == pA
    @test c5.attrs == (id=1,)
end

# ─────────────────────────────────────────────────────────────────────────────
# MassScan
# ─────────────────────────────────────────────────────────────────────────────

@testset "MassScan – inner constructor" begin
    # Valid construction (all unitless vectors)
    ms = MassScan{Float64, typeof(u"s"),
                  Vector{Float64}, Nothing,
                  Vector{Float64}, Nothing,
                  Int, NamedTuple}(
        1.0, u"s",
        [100.0, 150.0], nothing,
        [10.0, 20.0],  nothing,
        1, NamedTuple())
    @test ms.retention === 1.0
    @test ms.retention_unit == u"s"
    @test ms.mz_values == [100.0, 150.0]
    @test ms.mz_unit === nothing
    @test ms.intensities == [10.0, 20.0]
    @test ms.intensity_unit === nothing
    @test ms.level === 1
    @test ms.attrs == NamedTuple()

    # Attrs preserved
    ms2 = MassScan{Float64, Nothing,
                   Vector{Float64}, Nothing,
                   Vector{Float64}, Nothing,
                   Int, typeof(ATTRS_TUPLE)}(
        0.5, nothing, [50.0], nothing, [5.0], nothing, 2, ATTRS_TUPLE)
    @test ms2.attrs == ATTRS_TUPLE

    # Rejections (constructor invariants)
    @test_throws ArgumentError MassScan{Float64, Nothing,
                                        Vector{Float64}, Nothing,
                                        Vector{Float64}, Nothing,
                                        Int, NamedTuple}(
        NaN, nothing, [100.0], nothing, [1.0], nothing, 1, NamedTuple())
    @test_throws ArgumentError MassScan{Float64, Nothing,
                                        Vector{Float64}, Nothing,
                                        Vector{Float64}, Nothing,
                                        Int, NamedTuple}(
        0.0, nothing, Float64[], nothing, [1.0], nothing, 1, NamedTuple()) # empty mz
    @test_throws ArgumentError MassScan{Float64, Nothing,
                                        Vector{Float64}, Nothing,
                                        Vector{Float64}, Nothing,
                                        Int, NamedTuple}(
        0.0, nothing, [100.0], nothing, Float64[], nothing, 1, NamedTuple()) # empty I
    @test_throws DimensionMismatch MassScan{Float64, Nothing,
                                            Vector{Float64}, Nothing,
                                            Vector{Float64}, Nothing,
                                            Int, NamedTuple}(
        0.0, nothing, [100.0, 101.0], nothing, [1.0], nothing, 1, NamedTuple())
    @test_throws ArgumentError MassScan{Float64, Nothing,
                                        Vector{Float64}, Nothing,
                                        Vector{Float64}, Nothing,
                                        Int, NamedTuple}(
        0.0, nothing, [100.0, 100.0], nothing, [1.0, 2.0], nothing, 1, NamedTuple()) # not strictly inc
    @test_throws ArgumentError MassScan{Float64, Nothing,
                                        Vector{Float64}, Nothing,
                                        Vector{Float64}, Nothing,
                                        Int, NamedTuple}(
        0.0, nothing, [101.0, 100.0], nothing, [1.0, 2.0], nothing, 1, NamedTuple()) # decreasing
    @test_throws ArgumentError MassScan{Float64, Nothing,
                                        Vector{Float64}, Nothing,
                                        Vector{Float64}, Nothing,
                                        Int, NamedTuple}(
        0.0, nothing, [-1.0, 2.0], nothing, [1.0, 2.0], nothing, 1, NamedTuple()) # nonpositive m/z
    @test_throws ArgumentError MassScan{Float64, Nothing,
                                        Vector{Float64}, Nothing,
                                        Vector{Float64}, Nothing,
                                        Int, NamedTuple}(
        0.0, nothing, [100.0, 150.0], nothing, [1.0, NaN], nothing, 1, NamedTuple()) # non-finite I
    @test_throws ArgumentError MassScan{Float64, Nothing,
                                        Vector{Float64}, Nothing,
                                        Vector{Float64}, Nothing,
                                        Int, NamedTuple}(
        0.0, nothing, [NaN, 150.0], nothing, [1.0, 2.0], nothing, 1, NamedTuple()) # non-finite m/z
    @test_throws ArgumentError MassScan{Float64, Nothing,
                                        Vector{Float64}, Nothing,
                                        Vector{Float64}, Nothing,
                                        Int, NamedTuple}(
        0.0, nothing, [100.0], nothing, [1.0], nothing, 0, NamedTuple()) # level < 1
end

@testset "MassScan – outer constructors" begin
    # Defaults: level=1, attrs=()
    m1 = MassScan(2.0u"s", [100.0, 150.0], [10.0, 20.0])
    @test m1.retention == 2.0
    @test m1.retention_unit == u"s"
    @test m1.mz_values == [100.0, 150.0]
    @test m1.mz_unit === nothing
    @test m1.intensities == [10.0, 20.0]
    @test m1.intensity_unit === nothing
    @test m1.level === 1
    @test m1.attrs == NamedTuple()

    # Unitful intensities vector → strip values, keep unit
    m2 = MassScan(1.0u"s", [120.0, 160.0], [600.0, 1100.0]u"pA"; level=2, attrs=ATTRS_TUPLE)
    @test m2.retention == 1.0
    @test m2.retention_unit == u"s"
    @test m2.intensities == [600.0, 1100.0]
    @test m2.intensity_unit == pA
    @test m2.level === 2
    @test m2.attrs == ATTRS_TUPLE

    # Mixed unitful/unitless within a vector should be rejected (outer constructor validation)
    @test_throws ArgumentError MassScan(1.0u"s", [100.0, 150.0], [1.0, 2.0u"pA"])
    @test_throws ArgumentError MassScan(1.0u"s", [100.0, 150.0u"s"], [1.0, 2.0])

    # Large-ish vectors (basic shape checks)
    ions = collect(100.0:0.01:100.3)
    ints = fill(1.0, length(ions))
    m3 = MassScan(0.5u"s", ions, ints)
    @test length(m3.mz_values) == length(m3.intensities)

    # Unitful m/z and intensity vectors with unitless retention
    m4 = MassScan(2.5, [100.0, 150.0]u"Th", [10.0, 20.0]u"pA"; level=3, attrs=ATTRS_TUPLE)
    @test m4.retention == 2.5
    @test m4.retention_unit === nothing
    @test m4.mz_values == [100.0, 150.0]
    @test m4.mz_unit == u"Th"
    @test m4.intensities == [10.0, 20.0]
    @test m4.intensity_unit == pA
    @test m4.level === 3
    @test m4.attrs == ATTRS_TUPLE
end

end # module
