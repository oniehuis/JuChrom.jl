module TestMassSpectra

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core/containers/massspectra.jl  (containers only; no getters)
# ─────────────────────────────────────────────────────────────────────────────

using Test
using Unitful
using Unitful: pA

using JuChrom

const ATTRS_TUPLE = (instrument="Test",)

# ─────────────────────────────────────────────────────────────────────────────
# MassSpectrum
# ─────────────────────────────────────────────────────────────────────────────

@testset "MassSpectrum – inner constructor" begin
    mz = [100.0, 150.0]
    ints = [10.0, 20.0]

    ms = MassSpectrum{Vector{Float64}, typeof(u"Th"), Vector{Float64}, typeof(pA), NamedTuple}(
        mz, u"Th", ints, pA, NamedTuple())
    @test ms.mzvalues === mz
    @test ms.mzunit == u"Th"
    @test ms.intensities === ints
    @test ms.intensityunit == pA
    @test ms.attrs == NamedTuple()
    @test ms isa AbstractMassSpectrum{typeof(u"Th"), typeof(pA)}

    ms2 = MassSpectrum{Vector{Float64}, Nothing, Vector{Float64}, Nothing, typeof(ATTRS_TUPLE)}(
        mz, nothing, ints, nothing, ATTRS_TUPLE)
    @test ms2.attrs == ATTRS_TUPLE

    # Rejections: invalid m/z values
    @test_throws ArgumentError MassSpectrum{Vector{Float64}, Nothing, Vector{Float64}, Nothing, NamedTuple}(
        Float64[], nothing, [1.0], nothing, NamedTuple())

    @test_throws ArgumentError MassSpectrum{Vector{Float64}, Nothing, Vector{Float64}, Nothing, NamedTuple}(
        [100.0, NaN], nothing, [1.0, 2.0], nothing, NamedTuple())

    @test_throws ArgumentError MassSpectrum{Vector{Float64}, Nothing, Vector{Float64}, Nothing, NamedTuple}(
        [-1.0, 100.0], nothing, [1.0, 2.0], nothing, NamedTuple())

    @test_throws ArgumentError MassSpectrum{Vector{Float64}, Nothing, Vector{Float64}, Nothing, NamedTuple}(
        [100.0, 100.0], nothing, [1.0, 2.0], nothing, NamedTuple())

    # Rejections: invalid intensities
    @test_throws ArgumentError MassSpectrum{Vector{Float64}, Nothing, Vector{Float64}, Nothing, NamedTuple}(
        [100.0], nothing, Float64[], nothing, NamedTuple())

    @test_throws ArgumentError MassSpectrum{Vector{Float64}, Nothing, Vector{Float64}, Nothing, NamedTuple}(
        [100.0, 150.0], nothing, [1.0, NaN], nothing, NamedTuple())

    @test_throws DimensionMismatch MassSpectrum{Vector{Float64}, Nothing, Vector{Float64}, Nothing, NamedTuple}(
        [100.0, 150.0], nothing, [1.0], nothing, NamedTuple())
end

@testset "MassSpectrum – outer constructors" begin
    # Unitful m/z and intensity
    ms1 = MassSpectrum([100.0, 150.0]u"Th", [10.0, 20.0]u"pA"; attrs=ATTRS_TUPLE)
    @test ms1.mzvalues == [100.0, 150.0]
    @test ms1.mzunit == u"Th"
    @test ms1.intensities == [10.0, 20.0]
    @test ms1.intensityunit == pA
    @test ms1.attrs == ATTRS_TUPLE

    # Unitful m/z, unitless intensities
    ms2 = MassSpectrum([100.0, 150.0]u"Th", [10.0, 20.0])
    @test ms2.mzunit == u"Th"
    @test ms2.intensityunit === nothing

    # Unitless m/z, unitful intensities
    ms3 = MassSpectrum([100.0, 150.0], [10.0, 20.0]u"pA")
    @test ms3.mzunit === nothing
    @test ms3.intensityunit == pA

    # Unitless m/z and intensity
    ms4 = MassSpectrum([100.0, 150.0], [10.0, 20.0])
    @test ms4.mzunit === nothing
    @test ms4.intensityunit === nothing

    # Explicit unit metadata constructor
    ms5 = MassSpectrum([100.0, 150.0], u"Th", [10.0, 20.0], pA; attrs=(scan_id=7,))
    @test ms5.mzunit == u"Th"
    @test ms5.intensityunit == pA
    @test ms5.attrs == (scan_id=7,)

    # Rejections: inconsistent units
    bad_mz = Unitful.AbstractQuantity{Float64}[100.0u"Th", 1.0u"s"]
    bad_ints = Unitful.AbstractQuantity{Float64}[10.0u"pA", 1.0u"s"]
    @test_throws ArgumentError MassSpectrum(bad_mz, [10.0, 20.0])
    @test_throws ArgumentError MassSpectrum([100.0, 150.0], bad_ints)
end

end  # module TestMassSpectra
