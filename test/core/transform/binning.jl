module TestBinning

using Test
using Unitful
using JuChrom

# ── Helpers ──────────────────────────────────────────────────────────────────
_make_ms(mz::AbstractVector, ints::AbstractVector; rt=1.0u"s", lvl::Int=1) =
    MassScan(rt, mz, ints; level=lvl)

_make_mss(scans::Vector{<:MassScan}; meta...) = MassScanSeries(scans; meta...)

struct _UncheckedVarianceMassScanMatrix{
    R,
    M,
    I,
    T<:AbstractMassScanMatrix{R,M,I},
    A<:AbstractMatrix{<:Real},
    V<:Union{Nothing, Unitful.Units},
} <: AbstractVarianceMassScanMatrix{R,M,I,V}
    msm::T
    variances::A
    varianceunit::V
end

Base.parent(vmsm::_UncheckedVarianceMassScanMatrix) = vmsm.msm

@testset "binretentions(vmsm::AbstractVarianceMassScanMatrix, rgrid, rho_lag1)" begin
    rets = [0.0, 0.5, 1.0] .* u"s"
    mzs = [100.0, 101.0]
    ints = [2 4; 6 8; 10 12] .* u"pA"
    msm = MassScanMatrix(rets, mzs, ints; extras=Dict("kind" => "signal"))
    rgrid = RetentionGrid([0.0, 1.0, 2.0], 1.0, u"s", 1e-8, 0.0, 2.0)
    vars = fill(0.2, size(ints))
    rho = fill(0.0, length(mzs))
    vmsm = VarianceMassScanMatrix(msm, vars)

    binned = binretentions(vmsm, rgrid, rho)
    @test binned isa VarianceMassScanMatrix
    @test parent(binned) isa MassScanMatrix
    @test retentionunit(binned) == u"s"
    @test rawretentions(binned) ≈ [0.5, 1.5]
    @test intensityunit(binned) == u"pA"
    @test varianceunit(binned) == u"pA"^2
    @test JuChrom.rawintensities(binned) ≈ [4 6; 10 12]
    @test rawvariances(binned) ≈ [0.1 0.1; 0.2 0.2]
    @test extras(binned)["kind"] == "signal"

    binned_default = binretentions(vmsm, rgrid)
    @test JuChrom.rawintensities(binned_default) ≈ JuChrom.rawintensities(binned)
    @test rawvariances(binned_default) ≈ rawvariances(binned)

    rgrid_ms = RetentionGrid([0.0, 1000.0, 2000.0], 1000.0, u"ms", 1e-5, 0.0, 2000.0)
    binned_ms = binretentions(vmsm, rgrid_ms, rho)
    @test retentionunit(binned_ms) == u"ms"
    @test rawretentions(binned_ms) ≈ [500.0, 1500.0]
    @test JuChrom.rawintensities(binned_ms) ≈ JuChrom.rawintensities(binned)
    @test rawvariances(binned_ms) ≈ rawvariances(binned)

    vmsm_unitful = VarianceMassScanMatrix(msm, vars .* u"pA"^2)
    binned_zth = binretentions(vmsm_unitful, rgrid, rho;
        zero_threshold=1e-6u"pA"^2)
    @test JuChrom.rawintensities(binned_zth) ≈ JuChrom.rawintensities(binned)
    @test rawvariances(binned_zth) ≈ rawvariances(binned)

    unchecked_vmsm = _UncheckedVarianceMassScanMatrix(msm, vars, nothing)
    binned_attached_vars = binretentions(unchecked_vmsm, rgrid, rho)
    @test intensityunit(binned_attached_vars) == u"pA"
    @test varianceunit(binned_attached_vars) == u"pA"^2
    @test JuChrom.rawintensities(binned_attached_vars) ≈ JuChrom.rawintensities(binned)
    @test rawvariances(binned_attached_vars) ≈ rawvariances(binned)

    _, unitful_empty_avg, unitful_empty_vars = JuChrom._binretention_trace(
        [0.25]u"s",
        [2.0]u"pA",
        [0.0, 1.0, 2.0]u"s",
        [0.5]u"pA"^2;
        zero_threshold=1e-6,
    )
    @test unitful_empty_avg[2] == 0.0u"pA"
    @test unitful_empty_vars[2] == Inf * u"pA"^2

    msm_unitless = MassScanMatrix([0.0, 0.5, 1.0], [100.0],
        reshape([2.0, 6.0, 10.0], 3, 1))
    vmsm_unitless = VarianceMassScanMatrix(msm_unitless, fill(0.2, 3, 1))
    rgrid_unitless = RetentionGrid([0.0, 1.0, 2.0], 1.0, nothing, 1e-8, 0.0, 2.0)
    binned_unitless = binretentions(vmsm_unitless, rgrid_unitless)
    @test retentionunit(binned_unitless) === nothing
    @test rawretentions(binned_unitless) ≈ [0.5, 1.5]
    @test rawintensities(binned_unitless) ≈ reshape([4.0, 10.0], 2, 1)
    @test rawvariances(binned_unitless) ≈ reshape([0.1, 0.2], 2, 1)

    bad_grid = RetentionGrid([0.0, 1.0, 2.0], 1.0, u"Th", 1e-8, 0.0, 2.0)
    @test_throws ArgumentError binretentions(vmsm, rgrid, [0.0])
    @test_throws ArgumentError binretentions(vmsm_unitless, rgrid)
    @test_throws ArgumentError binretentions(vmsm, rgrid_unitless)
    @test_throws Unitful.DimensionError binretentions(vmsm, bad_grid)
    @test_throws Unitful.DimensionError JuChrom._binretention_trace(
        [0.0, 0.5]u"s",
        [1.0, 2.0],
        [0.0, 1.0]u"kg",
        [0.2, 0.2],
    )
    @test_throws ArgumentError JuChrom._binretention_trace(
        [0.25]u"s",
        [2.0]u"pA",
        [0.0, 1.0]u"s",
        [0.5],
    )
    @test_throws MethodError binretentions(vmsm, rawbinedges(rgrid), rho)
    @test_throws MethodError binretentions(vmsm, rgrid, rho;
        jacobian_scale=fill(2.0, length(rets)))
    @test_throws MethodError binretentions(msm, rgrid, rho)
    @test_throws MethodError binretentions([0.0, 0.5], [1.0, 2.0], [0.0, 1.0], [0.1, 0.1])
    @test_throws MethodError binretentions(
        msm,
        rgrid,
        LinearObservedIntensityVarianceModel(0.1, 0.0, 0.0, 12.0, 0.0),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# binmzvalues
# ─────────────────────────────────────────────────────────────────────────────

@testset "binmzvalues(series::AbstractMassScanSeries; ionbin=integer)" begin
    @testset "basic integer binning" begin
        s = _make_ms([100.1, 100.9, 101.2], [10.0, 5.0, 20.0])
        mss = _make_mss([s])
        b = binmzvalues(mss)                 # default integer binning
        sb = scan(b, 1)
        @test JuChrom.mzvalues(sb) == [100.0, 101.0]
        @test JuChrom.intensities(sb) ≈ [10.0, 25.0]
    end

    @testset "custom binning function" begin
        s = _make_ms([100.11, 100.13, 100.85], [3.0, 2.0, 5.0])
        mss = _make_mss([s])
        b = binmzvalues(mss, x -> round(x, digits=1))
        sb = scan(b, 1)
        @test JuChrom.mzvalues(sb) == [100.1, 100.8]
        @test JuChrom.intensities(sb) ≈ [5.0, 5.0]
    end

    @testset "all ions in same bin" begin
        s = _make_ms([101.1, 101.2, 101.3], [1.0, 2.0, 3.0])
        mss = _make_mss([s])
        b = binmzvalues(mss, _ -> 101.0)
        sb = scan(b, 1)
        @test JuChrom.mzvalues(sb) == [101.0]
        @test JuChrom.intensities(sb) ≈ [6.0]
    end

    @testset "output respects intensity eltype" begin
        s = _make_ms([99.5, 100.5], Float32[1.0, 2.0])
        mss = _make_mss([s])
        b = binmzvalues(mss)
        sb = scan(b, 1)
        @test eltype(JuChrom.intensities(sb)) == Float32
    end

    @testset "multi-scan series" begin
        s1 = _make_ms([100.1, 100.9], [5.0, 10.0])
        s2 = _make_ms([101.7, 102.1], [20.0, 30.0])
        mss = _make_mss([s1, s2])
        b = binmzvalues(mss)
        @test scancount(b) == 2

        b1 = scan(b, 1)
        @test JuChrom.mzvalues(b1) == [100.0, 101.0]
        @test JuChrom.intensities(b1) ≈ [5.0, 10.0]

        b2 = scan(b, 2)
        @test JuChrom.mzvalues(b2) == [102.0]
        @test JuChrom.intensities(b2) ≈ [50.0]
    end

    @testset "validmzvalues filtering" begin
        s = _make_ms([28.1, 29.1, 30.2], [1.0, 2.0, 3.0])
        mss = _make_mss([s])
        b = binmzvalues(mss; validmzvalues=29:1:30)
        sb = scan(b, 1)
        @test JuChrom.mzvalues(sb) == [29.0, 30.0]
        @test JuChrom.intensities(sb) ≈ [2.0, 3.0]
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# integer
# ─────────────────────────────────────────────────────────────────────────────

@testset "integer() binning function" begin
    @test integer(200.0) == 200
    @test integer(200.3) == 200
    @test integer(200.69) == 200
    @test integer(200.7) == 201
    @test integer(200.9999999999999) == 201
    @test integer(201.0) == 201
    @test integer(201.0000000000001) == 201
    @test integer(201.6999999999999) == 201
    @test integer(201.7) == 202

    @test integer(199.75, 0.25) == 200
    @test integer(199.74, 0.25) == 199
    @test integer(201.0, 0.0) == 201
    @test integer(201.0, 0.999999999999) == 201

    @test integer(prevfloat(201.7)) == 201
    @test integer(201.7) == 202
    @test integer(prevfloat(201.0)) == 201
    @test integer(201.0) == 201
    @test integer(201.7 - 1e-13) == 201

    @test_throws ArgumentError integer(200.0, -0.1)
    @test_throws ArgumentError integer(200.0, 1.0)
end

end # module TestBinning
