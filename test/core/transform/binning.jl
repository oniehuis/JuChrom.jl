module TestBinning

using Statistics: mean
using Test
using Unitful
using JuChrom

# ── Helpers ──────────────────────────────────────────────────────────────────
_make_ms(mz::AbstractVector, ints::AbstractVector; rt=1.0u"s", lvl::Int=1) =
    MassScan(rt, mz, ints; level=lvl)

_make_mss(scans::Vector{<:MassScan}; meta...) = MassScanSeries(scans; meta...)

# ─────────────────────────────────────────────────────────────────────────────
# binretentions(retentions, intensities, bin_edges, variances)
# ─────────────────────────────────────────────────────────────────────────────

@testset "binretentions(retentions, intensities, bin_edges, variances)" begin
    @testset "returns correct centers, means, and baseline variances" begin
        rets = collect(0.0:0.5:2.5)
        ints = [10.0, 12.0, 8.0, 4.0, 6.0, 0.0]
        edges = collect(0.0:1.0:3.0)
        vars = fill(0.25, length(ints))

        centers, means, varŝ = binretentions(rets, ints, edges, vars)
        @test centers == [0.5, 1.5, 2.5]
        @test means ≈ [11.0, 6.0, 3.0]
        @test varŝ ≈ fill(0.125, 3)
    end

    @testset "inverse-variance weighting and zero-threshold clamping" begin
        rets = [0.0, 0.4, 0.8]
        ints = [1.0, 9.0, 5.0]
        edges = [0.0, 0.5, 1.0]
        vars = [0.01, 1.0, 0.25]

        centers, means, varŝ = binretentions(rets, ints, edges, vars; zero_threshold=0.02)
        w1 = 1 / 0.02
        w2 = 1 / 1.0
        expected_bin1 = (w1 * ints[1] + w2 * ints[2]) / (w1 + w2)

        @test centers == [0.25, 0.75]
        @test means[1] ≈ expected_bin1
        @test means[2] == ints[3]
        @test varŝ[1] ≈ 1 / (w1 + w2)
        @test varŝ[2] ≈ 0.25
    end

    @testset "rho inflation and negative intensities" begin
        rets = [0.0, 0.1, 0.2]
        ints = [2.0, -1.0, 3.0]
        edges = [0.0, 0.5, 1.0]
        vars = fill(0.5, length(ints))

        centers, means, varŝ = binretentions(
            rets, ints, edges, vars;
            rho_lag1=0.5, rho_max=0.8)

        wsum = sum(fill(1 / 0.5, 3))
        vif = JuChrom.vif(0.5, 3; ρmax=0.8, nonnegative=true, nmin=1)

        @test centers == [0.25, 0.75]
        @test means[1] ≈ mean(ints)
        @test varŝ[1] ≈ (1 / wsum) * vif
        @test isinf(varŝ[2])
    end

    @testset "unit compatibility" begin
        rets = [0.0, 0.5, 1.0] .* u"s"
        edges = (0.0:1.0:2.0) .* u"s"
        ints = [5.0, 7.0, 9.0]
        vars = fill(0.25, length(ints))

        centers, means, _ = binretentions(rets, ints, edges, vars)
        @test all(JuChrom.isunitful, centers)
        @test means[1] ≈ 6.0
        @test_throws ArgumentError binretentions(rets, ints, collect(0.0:1.0:2.0), vars)
    end

    @testset "unitful variances and empty bins" begin
        rets = [0.10, 0.20, 0.35] .* u"s"
        edges = [0.0, 0.3, 0.5, 0.8] .* u"s"
        ints = [10.0, 20.0, 5.0] .* u"pA"
        vars = [0.04, 0.16, 0.25] .* (u"pA"^2)

        centers, means, varŝ = binretentions(
            rets, ints, edges, vars; zero_threshold=1e-9)

        w1 = inv(vars[1])
        w2 = inv(vars[2])
        expected_mean = (w1 * ints[1] + w2 * ints[2]) / (w1 + w2)
        expected_var = inv(w1 + w2)

        @test all(JuChrom.isunitful, means)
        @test all(JuChrom.isunitful, varŝ)
        @test means[1] ≈ expected_mean
        @test varŝ[1] ≈ expected_var
        @test means[2] == ints[3]
        @test varŝ[2] == vars[3]
        @test isinf(varŝ[3])
        @test Unitful.unit(varŝ[3]) == Unitful.unit(vars[1])
    end

    @testset "input validation" begin
        rets = [0.0, 1.0]
        ints = [1.0, 2.0]
        vars = [0.1, 0.1]

        @test_throws ArgumentError binretentions(Float64[], ints, [0.0, 1.0], vars)
        @test_throws ArgumentError binretentions(rets, [1.0], [0.0, 1.0], [0.1])
        @test_throws ArgumentError binretentions(rets, ints, [0.0], vars)
        @test_throws ArgumentError binretentions(rets, ints, [0.0, 1.0], [-0.1, 0.2])

        ints_u = [1.0, 2.0] .* u"pA"
        vars_bad_dim = fill(0.1, 2) .* u"pA"
        @test_throws Unitful.DimensionError binretentions(rets .* u"s", ints_u, 
            [0.0, 1.0, 2.0] .* u"s", vars_bad_dim)

        vars_unitful = fill(0.1, 2) .* (u"pA"^2)
        @test_throws ArgumentError binretentions(rets, ints, [0.0, 1.0, 2.0], vars_unitful)

        rets_u = rets .* u"s"
        ints_u2 = ints .* u"pA"
        vars_unitless = [0.1, 0.2]
        @test_throws ArgumentError binretentions(rets_u, ints_u2, [0.0, 1.0, 2.0] .* u"s", vars_unitless)
    end
end

@testset "binretentions(msm::MassScanMatrix, bin_edges, quadvar_params, rho_lag1)" begin
    rets = [0.0, 0.5, 1.0] .* u"s"
    mzs = [100.0, 101.0]
    ints = [2 4; 6 8; 10 12] .* u"pA"
    msm = MassScanMatrix(rets, mzs, ints)
    edges = [0.0, 1.0, 2.0] .* u"s"
    params = fill(QuadVarParams(0.1, 0.0, 0.0), length(mzs))
    rho = fill(0.0, length(mzs))

    msm_binned, vars = binretentions(msm, edges, params, rho)
    @test intensityunit(msm_binned) == u"pA"
    @test JuChrom.rawintensities(msm_binned) ≈ [4 6; 10 12]
    @test all(JuChrom.isunitful, vars[:, 1])
    @test Unitful.unit(vars[1, 1]) == u"pA"^2
    @test vars[1, 1] ≈ 0.05u"pA"^2
    @test vars[2, 2] ≈ 0.1u"pA"^2

    # Keep the linear term explicit in the intensity units to avoid DimensionError.
    params_unitful = fill(QuadVarParams(0.2u"pA^2", 0.0u"pA", 0.0), length(mzs))
    msm_binned_u, vars_u = binretentions(msm, edges, params_unitful, rho)
    @test JuChrom.rawintensities(msm_binned_u) ≈ JuChrom.rawintensities(msm_binned)
    @test all(Unitful.unit.(vars_u[:, 1]) .== u"pA"^2)
    @test vars_u[1, 1] ≈ 0.1u"pA"^2

    msm_unitless = MassScanMatrix([0.0, 0.5], [100.0], reshape([1.0, 2.0], 2, 1))
    @test_throws ArgumentError binretentions(msm_unitless, [0.0, 1.0], 
        QuadVarParams(0.1, 0.0, 0.0), 0.0; zero_threshold=1e-6u"pA"^2)
    msm_unitless_ok = MassScanMatrix([0.0, 0.5, 1.0], [100.0],
        reshape([1.0, 2.0, 3.0], 3, 1))
    msm_binned_unitless, vars_unitless = binretentions(
        msm_unitless_ok,
        [0.0, 0.5, 1.0],
        QuadVarParams(0.1, 0.0, 0.0),
        0.0;
        jacobian_scale=_ -> 2.0,
    )
    @test intensityunit(msm_binned_unitless) === nothing
    @test all(!JuChrom.isunitful, vars_unitless[:, 1])

    @test_throws ArgumentError binretentions(msm, edges, params, rho;
        jacobian_scale=[1.0, 2.0])
    @test_throws ArgumentError binretentions(msm, edges, params, rho;
        jacobian_scale=[1.0, 0.0, 2.0])
    @test_throws ArgumentError binretentions(msm, edges, params, rho;
        jacobian_scale=t -> (t == 0.5u"s" ? 0.0 : 1.0))

    @test_throws Unitful.DimensionError QuadVarParams(0.1u"pA^2", 0.0, 0.0)

    @test_throws Unitful.DimensionError binretentions(msm, edges, params, rho;
        zero_threshold=1e-6u"s")
end

@testset "binretentions(msm::MassScanMatrix, bin_edges, variances, rho_lag1)" begin
    rets = [0.0, 0.5, 1.0] .* u"s"
    mzs = [100.0, 101.0]
    ints = [2 4; 6 8; 10 12] .* u"pA"
    msm = MassScanMatrix(rets, mzs, ints)
    edges = [0.0, 1.0, 2.0] .* u"s"
    vars = fill(0.2, size(ints))
    rho = fill(0.0, length(mzs))

    msm_binned, vars_b = binretentions(msm, edges, vars, rho)
    @test intensityunit(msm_binned) == u"pA"
    @test JuChrom.rawintensities(msm_binned) ≈ [4 6; 10 12]
    @test Unitful.unit(vars_b[1, 1]) == u"pA"^2
    @test vars_b[1, 1] ≈ 0.1u"pA"^2
    @test vars_b[2, 2] ≈ 0.2u"pA"^2

    msm_unitless = MassScanMatrix([0.0, 0.5], [100.0], reshape([1.0, 2.0], 2, 1))
    vars_unitful = fill(0.2u"pA"^2, 2, 1)
    @test_throws ArgumentError binretentions(msm_unitless, [0.0, 1.0], vars_unitful, 0.0)
    @test_throws ArgumentError binretentions(msm, edges, vars[1:2, :], rho)
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
