module TestQuadVarFitMakie

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./ext/makie_extensions/quadvarfit.jl
# ─────────────────────────────────────────────────────────────────────────────

using Test
using CairoMakie
using Makie
using JuChrom
using JuChrom: QuadVarFit, QuadVarParams
using Unitful
# Force-load the Makie extension so QuadVarFit recipes are available
let ext = Base.get_extension(JuChrom, :MakieExtension)
    if ext === nothing
        Base.require_extension(JuChrom, :MakieExtension)
        ext = Base.get_extension(JuChrom, :MakieExtension)
    end
    ext === nothing && error("JuChrom MakieExtension could not be loaded; ensure Makie is available.")
end

CairoMakie.activate!()  # headless backend for tests

function make_dummy_qvf()
    mz_ref = [100.0]
    mz_unit = nothing
    mz_idx = [1]
    mz_values = [100.0]
    batchcount = 1
    n_reps_per_batch = [1]
    n_scans_per_batch = [3]
    params = [QuadVarParams(0.1, 0.01, 0.0)]
    intensity_unit = nothing
    signal = [[ [1.0, 2.0, 3.0] ]]
    offsets = [[ fill(0.1, 3, 1) ]]
    gains = [[ [1.0] ]]
    scale_c = [1.0]
    acf = [0.0]
    acf_lag = [1]
    n_pairs = [3]
    qc_z_rms = [0.0]
    qc_cov68 = [0.0]
    qc_cov95 = [0.0]
    qc_s_min = [0.0]
    qc_s_max = [0.0]
    qc_nz = [1]
    observed = [[ [1.0; 2.0; 3.0] ] .|> x -> reshape(x, 3, 1)]

    QuadVarFit(
        mz_unit,
        mz_ref,
        mz_idx,
        mz_values,
        batchcount,
        n_reps_per_batch,
        n_scans_per_batch,
        params,
        intensity_unit,
        signal,
        offsets,
        gains,
        scale_c,
        acf,
        acf_lag,
        n_pairs,
        qc_z_rms,
        qc_cov68,
        qc_cov95,
        qc_s_min,
        qc_s_max,
        qc_nz,
        observed
    )
end

@testset "QuadVarFit Makie plot" begin
    qvf = make_dummy_qvf()
    fig = Makie.plot(qvf; mzi=1, batch=1, mode=:align, figsize=(300, 200), legend=false)
    @test fig isa Makie.Figure
    axes = [c for c in fig.content if c isa Makie.Axis]
    @test length(axes) == 1
    @test length([p for p in axes[1].scene.plots if p isa Makie.Lines]) >= 1

    fig_offsets = Makie.plot(qvf; mzi=1, batch=1, mode=:offsets, figsize=(300, 200), legend=false)
    axes_off = [c for c in fig_offsets.content if c isa Makie.Axis]
    @test length(axes_off) == 1
    @test length([p for p in axes_off[1].scene.plots if p isa Makie.Lines]) >= 1

    # unitful m/z title renders unit
    mz_unit = u"Th"
    qvf_unit = QuadVarFit(
        mz_unit,
        [100.0, 101.0],
        [1],
        [100.0],
        1,
        [1],
        [3],
        [QuadVarParams(0.1, 0.01, 0.0)],
        nothing,
        [[ [1.0, 2.0, 3.0] ]],
        [[ fill(0.1, 3, 1) ]],
        [[ [1.0] ]],
        [1.0],
        [0.0],
        [1],
        [3],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [3],
        [[ reshape([1.0, 2.0, 3.0], 3, 1) ]]
    )

    fig_unit = Makie.plot(qvf_unit; mzi=1, batch=1, mode=:align, figsize=(300, 200), legend=false)
    ax_unit = only([c for c in fig_unit.content if c isa Makie.Axis])
    @test occursin("Th", ax_unit.title[])

    # observed uses full grid columns
    qvf_full_obs = QuadVarFit(
        mz_unit,
        [100.0, 101.0],
        [2],
        [101.0],
        1,
        [2],
        [3],
        [QuadVarParams(0.1, 0.01, 0.0)],
        nothing,
        [[ [1.0, 2.0, 3.0] ]],
        [[ fill(0.1, 3, 2) ]],
        [[ [1.0, 1.1] ]],
        [1.0],
        [0.0],
        [1],
        [3],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [3],
        [[ reshape([1.0, 2.0, 3.0, 1.1, 2.1, 3.1], 3, 2),
           reshape([1.5, 2.5, 3.5, 1.6, 2.6, 3.6], 3, 2) ]]
    )
    @test Makie.plot(qvf_full_obs; mzi=1, batch=1, mode=:align, figsize=(300, 200), legend=false) isa Makie.Figure

    # guardrails: invalid indices and observed column mismatch
    @test_throws ArgumentError Makie.plot(qvf; mzi=2, batch=1, mode=:align, figsize=(300, 200), legend=false)
    @test_throws ArgumentError Makie.plot(qvf; mzi=1, batch=2, mode=:align, figsize=(300, 200), legend=false)

    qvf_bad_obs = make_dummy_qvf()
    qvf_bad_obs = QuadVarFit(
        qvf_bad_obs.mz_unit,
        qvf_bad_obs.mz_ref,
        qvf_bad_obs.mz_idx,
        qvf_bad_obs.mz_values,
        qvf_bad_obs.batchcount,
        qvf_bad_obs.n_reps_per_batch,
        qvf_bad_obs.n_scans_per_batch,
        qvf_bad_obs.params,
        qvf_bad_obs.intensity_unit,
        qvf_bad_obs.signal,
        qvf_bad_obs.offsets,
        qvf_bad_obs.gains,
        qvf_bad_obs.scale_c,
        qvf_bad_obs.acf,
        qvf_bad_obs.acf_lag,
        qvf_bad_obs.n_pairs,
        qvf_bad_obs.qc_z_rms,
        qvf_bad_obs.qc_cov68,
        qvf_bad_obs.qc_cov95,
        qvf_bad_obs.qc_s_min,
        qvf_bad_obs.qc_s_max,
        qvf_bad_obs.qc_nz,
        [[fill(0.0, 3, 2)]]  # wrong column count vs selection
    )
    @test_throws ArgumentError Makie.plot(qvf_bad_obs; mzi=1, batch=1, mode=:align, figsize=(300, 200), legend=false)
end

end  # module
