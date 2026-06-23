module TestMscanmatrix_ops

using Test
using Unitful
using JuChrom

# ── Helpers ──────────────────────────────────────────────────────────────────

make_msm(; ret, mz, ints, level=1, instrument=NamedTuple(), acquisition=NamedTuple(),
         user=NamedTuple(), sample=NamedTuple(), extras=Dict{String, Any}()) =
    MassScanMatrix(
        ret,
        mz,
        ints;
        level=level,
        instrument=instrument,
        acquisition=acquisition,
        user=user,
        sample=sample,
        extras=extras,
    )

# ─────────────────────────────────────────────────────────────────────────────
# subtractbaseline
# ─────────────────────────────────────────────────────────────────────────────

@testset "subtractbaseline" begin
    ret = [0.0, 1.0]
    mz = [100.0, 200.0, 300.0]
    ints_a = [5.0 2.0 1.0; 3.0 0.2 4.0]
    baseline_ints = [1.0 1.5 2.0; 4.0 0.1 10.0]
    source_meta = (instrument=(model="X",), acquisition=(mode="scan",),
                   user=(name="tester",), sample=(id="S1",),
                   extras=Dict("note" => "source"))
    baseline_meta = (instrument=(model="Y",), acquisition=(mode="baseline",),
                     user=(name="other",), sample=(id="baseline",),
                     extras=Dict("note" => "baseline"))

    msm_a = make_msm(; ret=ret, mz=mz, ints=ints_a, level=2,
                     instrument=source_meta.instrument, acquisition=source_meta.acquisition,
                     user=source_meta.user, sample=source_meta.sample,
                     extras=source_meta.extras)
    baseline = make_msm(; ret=ret, mz=mz, ints=baseline_ints, level=2,
                        instrument=baseline_meta.instrument,
                        acquisition=baseline_meta.acquisition,
                        user=baseline_meta.user, sample=baseline_meta.sample,
                        extras=baseline_meta.extras)

    out = subtractbaseline(msm_a, baseline)
    @test retentions(out) == ret
    @test retentionunit(out) ≡ nothing
    @test mzvalues(out) == mz
    @test mzunit(out) ≡ nothing
    @test intensityunit(out) ≡ nothing
    @test level(out) == 2
    @test instrument(out) == source_meta.instrument
    @test acquisition(out) == source_meta.acquisition
    @test user(out) == source_meta.user
    @test sample(out) == source_meta.sample
    @test extras(out) == source_meta.extras

    expected = ints_a .- baseline_ints
    @test rawintensities(out) == expected

    extras(msm_a)["note"] = "changed"
    @test extras(out)["note"] == "source"

    variances = [0.2 0.3 0.4; 0.5 0.6 0.7]
    vmsm = VarianceMassScanMatrix(msm_a, variances)
    vout = subtractbaseline(vmsm, baseline)
    @test vout isa VarianceMassScanMatrix
    @test rawintensities(parent(vout)) == expected
    @test rawvariances(vout) == variances
    @test rawvariances(vout) !== rawvariances(vmsm)
    @test varianceunit(vout) == varianceunit(vmsm)
    @test extras(parent(vout)) == extras(msm_a)

    baseline_bad_ret = make_msm(; ret=[0.0, 2.0], mz=mz, ints=baseline_ints, level=2)
    @test_throws DimensionMismatch subtractbaseline(msm_a, baseline_bad_ret)

    baseline_bad_mz = make_msm(; ret=ret, mz=[100.0, 250.0, 300.0],
                               ints=baseline_ints, level=2)
    @test_throws DimensionMismatch subtractbaseline(msm_a, baseline_bad_mz)

    baseline_bad_level = make_msm(; ret=ret, mz=mz, ints=baseline_ints, level=1)
    @test_throws DimensionMismatch subtractbaseline(msm_a, baseline_bad_level)
end

end  # module TestMscanmatrix_ops
