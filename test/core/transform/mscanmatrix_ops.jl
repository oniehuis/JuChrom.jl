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
# subtract
# ─────────────────────────────────────────────────────────────────────────────

@testset "subtract" begin
    ret = [0.0, 1.0]
    mz = [100.0, 200.0, 300.0]
    ints_a = [5.0 2.0 1.0; 3.0 0.2 4.0]
    ints_b = [1.0 1.5 2.0; 4.0 0.1 10.0]
    meta = (instrument=(model="X",), acquisition=(mode="scan",),
            user=(name="tester",), sample=(id="S1",), extras=Dict("note" => "ok"))

    msm_a = make_msm(; ret=ret, mz=mz, ints=ints_a, level=2,
                     instrument=meta.instrument, acquisition=meta.acquisition,
                     user=meta.user, sample=meta.sample, extras=meta.extras)
    msm_b = make_msm(; ret=ret, mz=mz, ints=ints_b, level=2,
                     instrument=meta.instrument, acquisition=meta.acquisition,
                     user=meta.user, sample=meta.sample, extras=meta.extras)

    out = JuChrom.subtract(msm_a, msm_b, 0.5)
    @test retentions(out) == ret
    @test retentionunit(out) === nothing
    @test mzvalues(out) == mz
    @test mzunit(out) === nothing
    @test intensityunit(out) === nothing
    @test level(out) == 2
    @test instrument(out) == meta.instrument
    @test acquisition(out) == meta.acquisition
    @test user(out) == meta.user
    @test sample(out) == meta.sample
    @test extras(out) == meta.extras

    expected = ints_a .- ints_b
    expected[expected .< 0.5] .= 0.5
    @test rawintensities(out) == expected

    extras(msm_a)["note"] = "changed"
    @test extras(out)["note"] == "ok"

    msm_bad = make_msm(; ret=[0.0, 2.0], mz=mz, ints=ints_b, level=2)
    @test_throws ArgumentError JuChrom.subtract(msm_a, msm_bad, 0.0)
end

end  # module TestMscanmatrix_ops
