module TestRetentionMapper

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./ext/jld2_extensions/retentionmapper.jl
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

using Test
using JLD2
using Unitful

# ── Internal Project Imports ─────────────────────────────────────────────────

using JuChrom
using JuChrom: retentions_A, retentions_B, rawmapmin, rawmapmax

# ── Roundtrip serialization with JLD2 ────────────────────────────────────────

function roundtrip(rmapper)
    path, io = mktemp()
    close(io)
    path_jld2 = string(path, ".jld2")
    try
        JLD2.jldsave(path_jld2; rmapper)
        JLD2.load(path_jld2, "rmapper")
    finally
        isfile(path_jld2) && rm(path_jld2; force=true)
    end
end

@testset "RetentionMapper JLD2 roundtrip (unitless)" begin
    rmapper = fitmap([1.0, 2.0, 4.0], [10.0, 20.0, 40.0])
    rmapper2 = roundtrip(rmapper)

    @test rmapper2 isa JuChrom.RetentionMapper
    @test retentions_A(rmapper2) ≈ retentions_A(rmapper)
    @test retentions_B(rmapper2) ≈ retentions_B(rmapper)
    @test rawmapmin(rmapper2) ≈ rawmapmin(rmapper)
    @test rawmapmax(rmapper2) ≈ rawmapmax(rmapper)
    @test rmapper2.lambda ≈ rmapper.lambda
end

@testset "RetentionMapper JLD2 roundtrip (unitful)" begin
    rA = [1.0, 2.0, 4.0]u"s"
    rB = [10.0, 20.0, 40.0]u"s"
    rmapper = fitmap(rA, rB)
    rmapper2 = roundtrip(rmapper)

    @test rmapper2 isa JuChrom.RetentionMapper
    @test retentions_A(rmapper2) ≈ retentions_A(rmapper)
    @test retentions_B(rmapper2) ≈ retentions_B(rmapper)
    @test rawmapmin(rmapper2; unit=u"s") ≈ rawmapmin(rmapper; unit=u"s")
    @test rawmapmax(rmapper2; unit=u"s") ≈ rawmapmax(rmapper; unit=u"s")
    @test rmapper2.lambda ≈ rmapper.lambda
end

@testset "RetentionMapper JLD2 rconvert V1" begin
    ext = Base.get_extension(JuChrom, :JLD2Extension)
    @test ext !== nothing

    rmapper = fitmap([1.0, 2.0, 4.0], [10.0, 20.0, 40.0])
    T1 = typeof(rmapper.rA)
    T2 = Nothing
    T3 = typeof(rmapper.rA_min)
    T4 = typeof(rmapper.rA_max)
    T5 = typeof(rmapper.rA_norm_min)
    T6 = typeof(rmapper.rA_norm_max)
    T7 = typeof(rmapper.rB)
    T8 = Nothing
    T9 = typeof(rmapper.rB_min)
    T10 = typeof(rmapper.rB_max)
    T11 = typeof(rmapper.rB_pred_min)
    T12 = typeof(rmapper.rB_pred_max)
    T13 = typeof(rmapper.rB_pred_norm_min)
    T14 = typeof(rmapper.rB_pred_norm_max)
    T15 = typeof(rmapper.knots)
    T16 = typeof(rmapper.coefs)

    rmapper_v1 = ext.RetentionMapperJLD2V1{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16}(
        rmapper.rA,
        rmapper.rA_unit,
        rmapper.rA_min,
        rmapper.rA_max,
        rmapper.rA_norm_min,
        rmapper.rA_norm_max,
        rmapper.rB,
        rmapper.rB_unit,
        rmapper.rB_min,
        rmapper.rB_max,
        rmapper.rB_pred_min,
        rmapper.rB_pred_max,
        rmapper.rB_pred_norm_min,
        rmapper.rB_pred_norm_max,
        rmapper.knots,
        rmapper.coefs,
        rmapper.extras
    )

    rmapper2 = JLD2.rconvert(typeof(rmapper), rmapper_v1)

    @test rmapper2 isa JuChrom.RetentionMapper
    @test retentions_A(rmapper2) ≈ retentions_A(rmapper)
    @test retentions_B(rmapper2) ≈ retentions_B(rmapper)
    @test rawmapmin(rmapper2) ≈ rawmapmin(rmapper)
    @test rawmapmax(rmapper2) ≈ rawmapmax(rmapper)
    @test rmapper2.knots == rmapper.knots
    @test rmapper2.coefs == rmapper.coefs
    @test rmapper2.lambda === nothing
    @test rmapper2.extras == rmapper.extras
end

end  # module TestRetentionMapper
