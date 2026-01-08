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
end

end  # module TestRetentionMapper
