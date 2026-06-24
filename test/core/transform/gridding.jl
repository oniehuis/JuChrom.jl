module TestGridding

using Test
using Unitful
using JuChrom

# ─────────────────────────────────────────────────────────────────────────────
# densestgrid
# ─────────────────────────────────────────────────────────────────────────────

@testset "densestgrid" begin
    @testset "numeric core" begin
        runs = [[0.0, 0.5, 1.0],
                [0.0, 0.5, 1.0]]

        @test_throws ArgumentError densestgrid(runs)

        grid = densestgrid(runs; retentionunit=nothing)
        @test grid isa RetentionGrid
        @test retentionunit(grid) ≡ nothing
        @test rawbinedges(grid) ≈ [0.0, 0.5, 1.0]
        @test rawbinwidth(grid) ≈ 0.5
        @test rawtolerance(grid) ≈ 1e-8
        @test rawoverlapmin(grid) ≈ 0.0
        @test rawoverlapmax(grid) ≈ 1.0
        @test binedges(grid) ≈ [0.0, 0.5, 1.0]
        @test binwidth(grid) ≈ 0.5
        @test overlapmin(grid) ≈ 0.0
        @test overlapmax(grid) ≈ 1.0
        @test_throws ArgumentError binedges(grid; unit=u"s")
        @test identity.(grid) === grid

        edges, width = grid
        @test isapprox(width, 0.5; atol=1e-10)
        @test edges ≈ [0.0, 0.5, 1.0]

        grid_s = densestgrid(runs; retentionunit=u"s")
        @test retentionunit(grid_s) == u"s"
        @test rawbinedges(grid_s) ≈ [0.0, 0.5, 1.0]
        @test binedges(grid_s) ≈ [0.0, 0.5, 1.0] .* u"s"
        @test binwidth(grid_s; unit=u"ms") ≈ 500.0u"ms"
        @test occursin("RetentionGrid(2 bins", sprint(show, grid_s))
        @test occursin("width=0.5 s", sprint(show, grid_s))
        @test occursin("unit=s", sprint(show, grid_s))

        display_text = sprint(show, MIME"text/plain"(), grid_s)
        @test occursin("RetentionGrid with 2 bins", display_text)
        @test occursin("Retention unit: s", display_text)
        @test occursin("Bin edges: 0.0 s to 1.0 s (3 edges)", display_text)
        @test occursin("Tolerance: 1.0e-8 s", display_text)
        @test occursin("Raw edge type: Float64", display_text)
        @test occursin("Extras: none", display_text)

        grid_extra = RetentionGrid(
            [0.0, 0.5, 1.0],
            0.5,
            nothing,
            1e-8,
            0.0,
            1.0;
            extras=Dict("source" => "test"),
        )
        @test occursin("unit=unitless", sprint(show, grid_extra))
        @test occursin("Extras: 1 entry", sprint(show, MIME"text/plain"(), grid_extra))

        shifted = [[0.0, 1.0, 2.0, 3.0],
                   [0.5, 1.5, 2.5, 3.5]]
        edges_shift, width_shift = densestgrid(shifted; retentionunit=nothing)
        @test isapprox(width_shift, 1.0; atol=1e-8)
        @test first(edges_shift) ≈ 0.5 atol=1e-8

        _, constrained = densestgrid(runs; retentionunit=nothing, minwidth=0.75, maxwidth=1.5)
        @test constrained ≥ 0.75
        @test constrained ≤ 1.5

        _, constrained_s = densestgrid(
            runs;
            retentionunit=u"s",
            minwidth=750u"ms",
            maxwidth=1500u"ms",
        )
        @test constrained_s ≥ 0.75u"s"
        @test constrained_s ≤ 1.5u"s"

        @test_throws ArgumentError densestgrid(runs; retentionunit=nothing, maxwidth=0.1)
        @test_throws ArgumentError densestgrid(Vector{Vector{Float64}}(); retentionunit=nothing)
        @test_throws ArgumentError densestgrid([[0.0, 0.1], [1.0, 1.1]]; retentionunit=nothing)
        @test_throws ArgumentError densestgrid(runs; retentionunit=nothing, minwidth=0.1u"s")

        runs_refine = [[0.0, 2.0, 4.0],
                       [1.0, 3.0, 5.0]]
        edges_refine, width_refine = densestgrid(
            runs_refine;
            retentionunit=nothing,
            coarse_inflation=2.0,
            primary_refine_iters=2,
            secondary_refine_iters=0,
        )
        @test isapprox(width_refine, 2.0; atol=1e-12)
        @test edges_refine == [1.0, 3.0, 5.0]

        edges_refine_ok, width_refine_ok = densestgrid(
            runs_refine;
            retentionunit=nothing,
            coarse_inflation=10.0,
            primary_refine_iters=1,
            secondary_refine_iters=0,
        )
        @test width_refine_ok < 10.0
        @test width_refine_ok > 1.0
        @test length(edges_refine_ok) == 3
        @test edges_refine_ok[1] == 1.0
    end

    @testset "unitful vectors" begin
        runs_u = [[0.0, 0.5, 1.0] .* u"s",
                  [0.0, 0.5, 1.0] .* u"s"]

        grid_u = densestgrid(runs_u)
        @test grid_u isa RetentionGrid
        @test retentionunit(grid_u) == u"s"
        @test rawbinedges(grid_u) ≈ [0.0, 0.5, 1.0]
        @test rawbinwidth(grid_u) ≈ 0.5
        @test binedges(grid_u) ≈ [0.0, 0.5, 1.0] .* u"s"
        @test binedges(grid_u; unit=u"ms") ≈ [0.0, 500.0, 1000.0] .* u"ms"
        @test binwidth(grid_u; unit=u"ms") ≈ 500.0u"ms"
        @test overlapmax(grid_u; unit=u"ms") ≈ 1000.0u"ms"

        grid_u_ms = densestgrid(runs_u; retentionunit=u"ms")
        @test retentionunit(grid_u_ms) == u"ms"
        @test rawbinedges(grid_u_ms) ≈ [0.0, 500.0, 1000.0]
        @test binwidth(grid_u_ms) ≈ 500.0u"ms"
        @test_throws ArgumentError densestgrid(runs_u; retentionunit=nothing)

        edges_u, width_u = grid_u
        expected_edges = [0.0, 0.5, 1.0] .* u"s"
        @test all(isapprox.(edges_u, expected_edges; atol=1e-10u"s"))
        @test isapprox(width_u, 0.5u"s"; atol=1e-10u"s")
    end

    @testset "MassScanMatrix wrappers" begin
        rets = [0.0, 0.5, 1.0]
        mzs = [100.0]
        ints1 = reshape([1.0, 2.0, 3.0], 3, 1)
        ints2 = reshape([2.0, 3.0, 4.0], 3, 1)
        msm1 = MassScanMatrix(rets, mzs, ints1)
        msm2 = MassScanMatrix(rets, mzs, ints2)

        edges_dict, width_dict = densestgrid(Dict(:a => msm1, :b => msm2))
        @test edges_dict ≈ [0.0, 0.5, 1.0]
        @test isapprox(width_dict, 0.5; atol=1e-10)

        @test_throws ArgumentError densestgrid(Dict(:a => msm1, :b => msm2); minwidth=0.1u"s")
        @test_throws ArgumentError densestgrid(Dict(:a => msm1, :b => msm2); maxwidth=0.1u"s")
        @test_throws ArgumentError densestgrid(Dict(:a => msm1, :b => msm2); tolerance=1e-8u"s")

        rets_u = [0.0, 0.5, 1.0] .* u"s"
        expected_edges = [0.0, 0.5, 1.0] .* u"s"
        mzs2 = [100.0, 101.0]
        ints_u1 = [1 2;
                   3 4;
                   5 6] .* u"pA"
        ints_u2 = [2 3;
                   4 5;
                   6 7] .* u"pA"
        msmu1 = MassScanMatrix(rets_u, mzs2, ints_u1)
        msmu2 = MassScanMatrix(rets_u, mzs2, ints_u2)

        @test_throws ArgumentError densestgrid(Dict(:a => msm1, :b => msmu1))

        edges_vec, width_vec = densestgrid([msmu1, msmu2]; tolerance=1e-8u"s")
        @test all(isapprox.(edges_vec, expected_edges; atol=1e-10u"s"))
        @test isapprox(width_vec, 0.5u"s"; atol=1e-10u"s")

        edges_dict_u, width_dict_u = densestgrid(
            Dict(:a => msmu1, :b => msmu2);
            minwidth=0.1u"s",
            maxwidth=2.0u"s",
        )
        @test width_dict_u ≥ 0.1u"s"
        @test width_dict_u ≤ 2.0u"s"
        @test edges_dict_u[1] ≈ 0.0u"s"
        @test edges_dict_u[end] ≥ 1.0u"s"
        @test issorted(edges_dict_u)

        @test_throws ArgumentError densestgrid(Dict(:a => msmu1, :b => msmu2); minwidth=0.1)
        @test_throws ArgumentError densestgrid(Dict(:a => msmu1, :b => msmu2); maxwidth=0.1)
        @test_throws ArgumentError densestgrid(Dict(:a => msmu1, :b => msmu2); tolerance=1e-8)

        msmu_bad = MassScanMatrix([0.0, 0.5, 1.0] .* u"Th", mzs2, ints_u1)
        @test_throws ArgumentError densestgrid(Dict(:a => msmu1, :b => msmu_bad))

        edges_vec_unitless, width_vec_unitless = densestgrid(
            [msm1, msm2];
            minwidth=0.1,
            maxwidth=2.0,
            tolerance=1e-9,
        )
        @test width_vec_unitless ≥ 0.1
        @test width_vec_unitless ≤ 2.0
        @test edges_vec_unitless[1] ≈ 0.0

        @test_throws ArgumentError densestgrid([msm1, msm2]; minwidth=0.1u"s")
        @test_throws ArgumentError densestgrid([msm1, msm2]; maxwidth=0.1u"s")
        @test_throws ArgumentError densestgrid([msm1, msm2]; tolerance=1e-8u"s")
        @test_throws ArgumentError densestgrid([msm1, msmu1])
    end
end

end  # module TestGridding
