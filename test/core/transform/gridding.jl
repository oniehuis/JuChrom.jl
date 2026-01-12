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

        edges, width = densestgrid(runs)
        @test isapprox(width, 0.5; atol=1e-10)
        @test edges ≈ [0.0, 0.5, 1.0]

        shifted = [[0.0, 1.0, 2.0, 3.0],
                   [0.5, 1.5, 2.5, 3.5]]
        edges_shift, width_shift = densestgrid(shifted)
        @test isapprox(width_shift, 1.0; atol=1e-8)
        @test first(edges_shift) ≈ 0.5 atol=1e-8

        _, constrained = densestgrid(runs; minwidth=0.75, maxwidth=1.5)
        @test constrained ≥ 0.75
        @test constrained ≤ 1.5

        @test_throws ArgumentError densestgrid(runs; maxwidth=0.1)
        @test_throws ArgumentError densestgrid(Vector{Vector{Float64}}())
        @test_throws ArgumentError densestgrid([[0.0, 0.1], [1.0, 1.1]])

        runs_refine = [[0.0, 2.0, 4.0],
                       [1.0, 3.0, 5.0]]
        edges_refine, width_refine = densestgrid(
            runs_refine;
            coarse_inflation=2.0,
            primary_refine_iters=2,
            secondary_refine_iters=0,
        )
        @test isapprox(width_refine, 2.0; atol=1e-12)
        @test edges_refine == [1.0, 3.0, 5.0]

        edges_refine_ok, width_refine_ok = densestgrid(
            runs_refine;
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

        edges_u, width_u = densestgrid(runs_u)
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
