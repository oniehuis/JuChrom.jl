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

        edges_vec, width_vec = densestgrid([msmu1, msmu2])
        @test all(isapprox.(edges_vec, expected_edges; atol=1e-10u"s"))
        @test isapprox(width_vec, 0.5u"s"; atol=1e-10u"s")
    end
end

end  # module TestGridding