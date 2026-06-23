module TestCensoring

using Test
using Unitful
using JuChrom

@testset "replacecensored(vmsm)" begin
    rets = [1.0, 2.0, 3.0]u"s"
    mzs = [100.0, 101.0]
    x = [0.0 10.0;
         0.4  0.03;
         2.0  5.0]
    v = [0.0025 0.01;
         0.0001 0.0001;
         0.01   0.04]
    msm = MassScanMatrix(rets, mzs, x .* u"pA"; extras=Dict("run" => "A"))
    vmsm = VarianceMassScanMatrix(msm, v .* u"pA"^2)

    replaced, info = replacecensored(vmsm; floor=0.2, k=2.0, q=0.0, returninfo=true)

    @test replaced isa VarianceMassScanMatrix
    @test info isa CensoredReplacementInfo
    @test parent(replaced) isa MassScanMatrix
    @test retentionunit(replaced) == u"s"
    @test intensityunit(replaced) == u"pA"
    @test varianceunit(replaced) == u"pA"^2
    @test rawretentions(replaced) == rawretentions(vmsm)
    @test rawmzvalues(replaced) == rawmzvalues(vmsm)
    @test extras(replaced)["run"] == "A"

    expected_x = copy(x)
    expected_x[1, 1] = 0.1
    expected_x[2, 2] = 0.1
    @test rawintensities(replaced) ≈ expected_x

    expected_v = copy(v)
    expected_v[1, 1] = 0.2^2 / 12
    expected_v[2, 2] = 0.2^2 / 12
    @test rawvariances(replaced) ≈ expected_v

    @test rawintensities(vmsm) == x
    @test rawvariances(vmsm) ≈ v

    @test info.k == 2.0
    @test info.q == 0.0
    @test info.floor_rule == :scalar
    @test info.global_floor ≈ 0.03
    @test info.column_floors == [0.2, 0.2]
    @test info.n_replaced == 2
    @test info.fraction_replaced ≈ 2 / length(x)
    @test info.n_nonpositive == 1
    @test info.n_below_limit_positive == 1
    @test info.replaced_by_column == [1, 1]
    @test sprint(show, info) ==
        "CensoredReplacementInfo(2 replaced, fraction=0.3333, floor_rule=:scalar, k=2.0, q=0.0)"
    @test occursin("Replaced cells: 2", sprint(show, MIME"text/plain"(), info))

    @test !haskey(extras(replaced), "replacecensored")

    @test replacecensored(vmsm; floor=0.2, k=2.0) isa VarianceMassScanMatrix
end

@testset "replacecensored floor rules and validation" begin
    rets = [1.0, 2.0, 3.0]u"s"
    mzs = [100.0, 101.0, 102.0]
    x = [0.0  0.0 7.0;
         1.0 10.0 0.0;
         3.0 30.0 0.0]
    v = fill(0.01, size(x))
    msm = MassScanMatrix(rets, mzs, x)
    vmsm = VarianceMassScanMatrix(msm, v)

    _, info = replacecensored(vmsm; k=0.0, q=0.5, minpositives=2, returninfo=true)
    @test info.floor_rule == :column_quantile
    @test info.global_floor ≈ 7.0
    @test info.column_floors ≈ [2.0, 20.0, 7.0]

    _, global_info = replacecensored(
        vmsm; floor=:global_quantile, k=0.0, q=0.5, returninfo=true)
    @test global_info.floor_rule == :global_quantile
    @test global_info.column_floors ≈ fill(7.0, 3)

    _, vector_info = replacecensored(
        vmsm; floor=[0.5, 2.0, 4.0], k=0.0, returninfo=true)
    @test vector_info.floor_rule == :vector
    @test vector_info.column_floors == [0.5, 2.0, 4.0]

    @test_throws ArgumentError replacecensored(vmsm; k=-1.0)
    @test_throws ArgumentError replacecensored(vmsm; q=-0.1)
    @test_throws ArgumentError replacecensored(vmsm; q=1.1)
    @test_throws ArgumentError replacecensored(vmsm; minpositives=0)
    @test_throws ArgumentError replacecensored(vmsm; floor=0.0)
    @test_throws ArgumentError replacecensored(vmsm; floor=[1.0, 2.0])
    @test_throws ArgumentError replacecensored(vmsm; floor=[1.0, 0.0, 2.0])
    @test_throws ArgumentError replacecensored(vmsm; floor=:unknown)

    x_no_positive = fill(0.0, 3, 2)
    v_no_positive = fill(0.01, 3, 2)
    vmsm_no_positive = VarianceMassScanMatrix(
        MassScanMatrix(rets, [100.0, 101.0], x_no_positive),
        v_no_positive,
    )
    @test_throws ArgumentError replacecensored(vmsm_no_positive)
end

end # module TestCensoring
