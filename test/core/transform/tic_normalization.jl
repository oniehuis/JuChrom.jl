module TestTicNormalization

using Test
using Unitful
using JuChrom

@testset "ticnormalize(msm::MassScanMatrix)" begin
    rets = [1.0, 2.0, 3.0]u"s"
    mzs = [100.0, 101.0]u"Th"
    x = [1.0 2.0;
         -1.0 3.0;
         -5.0 1.0]
    msm = MassScanMatrix(
        rets,
        mzs,
        x .* u"pA";
        level=2,
        instrument=(name="test instrument",),
        acquisition=(mode="scan",),
        user=(analyst="tester",),
        sample=(id="sample-a",),
        extras=Dict("run" => "A")
    )

    normalized = ticnormalize(msm)
    scale = 5.0

    @test normalized isa MassScanMatrix
    @test retentionunit(normalized) == u"s"
    @test mzunit(normalized) == u"Th"
    @test intensityunit(normalized) === nothing
    @test level(normalized) == 2
    @test instrument(normalized) == instrument(msm)
    @test acquisition(normalized) == acquisition(msm)
    @test user(normalized) == user(msm)
    @test sample(normalized) == sample(msm)
    @test extras(normalized) == extras(msm)
    @test rawretentions(normalized) == rawretentions(msm)
    @test rawmzvalues(normalized) == rawmzvalues(msm)
    @test rawintensities(normalized) ≈ x ./ scale
    @test sum(filter(>(0), vec(sum(rawintensities(normalized); dims=2)))) ≈ 1.0

    @test rawintensities(msm) == x
end

@testset "ticnormalize(vmsm::AbstractVarianceMassScanMatrix)" begin
    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]u"Th"
    x = [1.0 2.0;
         4.0 8.0]
    v = [0.01 0.04;
         0.16 0.64]
    msm = MassScanMatrix(
        rets,
        mzs,
        x .* u"pA";
        level=2,
        instrument=(name="test instrument",),
        acquisition=(mode="scan",),
        user=(analyst="tester",),
        sample=(id="sample-a",),
        extras=Dict("run" => "A")
    )
    vmsm = VarianceMassScanMatrix(msm, uconvert.(u"nA"^2, v .* u"pA"^2))

    normalized = ticnormalize(vmsm)
    scale = 15.0

    @test normalized isa VarianceMassScanMatrix
    @test parent(normalized) isa MassScanMatrix
    @test retentionunit(normalized) == u"s"
    @test mzunit(normalized) == u"Th"
    @test intensityunit(normalized) === nothing
    @test varianceunit(normalized) === nothing
    @test level(normalized) == 2
    @test instrument(normalized) == instrument(vmsm)
    @test acquisition(normalized) == acquisition(vmsm)
    @test user(normalized) == user(vmsm)
    @test sample(normalized) == sample(vmsm)
    @test extras(normalized) == extras(vmsm)
    @test rawretentions(normalized) == rawretentions(vmsm)
    @test rawmzvalues(normalized) == rawmzvalues(vmsm)
    @test rawintensities(normalized) ≈ x ./ scale
    @test rawvariances(normalized) ≈ v ./ scale^2

    @test rawintensities(vmsm) == x
    @test rawvariances(vmsm; unit=u"pA"^2) ≈ v
end

@testset "ticnormalize preserves CLR results for positive VMSM data" begin
    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]
    x = [1.0 2.0;
         4.0 8.0]
    v = [0.01 0.04;
         0.16 0.64]
    vmsm = VarianceMassScanMatrix(
        MassScanMatrix(rets, mzs, x .* u"pA"),
        v .* u"pA"^2
    )

    direct = clr(vmsm)
    normalized = clr(ticnormalize(vmsm))

    @test rawintensities(normalized) ≈ rawintensities(direct)
    @test rawvariances(normalized) ≈ rawvariances(direct)
end

@testset "ticnormalize validation and scope" begin
    @test_throws ArgumentError ticnormalize(MassScanMatrix(
        [1.0, 2.0],
        [100.0, 101.0],
        [-1.0 0.0; 2.0 -3.0]
    ))

    vmsm = VarianceMassScanMatrix(
        MassScanMatrix([1.0, 2.0], [100.0], reshape([-1.0, -2.0], 2, 1)),
        fill(0.1, 2, 1)
    )
    @test_throws ArgumentError ticnormalize(vmsm)

    mss = MassScanSeries([
        MassScan(1.0u"s", [100.0], [1.0]u"pA"),
        MassScan(2.0u"s", [100.0], [2.0]u"pA")
    ])
    @test_throws MethodError ticnormalize(mss)
end

end # module TestTicNormalization
