module TestIntensityUnits

using Test
using Unitful
using JuChrom

@testset "dwellnormalize MassScanMatrix" begin
    ret = [1.0, 2.0]u"s"
    mz = [100.0, 200.0, 300.0]
    ints = [10.0 20.0 30.0; 40.0 50.0 60.0]
    msm = MassScanMatrix(
        ret,
        mz,
        ints;
        level=2,
        instrument=(model="MS",),
        acquisition=(mode="scan",),
        user=(name="tester",),
        sample=(id="S1",),
        extras=Dict("note" => "raw counts"),
    )

    normalized = dwellnormalize(msm, [1.0, 2.0, 5.0], u"s")

    @test rawintensities(normalized) == [10.0 10.0 6.0; 40.0 25.0 12.0]
    @test intensityunit(normalized) == JuChrom.inverse(u"s")
    @test rawretentions(normalized) == rawretentions(msm)
    @test retentionunit(normalized) == retentionunit(msm)
    @test rawmzvalues(normalized) == rawmzvalues(msm)
    @test mzunit(normalized) == mzunit(msm)
    @test level(normalized) == level(msm)
    @test instrument(normalized) == instrument(msm)
    @test acquisition(normalized) == acquisition(msm)
    @test user(normalized) == user(msm)
    @test sample(normalized) == sample(msm)
    @test extras(normalized) == extras(msm)

    extras(msm)["note"] = "changed"
    @test extras(normalized)["note"] == "raw counts"

    unitful = dwellnormalize(msm, [1.0, 2.0, 5.0]u"ms")
    @test rawintensities(unitful) == rawintensities(normalized)
    @test intensityunit(unitful) == JuChrom.inverse(u"ms")

    count_msm = withintensityunit(msm, u"count")
    count_normalized = dwellnormalize(count_msm, [1.0, 2.0, 5.0], u"s")
    @test rawintensities(count_normalized) == rawintensities(normalized)
    @test intensityunit(count_normalized) == u"count" * JuChrom.inverse(u"s")

    count_unitful = dwellnormalize(count_msm, [1.0, 2.0, 5.0]u"ms")
    @test rawintensities(count_unitful) == rawintensities(normalized)
    @test intensityunit(count_unitful) == u"count" * JuChrom.inverse(u"ms")

    inferred = dwellnormalize(msm)
    @test rawintensities(inferred) == ints .* 3.0
    @test intensityunit(inferred) == JuChrom.inverse(u"s")

    count_inferred = dwellnormalize(count_msm)
    @test rawintensities(count_inferred) == ints .* 3.0
    @test intensityunit(count_inferred) == u"count" * JuChrom.inverse(u"s")

    @test_throws DimensionMismatch dwellnormalize(msm, [1.0, 2.0], u"s")
    @test_throws ArgumentError dwellnormalize(msm, [1.0, Inf, 5.0], u"s")
    @test_throws ArgumentError dwellnormalize(msm, [1.0, 0.0, 5.0], u"s")

    unitful_msm = MassScanMatrix(ret, mz, ints .* u"pA")
    @test_throws ArgumentError dwellnormalize(unitful_msm, [1.0, 2.0, 5.0], u"s")

    rate_msm = withintensityunit(msm, u"count" * JuChrom.inverse(u"s"))
    @test_throws ArgumentError dwellnormalize(rate_msm, [1.0, 2.0, 5.0], u"s")

    unitless_retention = MassScanMatrix([1.0, 2.0], mz, ints)
    @test_throws ArgumentError dwellnormalize(unitless_retention)

    single_scan = MassScanMatrix([1.0]u"s", mz, reshape(ints[1, :], 1, :))
    @test_throws ArgumentError dwellnormalize(single_scan)
end

@testset "withintensityunit MassScanMatrix" begin
    ret = [1.0, 2.0]u"s"
    mz = [100.0, 200.0, 300.0]
    ints = [10.0 20.0 30.0; 40.0 50.0 60.0]
    msm = MassScanMatrix(
        ret,
        mz,
        ints;
        level=2,
        instrument=(model="MS",),
        acquisition=(mode="scan",),
        user=(name="tester",),
        sample=(id="S1",),
        extras=Dict("note" => "raw counts"),
    )

    annotated = withintensityunit(msm, u"pA")

    @test rawintensities(annotated) == rawintensities(msm)
    @test intensityunit(annotated) == u"pA"
    @test rawretentions(annotated) == rawretentions(msm)
    @test retentionunit(annotated) == retentionunit(msm)
    @test rawmzvalues(annotated) == rawmzvalues(msm)
    @test mzunit(annotated) == mzunit(msm)
    @test level(annotated) == level(msm)
    @test instrument(annotated) == instrument(msm)
    @test acquisition(annotated) == acquisition(msm)
    @test user(annotated) == user(msm)
    @test sample(annotated) == sample(msm)
    @test extras(annotated) == extras(msm)

    extras(msm)["note"] = "changed"
    @test extras(annotated)["note"] == "raw counts"

    unitful_msm = MassScanMatrix(ret, mz, ints .* u"pA")
    @test_throws ArgumentError withintensityunit(unitful_msm, u"nA")
end

end
