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

    normalized = dwellnormalize(msm, [0.2, 0.3, 0.5], u"s")

    @test rawintensities(normalized) ≈ [50.0 66.66666666666667 60.0; 200.0 166.66666666666669 120.0]
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

    unitful = dwellnormalize(msm, [200.0, 300.0, 500.0]u"ms")
    @test rawintensities(unitful) ≈ [0.05 0.06666666666666667 0.06; 0.2 0.16666666666666666 0.12]
    @test intensityunit(unitful) == JuChrom.inverse(u"ms")

    inferred = dwellnormalize(msm)
    @test rawintensities(inferred) == ints .* 3.0
    @test intensityunit(inferred) == JuChrom.inverse(u"s")

    simultaneous = dwellnormalize(msm; acquisition=:simultaneous)
    @test rawintensities(simultaneous) == ints
    @test intensityunit(simultaneous) == JuChrom.inverse(u"s")

    nonuniform = MassScanMatrix(
        [1.0, 2.1, 3.1]u"s",
        mz,
        [10.0 20.0 30.0; 40.0 50.0 60.0; 70.0 80.0 90.0]
    )
    nonuniforminferred = dwellnormalize(nonuniform)
    @test rawintensities(nonuniforminferred) ≈
        [10.0 20.0 30.0; 40.0 50.0 60.0; 70.0 80.0 90.0] .* 3.0

    many_mz = collect(50.0:1049.0)
    float32retention = MassScanMatrix(
        Float32[0.0, 100.0, 200.0]u"ms",
        many_mz,
        ones(Float32, 3, length(many_mz))
    )
    float32inferred = dwellnormalize(float32retention)
    @test rawintensities(float32inferred) ≈ fill(10.0, 3, length(many_mz))

    @test_throws DimensionMismatch dwellnormalize(msm, [1.0, 2.0], u"s")
    @test_throws ArgumentError dwellnormalize(msm, [1.0, Inf, 5.0], u"s")
    @test_throws ArgumentError dwellnormalize(msm, [1.0, 0.0, 5.0], u"s")
    @test_throws ArgumentError dwellnormalize(msm, [0.5, 0.5, 0.1], u"s")
    @test_throws ArgumentError dwellnormalize(msm; acquisition=:unknown)

    unitful_msm = MassScanMatrix(ret, mz, ints .* u"pA")
    @test_throws ArgumentError dwellnormalize(unitful_msm, [1.0, 2.0, 5.0], u"s")

    annotated_msm = withintensityunit(msm, u"pA")
    @test_throws ArgumentError dwellnormalize(annotated_msm, [1.0, 2.0, 5.0], u"s")

    unitless_retention = MassScanMatrix([1.0, 2.0], mz, ints)
    @test_throws ArgumentError dwellnormalize(unitless_retention)

    single_scan = MassScanMatrix([1.0]u"s", mz, reshape(ints[1, :], 1, :))
    @test_throws ArgumentError dwellnormalize(single_scan)
    @test_throws ArgumentError dwellnormalize(single_scan, [0.2, 0.3, 0.5], u"s")
end

@testset "dwellnormalize VarianceMassScanMatrix" begin
    ret = [1.0, 2.0]u"s"
    mz = [100.0, 200.0, 300.0]
    ints = [10.0 20.0 30.0; 40.0 50.0 60.0]
    vars = [1.0 4.0 9.0; 16.0 25.0 36.0]
    msm = MassScanMatrix(
        ret,
        mz,
        ints;
        level=2,
        instrument=(model="MS",),
        acquisition=(mode="scan",),
        user=(name="tester",),
        sample=(id="S1",),
        extras=Dict("note" => "raw counts")
    )
    vmsm = VarianceMassScanMatrix(msm, vars)

    normalized = dwellnormalize(vmsm, [0.2, 0.3, 0.5], u"s")

    @test normalized isa VarianceMassScanMatrix
    @test rawintensities(normalized) ≈ [50.0 66.66666666666667 60.0; 200.0 166.66666666666669 120.0]
    @test rawvariances(normalized) ≈ [25.0 44.44444444444445 36.0; 400.0 277.77777777777777 144.0]
    @test intensityunit(normalized) == JuChrom.inverse(u"s")
    @test varianceunit(normalized) == JuChrom.inverse(u"s")^2
    @test rawretentions(normalized) == rawretentions(vmsm)
    @test retentionunit(normalized) == retentionunit(vmsm)
    @test rawmzvalues(normalized) == rawmzvalues(vmsm)
    @test mzunit(normalized) == mzunit(vmsm)
    @test level(normalized) == level(vmsm)
    @test instrument(normalized) == instrument(vmsm)
    @test acquisition(normalized) == acquisition(vmsm)
    @test user(normalized) == user(vmsm)
    @test sample(normalized) == sample(vmsm)
    @test extras(normalized) == extras(vmsm)

    extras(msm)["note"] = "changed"
    @test extras(normalized)["note"] == "raw counts"

    unitful = dwellnormalize(vmsm, [200.0, 300.0, 500.0]u"ms")
    @test rawintensities(unitful) ≈ [0.05 0.06666666666666667 0.06; 0.2 0.16666666666666666 0.12]
    @test rawvariances(unitful) ≈ [2.5e-5 4.4444444444444447e-5 3.6e-5; 0.0004 0.0002777777777777778 0.000144]
    @test intensityunit(unitful) == JuChrom.inverse(u"ms")
    @test varianceunit(unitful) == JuChrom.inverse(u"ms")^2

    inferred = dwellnormalize(vmsm)
    @test rawintensities(inferred) == ints .* 3.0
    @test rawvariances(inferred) == vars .* 9.0
    @test intensityunit(inferred) == JuChrom.inverse(u"s")
    @test varianceunit(inferred) == JuChrom.inverse(u"s")^2

    simultaneous = dwellnormalize(vmsm; acquisition=:simultaneous)
    @test rawintensities(simultaneous) == ints
    @test rawvariances(simultaneous) == vars
    @test intensityunit(simultaneous) == JuChrom.inverse(u"s")
    @test varianceunit(simultaneous) == JuChrom.inverse(u"s")^2

    nonuniform = MassScanMatrix(
        [1.0, 2.1, 3.1]u"s",
        mz,
        [10.0 20.0 30.0; 40.0 50.0 60.0; 70.0 80.0 90.0]
    )
    nonuniformvmsm = VarianceMassScanMatrix(
        nonuniform,
        [1.0 4.0 9.0; 16.0 25.0 36.0; 49.0 64.0 81.0]
    )
    nonuniforminferred = dwellnormalize(nonuniformvmsm)
    @test rawintensities(nonuniforminferred) ≈
        [10.0 20.0 30.0; 40.0 50.0 60.0; 70.0 80.0 90.0] .* 3.0
    @test rawvariances(nonuniforminferred) ≈
        [1.0 4.0 9.0; 16.0 25.0 36.0; 49.0 64.0 81.0] .* 9.0

    @test_throws ArgumentError dwellnormalize(vmsm, [0.5, 0.5, 0.1], u"s")
    @test_throws ArgumentError dwellnormalize(vmsm; acquisition=:unknown)
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
