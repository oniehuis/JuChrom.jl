module TestWhitening

using Test
using Statistics: quantile
using Unitful
using JuChrom

@testset "whiten(vmsm::AbstractVarianceMassScanMatrix, sigmafloor)" begin
    @test JuChrom.isdimensionlessvarianceunit(nothing)

    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]u"Th"
    x = [1.0 2.0;
         3.0 4.0]
    v = [1.0 0.25;
         1e-12 9.0]
    msm = MassScanMatrix(
        rets,
        mzs,
        x;
        level=2,
        instrument=(name="test instrument",),
        acquisition=(mode="scan",),
        user=(analyst="tester",),
        sample=(id="sample-a",),
        extras=Dict("run" => "A"),
    )
    vmsm = VarianceMassScanMatrix(msm, v)

    sigmafloor = 1e-3
    transformed = whiten(vmsm, sigmafloor)

    denominator_variances = max.(v, sigmafloor^2)
    expected_x = x ./ sqrt.(denominator_variances)
    expected_v = v ./ denominator_variances

    @test transformed isa VarianceMassScanMatrix
    @test parent(transformed) isa MassScanMatrix
    @test retentionunit(transformed) == u"s"
    @test mzunit(transformed) == u"Th"
    @test intensityunit(transformed) === nothing
    @test varianceunit(transformed) === nothing
    @test level(transformed) == 2
    @test instrument(transformed) == instrument(vmsm)
    @test acquisition(transformed) == acquisition(vmsm)
    @test user(transformed) == user(vmsm)
    @test sample(transformed) == sample(vmsm)
    @test extras(transformed) == extras(vmsm)
    @test rawretentions(transformed) == rawretentions(vmsm)
    @test rawmzvalues(transformed) == rawmzvalues(vmsm)
    @test rawintensities(transformed) ≈ expected_x
    @test rawvariances(transformed) ≈ expected_v
    @test rawvariances(transformed)[1, 1] == 1.0
    @test rawvariances(transformed)[2, 1] ≈ 1e-6

    @test rawintensities(vmsm) == x
    @test rawvariances(vmsm) == v
end

@testset "whiten infers sigmafloor automatically" begin
    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]u"Th"
    x = [1.0 2.0;
         3.0 4.0]
    v = [1.0 0.25;
         1e-12 9.0]
    vmsm = VarianceMassScanMatrix(MassScanMatrix(rets, mzs, x), v)

    floorquantile = 0.5
    sigmafloor = quantile(sqrt.(filter(>(0), vec(v))), floorquantile)
    denominator_variances = max.(v, sigmafloor^2)

    transformed = whiten(vmsm; floorquantile=floorquantile)
    positional = whiten(vmsm, :auto; floorquantile=floorquantile)

    @test rawintensities(transformed) ≈ x ./ sqrt.(denominator_variances)
    @test rawvariances(transformed) ≈ v ./ denominator_variances
    @test rawintensities(positional) ≈ rawintensities(transformed)
    @test rawvariances(positional) ≈ rawvariances(transformed)
    @test intensityunit(transformed) === nothing
    @test varianceunit(transformed) === nothing
end

@testset "whiten accepts dimensionless variance units" begin
    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]
    x = [1.0 2.0;
         3.0 4.0]
    v = [1.0 4.0;
         9.0 16.0]
    msm = MassScanMatrix(rets, mzs, x .* u"rad")
    vmsm = VarianceMassScanMatrix(msm, v .* u"rad"^2)

    transformed = whiten(vmsm, 0.1)

    @test intensityunit(transformed) === nothing
    @test varianceunit(transformed) === nothing
    @test rawintensities(transformed) ≈ x ./ sqrt.(v)
    @test rawvariances(transformed) ≈ ones(size(v))
end

@testset "whiten accepts unitful floors for physical variances" begin
    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]
    x = [1.0 2.0;
         3.0 4.0]
    v = [1.0 0.25;
         1e-12 9.0]
    vmsm = VarianceMassScanMatrix(
        MassScanMatrix(rets, mzs, x .* u"pA"),
        v .* u"pA"^2,
    )

    sigmafloor = 0.5
    denominator_variances = max.(v, sigmafloor^2)
    transformed = whiten(vmsm, sigmafloor * u"pA")
    converted = whiten(vmsm, 500u"fA")

    @test rawintensities(transformed) ≈ x ./ sqrt.(denominator_variances)
    @test rawvariances(transformed) ≈ v ./ denominator_variances
    @test rawintensities(converted) ≈ rawintensities(transformed)
    @test rawvariances(converted) ≈ rawvariances(transformed)
    @test intensityunit(transformed) === nothing
    @test varianceunit(transformed) === nothing

    floorquantile = 0.5
    autofloor = quantile(sqrt.(filter(>(0), vec(v))), floorquantile)
    auto_denominator_variances = max.(v, autofloor^2)
    automatic = whiten(vmsm; floorquantile=floorquantile)

    @test rawintensities(automatic) ≈ x ./ sqrt.(auto_denominator_variances)
    @test rawvariances(automatic) ≈ v ./ auto_denominator_variances
    @test intensityunit(automatic) === nothing
    @test varianceunit(automatic) === nothing
end

@testset "whiten rejects invalid inputs" begin
    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]
    x = [1.0 2.0;
         3.0 4.0]
    vmsm = VarianceMassScanMatrix(MassScanMatrix(rets, mzs, x), ones(size(x)))

    @test_throws ArgumentError whiten(vmsm, 0.0)
    @test_throws ArgumentError whiten(vmsm, -1.0)
    @test_throws ArgumentError whiten(vmsm; floorquantile=0.0)
    @test_throws ArgumentError whiten(vmsm; floorquantile=1.1)
    @test_throws ArgumentError whiten(vmsm, :median)
    @test_throws ArgumentError whiten(vmsm, 0.1u"pA")

    zero_vmsm = VarianceMassScanMatrix(MassScanMatrix(rets, mzs, x), zeros(size(x)))
    @test_throws ArgumentError whiten(zero_vmsm)

    physical_msm = MassScanMatrix(rets, mzs, x .* u"pA")
    physical_vmsm = VarianceMassScanMatrix(physical_msm, ones(size(x)) .* u"pA"^2)
    @test_throws ArgumentError whiten(physical_vmsm, 0.1)
    @test_throws ArgumentError whiten(physical_vmsm, 0.0u"pA")
    @test_throws ArgumentError whiten(physical_vmsm, -1.0u"pA")
    @test_throws Unitful.DimensionError whiten(physical_vmsm, 0.1u"s")

    @test_throws MethodError whiten(x, ones(size(x)), 0.1)
end

end # module TestWhitening
