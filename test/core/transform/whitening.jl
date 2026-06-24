module TestWhitening

using Test
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

@testset "whiten rejects invalid inputs" begin
    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]
    x = [1.0 2.0;
         3.0 4.0]
    vmsm = VarianceMassScanMatrix(MassScanMatrix(rets, mzs, x), ones(size(x)))

    @test_throws ArgumentError whiten(vmsm, 0.0)
    @test_throws ArgumentError whiten(vmsm, -1.0)

    physical_msm = MassScanMatrix(rets, mzs, x .* u"pA")
    physical_vmsm = VarianceMassScanMatrix(physical_msm, ones(size(x)) .* u"pA"^2)
    @test_throws ArgumentError whiten(physical_vmsm, 0.1)

    @test_throws MethodError whiten(x, ones(size(x)), 0.1)
end

end # module TestWhitening
