module TestClr

using Statistics: mean
using Test
using Unitful
using JuChrom

@testset "clr(vmsm::AbstractVarianceMassScanMatrix)" begin
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
        extras=Dict("run" => "A"),
    )
    vmsm = VarianceMassScanMatrix(msm, v .* u"pA"^2)

    transformed = clr(vmsm)

    expected_x = log.(x) .- mean(log.(x))
    sigma2_log = v ./ abs2.(x)
    N = length(x)
    total_sigma2_log = sum(sigma2_log)
    expected_v = @. sigma2_log * (1 - 2 / N) + total_sigma2_log / N^2

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
    @test isapprox(sum(rawintensities(transformed)), 0.0; atol=1e-12)
    @test all(≥(0), rawvariances(transformed))

    @test rawintensities(vmsm) == x
    @test rawvariances(vmsm) ≈ v
end

@testset "clr converts variance units to the intensity scale" begin
    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]
    x = [1000.0 2000.0;
         4000.0 8000.0]
    v_pA2 = [1.0 4.0;
             16.0 64.0]
    msm = MassScanMatrix(rets, mzs, x .* u"pA")
    vmsm = VarianceMassScanMatrix(msm, uconvert.(u"nA"^2, v_pA2 .* u"pA"^2))

    transformed = clr(vmsm)

    sigma2_log = v_pA2 ./ abs2.(x)
    N = length(x)
    total_sigma2_log = sum(sigma2_log)
    expected_v = @. sigma2_log * (1 - 2 / N) + total_sigma2_log / N^2

    @test rawvariances(transformed) ≈ expected_v
end

@testset "clr requires strictly positive intensities" begin
    rets = [1.0, 2.0]u"s"
    mzs = [100.0, 101.0]
    x = [1.0 0.0;
         2.0 3.0]
    msm = MassScanMatrix(rets, mzs, x .* u"pA")
    vmsm = VarianceMassScanMatrix(msm, ones(size(x)) .* u"pA"^2)

    @test_throws DomainError clr(vmsm)
end

@testset "array clr methods are not public API" begin
    @test_throws MethodError clr([1.0, 2.0, 3.0])
    @test_throws MethodError clr([1.0, 2.0, 3.0], [0.1, 0.1, 0.1])
end

end # module TestClr
