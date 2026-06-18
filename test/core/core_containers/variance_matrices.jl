module TestVarianceMatrices

using Test
using Unitful
using Unitful: nA, pA, s
using JuChrom

@testset "VarianceMassScanMatrix construction and getter forwarding" begin
    ret = [1.0, 2.0]u"s"
    mzs = [100.0, 200.0]
    ints = [1.0 2.0; 3.0 4.0]u"pA"
    vars = [0.1 0.2; 0.3 0.4]

    msm = MassScanMatrix(ret, mzs, ints;
        level=2,
        instrument=(vendor="Acme",),
        acquisition=(mode="MS",),
        user=(name="test",),
        sample=(id="sample",),
        extras=Dict("key" => "value"))

    vmsm = VarianceMassScanMatrix(msm, vars)

    @test vmsm isa AbstractMassScanMatrix{typeof(u"s"), Nothing, typeof(u"pA")}
    @test vmsm isa AbstractVarianceMassScanMatrix{
        typeof(u"s"), Nothing, typeof(u"pA"), typeof(u"pA"^2)}
    @test parent(vmsm) ≡ msm
    @test Base.broadcastable(vmsm)[] ≡ vmsm
    @test varianceunit(vmsm) == u"pA"^2
    @test rawvariances(vmsm) == vars
    @test variances(vmsm) == vars .* (u"pA"^2)
    @test rawvariances(vmsm; unit=u"nA"^2) == vars .* 1e-6

    @test retentions(vmsm) == retentions(msm)
    @test rawretentions(vmsm) == rawretentions(msm)
    @test mzvalues(vmsm) == mzvalues(msm)
    @test rawmzvalues(vmsm) == rawmzvalues(msm)
    @test intensities(vmsm) == intensities(msm)
    @test rawintensities(vmsm) == rawintensities(msm)
    @test retentionunit(vmsm) == retentionunit(msm)
    @test mzunit(vmsm) == mzunit(msm)
    @test intensityunit(vmsm) == intensityunit(msm)
    @test scancount(vmsm) == scancount(msm)
    @test mzcount(vmsm) == mzcount(msm)
    @test level(vmsm) == level(msm)
    @test instrument(vmsm) == instrument(msm)
    @test acquisition(vmsm) == acquisition(msm)
    @test user(vmsm) == user(msm)
    @test sample(vmsm) == sample(msm)
    @test extras(vmsm) == extras(msm)

    qvars = vars .* (u"pA"^2)
    vmsm_unitful_vars = VarianceMassScanMatrix(msm, qvars)
    @test varianceunit(vmsm_unitful_vars) == u"pA"^2
    @test rawvariances(vmsm_unitful_vars) == vars
end

@testset "VarianceMassScanMatrix validation" begin
    msm_unitless = MassScanMatrix([1.0, 2.0], [100.0, 200.0], [1.0 2.0; 3.0 4.0])
    vars = [0.1 0.2; 0.3 0.4]
    vmsm_unitless = VarianceMassScanMatrix(msm_unitless, vars)

    @test isnothing(varianceunit(vmsm_unitless))
    @test variances(vmsm_unitless) ≡ vars
    @test rawvariances(vmsm_unitless) ≡ vars
    @test_throws ArgumentError variances(vmsm_unitless; unit=u"pA"^2)

    @test_throws DimensionMismatch VarianceMassScanMatrix(msm_unitless, vars[1:1, :])
    @test_throws ArgumentError VarianceMassScanMatrix(msm_unitless, [0.1 -0.2; 0.3 0.4])
    @test_throws ArgumentError VarianceMassScanMatrix(msm_unitless, [0.1 NaN; 0.3 0.4])
    @test_throws ArgumentError VarianceMassScanMatrix(msm_unitless, vars, u"pA"^2)
    @test_throws ArgumentError VarianceMassScanMatrix(msm_unitless, vars .* (u"pA"^2))

    msm_unitful = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0],
        [1.0 2.0; 3.0 4.0]u"pA")

    @test_throws ArgumentError VarianceMassScanMatrix(msm_unitful, vars, nothing)
    @test_throws Unitful.DimensionError VarianceMassScanMatrix(msm_unitful, vars .* u"pA")
end

@testset "VarianceMassScanMatrix show methods" begin
    msm_plain = MassScanMatrix([1.0, 2.0], [100.0, 200.0], [1.0 2.0; 3.0 4.0])
    vmsm_plain = VarianceMassScanMatrix(msm_plain, [0.1 0.2; 0.3 0.4])
    out_plain = sprint(show, vmsm_plain)

    @test occursin("VarianceMassScanMatrix with 2 scans", out_plain)
    @test occursin("├─ Retention:", out_plain)
    @test occursin("Range: 1.0 to 2.0 (unitless)", out_plain)
    @test occursin("├─ M/Z values:", out_plain)
    @test occursin("Total data points: 2", out_plain)
    @test occursin("├─ Intensity:", out_plain)
    @test occursin("MS Level: 1", out_plain)
    @test occursin("└─ Variance:", out_plain)
    @test occursin("Range: 0.1 to 0.4 (unitless)", out_plain)
    @test !occursin("Scan type", out_plain)

    msm_meta = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0]u"Th",
                              [1.0 2.0; 3.0 4.0]u"pA";
                              instrument=(vendor="A",))
    vmsm_meta = VarianceMassScanMatrix(msm_meta, [0.1 0.2; 0.3 0.4])
    out_meta = sprint(io -> show(io, MIME"text/plain"(), vmsm_meta))

    @test occursin("VarianceMassScanMatrix with 2 scans", out_meta)
    @test occursin("Range: 1.0 to 2.0 (s)", out_meta)
    @test occursin("Range: 100.0 to 200.0 (Th)", out_meta)
    @test occursin("Range: 1.0 to 4.0 (pA)", out_meta)
    @test occursin("├─ Variance:", out_meta)
    @test occursin("Range: 0.1 to 0.4 ($(string(varianceunit(vmsm_meta))))", out_meta)
    @test occursin("Instrument", out_meta)
    @test occursin("vendor = A", out_meta)
    @test !occursin("Scan type", out_meta)
end

end
