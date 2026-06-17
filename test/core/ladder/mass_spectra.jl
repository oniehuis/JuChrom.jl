using Test
using JuChrom

function test_ladder_mass_spectrum_inputs(; baseline=false)
    retentions = collect(1.0:7.0)
    mzs = [100.0, 101.0, 102.0, 103.0]
    apexretention = 3.25
    profile = exp.(-0.5 .* abs2.((retentions .- apexretention) ./ 1.0))
    heights = [100.0, 60.0, 30.0, -20.0]
    signal_intensities = profile * heights'
    signal = MassScanMatrix(retentions, mzs, signal_intensities)
    variances = ones(size(signal_intensities))
    abundance = profile .* sum(heights[1:3])
    abundanceinfo = (
        abundances=Dict(8 => abundance),
        abundancevariances=Dict(8 => ones(length(retentions)))
    )
    candidate = (
        ladderstep=8,
        scanindex=3,
        window=(leftindex=1, rightindex=length(retentions))
    )
    mzkwargs = (
        retentionref=:middle,
        scaninterval=1e-9,
        mzcount=length(mzs),
        order=:ascending,
        dwellref=:middle,
        dwell=:homogeneous
    )
    apex = alkaneladderapex(
        signal,
        variances,
        abundanceinfo,
        candidate;
        scanwindow=2,
        apexionmzvalues=[100.0, 101.0, 102.0],
        minioncount=3,
        mzretentionkwargs=mzkwargs,
        variancefloor=1.0,
        logfloorfraction=1e-9
    )

    raw = signal
    baselineinfo = nothing
    if baseline
        baseline_intensities = fill(10.0, size(signal_intensities))
        baselines = MassScanMatrix(retentions, mzs, baseline_intensities)
        raw = MassScanMatrix(retentions, mzs, signal_intensities .+ baseline_intensities)
        signal = raw - baselines
        baselineinfo = (baselines=baselines, estimator=:test)
    end

    result = AlkaneSeriesResult(
        nothing,
        variances,
        nothing,
        baselineinfo,
        nothing,
        abundanceinfo,
        nothing,
        (path=NamedTuple[],),
        (apexes=[apex],),
        (gapfilled=NamedTuple[], leftextended=NamedTuple[], rightextended=NamedTuple[]),
        JuChrom.alkane_series_datainfo(raw, signal, variances)
    )

    raw, signal, result
end

@testset "alkaneladdermassspectrum extracts fitted peak-model spectrum" begin
    msm, _, result = test_ladder_mass_spectrum_inputs()

    spectrum = alkaneladdermassspectrum(
        msm,
        result,
        8;
        threaded=false
    )

    ints = intensities(spectrum)
    @test attrs(spectrum).ladderstep == 8
    @test attrs(spectrum).nonnegative
    @test attrs(spectrum).fit_success == trues(4)
    @test ints[1] ≈ 100.0 rtol = 1e-3
    @test ints[2] ≈ 60.0 rtol = 1e-3
    @test ints[3] ≈ 30.0 rtol = 1e-3
    @test ints[4] == 0.0

    signed = alkaneladdermassspectrum(
        msm,
        result,
        8;
        nonnegative=false,
        threaded=false
    )
    @test attrs(signed).allownegative
    @test intensities(signed)[4] ≈ -20.0 rtol = 1e-3
end

@testset "alkaneladdermassspectra validates checksums and reconstructs baseline signal" begin
    msm, signal, result = test_ladder_mass_spectrum_inputs(; baseline=true)

    extracted = alkaneladdermassspectra(
        msm,
        result;
        threaded=false
    )
    @test isempty(extracted.failures)
    @test haskey(extracted.spectra, 8)
    @test intensities(extracted.spectra[8])[1] ≈ 100.0 rtol = 1e-3

    @test_throws ArgumentError alkaneladdermassspectra(
        signal,
        result;
        threaded=false
    )

    bad = MassScanMatrix(
        rawretentions(msm),
        rawmzvalues(msm),
        rawintensities(msm) .+ 1.0
    )
    @test_throws ArgumentError alkaneladdermassspectra(
        bad,
        result;
        threaded=false
    )
end
