using Test
using Unitful
using JuChrom

function test_ladder_mass_spectrum_inputs(; baseline=false)
    retentions = collect(1.0:7.0)
    mzs = [100.0, 101.0, 102.0, 103.0]
    apexretention = 3.25
    profile = exp.(-0.5 .* abs2.((retentions .- apexretention) ./ 1.0))
    heights = [100.0, 60.0, 30.0, -20.0]
    signal_intensities = profile * heights'
    signal = MassScanMatrix(retentions, u"s", mzs, nothing, signal_intensities, nothing)
    variances = ones(size(signal_intensities))
    abundance = profile .* sum(heights[1:3])
    abundanceinfo = AlkaneAbundanceInfo(
        Dict(8 => abundance),
        Dict(8 => ones(length(retentions))),
        Dict{Int, Vector{AlkaneAbundanceWindow}}(),
        NamedTuple()
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
    settings = JuChrom.AlkaneLadderApexSettings(
        defaultalkanestandard(),
        2,
        1.0,
        1e-9,
        mzkwargs,
        :ascending,
        (),
        [100.0, 101.0, 102.0],
        0.1,
        3,
        3.0,
        0.25,
        3.0,
        6,
        3,
        3,
        1.25,
        3,
        14,
        5,
        8
    )
    apex = alkaneladderapex(
        signal,
        variances,
        abundanceinfo,
        candidate,
        settings
    )

    raw = signal
    baselineinfo = nothing
    if baseline
        baseline_intensities = fill(10.0, size(signal_intensities))
        baselines = MassScanMatrix(retentions, u"s", mzs, nothing, baseline_intensities, nothing)
        raw = MassScanMatrix(
            retentions,
            u"s",
            mzs,
            nothing,
            signal_intensities .+ baseline_intensities,
            nothing
        )
        signal = raw - baselines
        baselineinfo = JuChrom.AlkaneBaselineInfo(
            baselines,
            :test,
            1.0,
            true,
            10.0,
            0.5,
            1.0,
            0.2
        )
    end

    apexinfo = JuChrom.AlkaneLadderApexInfo(
        :success,
        :success,
        [apex],
        Dict(apex.ladderstep => apex),
        settings,
        JuChrom.alkane_ladder_no_scan_order_info(:test),
        [apex.calibration_excluded],
        [apex.good_for_calibration],
        [apex.apex_fit_quality_score],
        [apex.apex_fit_quality_zscore]
    )
    additionsettings = JuChrom.AlkaneLadderAdditionSettings(
        5.0,
        0.15,
        0.05,
        100,
        false,
        0.85,
        0.03,
        6,
        5.0,
        0.2,
        0.9,
        0.03,
        3,
        0.1,
        1.0,
        2,
        1.0,
        1e-9,
        JuChrom.DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
        [100.0, 101.0, 102.0],
        0.1,
        3,
        mzkwargs,
        3.0,
        0.25,
        3.0,
        6,
        [8]
    )
    additioninfo = JuChrom.AlkaneLadderAdditionInfo(
        :empty,
        JuChrom.AlkaneLadderAddition[],
        JuChrom.AlkaneLadderAddition[],
        JuChrom.AlkaneLadderAddition[],
        JuChrom.AlkaneLadderAddition[],
        JuChrom.AlkaneLadderAdditionDiagnostics(
            JuChrom.AlkaneLadderAdditionDiagnostic[],
            JuChrom.AlkaneLadderAdditionDiagnostic[],
            JuChrom.AlkaneLadderAdditionDiagnostic[]
        ),
        additionsettings
    )

    result = AlkaneSeriesResult(
        nothing,
        variances,
        nothing,
        baselineinfo,
        nothing,
        abundanceinfo,
        nothing,
        nothing,
        apexinfo,
        additioninfo,
        JuChrom.alkane_series_datainfo(raw, signal, variances),
        retentionunit(signal)
    )

    @test result.retentionunit == u"s"
    @test retentionunit(result) == retentionunit(signal)

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
    @test result.datainfo isa JuChrom.AlkaneSeriesDataInfo

    extracted = alkaneladdermassspectra(
        msm,
        result;
        threaded=false
    )
    @test extracted isa JuChrom.AlkaneLadderMassSpectrumExtraction
    @test extracted.settings.threaded == false
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
