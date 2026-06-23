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
        signal = subtractbaseline(raw, baselines)
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
        false,
        :missing_standard,
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

@testset "alkane ladder mass spectrum extraction records per-step failures" begin
    msm, _, result = test_ladder_mass_spectrum_inputs()
    step = only(alkaneladdersteps(result))
    settings = JuChrom.alkane_ladder_mass_spectrum_settings(
        true,
        1.0,
        false,
        true,
        false,
        true,
        true,
        true
    )

    gatedstep = AlkaneLadderStep(
        step.ladderstep,
        step.apexscanindex,
        step.apexretention,
        step.source,
        step.massspectrumcosine,
        step.requiredcosine,
        false,
        step.apex
    )
    gated = JuChrom.alkane_ladder_extract_step_mass_spectrum(
        msm,
        gatedstep,
        result.variances,
        settings,
        false
    )
    @test gated.failure == "ladder step did not pass apex fit quality gate"

    apexfields = fieldnames(typeof(step.apex))
    badapex = JuChrom.AlkaneLadderApex(
        (field ≡ :fit ? nothing : getfield(step.apex, field) for field in apexfields)...
    )
    badapexinfo = JuChrom.AlkaneLadderApexInfo(
        result.apexinfo.status,
        result.apexinfo.reason,
        [badapex],
        Dict(badapex.ladderstep => badapex),
        result.apexinfo.settings,
        result.apexinfo.scanorderinfo,
        result.apexinfo.calibrationexcluded,
        result.apexinfo.goodforcalibration,
        result.apexinfo.apexfitqualityscores,
        result.apexinfo.apexfitqualityzscores
    )
    badresult = AlkaneSeriesResult(
        result.success,
        result.status,
        result.standard,
        result.variances,
        result.varianceinfo,
        result.baselineinfo,
        result.channelinfo,
        result.abundanceinfo,
        result.molecularioninfo,
        result.pathinfo,
        badapexinfo,
        result.additioninfo,
        result.datainfo,
        result.retentionunit
    )

    extraction = alkaneladdermassspectra(
        msm,
        badresult;
        threaded=false,
        validatechecksum=false
    )

    @test isempty(extraction.spectra)
    @test occursin(
        "does not contain a usable peak model apex",
        extraction.failures[8]
    )
end

@testset "alkane ladder mass spectrum extraction supports ion-threaded fitting" begin
    msm, _, result = test_ladder_mass_spectrum_inputs()
    step = only(alkaneladdersteps(result))
    settings = JuChrom.alkane_ladder_mass_spectrum_settings(
        true,
        1.0,
        true,
        true,
        false,
        true,
        true,
        true
    )

    spectrum = JuChrom.alkane_ladder_step_mass_spectrum(
        msm,
        step,
        result.variances,
        settings,
        true
    )

    @test attrs(spectrum).ion_threaded
    @test intensities(spectrum)[1] ≈ 100.0 rtol = 1e-3
end

@testset "alkane ladder mass spectrum extraction supports step-threaded fitting" begin
    msm, _, result = test_ladder_mass_spectrum_inputs()
    apex = only(result.apexinfo.apexes)
    apexfields = fieldnames(typeof(apex))
    secondapex = JuChrom.AlkaneLadderApex(
        (
            field ≡ :ladderstep ? 9 :
            field ≡ :source ? :molecularion :
            getfield(apex, field)
            for field in apexfields
        )...
    )
    apexinfo = JuChrom.AlkaneLadderApexInfo(
        result.apexinfo.status,
        result.apexinfo.reason,
        [apex, secondapex],
        Dict(8 => apex, 9 => secondapex),
        result.apexinfo.settings,
        result.apexinfo.scanorderinfo,
        [false, false],
        [true, true],
        [apex.apex_fit_quality_score, secondapex.apex_fit_quality_score],
        [apex.apex_fit_quality_zscore, secondapex.apex_fit_quality_zscore]
    )
    multiresult = AlkaneSeriesResult(
        result.success,
        result.status,
        result.standard,
        result.variances,
        result.varianceinfo,
        result.baselineinfo,
        result.channelinfo,
        result.abundanceinfo,
        result.molecularioninfo,
        result.pathinfo,
        apexinfo,
        result.additioninfo,
        result.datainfo,
        result.retentionunit
    )

    extraction = alkaneladdermassspectra(
        msm,
        multiresult;
        threaded=true,
        validatechecksum=false
    )

    @test sort(collect(keys(extraction.spectra))) == [8, 9]
    @test isempty(extraction.failures)
    @test !attrs(extraction.spectra[8]).ion_threaded
    @test !attrs(extraction.spectra[9]).ion_threaded
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
