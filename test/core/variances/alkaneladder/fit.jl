using Test
using Unitful
using JuChrom

function _alkane_variance_entrypoint_inputs()
    intensities = reshape([10.0, 20.0], 2, 1)
    msm = MassScanMatrix([1.0, 2.0], [57.0], intensities)
    variances = ones(2, 1)
    result = AlkaneSeriesResult(
        nothing,
        variances,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        JuChrom.alkane_series_datainfo(msm, msm, variances),
        retentionunit(msm),
    )

    msm, result
end

function _alkane_variance_mass_spectrum_inputs(; baseline=false)
    retentions = collect(1.0:41.0)
    mzs = [100.0, 101.0, 102.0, 103.0]
    truepeakretention = 21.25
    profile = exp.(-0.5 .* abs2.((retentions .- truepeakretention) ./ 1.0))
    heights = [100.0, 60.0, 30.0, -20.0]
    signal_intensities = profile * heights'
    signal = MassScanMatrix(retentions, u"s", mzs, nothing, signal_intensities, nothing)
    variances = ones(size(signal_intensities))
    abundance = profile .* sum(heights[1:3])
    window = AlkaneAbundanceWindow(
        8,
        18,
        21,
        24,
        abundance[18],
        abundance[21],
        abundance[24],
        0.0,
        :boundary,
        :boundary,
    )
    abundanceinfo = AlkaneAbundanceInfo(
        Dict(8 => abundance),
        Dict(8 => ones(length(retentions))),
        Dict(8 => [window]),
        NamedTuple(),
    )
    candidate = (
        ladderstep=8,
        scanindex=21,
        window=(leftindex=18, rightindex=24),
    )
    mzkwargs = (
        retentionref=:middle,
        scaninterval=1e-9,
        mzcount=length(mzs),
        order=:ascending,
        dwellref=:middle,
        dwell=:homogeneous,
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
        8,
    )
    apex = alkaneladderapex(
        signal,
        variances,
        abundanceinfo,
        candidate,
        settings,
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
            nothing,
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
            0.2,
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
        [apex.apex_fit_quality_zscore],
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
        [8],
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
            JuChrom.AlkaneLadderAdditionDiagnostic[],
        ),
        additionsettings,
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
        retentionunit(signal),
    )

    raw, signal, result
end

function _alkane_variance_result_with_abundance_window(result, window)
    abundanceinfo = result.abundanceinfo
    replacement = AlkaneAbundanceInfo(
        copy(abundanceinfo.abundances),
        copy(abundanceinfo.abundancevariances),
        Dict(Int(window.ladderstep) => [window]),
        abundanceinfo.settings,
    )

    AlkaneSeriesResult(
        result.standard,
        result.variances,
        result.varianceinfo,
        result.baselineinfo,
        result.channelinfo,
        replacement,
        result.molecularioninfo,
        result.pathinfo,
        result.apexinfo,
        result.additioninfo,
        result.datainfo,
        result.retentionunit,
    )
end

@testset "fitalkanevariancemodel entry point" begin
    oldmodel = JuChrom.LinearObservedIntensityVarianceModel(10.0, 2.0, 0.0, 100.0, 0.0)
    @test oldmodel.intensity_offset == 0.0
    @test varpred(5.0, oldmodel) == 20.0

    offsetmodel = JuChrom.LinearObservedIntensityVarianceModel(
        10.0,
        2.0,
        5.0,
        0.0,
        100.0,
        0.0,
    )
    @test offsetmodel.intensity_offset == 5.0
    @test varpred(4.0, offsetmodel) == 10.0
    @test varpred(5.0, offsetmodel) == 10.0
    @test varpred(8.0, offsetmodel) == 16.0
    @test varpred([-1.0, 5.0, 8.0], offsetmodel) == [10.0, 10.0, 16.0]
    @test varpred(101.0, offsetmodel; extrapolation=:allow) == 202.0
    @test_logs (:warn, r"outside the calibrated range") varpred(
        101.0,
        offsetmodel;
        extrapolation=:warn,
    )
    @test_throws ArgumentError varpred(101.0, offsetmodel; extrapolation=:throw)
    @test_throws ArgumentError varpred(
        8.0,
        offsetmodel;
        extrapolation=:unknown,
    )

    unitmodel = JuChrom.LinearObservedIntensityVarianceModel(
        10.0u"pA^2",
        2.0u"pA",
        5.0u"pA",
        0.0,
        100.0,
        0.0,
    )
    @test varpred(8.0u"pA", unitmodel) == 16.0u"pA^2"
    @test_throws ArgumentError varpred(101.0u"pA", unitmodel; extrapolation=:throw)

    raw, signal, result = _alkane_variance_mass_spectrum_inputs(; baseline=true)

    fit = fitalkanevariancemodel(raw, result)
    @test fit isa AlkaneVarianceFit
    @test fit.success
    @test fit.status === :ok
    @test fit.model isa JuChrom.LinearObservedIntensityVarianceModel
    modeldisplay = sprint(io -> show(io, MIME"text/plain"(), fit.model))
    @test occursin("LinearObservedIntensityVarianceModel", modeldisplay)
    @test occursin("variance(I)", modeldisplay)
    fitdisplay = sprint(io -> show(io, MIME"text/plain"(), fit))
    @test occursin("AlkaneVarianceFit", fitdisplay)
    @test occursin("success: true", fitdisplay)
    @test !occursin("status: ok", fitdisplay)
    @test occursin("model:", fitdisplay)
    @test occursin("range:", fitdisplay)
    @test occursin("data:", fitdisplay)
    @test occursin("diagnostics:", fitdisplay)
    @test occursin("lag1 rho", fitdisplay)
    @test !occursin("residualrecords", fitdisplay)
    @test fit.model.intercept ≥ 0
    @test fit.model.slope ≥ 0
    @test isfinite(fit.model.intensity_offset)
    @test fit.model.intensity_min == 0.0
    @test fit.model.intensity_max > fit.model.intensity_min
    @test abs(fit.model.rho_lag1) < 1
    predicted_variances = varpred(raw, fit.model)
    @test predicted_variances isa VarianceMassScanMatrix
    @test parent(predicted_variances) === raw
    @test rawvariances(predicted_variances) ≈ varpred(rawintensities(raw), fit.model)
    @test rawvariances(varpred(raw, fit)) ≈ rawvariances(predicted_variances)
    @test varpred(8.0, fit) == varpred(8.0, fit.model)
    @test size(rawvariances(predicted_variances)) == size(rawintensities(raw))
    @test fit.qc.status === :ok
    @test fit.qc.accept
    @test fit.qc.flatrecordcount == length(fit.flatrecords)
    @test fit.qc.peakrecordcount == length(fit.residualrecords)
    @test fit.variancefit.model == fit.model
    @test fit.variancefit.flatvarianceanchor ≈ fit.model.intercept
    @test fit.variancefit.intensity_offset ≈ fit.model.intensity_offset
    @test fit.variancefit.robustslope.slope ≈ fit.model.slope
    @test fit.variancefit.lag1.rho ≈ fit.model.rho_lag1
    @test fit.variancefit.lag1.source === :flat_nonpeak
    @test fit.variancefit.lag1.paircount == fit.qc.lag1paircount
    @test fit.variancefit.lag1.paircount ≤ length(fit.flatrecords) - 1
    @test !isempty(fit.flatrecords)
    @test all(record -> record.ladderstep == 0, fit.flatrecords)
    @test all(record -> !(13 ≤ record.scanindex ≤ 29), fit.flatrecords)
    @test fit.excludeladdersteps == Int[]
    @test fit.includedladdersteps == [8]
    @test rawintensities(fit.signal) ≈ rawintensities(signal)
    @test Set(keys(fit.spectra)) == Set([8])
    @test isempty(fit.failures)
    @test attrs(fit.spectra[8]).ladderstep == 8
    @test Set(keys(fit.peakinputs)) == Set([8])
    @test isempty(fit.peakinputfailures)
    peakinput = fit.peakinputs[8]
    @test peakinput.ladderstep == 8
    @test Set(propertynames(peakinput)) == Set([
        :ladderstep,
        :scanindices,
        :abundancewindow,
        :mzindices,
        :mzvalues,
        :spectrumintensities,
        :spectrumweights,
        :fitionmask,
        :fitionpositions,
        :observationretentions,
        :normalizedretentions,
        :retentioncenter,
        :retentionscale,
        :observed,
        :observedraw,
        :baseline,
        :normalizedobserved,
        :normalizedobservedraw,
        :normalizedbaseline,
        :intensityscale,
    ])
    @test peakinput.scanindices == collect(18:24)
    @test peakinput.abundancewindow === result.abundanceinfo.windows[8][1]
    @test peakinput.fitionpositions == [1, 2, 3]
    @test peakinput.fitionmask == trues(3)
    @test peakinput.spectrumweights ≈ [1.0, 0.6, 0.3]
    @test peakinput.observedraw ≈ peakinput.observed .+ peakinput.baseline
    @test maximum(abs, peakinput.normalizedobserved) ≤ 1.0 + eps(Float64)
    finite_retentions = Float64[
        value for value in vec(peakinput.observationretentions)
        if isfinite(value)
    ]
    retention_left, retention_right = extrema(finite_retentions)
    @test peakinput.retentioncenter ≈ (retention_left + retention_right) / 2
    @test peakinput.retentionscale ≈
        maximum(abs.(finite_retentions .- peakinput.retentioncenter))
    @test maximum(abs, peakinput.normalizedretentions) ≤ 1.0 + eps(Float64)
    @test !hasproperty(peakinput, :knots)
    @test Set(keys(fit.peakfits)) == Set([8])
    @test isempty(fit.peakfitfailures)
    peakfit = fit.peakfits[8]
    @test peakfit.ladderstep == 8
    @test Set(propertynames(peakfit)) == Set([
        :ladderstep,
        :envelopegrid,
        :envelopevalues,
        :envelopeknots,
        :envelopecoefficients,
        :envelopemethod,
        :smoothnessfactor,
        :baselinemodel,
        :normalizedbaseline,
        :baseline,
        :normalizedfitted,
        :normalizedfittedraw,
        :normalizedresiduals,
        :fitted,
        :fittedraw,
        :residuals,
        :rmse,
        :normalizedrmse,
        :fitrmse,
        :normalizedfitrmse,
        :normalizedfitrss,
        :observationcount,
        :fitobservationcount,
        :fitionpositions,
        :fitionpositions_initial,
        :excludedfitionpositions,
        :replacementfitionpositions,
    ])
    @test peakfit.baselinemodel === :linear
    @test size(peakfit.fitted) == size(peakinput.observed)
    @test size(peakfit.residuals) == size(peakinput.observed)
    @test size(peakfit.normalizedfitted) == size(peakinput.observed)
    @test size(peakfit.baseline) == size(peakinput.observedraw)
    @test size(peakfit.fittedraw) == size(peakinput.observedraw)
    @test all(isfinite, peakfit.fitted)
    @test all(isfinite, peakfit.baseline)
    @test all(isfinite, peakfit.fittedraw)
    @test all(≥(0), peakfit.normalizedfitted)
    @test peakfit.fitted ≈ peakfit.normalizedfitted .* peakinput.intensityscale
    @test peakfit.fittedraw ≈ peakfit.fitted .+ peakfit.baseline
    @test peakfit.residuals ≈ peakinput.observedraw .- peakfit.fittedraw
    @test maximum(abs, peakfit.baseline .- 10.0) < 3.0
    @test peakfit.fitobservationcount == count(peakinput.fitionmask) * length(peakinput.scanindices)
    @test peakfit.observationcount == length(peakinput.observed)
    @test all(isfinite, peakfit.envelopevalues)
    @test all(≥(0), peakfit.envelopevalues)
    @test peakfit.envelopemethod === :curvature_adaptive_penalized_kernel
    @test peakfit.smoothnessfactor == JuChrom.ALKANE_VARIANCE_PEAK_SMOOTHNESS_FACTOR
    @test first(peakfit.envelopeknots) ≈ first(peakfit.envelopegrid)
    @test last(peakfit.envelopeknots) ≈ last(peakfit.envelopegrid)
    @test length(peakfit.envelopegrid) > length(peakfit.envelopeknots)
    @test all(>(0), diff(peakfit.envelopeknots))
    @test isfinite(peakfit.rmse)
    @test isfinite(peakfit.fitrmse)
    @test isfinite(peakfit.normalizedfitrss)
    @test peakfit.fitionpositions_initial == [1, 2, 3]
    @test peakfit.fitionpositions == [1, 2, 3]
    @test peakfit.excludedfitionpositions == Int[]
    @test peakfit.replacementfitionpositions == Int[]
    @test fit.peakqc.status === :too_few_peaks_for_robust_qc
    @test fit.peakqc.includedladdersteps == [8]
    @test fit.peakqc.excludedladdersteps == Int[]
    @test all(row -> Set(propertynames(row)) == Set([
        :runindex,
        :ladderstep,
        :scan_count,
        :ion_count,
        :initial_fition_count,
        :final_fition_count,
        :replacement_count,
        :replacement_fraction,
        :rmse,
        :normalizedrmse,
        :fitrmse,
        :normalizedfitrmse,
        :robustz,
        :exclude,
        :reason,
        :reasons,
    ]), fit.peakqc.rows)
    @test isfinite(fit.peakqc.robustcenter)
    @test isfinite(fit.peakqc.robustscale)
    peakqcrow = only(fit.peakqc.rows)
    @test peakqcrow.ladderstep == 8
    @test !peakqcrow.exclude
    @test peakqcrow.reason === :none
    @test peakqcrow.initial_fition_count == 3
    @test peakqcrow.final_fition_count == 3
    @test peakqcrow.replacement_fraction == 0.0
    @test length(fit.residualrecords) == length(peakinput.observed)
    @test all(record -> Set(propertynames(record)) == Set([
        :runindex,
        :ladderstep,
        :ionposition,
        :mzindex,
        :mzvalue,
        :scanposition,
        :scanindex,
        :retention,
        :normalizedretention,
        :spectrumweight,
        :isfition,
        :isinitialfition,
        :fittedsignal,
        :observedsignal,
        :baseline,
        :fittedintensity,
        :observedintensity,
        :residual,
        :residual2,
        :normalizedfitted,
        :normalizedresidual,
        :fittedsignalslope,
        :fittedsignalslope2,
    ]), fit.residualrecords)
    @test all(record -> record.ladderstep == 8, fit.residualrecords)
    @test all(record -> record.isfition, fit.residualrecords)
    @test all(record -> record.isinitialfition, fit.residualrecords)
    @test all(record -> record.fittedintensity ≈
        record.fittedsignal + record.baseline, fit.residualrecords)
    @test all(record -> record.observedintensity ≈
        record.observedsignal + record.baseline, fit.residualrecords)
    @test all(record -> record.residual ≈
        record.observedintensity - record.fittedintensity, fit.residualrecords)
    @test all(record -> record.residual2 ≈ abs2(record.residual), fit.residualrecords)
    @test all(record -> isfinite(record.fittedsignalslope), fit.residualrecords)

    full_window = AlkaneAbundanceWindow(
        8,
        firstindex(rawretentions(raw)),
        result.abundanceinfo.windows[8][1].apexindex,
        lastindex(rawretentions(raw)),
        0.0,
        maximum(result.abundanceinfo.abundances[8]),
        0.0,
        0.0,
        :test,
        :test,
    )
    noflatresult = _alkane_variance_result_with_abundance_window(result, full_window)
    noflatfit = fitalkanevariancemodel(raw, noflatresult)
    @test noflatfit isa AlkaneVarianceFit
    @test !noflatfit.success
    @test noflatfit.status === :failed
    @test isnothing(noflatfit.model)
    @test noflatfit.qc.accept === false
    @test :variance_model_fit_failed in noflatfit.qc.reasons
    @test_throws ArgumentError varpred(8.0, noflatfit)
    noflatdisplay = sprint(io -> show(io, MIME"text/plain"(), noflatfit))
    @test occursin("success: false", noflatdisplay)
    @test occursin("status: failed", noflatdisplay)
    @test occursin("reasons:", noflatdisplay)

    narrow_bandwidth = JuChrom.alkane_variance_peak_local_quadratic_bandwidth(
        [0.0, 1.0, 2.0],
        [1, 2, 3],
    )
    broad_bandwidth = JuChrom.alkane_variance_peak_local_quadratic_bandwidth(
        collect(0.0:1.0:8.0),
        collect(1:9),
    )
    @test broad_bandwidth > narrow_bandwidth

    twoionfit = fitalkanevariancemodel(raw, result; fitioncount=2)
    @test twoionfit.success
    twoioninput = twoionfit.peakinputs[8]
    @test twoioninput.fitionpositions == [1, 2]
    @test twoioninput.fitionmask == [true, true, false]
    @test twoionfit.peakfits[8].fitionpositions == [1, 2]
    @test twoionfit.peakfits[8].fitionpositions_initial == [1, 2]

    fakefit = (normalizedresiduals=[
        0.1 0.1 0.1
        0.1 0.1 0.1
        10.0 10.0 10.0
    ],)
    @test JuChrom.alkane_variance_peak_fit_outlier_rows(fakefit, [1, 2, 3]) == [3]
    fakepeakinput = (spectrumweights=[1.0, 0.8, 0.6, 0.4],)
    finalrows, excludedrows, replacementrows =
        JuChrom.alkane_variance_peak_replacement_fitrows(
            fakepeakinput,
            fakefit,
            [1, 2, 3],
        )
    @test finalrows == [1, 2, 4]
    @test excludedrows == [3]
    @test replacementrows == [4]

    badmsm = MassScanMatrix(
        rawretentions(raw),
        retentionunit(raw),
        rawmzvalues(raw),
        mzunit(raw),
        rawintensities(raw) .+ 1.0,
        intensityunit(raw),
    )
    @test_throws ArgumentError fitalkanevariancemodel(badmsm, result)

    msm, emptyresult = _alkane_variance_entrypoint_inputs()
    @test_throws ArgumentError fitalkanevariancemodel(msm, emptyresult)
    @test_throws MethodError fitalkanevariancemodel([(raw, raw)])
    @test_throws MethodError fitalkanevariancemodel(
        Tuple{MassScanMatrix, AlkaneSeriesResult}[],
    )
    @test_throws ArgumentError fitalkanevariancemodel(raw, result; laddersteps=[8])
    @test_throws ArgumentError fitalkanevariancemodel(raw, result; excludeladdersteps=[8])
    @test_throws ArgumentError fitalkanevariancemodel(raw, result; excludeladdersteps=[9])
    @test_throws ArgumentError fitalkanevariancemodel(raw, result; excludeladdersteps=[8.0])
    @test_throws ArgumentError fitalkanevariancemodel(raw, result; fitioncount=0)
    @test_throws ArgumentError fitalkanevariancemodel(
        raw,
        result;
        peakbaselinemodel=:linear,
    )
    @test_throws MethodError fitalkanevariancemodel([
        (raw, result, (; laddersteps=[8])),
    ])
    badpeakinput = merge(peakinput, (fitionmask=falses(length(peakinput.fitionmask)),))
    @test_throws ArgumentError JuChrom.alkane_variance_fit_peak_envelope(badpeakinput)
end
