"""
    fitalkanevariancemodel(msm, result; excludeladdersteps=(), fitioncount=5)
    fitalkanevariancemodel(runs; excludeladdersteps=(), fitioncount=5)

Entry point for fitting an alkane-ladder-based variance model.

`result` must be the [`AlkaneSeriesResult`](@ref) produced from `msm`. The function
validates that the mass-scan matrix and alkane-ladder result belong together before
fitting proceeds.

The one-argument form accepts any iterable whose elements are `(msm, result)` tuples or
`msm => result` pairs. Run-specific exclusions can be supplied as
`(msm, result, (; excludeladdersteps=[23]))`. Global `excludeladdersteps` are combined
with run-specific exclusions. Ladder detection is never run implicitly.

The implementation fits one curvature-adaptive smooth nonnegative peak envelope per
ladder step using the selected fit ions, allows a small internal apex readjustment, and
flags poor peak fits by within-run robust QC. Peaks are fit on raw intensities with one
ion-specific linear baseline over the peak window. The baseline stored in `result` is used
only as the source of a soft baseline anchor, not as the final fitted baseline.

The returned object contains a fitted [`LinearObservedIntensityVarianceModel`](@ref) in
`model`, per-run details in `runs`, and minimal run-level quality control in `qc`. The
model is calibrated as

    variance(I) = a_flat + b * max(I - I_flat, 0)

where `a_flat` and `I_flat` are estimated from flat non-peak regions and `b` is estimated
robustly from accepted alkane peak residuals.
"""
function fitalkanevariancemodel(
    msm::MassScanMatrix,
    result::AlkaneSeriesResult;
    excludeladdersteps=(),
    fitioncount::Integer=5,
    kwargs...
)
    alkane_variance_reject_unimplemented_keywords(kwargs)
    alkane_variance_validate_fitioncount(fitioncount)
    global_exclusions = alkane_variance_ladder_step_vector(
        excludeladdersteps,
        "excludeladdersteps",
    )
    settings = (
        fitioncount=Int(fitioncount),
    )
    validated = [alkane_variance_validated_run(msm, result, global_exclusions, Int[])]
    validated = [merge(run, (settings=settings,)) for run in validated]
    alkane_variance_fit_validated_runs(validated)
end

function fitalkanevariancemodel(
    runs;
    excludeladdersteps=(),
    fitioncount::Integer=5,
    kwargs...
)
    alkane_variance_reject_unimplemented_keywords(kwargs)
    alkane_variance_validate_fitioncount(fitioncount)
    global_exclusions = alkane_variance_ladder_step_vector(
        excludeladdersteps,
        "excludeladdersteps",
    )
    validated = alkane_variance_validated_runs(runs, global_exclusions)
    settings = (
        fitioncount=Int(fitioncount),
    )
    validated = [merge(run, (settings=settings,)) for run in validated]
    alkane_variance_fit_validated_runs(validated)
end

function alkane_variance_reject_unimplemented_keywords(kwargs)
    isempty(kwargs) && return nothing
    names = join(string.(collect(keys(kwargs))), ", ")
    throw(ArgumentError(
        "fitalkanevariancemodel does not accept keyword arguments yet: $names"))
end

function alkane_variance_validate_fitioncount(fitioncount::Integer)
    fitioncount ≥ 1 || throw(ArgumentError("fitioncount must be at least 1"))
    nothing
end

function alkane_variance_ladder_step_vector(value, name::AbstractString)
    if isnothing(value)
        return Int[]
    elseif value isa Integer
        return [Int(value)]
    end

    steps = Int[]
    for step in value
        step isa Integer || throw(ArgumentError("$name must contain integer ladder steps"))
        push!(steps, Int(step))
    end
    unique!(steps)
    sort!(steps)
    steps
end

function alkane_variance_validated_runs(runs, global_exclusions::AbstractVector{<:Integer})
    validated = NamedTuple[]
    for run in runs
        msm, result, run_exclusions = alkane_variance_run_pair(run)
        push!(validated, alkane_variance_validated_run(
            msm,
            result,
            global_exclusions,
            run_exclusions,
        ))
    end
    !isempty(validated) || throw(ArgumentError(
        "at least one (msm, result) run pair is required"))

    validated
end

function alkane_variance_run_pair(run)
    if run isa Pair
        msm = first(run)
        result = last(run)
        run_exclusions = Int[]
    elseif run isa Tuple && length(run) == 2
        msm = run[1]
        result = run[2]
        run_exclusions = Int[]
    elseif run isa Tuple && length(run) == 3
        msm = run[1]
        result = run[2]
        run_exclusions = alkane_variance_run_exclusions(run[3])
    else
        throw(ArgumentError(
            "runs must yield (msm, result) tuples, msm => result pairs, or " *
            "(msm, result, (; excludeladdersteps=[...])) tuples"))
    end

    msm isa MassScanMatrix || throw(ArgumentError(
        "each variance calibration run must contain a MassScanMatrix"))
    result isa AlkaneSeriesResult || throw(ArgumentError(
        "each variance calibration run must contain an AlkaneSeriesResult"))

    msm, result, run_exclusions
end

function alkane_variance_run_exclusions(value)
    if value isa NamedTuple
        unknown = setdiff(collect(keys(value)), [:excludeladdersteps])
        isempty(unknown) || throw(ArgumentError(
            "unknown run metadata field(s): " * join(string.(unknown), ", ")))
        return alkane_variance_ladder_step_vector(
            get(value, :excludeladdersteps, ()),
            "run-specific excludeladdersteps",
        )
    end

    alkane_variance_ladder_step_vector(value, "run-specific excludeladdersteps")
end

function alkane_variance_validated_run(
    msm::MassScanMatrix,
    result::AlkaneSeriesResult,
    global_exclusions::AbstractVector{<:Integer},
    run_exclusions::AbstractVector{<:Integer},
)
    signal = alkane_ladder_extraction_signal(msm, result, true)
    excluded = alkane_variance_merge_ladder_step_exclusions(
        global_exclusions,
        run_exclusions,
    )
    available = alkane_variance_available_laddersteps(result)
    included = alkane_variance_included_laddersteps(available, excluded)

    (
        msm=msm,
        result=result,
        signal=signal,
        excludeladdersteps=excluded,
        includedladdersteps=included,
    )
end

function alkane_variance_merge_ladder_step_exclusions(global_exclusions, run_exclusions)
    excluded = Int[]
    append!(excluded, Int.(global_exclusions))
    append!(excluded, Int.(run_exclusions))
    unique!(excluded)
    sort!(excluded)
    excluded
end

function alkane_variance_available_laddersteps(result::AlkaneSeriesResult)
    if isnothing(result.apexinfo) && isnothing(result.additioninfo)
        return Int[]
    end

    sort!([step.ladderstep for step in alkaneladdersteps(result)])
end

function alkane_variance_included_laddersteps(available, excluded)
    isempty(available) && throw(ArgumentError(
        "alkane series result contains no refined ladder steps"))
    missing = sort!(collect(setdiff(Set(excluded), Set(available))))
    isempty(missing) || throw(ArgumentError(
        "excluded ladder step(s) are not present in result: " *
        join("C" .* string.(missing), ", ")))

    included = sort!(collect(setdiff(Set(available), Set(excluded))))
    if !isempty(available) && isempty(included)
        throw(ArgumentError("excludeladdersteps removes all available ladder steps"))
    end

    included
end

function alkane_variance_fit_validated_runs(validated_runs)
    prepared_runs = [
        alkane_variance_prepare_peak_inputs(run, runindex)
        for (runindex, run) in pairs(validated_runs)
    ]
    runs = alkane_variance_fit_run_variance_models(prepared_runs)
    accepted = [
        run for run in runs
        if !isnothing(run.model) && run.varianceqc.accept
    ]
    if isempty(accepted)
        reasons = [
            "run $(run.runindex): " *
            join(string.(run.varianceqc.reasons), ", ")
            for run in runs
        ]
        throw(ArgumentError(
            "no variance calibration run passed QC" *
            (isempty(reasons) ? "" : ": " * join(reasons, "; "))))
    end
    model = alkane_variance_aggregate_models(accepted)

    (
        status=length(accepted) == length(runs) ? :ok : :partial,
        model=model,
        runs=runs,
        qc=alkane_variance_fit_qc(runs),
    )
end

function alkane_variance_prepare_peak_inputs(run, runindex::Integer)
    with_spectra = alkane_variance_extract_mass_spectra(run, runindex)
    peakinputs, peakinputfailures = alkane_variance_peak_inputs(with_spectra)
    peakfits, peakfitfailures = alkane_variance_peak_envelope_fits(peakinputs)
    peakqc = alkane_variance_peak_qc(runindex, peakinputs, peakfits)
    residualrecords = alkane_variance_residual_records(
        runindex,
        peakinputs,
        peakfits,
        peakqc,
    )

    merge(with_spectra, (
        peakinputs=peakinputs,
        peakinputfailures=peakinputfailures,
        peakfits=peakfits,
        peakfitfailures=peakfitfailures,
        peakqc=peakqc,
        residualrecords=residualrecords,
    ))
end

function alkane_variance_extract_mass_spectra(run, runindex::Integer)
    settings = alkane_ladder_mass_spectrum_settings(
        true,
        1.0,
        false,
        true,
        false,
        true,
        true,
        true,
    )
    stepmap = Dict(step.ladderstep => step for step in alkaneladdersteps(run.result))
    spectra = Dict{Int, AbstractMassSpectrum}()
    failures = Dict{Int, String}()

    for ladderstep in run.includedladdersteps
        step = stepmap[ladderstep]
        extraction = alkane_ladder_extract_step_mass_spectrum(
            run.signal,
            step,
            run.result.variances,
            settings,
            false,
        )
        if isnothing(extraction.failure)
            spectra[ladderstep] = extraction.spectrum
        else
            failures[ladderstep] = extraction.failure
        end
    end

    merge(run, (
        runindex=Int(runindex),
        spectra=spectra,
        failures=failures,
    ))
end

function alkane_variance_peak_inputs(run)
    stepmap = Dict(step.ladderstep => step for step in alkaneladdersteps(run.result))
    inputs = Dict{Int, NamedTuple}()
    failures = Dict{Int, String}()

    for ladderstep in run.includedladdersteps
        haskey(run.spectra, ladderstep) || continue
        try
            inputs[ladderstep] = alkane_variance_peak_input(
                run,
                stepmap[ladderstep],
                run.spectra[ladderstep],
            )
        catch err
            failures[ladderstep] = sprint(showerror, err)
        end
    end

    inputs, failures
end

function alkane_variance_peak_input(run, step::AlkaneLadderStep, spectrum)
    scanindices, abundancewindow = alkane_variance_peak_input_scanindices(
        run.signal,
        run.result,
        step,
    )
    mzindices, mzvalues, spectrumintensities, spectrumweights, fitionmask =
        alkane_variance_peak_input_ions(run.signal, spectrum, run.settings.fitioncount)
    observationretentions = alkane_variance_peak_observation_retentions(
        run.signal,
        run.result,
        step,
        scanindices,
        mzindices,
    )
    normalizedretentions, retentioncenter, retentionscale =
        alkane_variance_normalized_retentions(observationretentions, step.apexretention)

    observed = permutedims(Float64.(rawintensities(run.signal)[scanindices, mzindices]))
    observedraw = permutedims(Float64.(rawintensities(run.msm)[scanindices, mzindices]))
    baseline = observedraw .- observed
    intensityscale = alkane_variance_intensity_scale(observed, observedraw)

    (
        ladderstep=step.ladderstep,
        scanindices=scanindices,
        abundancewindow=abundancewindow,
        apexscanindex=step.apexscanindex,
        apexretention=step.apexretention,
        mzindices=mzindices,
        mzvalues=mzvalues,
        spectrumintensities=spectrumintensities,
        spectrumweights=spectrumweights,
        fitionmask=fitionmask,
        fitionpositions=findall(fitionmask),
        observationretentions=observationretentions,
        normalizedretentions=normalizedretentions,
        retentioncenter=retentioncenter,
        retentionscale=retentionscale,
        observed=observed,
        observedraw=observedraw,
        baseline=baseline,
        normalizedobserved=observed ./ intensityscale,
        normalizedobservedraw=observedraw ./ intensityscale,
        normalizedbaseline=baseline ./ intensityscale,
        intensityscale=intensityscale,
    )
end

function alkane_variance_peak_input_scanindices(
    signal::MassScanMatrix,
    result::AlkaneSeriesResult,
    step::AlkaneLadderStep,
)
    result.abundanceinfo isa AlkaneAbundanceInfo || throw(ArgumentError(
        "alkane series result does not contain abundance windows"))
    haskey(result.abundanceinfo.windows, step.ladderstep) || throw(ArgumentError(
        "no abundance window is available for C$(step.ladderstep)"))
    windows = result.abundanceinfo.windows[step.ladderstep]
    !isempty(windows) || throw(ArgumentError(
        "no abundance window is available for C$(step.ladderstep)"))

    target = round(Int, step.apexscanindex)
    candidates = [
        index for index in eachindex(windows)
        if windows[index].leftindex ≤ target ≤ windows[index].rightindex
    ]
    isempty(candidates) && append!(candidates, eachindex(windows))
    distances = [abs(windows[index].apexindex - target) for index in candidates]
    window = windows[candidates[argmin(distances)]]
    left = clamp(window.leftindex, 1, scancount(signal))
    right = clamp(window.rightindex, 1, scancount(signal))
    left ≤ right || throw(ArgumentError(
        "abundance window has invalid scan boundaries for C$(step.ladderstep)"))

    collect(left:right), window
end

function alkane_variance_peak_input_ions(
    signal::MassScanMatrix,
    spectrum,
    fitioncount::Integer,
)
    msm_mzvalues = collect(rawmzvalues(signal))
    mzindices = Int[]
    mzvalues_selected = Float64[]
    spectrumintensities = Float64[]

    for (mz, intensity) in zip(spectrum.mzvalues, spectrum.intensities)
        mzindex = findfirst(==(mz), msm_mzvalues)
        isnothing(mzindex) && continue
        value = Float64(intensity)
        isfinite(value) && value > 0 || continue
        push!(mzindices, mzindex)
        push!(mzvalues_selected, Float64(msm_mzvalues[mzindex]))
        push!(spectrumintensities, value)
    end
    !isempty(mzindices) || throw(ArgumentError(
        "mass spectrum contains no positive ions on the signal m/z grid"))

    order = sortperm(spectrumintensities; rev=true)
    mzindices = mzindices[order]
    mzvalues_selected = mzvalues_selected[order]
    spectrumintensities = spectrumintensities[order]

    weights = spectrumintensities ./ maximum(spectrumintensities)
    fitionmask = falses(length(mzindices))
    fitionmask[1:min(fitioncount, length(fitionmask))] .= true

    mzindices, mzvalues_selected, spectrumintensities, weights, fitionmask
end

function alkane_variance_peak_observation_retentions(
    signal::MassScanMatrix,
    result::AlkaneSeriesResult,
    step::AlkaneLadderStep,
    scanindices::AbstractVector{<:Integer},
    mzindices::AbstractVector{<:Integer},
)
    raw_scan_retentions = collect(rawretentions(signal))
    mzkwargs = alkane_variance_peak_mzretention_kwargs(result, step)
    retentions = Matrix{Float64}(undef, length(mzindices), length(scanindices))
    for (i, mzindex) in pairs(mzindices)
        for (j, scanindex) in pairs(scanindices)
            retentions[i, j] = alkane_ladder_observation_retention(
                raw_scan_retentions,
                scanindex,
                mzindex,
                mzkwargs,
            )
        end
    end

    retentions
end

function alkane_variance_peak_mzretention_kwargs(
    result::AlkaneSeriesResult,
    step::AlkaneLadderStep,
)
    fitresult = step.apex.fit
    if !isnothing(fitresult) &&
            hasproperty(fitresult, :fit) &&
            hasproperty(fitresult.fit, :mzretentionkwargs)
        return fitresult.fit.mzretentionkwargs
    end
    if !isnothing(result.apexinfo) && !isnothing(result.apexinfo.settings.mzretentionkwargs)
        return result.apexinfo.settings.mzretentionkwargs
    end

    throw(ArgumentError(
        "m/z retention timing information is missing for C$(step.ladderstep)"))
end

function alkane_variance_normalized_retentions(
    retentions::AbstractMatrix{<:Real},
    center::Real,
)
    scale = maximum(abs.(retentions .- center))
    isfinite(scale) && scale > 0 || throw(ArgumentError(
        "retention normalization scale must be finite and positive"))

    (retentions .- center) ./ scale, Float64(center), Float64(scale)
end

function alkane_variance_intensity_scale(observed, observedraw)
    signal_values = abs.(vec(observed))
    finite_signal = signal_values[isfinite.(signal_values)]
    signal_scale = isempty(finite_signal) ? 0.0 : maximum(finite_signal)
    signal_scale > 0 && return signal_scale

    raw_values = abs.(vec(observedraw))
    finite_raw = raw_values[isfinite.(raw_values)]
    raw_scale = isempty(finite_raw) ? 0.0 : maximum(finite_raw)
    raw_scale > 0 && return raw_scale

    1.0
end

function alkane_variance_peak_envelope_fits(peakinputs)
    fits = Dict{Int, NamedTuple}()
    failures = Dict{Int, String}()

    for ladderstep in sort!(collect(keys(peakinputs)))
        try
            fits[ladderstep] = alkane_variance_fit_peak_envelope(peakinputs[ladderstep])
        catch err
            failures[ladderstep] = sprint(showerror, err)
        end
    end

    fits, failures
end

function alkane_variance_peak_qc(runindex::Integer, peakinputs, peakfits)
    laddersteps = sort!(collect(intersect(Set(keys(peakinputs)), Set(keys(peakfits)))))
    rows = NamedTuple[]
    for ladderstep in laddersteps
        peakinput = peakinputs[ladderstep]
        peakfit = peakfits[ladderstep]
        initialcount = hasproperty(peakfit, :fitionpositions_initial) ?
            length(peakfit.fitionpositions_initial) :
            count(peakinput.fitionmask)
        replacementcount = hasproperty(peakfit, :replacementfitionpositions) ?
            length(peakfit.replacementfitionpositions) :
            0
        replacementfraction = initialcount == 0 ? 0.0 : replacementcount / initialcount
        push!(rows, (
            runindex=Int(runindex),
            ladderstep=Int(ladderstep),
            scan_count=length(peakinput.scanindices),
            ion_count=length(peakinput.mzindices),
            initial_fition_count=Int(initialcount),
            final_fition_count=hasproperty(peakfit, :fitionpositions) ?
                length(peakfit.fitionpositions) :
                count(peakinput.fitionmask),
            replacement_count=Int(replacementcount),
            replacement_fraction=Float64(replacementfraction),
            rmse=Float64(peakfit.rmse),
            normalizedrmse=Float64(peakfit.normalizedrmse),
            fitrmse=Float64(peakfit.fitrmse),
            normalizedfitrmse=Float64(peakfit.normalizedfitrmse),
            apexscanshift=Float64(peakfit.apexscanshift),
        ))
    end

    fitrmses = [
        row.normalizedfitrmse for row in rows
        if isfinite(row.normalizedfitrmse)
    ]
    center, scale = alkane_variance_robust_center_scale(fitrmses)
    minpeakcount = 5
    robustzthreshold = 4.0
    robust_qc_active = length(fitrmses) ≥ minpeakcount
    qcrow = NamedTuple[]
    for row in rows
        z = isfinite(row.normalizedfitrmse) && isfinite(center) && scale > 0 ?
            (row.normalizedfitrmse - center) / scale :
            Inf
        reasons = Symbol[]
        isfinite(row.normalizedfitrmse) || push!(reasons, :nonfinite_fit_rmse)
        if robust_qc_active && z > robustzthreshold
            push!(reasons, :high_fit_rmse)
        end
        if row.replacement_fraction > 0.5
            push!(reasons, :many_fit_ion_replacements)
        end
        push!(qcrow, merge(row, (
            robustz=Float64(z),
            exclude=!isempty(reasons),
            reason=isempty(reasons) ? :none :
                length(reasons) == 1 ? first(reasons) : :multiple_qc_failures,
            reasons=Tuple(reasons),
        )))
    end

    excluded = sort!(Int[row.ladderstep for row in qcrow if row.exclude])
    included = sort!(Int[row.ladderstep for row in qcrow if !row.exclude])
    (
        status=isempty(qcrow) ? :empty :
            robust_qc_active ? :ok : :too_few_peaks_for_robust_qc,
        minpeakcount=minpeakcount,
        robustqcisactive=robust_qc_active,
        robustcenter=Float64(center),
        robustscale=Float64(scale),
        robustzthreshold=robustzthreshold,
        rows=qcrow,
        includedladdersteps=included,
        excludedladdersteps=excluded,
    )
end

function alkane_variance_robust_center_scale(values)
    finite = Float64[value for value in values if isfinite(value)]
    if isempty(finite)
        return NaN, NaN
    end

    center = median(finite)
    madscale = 1.4826 * median(abs.(finite .- center))
    iqrscale = if length(finite) ≥ 4
        (quantile(finite, 0.75) - quantile(finite, 0.25)) / 1.349
    else
        0.0
    end
    floor = sqrt(eps(Float64)) * max(abs(center), 1.0)
    scale = max(madscale, iqrscale, floor)

    Float64(center), Float64(scale)
end

function alkane_variance_residual_records(
    runindex::Integer,
    peakinputs,
    peakfits,
    peakqc,
)
    included = Set(Int.(peakqc.includedladdersteps))
    records = NamedTuple[]
    laddersteps = sort!(collect(intersect(Set(keys(peakinputs)), Set(keys(peakfits)))))
    for ladderstep in laddersteps
        ladderstep in included || continue
        append!(records, alkane_variance_peak_residual_records(
            runindex,
            peakinputs[ladderstep],
            peakfits[ladderstep],
        ))
    end

    records
end

function alkane_variance_peak_residual_records(runindex::Integer, peakinput, peakfit)
    records = NamedTuple[]
    fitionpositions = hasproperty(peakfit, :fitionpositions) ?
        Set(Int.(peakfit.fitionpositions)) :
        Set(Int.(findall(peakinput.fitionmask)))
    initialfitionpositions = hasproperty(peakfit, :fitionpositions_initial) ?
        Set(Int.(peakfit.fitionpositions_initial)) :
        Set(Int.(findall(peakinput.fitionmask)))

    for ionposition in axes(peakinput.observed, 1)
        spectrumweight = Float64(peakinput.spectrumweights[ionposition])
        for scanposition in axes(peakinput.observed, 2)
            fittedsignal = Float64(peakfit.fitted[ionposition, scanposition])
            baseline = Float64(
                hasproperty(peakfit, :baseline) ?
                    peakfit.baseline[ionposition, scanposition] :
                    peakinput.baseline[ionposition, scanposition]
            )
            fittedintensity = fittedsignal + baseline
            observedintensity = Float64(peakinput.observedraw[ionposition, scanposition])
            observedsignal = observedintensity - baseline
            residual = Float64(peakfit.residuals[ionposition, scanposition])
            normalizedfitted = Float64(peakfit.normalizedfitted[ionposition, scanposition])
            normalizedresidual =
                Float64(peakfit.normalizedresiduals[ionposition, scanposition])
            retention = Float64(peakinput.observationretentions[ionposition, scanposition])
            normalizedretention =
                Float64(peakinput.normalizedretentions[ionposition, scanposition])
            fittedsignalslope = alkane_variance_peak_fitted_signal_slope(
                peakinput,
                peakfit,
                ionposition,
                normalizedretention,
            )

            all(isfinite, (
                fittedsignal,
                baseline,
                fittedintensity,
                observedsignal,
                observedintensity,
                residual,
                normalizedfitted,
                normalizedresidual,
                retention,
                normalizedretention,
                spectrumweight,
            )) || continue

            push!(records, (
                runindex=Int(runindex),
                ladderstep=Int(peakinput.ladderstep),
                ionposition=Int(ionposition),
                mzindex=Int(peakinput.mzindices[ionposition]),
                mzvalue=Float64(peakinput.mzvalues[ionposition]),
                scanposition=Int(scanposition),
                scanindex=Int(peakinput.scanindices[scanposition]),
                retention=retention,
                normalizedretention=normalizedretention,
                adjustednormalizedretention=
                    normalizedretention - Float64(peakfit.normalizedapexshift),
                apexretention=Float64(peakinput.apexretention),
                fittedapexretention=Float64(peakfit.fittedapexretention),
                spectrumweight=spectrumweight,
                isfition=Int(ionposition) in fitionpositions,
                isinitialfition=Int(ionposition) in initialfitionpositions,
                fittedsignal=fittedsignal,
                observedsignal=observedsignal,
                baseline=baseline,
                fittedintensity=fittedintensity,
                observedintensity=observedintensity,
                residual=residual,
                residual2=abs2(residual),
                normalizedfitted=normalizedfitted,
                normalizedresidual=normalizedresidual,
                fittedsignalslope=fittedsignalslope,
                fittedsignalslope2=abs2(fittedsignalslope),
            ))
        end
    end

    records
end

function alkane_variance_fit_run_variance_models(runs)
    [
        alkane_variance_fit_run_variance_model(run)
        for run in runs
    ]
end

function alkane_variance_fit_run_variance_model(run)
    try
        flatrecords, flatqc = alkane_variance_flat_nonpeak_residual_records(
            run,
            run.residualrecords,
        )
        variancefit = alkane_variance_fit_linear_observed_model(
            run.residualrecords,
            flatrecords,
        )
        rho = alkane_variance_lag1_autocorrelation(
            flatrecords,
            variancefit.model,
        )
        model = LinearObservedIntensityVarianceModel(
            variancefit.model.intercept,
            variancefit.model.slope,
            variancefit.model.intensity_offset,
            variancefit.model.intensity_min,
            variancefit.model.intensity_max,
            rho.rho,
        )
        variancefit = merge(variancefit, (
            model=model,
            lag1=rho,
        ))
        varianceqc = alkane_variance_run_qc(run, flatqc, variancefit)

        merge(run, (
            flatrecords=flatrecords,
            flatqc=flatqc,
            variancefit=variancefit,
            model=model,
            varianceqc=varianceqc,
        ))
    catch err
        merge(run, (
            flatrecords=NamedTuple[],
            flatqc=alkane_variance_empty_flat_qc(),
            variancefit=nothing,
            model=nothing,
            varianceqc=(
                accept=false,
                status=:failed,
                reasons=(:variance_model_fit_failed,),
                error=sprint(showerror, err),
                flatrecordcount=0,
                peakrecordcount=length(run.residualrecords),
                keptpeakcount=length(run.peakqc.includedladdersteps),
                excludedpeakcount=length(run.peakqc.excludedladdersteps),
                excludedfraction=alkane_variance_excluded_peak_fraction(run.peakqc),
                intensity_min=NaN,
                intensity_max=NaN,
            ),
        ))
    end
end

function alkane_variance_flat_nonpeak_residual_records(
    run,
    peakrecords;
    topn::Integer=20,
    minrecordsperion::Integer=60,
    peakmargin_scans::Integer=5,
    flatwindowsize::Integer=17,
    flatwindowstride=nothing,
    flat_slope_quantile::Real=0.35,
    flat_curvature_quantile::Real=0.50,
    maxwindowsperion::Integer=250,
)
    X = Float64.(rawintensities(run.msm))
    nscans, _ = size(X)
    retentions = Float64.(rawretentions(run.msm))
    length(retentions) == nscans || throw(DimensionMismatch(
        "raw retention count must match the number of scans"))
    excluded = alkane_variance_flat_nonpeak_exclusion_mask(
        run.result,
        nscans;
        margin_scans=peakmargin_scans,
    )
    mzindices, mzvalues = alkane_variance_flat_nonpeak_mzindices(
        run.msm,
        peakrecords;
        topn=topn,
        minrecords=minrecordsperion,
    )
    !isempty(mzindices) || throw(ArgumentError(
        "no m/z channels with enough accepted peak records are available for " *
        "flat-region variance estimation"))

    maxusable = alkane_variance_max_contiguous_true(.!excluded)
    windowsize = min(max(5, Int(flatwindowsize)), maxusable)
    windowsize ≥ 5 || throw(ArgumentError(
        "not enough contiguous flat non-peak scans are available"))
    stride = isnothing(flatwindowstride) ? windowsize : max(1, Int(flatwindowstride))
    slopeq = clamp(Float64(flat_slope_quantile), 0.01, 0.95)
    curvatureq = clamp(Float64(flat_curvature_quantile), 0.01, 0.95)

    records = NamedTuple[]
    windowcount = 0
    usableioncount = 0
    for (mzindex, mz) in zip(mzindices, mzvalues)
        y = X[:, mzindex]
        usable = .!excluded .& isfinite.(y) .& isfinite.(retentions)
        sum(usable) ≥ windowsize || continue

        adjacent = usable[1:(end - 1)] .& usable[2:end]
        diffs = abs.(diff(y))
        diffvalues = diffs[adjacent]
        isempty(diffvalues) && continue
        slope_threshold = quantile(diffvalues, slopeq)

        curvature_threshold = Inf
        if nscans ≥ 3
            triples = usable[1:(end - 2)] .&
                usable[2:(end - 1)] .&
                usable[3:end]
            curvatures = abs.(y[3:end] .- 2 .* y[2:(end - 1)] .+ y[1:(end - 2)])
            curvvalues = curvatures[triples]
            isempty(curvvalues) || (curvature_threshold = quantile(curvvalues, curvatureq))
        end

        accepted_for_ion = 0
        for firstscan in 1:stride:(nscans - windowsize + 1)
            accepted_for_ion ≥ maxwindowsperion && break
            scanrange = firstscan:(firstscan + windowsize - 1)
            all(@view usable[scanrange]) || continue

            local_y = Float64.(@view y[scanrange])
            local_x = Float64.(@view retentions[scanrange])
            median(abs.(diff(local_y))) ≤ slope_threshold || continue
            if length(local_y) ≥ 3 && isfinite(curvature_threshold)
                local_curvature = abs.(
                    local_y[3:end] .-
                    2 .* local_y[2:(end - 1)] .+
                    local_y[1:(end - 2)]
                )
                median(local_curvature) ≤ curvature_threshold || continue
            end

            fitted, slope = alkane_variance_flat_linear_baseline(local_x, local_y)
            all(isfinite, fitted) && isfinite(slope) || continue
            for (localposition, scanindex) in pairs(scanrange)
                residual = local_y[localposition] - fitted[localposition]
                isfinite(residual) || continue
                fittedintensity = fitted[localposition]
                observedintensity = local_y[localposition]
                push!(records, (
                    runindex=Int(run.runindex),
                    ladderstep=0,
                    ionposition=0,
                    mzindex=Int(mzindex),
                    mzvalue=Float64(mz),
                    scanposition=Int(localposition),
                    scanindex=Int(scanindex),
                    retention=local_x[localposition],
                    normalizedretention=NaN,
                    adjustednormalizedretention=NaN,
                    apexretention=NaN,
                    fittedapexretention=NaN,
                    spectrumweight=NaN,
                    isfition=true,
                    isinitialfition=true,
                    fittedsignal=0.0,
                    observedsignal=observedintensity - fittedintensity,
                    baseline=fittedintensity,
                    fittedintensity=fittedintensity,
                    observedintensity=observedintensity,
                    residual=residual,
                    residual2=abs2(residual),
                    normalizedfitted=NaN,
                    normalizedresidual=NaN,
                    fittedsignalslope=slope,
                    fittedsignalslope2=abs2(slope),
                ))
            end
            accepted_for_ion += 1
            windowcount += 1
        end
        accepted_for_ion > 0 && (usableioncount += 1)
    end
    !isempty(records) || throw(ArgumentError(
        "flat-region variance estimation found no accepted flat windows"))

    qc = (
        status=:ok,
        recordcount=length(records),
        windowcount=windowcount,
        ioncount=usableioncount,
        requestedioncount=length(mzindices),
        peakmargin_scans=Int(peakmargin_scans),
        flatwindowsize=windowsize,
        flatwindowstride=stride,
        flat_slope_quantile=slopeq,
        flat_curvature_quantile=curvatureq,
    )

    records, qc
end

function alkane_variance_flat_nonpeak_exclusion_mask(
    result,
    nscans::Integer;
    margin_scans::Integer=5,
)
    excluded = falses(nscans)
    if isnothing(result.abundanceinfo) || !hasproperty(result.abundanceinfo, :windows)
        return excluded
    end

    margin = max(0, Int(margin_scans))
    for stepwindows in values(result.abundanceinfo.windows)
        for window in stepwindows
            left = clamp(Int(window.leftindex) - margin, 1, nscans)
            right = clamp(Int(window.rightindex) + margin, 1, nscans)
            left ≤ right || continue
            excluded[left:right] .= true
        end
    end

    excluded
end

function alkane_variance_flat_nonpeak_mzindices(
    msm,
    peakrecords;
    topn::Integer=20,
    minrecords::Integer=60,
)
    adaptive_minrecords = min(Int(minrecords), max(3, length(peakrecords) ÷ 20))
    topions = alkane_variance_top_mzvalues_by_peak_signal(
        peakrecords;
        n=topn,
        minrecords=adaptive_minrecords,
    )
    if isempty(topions.mzvalues) && adaptive_minrecords > 1
        topions = alkane_variance_top_mzvalues_by_peak_signal(
            peakrecords;
            n=topn,
            minrecords=1,
        )
    end
    msm_mzvalues = Float64.(collect(rawmzvalues(msm)))
    mzindices = Int[]
    mzvalues = Float64[]
    for mz in topions.mzvalues
        index = findfirst(==(mz), msm_mzvalues)
        isnothing(index) && continue
        push!(mzindices, index)
        push!(mzvalues, msm_mzvalues[index])
    end

    mzindices, mzvalues
end

function alkane_variance_top_mzvalues_by_peak_signal(
    records;
    n::Integer=20,
    minrecords::Integer=30,
)
    mzvalues = Float64[]
    totalsignal = Float64[]
    counts = Int[]
    for mz in sort(unique(record.mzvalue for record in records))
        rows = [record for record in records if record.mzvalue == mz]
        length(rows) ≥ minrecords || continue
        push!(mzvalues, Float64(mz))
        push!(totalsignal, sum(max(record.fittedsignal, 0.0) for record in rows))
        push!(counts, length(rows))
    end
    order = sortperm(totalsignal; rev=true)
    keep = order[1:min(Int(n), length(order))]

    (
        mzvalues=mzvalues[keep],
        totalsignal=totalsignal[keep],
        counts=counts[keep],
    )
end

function alkane_variance_max_contiguous_true(mask)
    best = 0
    current = 0
    for value in mask
        if value
            current += 1
            best = max(best, current)
        else
            current = 0
        end
    end

    best
end

function alkane_variance_flat_linear_baseline(x, y)
    center = mean(x)
    scale = maximum(abs.(x .- center))
    scale = isfinite(scale) && scale > 0 ? scale : 1.0
    xnorm = (Float64.(x) .- center) ./ scale
    design = hcat(ones(length(xnorm)), xnorm)
    coefficients = design \ Float64.(y)

    design * coefficients, coefficients[2] / scale
end

function alkane_variance_fit_linear_observed_model(peakrecords, flatrecords)
    length(peakrecords) ≥ 3 || throw(ArgumentError(
        "at least three accepted peak residual records are required"))
    length(flatrecords) ≥ 3 || throw(ArgumentError(
        "at least three flat-region residual records are required"))

    intercept = alkane_variance_flat_variance_anchor(flatrecords)
    isfinite(intercept) && intercept ≥ 0 || throw(ArgumentError(
        "flat-region residuals do not provide a finite variance anchor"))
    intensity_offset = alkane_variance_safe_median(
        [record.fittedintensity for record in flatrecords],
    )
    isfinite(intensity_offset) || throw(ArgumentError(
        "flat-region records do not provide a finite intensity offset"))
    robust = alkane_variance_fit_robust_intensity_slope(
        peakrecords,
        intercept;
        intensity_offset=intensity_offset,
    )
    intensities = [
        record.fittedintensity for record in Iterators.flatten((peakrecords, flatrecords))
        if isfinite(record.fittedintensity)
    ]
    !isempty(intensities) || throw(ArgumentError(
        "no finite fitted intensities are available for the variance model range"))
    model = LinearObservedIntensityVarianceModel(
        intercept,
        robust.slope,
        intensity_offset,
        0.0,
        maximum(intensities),
        0.0,
    )

    (
        status=:ok,
        model=model,
        flatvarianceanchor=intercept,
        intensity_offset=intensity_offset,
        robustslope=robust,
        intensity_min=0.0,
        intensity_max=maximum(intensities),
    )
end

function alkane_variance_flat_variance_anchor(
    flatrecords;
    trim_fraction::Real=0.2,
)
    values = Float64[
        record.residual2 for record in flatrecords
        if isfinite(record.residual2) && record.residual2 ≥ 0
    ]
    alkane_variance_trimmed_mean(values; trim_fraction=trim_fraction)
end

function alkane_variance_fit_robust_intensity_slope(
    records,
    intercept::Real;
    nbins::Integer=24,
    trim_fraction::Real=0.2,
    slope2_keep_quantile::Real=0.75,
    min_intensity_quantile::Real=0.05,
    minbins::Integer=5,
    intensity_offset::Real=0.0,
)
    fixed_intercept = max(Float64(intercept), 0.0)
    offset = Float64(intensity_offset)
    intensity = Float64[record.fittedintensity for record in records]
    shifted_intensity = max.(intensity .- offset, 0.0)
    slope2 = Float64[record.fittedsignalslope2 for record in records]
    residual2 = Float64[record.residual2 for record in records]
    finite = isfinite.(intensity) .&
        isfinite.(shifted_intensity) .&
        isfinite.(slope2) .&
        isfinite.(residual2) .&
        (shifted_intensity .> 0) .&
        (slope2 .≥ 0) .&
        (residual2 .≥ 0)
    finite_intensity = shifted_intensity[finite]
    finite_slope2 = slope2[finite]
    !isempty(finite_intensity) || throw(ArgumentError(
        "no finite positive-intensity peak residual records are available"))

    iq = clamp(Float64(min_intensity_quantile), 0.0, 0.90)
    sq = clamp(Float64(slope2_keep_quantile), 0.05, 1.0)
    intensity_threshold = quantile(finite_intensity, iq)
    slope2_threshold = sq ≥ 1.0 ? Inf : quantile(finite_slope2, sq)
    keep = finite .&
        (shifted_intensity .≥ intensity_threshold) .&
        (slope2 .≤ slope2_threshold)
    filteredrecordcount = count(keep)
    adaptive_minbins = min(Int(minbins), max(3, filteredrecordcount ÷ 20))
    filteredrecordcount ≥ adaptive_minbins || throw(ArgumentError(
        "too few residual records remain after robust slope filtering"))

    filtered = records[keep]
    bins = alkane_variance_binned_residual_records(
        filtered;
        nbins=min(nbins, max(adaptive_minbins, filteredrecordcount ÷ 20)),
        trim_fraction=trim_fraction,
    )
    length(bins.residual2) ≥ adaptive_minbins || throw(ArgumentError(
        "too few populated bins remain after robust slope filtering"))

    bin_shifted_intensity = max.(bins.intensity .- offset, 0.0)
    bin_slopes = max.(bins.residual2 .- fixed_intercept, 0.0) ./
        bin_shifted_intensity
    finite_bins = isfinite.(bin_slopes) .& (bin_slopes .≥ 0)
    sum(finite_bins) ≥ adaptive_minbins || throw(ArgumentError(
        "too few finite robust slope bins are available"))
    bin_slopes = bin_slopes[finite_bins]
    bin_shifted_intensity = bin_shifted_intensity[finite_bins]
    bin_residual2 = bins.residual2[finite_bins]
    bin_n = bins.n[finite_bins]

    slope = median(bin_slopes)
    prediction = fixed_intercept .+ slope .* bin_shifted_intensity
    sse = sum(Float64.(bin_n) .* abs2.(bin_residual2 .- prediction))
    center = median(bin_slopes)
    robust_cv = center > 0 ? median(abs.(bin_slopes .- center)) / center : NaN

    (
        status=:robust_bin_ratio,
        intercept=fixed_intercept,
        slope=slope,
        intensity_offset=offset,
        bins=bins,
        nbins=length(bin_slopes),
        trim_fraction=Float64(trim_fraction),
        filteredrecordcount=filteredrecordcount,
        originalrecordcount=length(records),
        slope2_keep_quantile=sq,
        slope2_threshold=slope2_threshold,
        min_intensity_quantile=iq,
        intensity_threshold=intensity_threshold,
        bin_shifted_intensity=bin_shifted_intensity,
        bin_slopes=bin_slopes,
        b_q25=quantile(bin_slopes, 0.25),
        b_median=slope,
        b_q75=quantile(bin_slopes, 0.75),
        b_robust_cv=robust_cv,
        weighted_sse=sse,
    )
end

function alkane_variance_binned_residual_records(
    records;
    nbins::Integer=40,
    trim_fraction::Real=0.2,
)
    indices = Int[
        index for index in eachindex(records)
        if isfinite(records[index].fittedintensity) &&
            isfinite(records[index].fittedsignalslope2) &&
            isfinite(records[index].residual2) &&
            records[index].fittedintensity ≥ 0 &&
            records[index].fittedsignalslope2 ≥ 0 &&
            records[index].residual2 ≥ 0
    ]
    isempty(indices) && return (
        intensity=Float64[],
        slope2=Float64[],
        residual2=Float64[],
        n=Int[],
    )

    sort!(indices; by=index -> records[index].fittedintensity)
    bin_count = min(max(1, Int(nbins)), length(indices))
    edges = round.(Int, range(1, length(indices) + 1; length=bin_count + 1))
    intensity_bins = Float64[]
    slope2_bins = Float64[]
    residual2_bins = Float64[]
    n_bins = Int[]

    for bin in 1:bin_count
        lo = edges[bin]
        hi = edges[bin + 1] - 1
        hi < lo && continue
        idx = indices[lo:hi]
        push!(intensity_bins, mean(record.fittedintensity for record in records[idx]))
        push!(slope2_bins, mean(record.fittedsignalslope2 for record in records[idx]))
        push!(residual2_bins, alkane_variance_trimmed_mean(
            [record.residual2 for record in records[idx]];
            trim_fraction=trim_fraction,
        ))
        push!(n_bins, length(idx))
    end

    (
        intensity=intensity_bins,
        slope2=slope2_bins,
        residual2=residual2_bins,
        n=n_bins,
    )
end

function alkane_variance_lag1_autocorrelation(
    records,
    model::LinearObservedIntensityVarianceModel;
    trim_fraction::Real=0.05,
)
    grouped = Dict{Tuple{Int, Int, Int}, Vector{NamedTuple}}()
    for record in records
        key = (Int(record.runindex), Int(record.ladderstep), Int(record.mzindex))
        push!(get!(grouped, key, NamedTuple[]), record)
    end

    left = Float64[]
    right = Float64[]
    tracecount = 0
    for rows in values(grouped)
        length(rows) ≥ 2 || continue
        sort!(rows; by=row -> row.scanindex)
        trace_has_pair = false
        for index in 1:(length(rows) - 1)
            rows[index + 1].scanindex == rows[index].scanindex + 1 || continue
            z1 = alkane_variance_standardized_residual(rows[index], model)
            z2 = alkane_variance_standardized_residual(rows[index + 1], model)
            isfinite(z1) && isfinite(z2) || continue
            push!(left, z1)
            push!(right, z2)
            trace_has_pair = true
        end
        trace_has_pair && (tracecount += 1)
    end

    paircount = length(left)
    if paircount < 3
        return (
            rho=0.0,
            status=:too_few_pairs,
            source=:flat_nonpeak,
            paircount=paircount,
            tracecount=tracecount,
            trim_fraction=Float64(trim_fraction),
        )
    end

    keep = trues(paircount)
    trim = clamp(Float64(trim_fraction), 0.0, 0.49)
    if trim > 0 && paircount ≥ 4
        score = max.(abs.(left), abs.(right))
        threshold = quantile(score, 1 - trim)
        keep .= score .≤ threshold
    end
    x = left[keep]
    y = right[keep]
    if length(x) < 3
        return (
            rho=0.0,
            status=:too_few_pairs_after_trimming,
            source=:flat_nonpeak,
            paircount=length(x),
            tracecount=tracecount,
            trim_fraction=trim,
        )
    end

    denominator = sqrt(sum(abs2, x) * sum(abs2, y))
    rho = denominator > 0 ? sum(x .* y) / denominator : 0.0
    rho = clamp(rho, -0.999999, 0.999999)

    (
        rho=Float64(rho),
        status=:ok,
        source=:flat_nonpeak,
        paircount=length(x),
        tracecount=tracecount,
        trim_fraction=trim,
    )
end

function alkane_variance_standardized_residual(
    record,
    model::LinearObservedIntensityVarianceModel,
)
    variance = varpred(record.fittedintensity, model)
    isfinite(variance) && variance > 0 || return NaN

    record.residual / sqrt(variance)
end

function alkane_variance_run_qc(run, flatqc, variancefit)
    reasons = Symbol[]
    keptpeakcount = length(run.peakqc.includedladdersteps)
    excludedpeakcount = length(run.peakqc.excludedladdersteps)
    peakrecordcount = length(run.residualrecords)
    flatrecordcount = flatqc.recordcount
    keptpeakcount ≥ 1 || push!(reasons, :too_few_kept_peaks)
    peakrecordcount ≥ 3 || push!(reasons, :too_few_peak_records)
    flatrecordcount ≥ 3 || push!(reasons, :too_few_flat_records)
    alkane_variance_excluded_peak_fraction(run.peakqc) ≤ 0.50 ||
        push!(reasons, :too_many_excluded_peaks)
    model = variancefit.model
    isfinite(model.intercept) && model.intercept ≥ 0 ||
        push!(reasons, :invalid_intercept)
    isfinite(model.slope) && model.slope ≥ 0 ||
        push!(reasons, :invalid_slope)
    isfinite(model.intensity_offset) ||
        push!(reasons, :invalid_intensity_offset)
    isfinite(model.rho_lag1) && abs(model.rho_lag1) < 1 ||
        push!(reasons, :invalid_lag1)

    (
        accept=isempty(reasons),
        status=isempty(reasons) ? :ok : :rejected,
        reasons=Tuple(reasons),
        error=nothing,
        flatrecordcount=flatrecordcount,
        flatwindowcount=flatqc.windowcount,
        flationcount=flatqc.ioncount,
        peakrecordcount=peakrecordcount,
        keptpeakcount=keptpeakcount,
        excludedpeakcount=excludedpeakcount,
        excludedfraction=alkane_variance_excluded_peak_fraction(run.peakqc),
        intensity_min=model.intensity_min,
        intensity_max=model.intensity_max,
        flatvarianceanchor=variancefit.flatvarianceanchor,
        intensity_offset=variancefit.intensity_offset,
        robustslope=model.slope,
        robustslope_q25=variancefit.robustslope.b_q25,
        robustslope_q75=variancefit.robustslope.b_q75,
        robustslope_cv=variancefit.robustslope.b_robust_cv,
        lag1=model.rho_lag1,
        lag1paircount=variancefit.lag1.paircount,
        lag1tracecount=variancefit.lag1.tracecount,
    )
end

function alkane_variance_empty_flat_qc()
    (
        status=:failed,
        recordcount=0,
        windowcount=0,
        ioncount=0,
        requestedioncount=0,
        peakmargin_scans=5,
        flatwindowsize=17,
        flatwindowstride=17,
        flat_slope_quantile=0.35,
        flat_curvature_quantile=0.50,
    )
end

function alkane_variance_excluded_peak_fraction(peakqc)
    total = length(peakqc.includedladdersteps) + length(peakqc.excludedladdersteps)
    total == 0 ? 1.0 : length(peakqc.excludedladdersteps) / total
end

function alkane_variance_aggregate_models(runs)
    models = [run.model for run in runs if !isnothing(run.model)]
    !isempty(models) || throw(ArgumentError(
        "at least one fitted run model is required for aggregation"))
    LinearObservedIntensityVarianceModel(
        median([model.intercept for model in models]),
        median([model.slope for model in models]),
        median([model.intensity_offset for model in models]),
        minimum([model.intensity_min for model in models]),
        maximum([model.intensity_max for model in models]),
        median([model.rho_lag1 for model in models]),
    )
end

function alkane_variance_fit_qc(runs)
    accepted = Int[
        run.runindex for run in runs
        if !isnothing(run.model) && run.varianceqc.accept
    ]
    rejected = Int[
        run.runindex for run in runs
        if isnothing(run.model) || !run.varianceqc.accept
    ]
    (
        status=isempty(rejected) ? :ok : :partial,
        acceptedrunindices=accepted,
        rejectedrunindices=rejected,
        runcount=length(runs),
        acceptedruncount=length(accepted),
        rejectedruncount=length(rejected),
        rows=[
            merge((
                runindex=run.runindex,
            ), run.varianceqc)
            for run in runs
        ],
    )
end

function alkane_variance_trimmed_mean(values; trim_fraction::Real=0.2)
    xs = sort(Float64[value for value in values if isfinite(value)])
    isempty(xs) && return NaN
    trim = clamp(Float64(trim_fraction), 0.0, 0.49)
    n = length(xs)
    k = floor(Int, trim * n)
    lo = min(k + 1, n)
    hi = max(n - k, lo)

    mean(@view xs[lo:hi])
end

function alkane_variance_safe_median(values)
    finite = Float64[value for value in values if isfinite(value)]
    isempty(finite) && return NaN

    median(finite)
end

function alkane_variance_peak_fitted_signal_slope(
    peakinput,
    peakfit,
    ionposition::Integer,
    normalizedretention::Real,
)
    retentionscale = Float64(peakinput.retentionscale)
    if !(isfinite(retentionscale) && retentionscale > 0)
        return NaN
    end
    spectrumweight = Float64(peakinput.spectrumweights[ionposition])
    derivative = alkane_variance_envelope_derivative(
        peakfit.envelopegrid,
        peakfit.envelopevalues,
        normalizedretention,
    )

    peakinput.intensityscale * spectrumweight * derivative / retentionscale
end

function alkane_variance_envelope_derivative(
    grid::AbstractVector{<:Real},
    values::AbstractVector{<:Real},
    retention::Real,
)
    length(grid) == length(values) || throw(DimensionMismatch(
        "envelope grid and values must have the same length"))
    length(grid) ≥ 2 || return NaN
    value = Float64(retention)
    right = if value ≤ first(grid)
        2
    elseif value ≥ last(grid)
        length(grid)
    else
        searchsortedfirst(grid, value)
    end
    right = clamp(right, 2, length(grid))
    left = right - 1
    dx = Float64(grid[right] - grid[left])
    dx > 0 || return NaN

    (Float64(values[right]) - Float64(values[left])) / dx
end

function alkane_variance_fit_peak_envelope(peakinput)
    initial_fitrows = findall(peakinput.fitionmask)
    !isempty(initial_fitrows) || throw(ArgumentError(
        "peak envelope fitting requires at least one fit ion"))
    size(peakinput.normalizedretentions) == size(peakinput.normalizedobserved) ||
        throw(DimensionMismatch(
            "normalized retention and intensity matrices must have the same size"))
    length(peakinput.spectrumweights) == size(peakinput.normalizedobserved, 1) ||
        throw(DimensionMismatch(
            "spectrumweights length must match the number of extracted ions"))

    firstfit = alkane_variance_fit_peak_envelope_with_rows(peakinput, initial_fitrows)
    final_fitrows, excluded_fitrows, replacement_fitrows =
        alkane_variance_peak_replacement_fitrows(peakinput, firstfit, initial_fitrows)

    finalfit = if final_fitrows == initial_fitrows
        firstfit
    else
        alkane_variance_fit_peak_envelope_with_rows(peakinput, final_fitrows)
    end

    alkane_variance_annotate_peak_fit(
        finalfit,
        initial_fitrows,
        excluded_fitrows,
        replacement_fitrows,
        final_fitrows,
    )
end

function alkane_variance_fit_peak_envelope_with_rows(peakinput, fitrows)
    normalizedapexshift = alkane_variance_peak_apex_shift_estimate(peakinput, fitrows)
    failures = String[]
    fit = try
        alkane_variance_fit_peak_envelope_at_shift(
            peakinput,
            fitrows,
            normalizedapexshift,
        )
    catch err
        push!(failures, sprint(showerror, err))
        fallbackshift = 0.0
        if !isapprox(normalizedapexshift, fallbackshift; atol=sqrt(eps(Float64)))
            try
                alkane_variance_fit_peak_envelope_at_shift(
                    peakinput,
                    fitrows,
                    fallbackshift,
                )
            catch fallbackerr
                push!(failures, sprint(showerror, fallbackerr))
                nothing
            end
        else
            nothing
        end
    end

    if isnothing(fit)
        message = isempty(failures) ? "" : ": " * join(unique(failures), "; ")
        throw(ArgumentError("peak envelope optimization failed" * message))
    end

    fit
end

function alkane_variance_annotate_peak_fit(
    fit,
    initial_fitrows,
    excluded_fitrows,
    replacement_fitrows,
    final_fitrows,
)
    merge(fit, (
        fitionpositions_initial=Int.(initial_fitrows),
        fitionpositions=Int.(final_fitrows),
        excludedfitionpositions=Int.(excluded_fitrows),
        replacementfitionpositions=Int.(replacement_fitrows),
    ))
end

function alkane_variance_peak_replacement_fitrows(peakinput, fit, initial_fitrows)
    outliers = alkane_variance_peak_candidate_outlier_rows(peakinput, fit, initial_fitrows)
    isempty(outliers) && return Int.(initial_fitrows), Int[], Int[]

    current = Int.(initial_fitrows)
    excluded = Int[]
    replacements = Int[]
    candidates = [
        index for index in eachindex(peakinput.spectrumweights)
        if !(index in current) &&
            isfinite(peakinput.spectrumweights[index]) &&
            peakinput.spectrumweights[index] > 0
    ]
    minfitcount = min(3, length(initial_fitrows))

    for outlier in outliers
        outlier in current || continue
        replacement = isempty(candidates) ? nothing : popfirst!(candidates)
        if !isnothing(replacement)
            deleteat!(current, findfirst(==(outlier), current))
            push!(current, replacement)
            push!(excluded, outlier)
            push!(replacements, replacement)
        elseif length(current) > minfitcount
            deleteat!(current, findfirst(==(outlier), current))
            push!(excluded, outlier)
        end
    end

    sort!(current)
    current, excluded, replacements
end

function alkane_variance_peak_candidate_outlier_rows(peakinput, fit, fitrows)
    outliers = Int[]
    if hasproperty(peakinput, :normalizedobserved) &&
            hasproperty(peakinput, :spectrumweights) &&
            hasproperty(fit, :normalizedapexshift)
        try
            append!(outliers, alkane_variance_peak_leave_one_ion_outlier_rows(
                peakinput,
                fitrows,
                fit.normalizedapexshift,
            ))
        catch
        end
    end

    append!(outliers, alkane_variance_peak_fit_outlier_rows(fit, fitrows))
    unique!(outliers)
    outlier_scores = alkane_variance_peak_fit_row_rmses(fit.normalizedresiduals, fitrows)
    sort!(outliers; by=row -> get(outlier_scores, row, -Inf), rev=true)
    Int.(outliers)
end

function alkane_variance_peak_leave_one_ion_outlier_rows(
    peakinput,
    fitrows,
    normalizedapexshift::Real,
)
    length(fitrows) ≥ 3 || return Int[]
    rowrmses = Dict{Int, Float64}()
    for heldout in fitrows
        trainrows = Int[row for row in fitrows if row != heldout]
        length(trainrows) ≥ 2 || continue
        fit = alkane_variance_fit_peak_envelope_at_shift(
            peakinput,
            trainrows,
            normalizedapexshift,
        )
        spectrumweight = Float64(peakinput.spectrumweights[heldout])
        isfinite(spectrumweight) && spectrumweight > 0 || continue
        values = Float64[
            value / spectrumweight for value in fit.normalizedresiduals[heldout, :]
            if isfinite(value)
        ]
        isempty(values) && continue
        rowrmses[Int(heldout)] = sqrt(mean(abs2, values))
    end

    alkane_variance_peak_rmse_outlier_rows(rowrmses, fitrows)
end

function alkane_variance_peak_fit_outlier_rows(fit, fitrows)
    length(fitrows) ≥ 3 || return Int[]
    rowrmses = alkane_variance_peak_fit_row_rmses(fit.normalizedresiduals, fitrows)
    alkane_variance_peak_rmse_outlier_rows(rowrmses, fitrows)
end

function alkane_variance_peak_rmse_outlier_rows(rowrmses, fitrows)
    finite = [value for value in values(rowrmses) if isfinite(value)]
    length(finite) ≥ 3 || return Int[]

    center = median(finite)
    scale = 1.4826 * median(abs.(finite .- center))
    threshold = max(center + 3 * scale, 1.5 * center, center + sqrt(eps(Float64)))
    outliers = [
        row for row in fitrows
        if haskey(rowrmses, row) && rowrmses[row] > threshold
    ]
    sort!(outliers; by=row -> rowrmses[row], rev=true)
    Int.(outliers)
end

function alkane_variance_peak_fit_row_rmses(residuals, fitrows)
    rmses = Dict{Int, Float64}()
    for row in fitrows
        values = Float64[
            value for value in residuals[row, :]
            if isfinite(value)
        ]
        isempty(values) && continue
        rmses[Int(row)] = sqrt(mean(abs2, values))
    end

    rmses
end

function alkane_variance_peak_apex_shift_estimate(peakinput, fitrows)
    scanstep = alkane_variance_peak_scan_interval(peakinput)
    if !(isfinite(scanstep) && scanstep > 0 && peakinput.retentionscale > 0)
        return 0.0
    end

    maxshift = 0.5 * scanstep / peakinput.retentionscale
    scanretentions = Float64[]
    scanprofile = Float64[]
    for scanindex in axes(peakinput.normalizedobserved, 2)
        values = Float64[]
        retentions = Float64[]
        for ionindex in fitrows
            spectrumweight = Float64(peakinput.spectrumweights[ionindex])
            spectrumweight > 0 && isfinite(spectrumweight) || continue
            intensity = Float64(peakinput.normalizedobserved[ionindex, scanindex])
            retention = Float64(peakinput.normalizedretentions[ionindex, scanindex])
            isfinite(intensity) && isfinite(retention) || continue
            intensity > 0 || continue
            push!(values, intensity / spectrumweight)
            push!(retentions, retention)
        end
        isempty(values) && continue
        push!(scanprofile, median(values))
        push!(scanretentions, mean(retentions))
    end

    length(scanprofile) ≥ 1 || return 0.0
    apexindex = argmax(scanprofile)
    estimate = scanretentions[apexindex]
    if 1 < apexindex < length(scanprofile)
        localestimate = alkane_variance_peak_quadratic_apex_estimate(
            scanretentions[(apexindex - 1):(apexindex + 1)],
            scanprofile[(apexindex - 1):(apexindex + 1)],
        )
        isfinite(localestimate) && (estimate = localestimate)
    end

    clamp(estimate, -maxshift, maxshift)
end

function alkane_variance_peak_quadratic_apex_estimate(retentions, profile)
    maximum(profile) > 0 || return NaN
    floorvalue = maximum(profile) * 1e-6
    y = log.(max.(Float64.(profile), floorvalue))
    X = hcat(Float64.(retentions).^2, Float64.(retentions), ones(3))
    coefs = X \ y
    a, b = coefs[1], coefs[2]
    a < 0 && isfinite(a) && isfinite(b) || return NaN
    estimate = -b / (2 * a)
    lo, hi = extrema(Float64.(retentions))
    lo ≤ estimate ≤ hi || return NaN

    estimate
end

function alkane_variance_fit_peak_envelope_at_shift(
    peakinput,
    fitrows,
    normalizedapexshift::Real,
    smoothnessfactor::Real=0.004,
)
    adjustedretentions = peakinput.normalizedretentions .- normalizedapexshift
    fitobservations = alkane_variance_peak_fit_observations(
        peakinput,
        fitrows,
        adjustedretentions,
    )
    grid, envelopevalues, envelopeknots, envelopecoefficients =
        alkane_variance_peak_direct_envelope_with_linear_baseline(
        peakinput,
        fitrows,
        adjustedretentions,
        fitobservations,
        smoothnessfactor,
    )
    normalizedfitted = alkane_variance_peak_envelope_predict(
        peakinput,
        adjustedretentions,
        grid,
        envelopevalues,
    )
    normalizedbaseline = alkane_variance_peak_fit_linear_baselines(
        peakinput,
        normalizedfitted,
    )
    normalizedfittedraw = normalizedfitted .+ normalizedbaseline
    normalizedresiduals = peakinput.normalizedobservedraw .- normalizedfittedraw
    fitted = normalizedfitted .* peakinput.intensityscale
    baseline = normalizedbaseline .* peakinput.intensityscale
    fittedraw = normalizedfittedraw .* peakinput.intensityscale
    residuals = peakinput.observedraw .- fittedraw
    fitrowmask = falses(size(residuals))
    fitrowmask[fitrows, :] .= true
    apexretentionshift = Float64(normalizedapexshift) * peakinput.retentionscale

    (
        ladderstep=peakinput.ladderstep,
        envelopegrid=grid .+ Float64(normalizedapexshift),
        envelopevalues=envelopevalues,
        envelopeknots=envelopeknots .+ Float64(normalizedapexshift),
        envelopecoefficients=envelopecoefficients,
        envelopemethod=:curvature_adaptive_penalized_kernel,
        smoothnessfactor=Float64(smoothnessfactor),
        normalizedapexshift=Float64(normalizedapexshift),
        apexretentionshift=apexretentionshift,
        fittedapexretention=peakinput.apexretention + apexretentionshift,
        apexscanshift=alkane_variance_peak_apex_scan_shift(
            peakinput,
            apexretentionshift,
        ),
        baselinemodel=:linear,
        normalizedbaseline=normalizedbaseline,
        baseline=baseline,
        normalizedfitted=normalizedfitted,
        normalizedfittedraw=normalizedfittedraw,
        normalizedresiduals=normalizedresiduals,
        fitted=fitted,
        fittedraw=fittedraw,
        residuals=residuals,
        rmse=alkane_variance_rmse(residuals),
        normalizedrmse=alkane_variance_rmse(normalizedresiduals),
        fitrmse=alkane_variance_rmse(residuals[fitrowmask]),
        normalizedfitrmse=alkane_variance_rmse(normalizedresiduals[fitrowmask]),
        normalizedfitrss=sum(abs2, normalizedresiduals[fitrowmask]),
        observationcount=alkane_variance_finite_count(residuals),
        fitobservationcount=length(fitobservations),
        fitionpositions=Int.(fitrows),
    )
end

function alkane_variance_peak_direct_envelope_with_linear_baseline(
    peakinput,
    fitrows,
    adjustedretentions::AbstractMatrix{<:Real},
    fitobservations,
    smoothnessfactor::Real=0.004,
)
    signalobservations = alkane_variance_peak_fit_observations_from_response(
        peakinput,
        fitrows,
        adjustedretentions,
        alkane_variance_peak_signal_fit_response(peakinput),
    )
    profile_retention, profile_intensity, profile_scanindices =
        alkane_variance_peak_observation_profile(signalobservations)
    left = min(minimum(profile_retention), 0.0)
    right = max(maximum(profile_retention), 0.0)
    left < right || throw(ArgumentError(
        "peak envelope retention range must be nonzero"))

    envelopeknots = alkane_variance_peak_adaptive_envelope_grid(
        left,
        right,
        profile_retention,
        profile_intensity,
        profile_scanindices,
    )
    envelopecoefficients =
        alkane_variance_peak_adaptive_penalized_coefficients_with_linear_baseline(
            peakinput,
            adjustedretentions,
            fitrows,
            fitobservations,
            envelopeknots,
            smoothnessfactor,
        )
    grid = alkane_variance_peak_envelope_dense_grid(
        left,
        right,
        size(adjustedretentions, 2),
    )
    envelopevalues = alkane_variance_peak_kernel_values(
        envelopeknots,
        envelopecoefficients,
        grid,
    )

    grid, envelopevalues, envelopeknots, envelopecoefficients
end

function alkane_variance_peak_adaptive_envelope_grid(
    left::Real,
    right::Real,
    retentions::AbstractVector{<:Real},
    intensities::AbstractVector{<:Real},
    scanindices::AbstractVector{<:Integer},
)
    scan_count = length(unique(Int.(scanindices)))
    pilotgrid = alkane_variance_peak_envelope_dense_grid(left, right, scan_count)
    pilotvalues = alkane_variance_peak_pilot_envelope(
        pilotgrid,
        retentions,
        intensities,
        scanindices,
    )
    curvature = alkane_variance_peak_abs_curvature(pilotgrid, pilotvalues)
    density = alkane_variance_peak_curvature_density(curvature)
    gridcount = alkane_variance_peak_adaptive_grid_count(scan_count)
    grid = alkane_variance_peak_density_quantile_grid(pilotgrid, density, gridcount)

    alkane_variance_peak_unique_grid(vcat(Float64(left), grid, 0.0, Float64(right)))
end

function alkane_variance_peak_pilot_envelope(
    grid::AbstractVector{<:Real},
    retentions::AbstractVector{<:Real},
    intensities::AbstractVector{<:Real},
    scanindices::AbstractVector{<:Integer},
)
    bandwidth = 1.75 * alkane_variance_peak_local_quadratic_bandwidth(
        retentions,
        scanindices,
    )

    [
        max(0.0, alkane_variance_peak_local_quadratic_value(
            Float64(point),
            retentions,
            intensities,
            bandwidth,
        ))
        for point in grid
    ]
end

function alkane_variance_peak_abs_curvature(grid, values)
    length(grid) == length(values) || throw(DimensionMismatch(
        "curvature grid and values must have the same length"))
    curvature = zeros(Float64, length(grid))
    length(grid) ≥ 3 || return curvature

    for index in 2:(length(grid) - 1)
        hleft = Float64(grid[index] - grid[index - 1])
        hright = Float64(grid[index + 1] - grid[index])
        hleft > 0 && hright > 0 || continue
        slopeleft = (Float64(values[index]) - Float64(values[index - 1])) / hleft
        sloperight = (Float64(values[index + 1]) - Float64(values[index])) / hright
        curvature[index] = abs(2 * (sloperight - slopeleft) / (hleft + hright))
    end
    curvature[1] = curvature[2]
    curvature[end] = curvature[end - 1]

    curvature
end

function alkane_variance_peak_curvature_density(curvature)
    finite = Float64[value for value in curvature if isfinite(value) && value ≥ 0]
    isempty(finite) && return ones(Float64, length(curvature))
    cap = quantile(finite, 0.9)
    scale = max(cap, median(finite), sqrt(eps(Float64)))
    adjusted = min.(max.(Float64.(curvature), 0.0), cap)

    1.0 .+ 4.0 .* sqrt.(adjusted ./ scale)
end

function alkane_variance_peak_adaptive_grid_count(scan_count::Integer)
    clamp(2 * Int(scan_count) + 5, 12, 45)
end

function alkane_variance_peak_density_quantile_grid(grid, density, gridcount::Integer)
    length(grid) == length(density) || throw(DimensionMismatch(
        "density grid and density values must have the same length"))
    length(grid) ≥ 2 || return Float64.(grid)
    cumulative = zeros(Float64, length(grid))
    for index in 2:length(grid)
        width = Float64(grid[index] - grid[index - 1])
        width > 0 || continue
        localdensity = 0.5 * (
            max(Float64(density[index - 1]), 0.0) +
            max(Float64(density[index]), 0.0)
        )
        cumulative[index] = cumulative[index - 1] + width * localdensity
    end
    total = last(cumulative)
    total > 0 || return collect(range(first(grid), last(grid); length=gridcount))

    targets = collect(range(0.0, total; length=gridcount))
    [
        alkane_variance_linear_interpolate(cumulative, grid, target)
        for target in targets
    ]
end

function alkane_variance_linear_interpolate(x, y, target::Real)
    length(x) == length(y) || throw(DimensionMismatch(
        "interpolation x and y vectors must have the same length"))
    value = Float64(target)
    value ≤ first(x) && return Float64(first(y))
    value ≥ last(x) && return Float64(last(y))

    right = searchsortedfirst(x, value)
    x[right] == value && return Float64(y[right])
    left = right - 1
    fraction = (value - x[left]) / (x[right] - x[left])

    (1 - fraction) * Float64(y[left]) + fraction * Float64(y[right])
end

function alkane_variance_peak_unique_grid(values)
    sorted = sort(Float64.(values))
    isempty(sorted) && return sorted
    tolerance = sqrt(eps(Float64)) * max(abs(first(sorted)), abs(last(sorted)), 1.0)
    grid = Float64[first(sorted)]
    for value in sorted[2:end]
        if value - last(grid) > tolerance
            push!(grid, value)
        end
    end

    grid
end

function alkane_variance_peak_adaptive_penalized_coefficients_with_linear_baseline(
    peakinput,
    adjustedretentions::AbstractMatrix{<:Real},
    fitrows,
    fitobservations,
    centers::AbstractVector{<:Real},
    smoothnessfactor::Real=0.004,
)
    length(centers) ≥ 3 || throw(ArgumentError(
        "peak envelope fitting requires at least three envelope centers"))
    length(fitobservations) ≥ 3 || throw(ArgumentError(
        "peak envelope fitting requires at least three fit observations"))

    envelope_design = alkane_variance_peak_kernel_matrix(
        [observation.retention for observation in fitobservations],
        centers,
    )
    envelope_design .*= reshape(
        [observation.spectrumweight for observation in fitobservations],
        :,
        1,
    )

    roughnessgrid = alkane_variance_peak_kernel_roughness_grid(centers)
    roughbasis = alkane_variance_peak_kernel_matrix(roughnessgrid, centers)
    roughness = alkane_variance_peak_spacing_aware_curvature_matrix(
        roughnessgrid,
    ) * roughbasis
    λ = alkane_variance_peak_scaled_smoothness_weight(
        envelope_design,
        roughness,
        smoothnessfactor,
    )

    design, y, envelope_column_count = alkane_variance_peak_raw_linear_baseline_design(
        peakinput,
        adjustedretentions,
        fitrows,
        fitobservations,
        centers,
    )
    roughness_with_baseline = hcat(
        roughness,
        zeros(Float64, size(roughness, 1), size(design, 2) - size(roughness, 2)),
    )
    coefficients = alkane_variance_peak_partly_nonnegative_penalized_solve(
        design,
        y,
        roughness_with_baseline,
        λ,
        envelope_column_count,
    )

    coefficients[1:envelope_column_count]
end

function alkane_variance_peak_raw_linear_baseline_design(
    peakinput,
    adjustedretentions::AbstractMatrix{<:Real},
    fitrows,
    fitobservations,
    centers::AbstractVector{<:Real},
)
    envelope_column_count = length(centers)
    baseline_column_count = 2 * length(fitrows)
    data_row_count = length(fitobservations)
    prior_row_count = 2 * length(fitrows)
    design = zeros(
        Float64,
        data_row_count + prior_row_count,
        envelope_column_count + baseline_column_count,
    )
    y = zeros(Float64, size(design, 1))
    fitrow_position = Dict(Int(row) => index for (index, row) in pairs(Int.(fitrows)))
    retention_points = [observation.retention for observation in fitobservations]
    envelope_design = alkane_variance_peak_kernel_matrix(retention_points, centers)

    for (rowindex, observation) in pairs(fitobservations)
        ionindex = Int(observation.ionindex)
        scanindex = Int(observation.scanindex)
        position = fitrow_position[ionindex]
        design[rowindex, 1:envelope_column_count] .=
            Float64(observation.spectrumweight) .* envelope_design[rowindex, :]
        baseline_col = envelope_column_count + 2 * (position - 1) + 1
        x = Float64(peakinput.normalizedretentions[ionindex, scanindex])
        design[rowindex, baseline_col] = 1.0
        design[rowindex, baseline_col + 1] = x
        y[rowindex] = Float64(observation.intensity)
    end

    priorweight = alkane_variance_peak_linear_baseline_prior_weight()
    prior_scale = sqrt(priorweight)
    rowindex = data_row_count
    for (position, ionindex) in pairs(Int.(fitrows))
        intercept, slope = alkane_variance_peak_arpls_linear_baseline(
            peakinput,
            ionindex,
        )
        baseline_col = envelope_column_count + 2 * (position - 1) + 1

        rowindex += 1
        design[rowindex, baseline_col] = prior_scale
        y[rowindex] = prior_scale * intercept

        rowindex += 1
        design[rowindex, baseline_col + 1] = prior_scale
        y[rowindex] = prior_scale * slope
    end

    design, y, envelope_column_count
end

function alkane_variance_peak_kernel_values(centers, coefficients, grid)
    design = alkane_variance_peak_kernel_matrix(grid, centers)
    max.(0.0, design * Float64.(coefficients))
end

function alkane_variance_peak_kernel_matrix(points, centers)
    center_values = Float64.(centers)
    widths = alkane_variance_peak_kernel_widths(center_values)
    matrix = Matrix{Float64}(undef, length(points), length(center_values))
    for (row, point) in pairs(points)
        x = Float64(point)
        for column in eachindex(center_values)
            z = (x - center_values[column]) / widths[column]
            matrix[row, column] = exp(-0.5 * z^2)
        end
    end

    matrix
end

function alkane_variance_peak_kernel_widths(centers)
    length(centers) == 1 && return [1.0]
    spacings = diff(centers)
    finite_spacings = Float64[
        spacing for spacing in spacings
        if isfinite(spacing) && spacing > 0
    ]
    fallback = isempty(finite_spacings) ? 1.0 : median(finite_spacings)
    widths = similar(Float64.(centers))
    for index in eachindex(centers)
        left = index == firstindex(centers) ? fallback : centers[index] - centers[index - 1]
        right = index == lastindex(centers) ? fallback : centers[index + 1] - centers[index]
        localspacing = max(left, right, fallback)
        widths[index] = max(1.25 * localspacing, 0.5 * fallback, sqrt(eps(Float64)))
    end

    widths
end

function alkane_variance_peak_kernel_roughness_grid(centers)
    count = max(4 * length(centers), 60)
    collect(range(first(centers), last(centers); length=count))
end

function alkane_variance_peak_spacing_aware_curvature_matrix(grid)
    rows = max(length(grid) - 2, 0)
    matrix = zeros(Float64, rows, length(grid))
    for row in 1:rows
        index = row + 1
        hleft = Float64(grid[index] - grid[index - 1])
        hright = Float64(grid[index + 1] - grid[index])
        hleft > 0 && hright > 0 || continue
        weight = sqrt(2 / (hleft + hright))
        matrix[row, index - 1] = weight / hleft
        matrix[row, index] = -weight * (1 / hleft + 1 / hright)
        matrix[row, index + 1] = weight / hright
    end

    matrix
end

function alkane_variance_peak_scaled_smoothness_weight(
    design,
    roughness,
    smoothnessfactor::Real=0.004,
)
    isempty(roughness) && return 0.0
    isfinite(smoothnessfactor) && smoothnessfactor ≥ 0 || throw(ArgumentError(
        "smoothness factor must be finite and nonnegative"))
    datascale = sum(abs2, design) / max(size(design, 2), 1)
    roughnessscale = sum(abs2, roughness) / max(size(roughness, 1), 1)
    roughnessscale > 0 || return 0.0

    Float64(smoothnessfactor) * datascale / roughnessscale
end

function alkane_variance_peak_partly_nonnegative_penalized_solve(
    design,
    y,
    roughness,
    λ::Real,
    nonnegativecount::Integer,
)
    nonnegativecount ≥ 0 || throw(ArgumentError(
        "nonnegative coefficient count must be nonnegative"))
    nonnegativecount ≤ size(design, 2) || throw(ArgumentError(
        "nonnegative coefficient count cannot exceed coefficient count"))
    active = trues(nonnegativecount)
    values = zeros(Float64, size(design, 2))
    freecolumns = collect((nonnegativecount + 1):size(design, 2))
    tolerance = sqrt(eps(Float64))

    for _ in 1:max(nonnegativecount, 1)
        columns = vcat(findall(active), freecolumns)
        isempty(columns) && return values
        reduced_design = design[:, columns]
        reduced_roughness = roughness[:, columns]
        normal_matrix =
            reduced_design' * reduced_design +
            Float64(λ) * (reduced_roughness' * reduced_roughness)
        rhs = reduced_design' * y
        ridge = sqrt(eps(Float64)) * max(1.0, maximum(abs, normal_matrix))
        normal_matrix[diagind(normal_matrix)] .+= ridge
        solution = normal_matrix \ rhs
        values .= 0.0
        values[columns] .= solution

        negative = Int[
            column for (position, column) in pairs(columns)
            if column ≤ nonnegativecount && solution[position] < -tolerance
        ]
        isempty(negative) && return values
        active[negative] .= false
    end

    values[1:nonnegativecount] .= max.(0.0, values[1:nonnegativecount])
    values
end

function alkane_variance_peak_fit_linear_baselines(peakinput, normalizedfitted)
    size(normalizedfitted) == size(peakinput.normalizedobservedraw) ||
        throw(DimensionMismatch(
            "fitted peak signal and raw observed matrices must have the same size"))

    fittedbaseline = Matrix{Float64}(undef, size(normalizedfitted))
    for ionindex in axes(normalizedfitted, 1)
        target = peakinput.normalizedobservedraw[ionindex, :] .-
            normalizedfitted[ionindex, :]
        prior = alkane_variance_peak_arpls_linear_baseline(peakinput, ionindex)
        intercept, slope = alkane_variance_linear_baseline_coefficients(
            peakinput.normalizedretentions[ionindex, :],
            target;
            prior=prior,
            priorweight=alkane_variance_peak_linear_baseline_prior_weight(),
        )
        for scanindex in axes(normalizedfitted, 2)
            x = Float64(peakinput.normalizedretentions[ionindex, scanindex])
            fittedbaseline[ionindex, scanindex] = intercept + slope * x
        end
    end

    fittedbaseline
end

function alkane_variance_peak_arpls_linear_baseline(peakinput, ionindex::Integer)
    alkane_variance_linear_baseline_coefficients(
        peakinput.normalizedretentions[ionindex, :],
        peakinput.normalizedbaseline[ionindex, :],
    )
end

alkane_variance_peak_linear_baseline_prior_weight() = 1.0

function alkane_variance_linear_baseline_coefficients(
    xvalues,
    yvalues;
    prior=nothing,
    priorweight::Real=0.0,
)
    rows = Tuple{Float64, Float64}[]
    for (x, y) in zip(xvalues, yvalues)
        xf = Float64(x)
        yf = Float64(y)
        isfinite(xf) && isfinite(yf) || continue
        push!(rows, (xf, yf))
    end

    priorweight = max(Float64(priorweight), 0.0)
    useprior = !isnothing(prior) && priorweight > 0
    rowcount = length(rows) + (useprior ? 2 : 0)
    if rowcount == 0
        return 0.0, 0.0
    end

    design = Matrix{Float64}(undef, rowcount, 2)
    y = Vector{Float64}(undef, rowcount)
    rowindex = 0
    for (x, value) in rows
        rowindex += 1
        design[rowindex, 1] = 1.0
        design[rowindex, 2] = x
        y[rowindex] = value
    end
    if useprior
        scale = sqrt(priorweight)
        rowindex += 1
        design[rowindex, 1] = scale
        design[rowindex, 2] = 0.0
        y[rowindex] = scale * Float64(prior[1])

        rowindex += 1
        design[rowindex, 1] = 0.0
        design[rowindex, 2] = scale
        y[rowindex] = scale * Float64(prior[2])
    end

    normal_matrix = design' * design
    rhs = design' * y
    ridge = sqrt(eps(Float64)) * max(1.0, maximum(abs, normal_matrix))
    normal_matrix[diagind(normal_matrix)] .+= ridge
    coefficients = normal_matrix \ rhs

    Float64(coefficients[1]), Float64(coefficients[2])
end

function alkane_variance_peak_observation_profile(fitobservations)
    retentions = Float64[]
    intensities = Float64[]
    scanindices = Int[]
    for observation in fitobservations
        retention = Float64(observation.retention)
        intensity = Float64(observation.intensity)
        spectrumweight = Float64(observation.spectrumweight)
        isfinite(retention) && isfinite(intensity) && isfinite(spectrumweight) ||
            continue
        spectrumweight > 0 || continue
        push!(retentions, retention)
        push!(intensities, intensity / spectrumweight)
        push!(scanindices, Int(observation.scanindex))
    end
    length(retentions) ≥ 3 || throw(ArgumentError(
        "peak envelope fitting requires at least three profile points"))

    order = sortperm(retentions)
    retentions[order], max.(0.0, intensities[order]), scanindices[order]
end

function alkane_variance_peak_local_quadratic_envelope(
    grid::AbstractVector{<:Real},
    retentions::AbstractVector{<:Real},
    intensities::AbstractVector{<:Real},
    scanindices::AbstractVector{<:Integer},
)
    length(retentions) == length(intensities) == length(scanindices) ||
        throw(DimensionMismatch(
            "peak profile retentions, intensities, and scanindices must have the same length"))
    length(retentions) ≥ 3 || throw(ArgumentError(
        "peak envelope fitting requires at least three profile points"))

    bandwidth = alkane_variance_peak_local_quadratic_bandwidth(retentions, scanindices)
    values = [
        max(0.0, alkane_variance_peak_local_quadratic_value(
            Float64(point),
            retentions,
            intensities,
            bandwidth,
        ))
        for point in grid
    ]
    values
end

function alkane_variance_peak_local_quadratic_bandwidth(retentions, scanindices)
    scanretentions = Dict{Int, Vector{Float64}}()
    for (retention, scanindex) in zip(retentions, scanindices)
        push!(get!(scanretentions, Int(scanindex), Float64[]), Float64(retention))
    end
    centers = sort!([
        mean(values) for values in values(scanretentions)
        if !isempty(values)
    ])
    if length(centers) ≥ 2
        spacings = [abs(value) for value in diff(centers) if isfinite(value) && value != 0]
        if !isempty(spacings)
            peakwidth = last(centers) - first(centers)
            scan_count = length(centers)
            if isfinite(peakwidth) && peakwidth > 0 && scan_count > 1
                width_scaled = 0.25 * peakwidth / sqrt(scan_count - 1)
                return max(width_scaled, 0.25 * median(spacings), sqrt(eps(Float64)))
            end

            return max(0.6 * median(spacings), sqrt(eps(Float64)))
        end
    end

    extent = maximum(retentions) - minimum(retentions)
    max(extent / 12, sqrt(eps(Float64)))
end

function alkane_variance_peak_local_quadratic_value(
    target::Float64,
    retentions,
    intensities,
    bandwidth::Float64,
)
    weighted_x0 = 0.0
    weighted_x1 = 0.0
    weighted_x2 = 0.0
    weighted_x3 = 0.0
    weighted_x4 = 0.0
    weighted_y0 = 0.0
    weighted_y1 = 0.0
    weighted_y2 = 0.0

    for (retention, intensity) in zip(retentions, intensities)
        x = (Float64(retention) - target) / bandwidth
        abs(x) ≤ 3.0 || continue
        weight = exp(-0.5 * x^2)
        y = Float64(intensity)
        weighted_x0 += weight
        weighted_x1 += weight * x
        weighted_x2 += weight * x^2
        weighted_x3 += weight * x^3
        weighted_x4 += weight * x^4
        weighted_y0 += weight * y
        weighted_y1 += weight * x * y
        weighted_y2 += weight * x^2 * y
    end

    weighted_x0 > 0 || return 0.0
    normal_matrix = [
        weighted_x0 weighted_x1 weighted_x2
        weighted_x1 weighted_x2 weighted_x3
        weighted_x2 weighted_x3 weighted_x4
    ]
    rhs = [weighted_y0, weighted_y1, weighted_y2]
    coefs = try
        normal_matrix \ rhs
    catch
        return weighted_y0 / weighted_x0
    end
    coefs[1]
end

function alkane_variance_peak_scan_interval(peakinput)
    intervals = Float64[]
    for row in axes(peakinput.observationretentions, 1)
        values = Float64[
            value for value in peakinput.observationretentions[row, :]
            if isfinite(value)
        ]
        length(values) ≥ 2 || continue
        append!(intervals, diff(values))
    end
    intervals = [abs(value) for value in intervals if isfinite(value) && value != 0]
    isempty(intervals) && return NaN

    median(intervals)
end

function alkane_variance_peak_apex_scan_shift(peakinput, apexretentionshift::Real)
    scanstep = alkane_variance_peak_scan_interval(peakinput)
    if !(isfinite(scanstep) && scanstep > 0)
        return NaN
    end

    Float64(apexretentionshift) / scanstep
end

function alkane_variance_peak_envelope_dense_grid(left::Real, right::Real, scan_count::Integer)
    gridcount = max(121, 8 * scan_count)
    leftwidth = abs(left)
    rightwidth = abs(right)
    totalwidth = leftwidth + rightwidth
    leftcount = max(2, round(Int, gridcount * leftwidth / totalwidth) + 1)
    rightcount = max(2, gridcount - leftcount + 2)

    uniform = collect(range(Float64(left), Float64(right); length=gridcount))
    leftgrid = if left < 0
        -reverse(alkane_variance_peak_apex_spaced_values(leftwidth, leftcount))
    else
        Float64[0.0]
    end
    rightgrid = if right > 0
        alkane_variance_peak_apex_spaced_values(rightwidth, rightcount)
    else
        Float64[0.0]
    end
    if left < 0
        pop!(leftgrid)
    end
    if right > 0
        popfirst!(rightgrid)
    end

    sort!(unique!(vcat(uniform, leftgrid, 0.0, rightgrid)))
end

function alkane_variance_peak_apex_spaced_values(width::Real, count::Integer)
    scale = Float64(width)
    scale > 0 || return zeros(Float64, count)
    collect(range(0.0, 1.0; length=count)).^1.5 .* scale
end

function alkane_variance_peak_fit_observations(
    peakinput,
    fitrows,
    normalizedretentions,
    scanpositions=axes(peakinput.normalizedobserved, 2);
    minobservations::Integer=3,
)
    alkane_variance_peak_fit_observations_from_response(
        peakinput,
        fitrows,
        normalizedretentions,
        alkane_variance_peak_fit_response(peakinput),
        scanpositions;
        minobservations=minobservations,
    )
end

function alkane_variance_peak_fit_response(peakinput)
    peakinput.normalizedobservedraw
end

function alkane_variance_peak_signal_fit_response(peakinput)
    peakinput.normalizedobserved
end

function alkane_variance_peak_fit_observations_from_response(
    peakinput,
    fitrows,
    normalizedretentions,
    response,
    scanpositions=axes(response, 2);
    minobservations::Integer=3,
)
    size(response) == size(normalizedretentions) || throw(DimensionMismatch(
        "peak fit response and normalized retention matrices must have the same size"))
    observations = NamedTuple[]
    for ionindex in fitrows
        spectrumweight = Float64(peakinput.spectrumweights[ionindex])
        isfinite(spectrumweight) && spectrumweight > 0 || continue
        for scanindex in scanpositions
            retention = Float64(normalizedretentions[ionindex, scanindex])
            intensity = Float64(response[ionindex, scanindex])
            isfinite(retention) && isfinite(intensity) || continue
            push!(observations, (
                ionindex=ionindex,
                scanindex=scanindex,
                retention=retention,
                intensity=intensity,
                spectrumweight=spectrumweight,
            ))
        end
    end
    length(observations) ≥ minobservations || throw(ArgumentError(
        "peak envelope fitting requires at least $minobservations finite " *
        "fit-ion observation(s)"))

    observations
end

function alkane_variance_peak_envelope_predict(
    peakinput,
    normalizedretentions,
    grid,
    envelopevalues,
)
    normalizedfitted = Matrix{Float64}(undef, size(peakinput.normalizedobserved))
    for ionindex in axes(normalizedfitted, 1)
        spectrumweight = Float64(peakinput.spectrumweights[ionindex])
        for scanindex in axes(normalizedfitted, 2)
            retention = Float64(normalizedretentions[ionindex, scanindex])
            if isfinite(retention) && isfinite(spectrumweight)
                normalizedfitted[ionindex, scanindex] =
                    spectrumweight * alkane_variance_envelope_interpolate(
                        grid,
                        envelopevalues,
                        retention,
                    )
            else
                normalizedfitted[ionindex, scanindex] = NaN
            end
        end
    end

    normalizedfitted
end

function alkane_variance_envelope_interpolate(
    grid::AbstractVector{<:Real},
    values::AbstractVector{<:Real},
    retention::Real,
)
    length(grid) == length(values) || throw(DimensionMismatch(
        "envelope grid and values must have the same length"))
    value = Float64(retention)
    value ≤ first(grid) && return Float64(first(values))
    value ≥ last(grid) && return Float64(last(values))

    right = searchsortedfirst(grid, value)
    grid[right] == value && return Float64(values[right])
    left = right - 1
    fraction = (value - grid[left]) / (grid[right] - grid[left])

    max(0.0, (1 - fraction) * Float64(values[left]) +
        fraction * Float64(values[right]))
end

function alkane_variance_rmse(values)
    finite = Float64[value for value in vec(values) if isfinite(value)]
    isempty(finite) && return NaN
    sqrt(mean(abs2, finite))
end

function alkane_variance_finite_count(values)
    count(isfinite, vec(values))
end
