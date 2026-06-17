const DEFAULT_ALKANE_APEX_ION_MIN_RELATIVE_INTENSITY = 0.1
const DEFAULT_ALKANE_APEX_FIT_QUALITY_OUTLIER_Z = 3.0
const DEFAULT_ALKANE_APEX_FIT_QUALITY_MIN_STEPS = 6
const DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES = (
    29,
    41,
    42,
    55,
    56,
    69,
    70,
    73,
    147,
    207,
    281,
    355
)

"""
    alkaneladderapexes(msm, variances, abundanceinfo, pathinfo; mzscanorder=:inferdirection, ...)

Refine alkane ladder apexes from selected ions.

`mzscanorder` controls the m/z timing model used for apex fitting. Use `:ascending` or
`:descending` when the sequential quadrupole scan direction is known. Use `:simultaneous`
for TOF-like or other full-spectrum acquisition where all m/z values in a scan are treated
as observed at the scan-level retention. The default `:inferdirection` means "infer the
sequential quadrupole direction" and currently compares only `:ascending` and
`:descending`; it does not test `:simultaneous`.
"""
function alkaneladderapexes(
    msm::MassScanMatrix,
    variances,
    abundanceinfo,
    pathinfo;
    scanwindow::Integer=2,
    standard=defaultalkanestandard(),
    apexionexcludemzvalues=DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
    apexionmzvalues=nothing,
    apexionminrelativeintensity::Real=DEFAULT_ALKANE_APEX_ION_MIN_RELATIVE_INTENSITY,
    minioncount::Integer=3,
    mzretentionkwargs=nothing,
    mzscanorder::Symbol=:inferdirection,
    variancefloor::Real=1.0,
    logfloorfraction::Real=1e-3,
    maxapexshiftfromguess::Real=3.0,
    apexfitqualityoutlierz=DEFAULT_ALKANE_APEX_FIT_QUALITY_OUTLIER_Z,
    mzscanordermaxpeaks=3,
    mzscanorderminpeaks::Integer=3,
    mzscanorderminapexvarianceratio::Real=1.25,
    mzscanordershapeioncount::Integer=3,
    mzscanordershapemzspacing::Integer=14,
    mzscanorderextremeioncount::Integer=5,
    mzscanorderminioncount::Integer=8
)
    validate_alkane_ladder_apex_settings(
        msm,
        variances,
        scanwindow,
        minioncount,
        apexionminrelativeintensity,
        variancefloor,
        logfloorfraction,
        maxapexshiftfromguess
    )
    hasproperty(abundanceinfo, :abundances) || throw(ArgumentError(
        "abundanceinfo must contain abundances"))
    hasproperty(pathinfo, :path) || throw(ArgumentError("pathinfo must contain path"))

    if isempty(pathinfo.path)
        base_mzkwargs = alkane_ladder_base_mzretention_kwargs(
            msm,
            mzretentionkwargs,
            mzscanorder
        )
        return alkane_ladder_empty_apexinfo(
            :no_path,
            base_mzkwargs,
            alkane_ladder_no_scan_order_info(:no_path);
            scanwindow=scanwindow,
            variancefloor=variancefloor,
            logfloorfraction=logfloorfraction,
            apexionexcludemzvalues=apexionexcludemzvalues,
            apexionmzvalues=apexionmzvalues,
            apexionminrelativeintensity=apexionminrelativeintensity,
            minioncount=minioncount,
            maxapexshiftfromguess=maxapexshiftfromguess,
            apexfitqualityoutlierz=apexfitqualityoutlierz,
            mzscanordermaxpeaks=mzscanordermaxpeaks,
            mzscanorderminpeaks=mzscanorderminpeaks,
            mzscanorderminapexvarianceratio=mzscanorderminapexvarianceratio,
            mzscanordershapeioncount=mzscanordershapeioncount,
            mzscanordershapemzspacing=mzscanordershapemzspacing,
            mzscanorderextremeioncount=mzscanorderextremeioncount,
            mzscanorderminioncount=mzscanorderminioncount
        )
    end

    scanorder = alkaneladderscanorder(
        msm,
        variances,
        pathinfo;
        standard=standard,
        mzretentionkwargs=mzretentionkwargs,
        mzscanorder=mzscanorder,
        scanwindow=scanwindow,
        apexionexcludemzvalues=apexionexcludemzvalues,
        apexionminrelativeintensity=apexionminrelativeintensity,
        variancefloor=variancefloor,
        logfloorfraction=logfloorfraction,
        maxapexshiftfromguess=maxapexshiftfromguess,
        maxpeaks=mzscanordermaxpeaks,
        minpeaks=mzscanorderminpeaks,
        minapexvarianceratio=mzscanorderminapexvarianceratio,
        shapeioncount=mzscanordershapeioncount,
        shapemzspacing=mzscanordershapemzspacing,
        extremeioncount=mzscanorderextremeioncount,
        minioncount=mzscanorderminioncount
    )
    resolved_mzretentionkwargs = scanorder.mzretentionkwargs
    if scanorder.info.selected_order == :unknown
        return alkane_ladder_empty_apexinfo(
            :scan_order_inference_failed,
            resolved_mzretentionkwargs,
            scanorder.info;
            scanwindow=scanwindow,
            variancefloor=variancefloor,
            logfloorfraction=logfloorfraction,
            apexionexcludemzvalues=apexionexcludemzvalues,
            apexionmzvalues=apexionmzvalues,
            apexionminrelativeintensity=apexionminrelativeintensity,
            minioncount=minioncount,
            maxapexshiftfromguess=maxapexshiftfromguess,
            apexfitqualityoutlierz=apexfitqualityoutlierz,
            mzscanordermaxpeaks=mzscanordermaxpeaks,
            mzscanorderminpeaks=mzscanorderminpeaks,
            mzscanorderminapexvarianceratio=mzscanorderminapexvarianceratio,
            mzscanordershapeioncount=mzscanordershapeioncount,
            mzscanordershapemzspacing=mzscanordershapemzspacing,
            mzscanorderextremeioncount=mzscanorderextremeioncount,
            mzscanorderminioncount=mzscanorderminioncount
        )
    end

    apexes = NamedTuple[]
    for candidate in pathinfo.path
        push!(
            apexes,
            alkaneladderapex(
                msm,
                variances,
                abundanceinfo,
                candidate;
                scanwindow=scanwindow,
                standard=standard,
                apexionexcludemzvalues=apexionexcludemzvalues,
                apexionmzvalues=apexionmzvalues,
                apexionminrelativeintensity=apexionminrelativeintensity,
                minioncount=minioncount,
                mzretentionkwargs=resolved_mzretentionkwargs,
                variancefloor=variancefloor,
                logfloorfraction=logfloorfraction,
                maxapexshiftfromguess=maxapexshiftfromguess
            )
        )
    end

    apexes = alkane_ladder_annotate_apex_fit_quality(
        apexes;
        apexfitqualityoutlierz=apexfitqualityoutlierz
    )
    bycarbon = Dict(apex.ladderstep => apex for apex in apexes)
    (
        status=isempty(apexes) ? :failed : :success,
        reason=isempty(apexes) ? :no_path : :success,
        apexes=apexes,
        bycarbon=bycarbon,
        settings=(
            scanwindow=Int(scanwindow),
            variancefloor=Float64(variancefloor),
            logfloorfraction=Float64(logfloorfraction),
            apexionexcludemzvalues=Float64.(collect(apexionexcludemzvalues)),
            apexionmzvalues=isnothing(apexionmzvalues) ?
                nothing :
                Float64.(collect(apexionmzvalues)),
            apexionminrelativeintensity=Float64(apexionminrelativeintensity),
            minioncount=Int(minioncount),
            maxapexshiftfromguess=Float64(maxapexshiftfromguess),
            apexfitqualityoutlierz=isnothing(apexfitqualityoutlierz) ?
                NaN :
                Float64(apexfitqualityoutlierz),
            mzretentionkwargs=resolved_mzretentionkwargs,
            mzscanordermaxpeaks=isnothing(mzscanordermaxpeaks) ?
                nothing :
                Int(mzscanordermaxpeaks),
            mzscanorderminpeaks=Int(mzscanorderminpeaks),
            mzscanorderminapexvarianceratio=Float64(mzscanorderminapexvarianceratio),
            mzscanordershapeioncount=Int(mzscanordershapeioncount),
            mzscanordershapemzspacing=Int(mzscanordershapemzspacing),
            mzscanorderextremeioncount=Int(mzscanorderextremeioncount),
            mzscanorderminioncount=Int(mzscanorderminioncount)
        ),
        scanorderinfo=scanorder.info,
        calibrationexcluded=[apex.calibration_excluded for apex in apexes],
        goodforcalibration=[apex.good_for_calibration for apex in apexes],
        apexfitqualityscores=[apex.apex_fit_quality_score for apex in apexes],
        apexfitqualityzscores=[apex.apex_fit_quality_zscore for apex in apexes]
    )
end

function alkane_ladder_empty_apexinfo(
    reason::Symbol,
    mzretentionkwargs,
    scanorderinfo;
    scanwindow,
    variancefloor,
    logfloorfraction,
    apexionexcludemzvalues,
    apexionmzvalues,
    apexionminrelativeintensity,
    minioncount,
    maxapexshiftfromguess,
    apexfitqualityoutlierz,
    mzscanordermaxpeaks,
    mzscanorderminpeaks,
    mzscanorderminapexvarianceratio,
    mzscanordershapeioncount,
    mzscanordershapemzspacing,
    mzscanorderextremeioncount,
    mzscanorderminioncount
)
    (
        status=:failed,
        reason=reason,
        apexes=NamedTuple[],
        bycarbon=Dict{Int,NamedTuple}(),
        settings=(
            scanwindow=Int(scanwindow),
            variancefloor=Float64(variancefloor),
            logfloorfraction=Float64(logfloorfraction),
            apexionexcludemzvalues=Float64.(collect(apexionexcludemzvalues)),
            apexionmzvalues=isnothing(apexionmzvalues) ?
                nothing :
                Float64.(collect(apexionmzvalues)),
            apexionminrelativeintensity=Float64(apexionminrelativeintensity),
            minioncount=Int(minioncount),
            maxapexshiftfromguess=Float64(maxapexshiftfromguess),
            apexfitqualityoutlierz=isnothing(apexfitqualityoutlierz) ?
                NaN :
                Float64(apexfitqualityoutlierz),
            mzretentionkwargs=mzretentionkwargs,
            mzscanordermaxpeaks=isnothing(mzscanordermaxpeaks) ?
                nothing :
                Int(mzscanordermaxpeaks),
            mzscanorderminpeaks=Int(mzscanorderminpeaks),
            mzscanorderminapexvarianceratio=Float64(mzscanorderminapexvarianceratio),
            mzscanordershapeioncount=Int(mzscanordershapeioncount),
            mzscanordershapemzspacing=Int(mzscanordershapemzspacing),
            mzscanorderextremeioncount=Int(mzscanorderextremeioncount),
            mzscanorderminioncount=Int(mzscanorderminioncount)
        ),
        scanorderinfo=scanorderinfo,
        calibrationexcluded=Bool[],
        goodforcalibration=Bool[],
        apexfitqualityscores=Float64[],
        apexfitqualityzscores=Float64[]
    )
end

function alkaneladderapexes(
    msm::MassScanMatrix,
    abundanceinfo,
    pathinfo;
    kwargs...
)
    alkaneladderapexes(
        msm,
        ones(size(rawintensities(msm))),
        abundanceinfo,
        pathinfo;
        kwargs...
    )
end

function validate_alkane_ladder_apex_settings(
    scanwindow,
    variancefloor,
    logfloorfraction
)
    scanwindow isa Integer || throw(ArgumentError("scanwindow must be an integer"))
    scanwindow >= 1 || throw(ArgumentError("scanwindow must be at least 1"))
    validate_alkane_abundance_variancefloor(variancefloor)
    isfinite(logfloorfraction) && logfloorfraction > 0 || throw(ArgumentError(
        "logfloorfraction must be finite and positive"))

    nothing
end

function validate_alkane_ladder_apex_settings(
    msm::MassScanMatrix,
    variances,
    scanwindow,
    minioncount,
    apexionminrelativeintensity,
    variancefloor,
    logfloorfraction,
    maxapexshiftfromguess
)
    validate_alkane_series_variances(msm, variances)
    validate_alkane_ladder_apex_settings(scanwindow, variancefloor, logfloorfraction)
    minioncount isa Integer || throw(ArgumentError("minioncount must be an integer"))
    minioncount >= 1 || throw(ArgumentError("minioncount must be at least 1"))
    isfinite(apexionminrelativeintensity) &&
        0 <= apexionminrelativeintensity < 1 ||
        throw(ArgumentError(
            "apexionminrelativeintensity must be finite and in [0, 1)"))
    isfinite(maxapexshiftfromguess) && maxapexshiftfromguess >= 0 ||
        throw(ArgumentError(
            "maxapexshiftfromguess must be finite and nonnegative"))

    nothing
end

function alkaneladderapex(
    msm::MassScanMatrix,
    variances,
    abundanceinfo,
    candidate;
    scanwindow::Integer=2,
    standard=defaultalkanestandard(),
    apexionexcludemzvalues=DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
    apexionmzvalues=nothing,
    apexionminrelativeintensity::Real=DEFAULT_ALKANE_APEX_ION_MIN_RELATIVE_INTENSITY,
    minioncount::Integer=3,
    mzretentionkwargs=nothing,
    variancefloor::Real=1.0,
    logfloorfraction::Real=1e-3,
    maxapexshiftfromguess::Real=3.0
)
    validate_alkane_ladder_apex_settings(
        msm,
        variances,
        scanwindow,
        minioncount,
        apexionminrelativeintensity,
        variancefloor,
        logfloorfraction,
        maxapexshiftfromguess
    )

    retentions = Float64.(rawretentions(msm))
    alkane_validate_retention_axis(retentions)
    step = alkane_ladder_candidate_step(candidate)
    abundance = alkane_ladder_candidate_abundance(abundanceinfo, step, retentions)
    inputscan = alkane_ladder_input_scan_index(candidate)

    try
        apex = alkane_ladder_ion_apex(
            msm,
            variances,
            candidate;
            scanwindow=scanwindow,
            standard=standard,
            ladderstep=step,
            apexionexcludemzvalues=apexionexcludemzvalues,
            apexionmzvalues=apexionmzvalues,
            apexionminrelativeintensity=apexionminrelativeintensity,
            minioncount=minioncount,
            mzretentionkwargs=mzretentionkwargs,
            variancefloor=variancefloor,
            logfloorfraction=logfloorfraction,
            maxapexshiftfromguess=maxapexshiftfromguess
        )
        return alkane_ladder_apex_public_result(
            apex,
            candidate,
            abundance,
            retentions
        )
    catch err
        (err isa ArgumentError || err isa DimensionMismatch) || rethrow()
        return alkane_ladder_discrete_apex(
            candidate,
            retentions,
            abundance;
            scanindices=alkane_ladder_fallback_scanindices(
                retentions,
                inputscan,
                Int(scanwindow),
                candidate
            ),
            reason=:apex_fit_failed,
            failurereason=sprint(showerror, err)
        )
    end
end

function alkaneladderapex(
    msm::MassScanMatrix,
    abundanceinfo,
    candidate;
    kwargs...
)
    alkaneladderapex(
        msm,
        ones(size(rawintensities(msm))),
        abundanceinfo,
        candidate;
        kwargs...
    )
end

function alkane_ladder_ion_apex(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidate;
    scanwindow::Integer,
    standard,
    ladderstep::Integer,
    apexionexcludemzvalues,
    apexionmzvalues,
    apexionminrelativeintensity::Real,
    minioncount::Integer,
    mzretentionkwargs,
    variancefloor::Real,
    logfloorfraction::Real,
    maxapexshiftfromguess::Real
)
    input_scanindex = alkane_ladder_input_scan_index(candidate)
    raw_scan_retentions = Float64.(rawretentions(msm))
    1 <= input_scanindex <= length(raw_scan_retentions) || throw(ArgumentError(
        "candidate scan index $(input_scanindex) is outside 1:$(length(raw_scan_retentions))"))

    attempts = NamedTuple[]
    function fit_at_center!(center_retention)
        apex = alkane_ladder_ion_apex_once(
            msm,
            variances,
            candidate;
            scanwindow=scanwindow,
            standard=standard,
            ladderstep=ladderstep,
            apexionexcludemzvalues=apexionexcludemzvalues,
            apexionmzvalues=apexionmzvalues,
            apexionminrelativeintensity=apexionminrelativeintensity,
            minioncount=minioncount,
            mzretentionkwargs=mzretentionkwargs,
            variancefloor=variancefloor,
            logfloorfraction=logfloorfraction,
            fitcenterretention=center_retention,
            maxapexshiftfromguess=maxapexshiftfromguess,
            attemptindex=length(attempts) + 1
        )
        push!(attempts, apex)

        apex
    end

    center_retention = raw_scan_retentions[input_scanindex]
    for _ in 1:alkane_ladder_max_center_passes(maxapexshiftfromguess)
        apex = fit_at_center!(center_retention)
        apex.success && alkane_ladder_apex_is_centered(apex) && break

        next_center = alkane_ladder_next_center_retention(
            apex,
            raw_scan_retentions,
            input_scanindex,
            maxapexshiftfromguess
        )
        isnothing(next_center) && break
        next_center_scan = alkane_ladder_fractional_scan_index(
            raw_scan_retentions,
            next_center
        )
        abs(next_center_scan - apex.fit_center_scan_index) < 0.25 && break
        center_retention = next_center
    end

    if !alkane_ladder_has_centered_success(attempts)
        for fallback_center in alkane_ladder_center_search_retentions(
            raw_scan_retentions,
            input_scanindex,
            maxapexshiftfromguess
        )
            fallback_scan = alkane_ladder_fractional_scan_index(
                raw_scan_retentions,
                fallback_center
            )
            alkane_ladder_center_already_attempted(attempts, fallback_scan) &&
                continue

            apex = fit_at_center!(fallback_center)
            apex.success && alkane_ladder_apex_is_centered(apex) && break
        end
    end

    final_index = alkane_ladder_best_apex_attempt_index(attempts)
    merge(
        attempts[final_index],
        (
            initial_apex=first(attempts),
            all_apex_attempts=attempts,
            recenter_attempts=length(attempts) - 1,
            recenter_used=final_index > 1,
            maxapexshiftfromguess=Float64(maxapexshiftfromguess)
        )
    )
end

function alkane_ladder_ion_apex_once(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidate;
    scanwindow::Integer,
    standard,
    ladderstep::Integer,
    apexionexcludemzvalues,
    apexionmzvalues,
    apexionminrelativeintensity::Real,
    minioncount::Integer,
    mzretentionkwargs,
    variancefloor::Real,
    logfloorfraction::Real,
    fitcenterretention,
    maxapexshiftfromguess::Real,
    attemptindex::Integer
)
    input_scanindex = alkane_ladder_input_scan_index(candidate)
    raw_scan_retentions = Float64.(rawretentions(msm))
    input_retention = raw_scan_retentions[input_scanindex]
    fit_center_retention = Float64(fitcenterretention)
    isfinite(fit_center_retention) || throw(ArgumentError(
        "fitcenterretention must be finite"))
    fit_center_scan_index = alkane_ladder_fractional_scan_index(
        raw_scan_retentions,
        fit_center_retention
    )
    localmargin = max(3, Int(scanwindow) + 2)

    nscans = min(scancount(msm), 2 * Int(scanwindow) + 1)
    model_scanindices = alkane_ladder_scan_indices_around_retention(
        raw_scan_retentions,
        fit_center_retention,
        nscans
    )
    length(model_scanindices) >= 3 || throw(ArgumentError(
        "at least three scans are needed for a log-quadratic apex fit"))

    ion_selection = isnothing(apexionmzvalues) ?
        alkane_ladder_reference_filtered_apex_mzvalues(
            msm,
            standard,
            ladderstep;
            excludemzvalues=apexionexcludemzvalues,
            minrelativeintensity=apexionminrelativeintensity,
            minioncount=minioncount
        ) :
        alkane_ladder_explicit_apex_mzvalue_selection(
            msm,
            apexionmzvalues;
            minioncount=minioncount
        )
    selected_mzindices = alkane_ladder_apex_mz_indices(
        msm,
        ion_selection.mzvalues;
        minioncount=minioncount
    )

    mzkwargs = isnothing(mzretentionkwargs) ?
        alkane_ladder_default_mzretention_kwargs(msm) :
        mzretentionkwargs
    scan_retentions = retentions(msm)
    nmz = length(selected_mzindices)
    scanindices_by_mz = alkane_ladder_scan_indices_by_mz(
        msm,
        scan_retentions,
        fit_center_retention,
        nscans,
        selected_mzindices,
        mzkwargs,
        fit_center_scan_index,
        localmargin
    )

    Xraw = rawintensities(msm)
    intensities_by_mz = Matrix{Float64}(undef, nscans, nmz)
    variances_by_mz = Matrix{Float64}(undef, nscans, nmz)
    observation_retentions = Matrix{Float64}(undef, nscans, nmz)
    raw_retentions_by_mz = Matrix{Float64}(undef, nscans, nmz)
    for mzcol in 1:nmz, scancol in 1:nscans
        selected_scan = scanindices_by_mz[scancol, mzcol]
        selected_mz = selected_mzindices[mzcol]
        intensities_by_mz[scancol, mzcol] =
            max(Float64(Xraw[selected_scan, selected_mz]), 0.0)
        variances_by_mz[scancol, mzcol] =
            Float64(variances[selected_scan, selected_mz])
        raw_retentions_by_mz[scancol, mzcol] =
            raw_scan_retentions[selected_scan]
        observation_retentions[scancol, mzcol] =
            alkane_ladder_observation_retention(
                msm,
                scan_retentions,
                selected_scan,
                selected_mz,
                mzkwargs
            )
    end

    ymax = maximum(intensities_by_mz)
    ymax > 0 || throw(ArgumentError("selected ions have no positive local signal"))
    logfloor = Float64(logfloorfraction) * ymax

    xvalues = observation_retentions .- fit_center_retention
    xscale = maximum(abs, xvalues)
    xscale > 0 || throw(ArgumentError("scan window has zero retention span"))

    nobs = nscans * nmz
    design = zeros(Float64, nobs, nmz + 2)
    logy = Vector{Float64}(undef, nobs)
    weights = Vector{Float64}(undef, nobs)
    row = 0
    for mzcol in 1:nmz
        for scancol in 1:nscans
            row += 1
            x = xvalues[scancol, mzcol] / xscale
            adjusted = intensities_by_mz[scancol, mzcol] + logfloor
            design[row, mzcol] = 1.0
            design[row, nmz + 1] = x
            design[row, nmz + 2] = abs2(x)
            logy[row] = log(adjusted)
            weights[row] = alkane_ladder_log_weight(
                adjusted,
                variances_by_mz[scancol, mzcol],
                variancefloor
            )
        end
    end
    any(>(0), weights) || throw(ArgumentError(
        "peak apex fit has no positive weights"))

    sqrtweights = sqrt.(weights)
    coefficients = (design .* sqrtweights) \ (logy .* sqrtweights)
    fitted_logy = design * coefficients
    beta = coefficients[nmz + 1]
    gamma = coefficients[nmz + 2]
    apexx = gamma < 0 ? -beta / (2 * gamma) : NaN
    apex_retention = isfinite(apexx) ? fit_center_retention + apexx * xscale : NaN
    apex_scan_index = isfinite(apex_retention) ?
        alkane_ladder_fractional_scan_index(raw_scan_retentions, apex_retention) :
        NaN
    apex_in_window = isfinite(apexx) &&
        minimum(xvalues) / xscale <= apexx <= maximum(xvalues) / xscale
    apex_shift_from_guess_scans = apex_scan_index - input_scanindex
    apex_within_allowed_shift = isfinite(apex_shift_from_guess_scans) &&
        abs(apex_shift_from_guess_scans) <= Float64(maxapexshiftfromguess) + 1e-9
    success = gamma < 0 &&
        apex_in_window &&
        isfinite(apex_scan_index) &&
        apex_within_allowed_shift

    logmean = sum(weights .* logy) / sum(weights)
    ssres = sum(weights .* abs2.(logy .- fitted_logy))
    sstot = sum(weights .* abs2.(logy .- logmean))
    r2 = sstot > 0 ? 1 - ssres / sstot : NaN
    fit_degrees_of_freedom = nobs - (nmz + 2)
    reduced_normalized_residual =
        fit_degrees_of_freedom > 0 ? ssres / fit_degrees_of_freedom : NaN

    normalized_log_residuals = Matrix{Float64}(undef, nscans, nmz)
    fitted_log_intensities = Matrix{Float64}(undef, nscans, nmz)
    row = 0
    for mzcol in 1:nmz
        for scancol in 1:nscans
            row += 1
            normalized_log_residuals[scancol, mzcol] =
                sqrt(weights[row]) * (logy[row] - fitted_logy[row])
            fitted_log_intensities[scancol, mzcol] = fitted_logy[row]
        end
    end

    model_x =
        (raw_scan_retentions[model_scanindices] .- fit_center_retention) ./ xscale
    peak_model = success ?
        exp.(gamma .* abs2.(model_x .- apexx)) :
        fill(Float64(NaN), length(model_x))
    observation_peak_model = success ?
        exp.(gamma .* abs2.(xvalues ./ xscale .- apexx)) :
        fill(Float64(NaN), size(xvalues))

    (
        success=success,
        apex_retention=apex_retention,
        apex_time=apex_retention,
        apex_scan_index=apex_scan_index,
        input_scan_index=input_scanindex,
        input_retention=input_retention,
        apex_offset_retention=apex_retention - input_retention,
        apex_offset_scans=apex_scan_index - input_scanindex,
        fit_center_retention=fit_center_retention,
        fit_center_scan_index=fit_center_scan_index,
        apex_offset_from_fit_center_retention=apex_retention - fit_center_retention,
        apex_offset_from_fit_center_scans=apex_scan_index - fit_center_scan_index,
        apex_shift_from_guess_scans=apex_shift_from_guess_scans,
        apex_within_allowed_shift=apex_within_allowed_shift,
        apex_in_window=apex_in_window,
        attemptindex=Int(attemptindex),
        beta=beta,
        gamma=gamma,
        apex_x=apexx,
        x_scale=xscale,
        r2=r2,
        weighted_log_residual_sum_squares=ssres,
        fit_degrees_of_freedom=fit_degrees_of_freedom,
        reduced_normalized_residual=reduced_normalized_residual,
        normalized_log_residuals=normalized_log_residuals,
        fitted_log_intensities=fitted_log_intensities,
        scanwindow=Int(scanwindow),
        scan_indices=model_scanindices,
        retentions=raw_scan_retentions[model_scanindices],
        peak_model=peak_model,
        observation_peak_model=observation_peak_model,
        scan_indices_by_mz=scanindices_by_mz,
        raw_retentions_by_mz=raw_retentions_by_mz,
        observation_retentions=observation_retentions,
        fit_intensities=intensities_by_mz,
        fit_variances=variances_by_mz,
        mz_indices=selected_mzindices,
        mz_values=collect(rawmzvalues(msm)[selected_mzindices]),
        apex_ion_selection=ion_selection.selection,
        apex_ion_excluded_mzvalues=ion_selection.excluded_mzvalues,
        apex_ion_reference_mzvalues=ion_selection.reference_mzvalues,
        apex_ion_reference_relative_intensities=
            ion_selection.reference_relative_intensities,
        apex_ion_min_relative_intensity=Float64(apexionminrelativeintensity),
        apex_ion_reference_key=ion_selection.reference_key,
        ion_intercepts=coefficients[1:nmz],
        log_floor=logfloor,
        logfloorfraction=Float64(logfloorfraction),
        variance_weighted=true,
        variancefloor=Float64(variancefloor),
        mzretentionkwargs=mzkwargs
    )
end

function alkane_ladder_apex_public_result(
    apex,
    candidate,
    abundance::AbstractVector{Float64},
    retentions::AbstractVector{Float64}
)
    step = alkane_ladder_candidate_step(candidate)
    scanindex = apex.input_scan_index
    reason = apex.success ? :success : Symbol(
        replace(alkane_ladder_apex_failure_reason(apex), ' ' => '_')
    )
    failurereason = apex.success ? nothing : alkane_ladder_apex_failure_reason(apex)
    merge(
        apex,
        (
            reason=reason,
            failurereason=failurereason,
            refined=Bool(apex.success),
            fallback=!Bool(apex.success),
            ladderstep=step,
            scanindex=scanindex,
            retention=retentions[scanindex],
            apexscanindex=apex.apex_scan_index,
            apexretention=apex.apex_retention,
            output_scan_index=apex.apex_scan_index,
            output_retention=apex.apex_retention,
            outputscanindex=apex.apex_scan_index,
            outputretention=apex.apex_retention,
            apexoffsetscans=apex.apex_offset_scans,
            apexoffsetretention=apex.apex_offset_retention,
            apexabundance=abundance[scanindex],
            scanindices=apex.scan_indices,
            abundance=abundance[apex.scan_indices],
            xscale=apex.x_scale,
            weightedlogresidualsumsquares=apex.weighted_log_residual_sum_squares,
            fit_coefficients=vcat(apex.ion_intercepts, apex.beta, apex.gamma),
            logfloor=apex.log_floor,
            source=hasproperty(candidate, :source) ? candidate.source : :molecular_ion_dp,
            gapfilled=hasproperty(candidate, :gapfilled) ? Bool(candidate.gapfilled) :
                hasproperty(candidate, :source) && candidate.source == :gapfilled,
            edgeextended=hasproperty(candidate, :edgeextended) ?
                Bool(candidate.edgeextended) :
                hasproperty(candidate, :source) &&
                    candidate.source in (:leftextended, :rightextended),
            mass_spectrum_cosine=hasproperty(candidate, :massspectrumcosine) ?
                Float64(candidate.massspectrumcosine) :
                NaN,
            required_cosine=hasproperty(candidate, :requiredcosine) ?
                Float64(candidate.requiredcosine) :
                NaN,
            apex_fit_quality_score=apex.reduced_normalized_residual,
            apex_fit_quality_log_score=(
                isfinite(apex.reduced_normalized_residual) &&
                    apex.reduced_normalized_residual > 0
            ) ? log(apex.reduced_normalized_residual) : NaN,
            apex_fit_quality_zscore=NaN,
            apex_fit_quality_zscore_excluded=false,
            apex_fit_quality_outlier_z=NaN,
            apex_fit_quality_baseline_median=NaN,
            apex_fit_quality_baseline_mad=NaN,
            apex_fit_quality_nsteps=0,
            apex_fit_quality_applied=false,
            calibration_excluded=false,
            calibration_exclusion_reason=nothing,
            good_for_calibration=Bool(apex.success),
            candidate=candidate
        )
    )
end

function alkane_ladder_annotate_apex_fit_quality(
    apexes::AbstractVector;
    apexfitqualityoutlierz=DEFAULT_ALKANE_APEX_FIT_QUALITY_OUTLIER_Z
)
    if !isnothing(apexfitqualityoutlierz)
        isfinite(apexfitqualityoutlierz) && apexfitqualityoutlierz > 0 ||
            throw(ArgumentError(
                "apexfitqualityoutlierz must be finite and positive, or nothing"))
    end

    qualityindices = Int[]
    logscores = Float64[]
    for (index, apex) in pairs(apexes)
        score = alkane_ladder_apex_fit_quality_score(apex)
        if hasproperty(apex, :refined) && apex.refined && isfinite(score) && score > 0
            push!(qualityindices, Int(index))
            push!(logscores, log(score))
        end
    end

    quality = alkane_ladder_apex_fit_quality_robust_zscores(logscores)
    zbyindex = Dict{Int, Float64}()
    for (index, zscore) in zip(qualityindices, quality.zscores)
        zbyindex[index] = zscore
    end

    updated = NamedTuple[]
    for (index, apex) in pairs(apexes)
        score = alkane_ladder_apex_fit_quality_score(apex)
        zscore = get(zbyindex, Int(index), NaN)
        refined = hasproperty(apex, :refined) ? Bool(apex.refined) : Bool(apex.success)
        zexcluded = !isnothing(apexfitqualityoutlierz) &&
            refined &&
            quality.applied &&
            isfinite(zscore) &&
            zscore > Float64(apexfitqualityoutlierz)
        calibrationexcluded = zexcluded
        push!(
            updated,
            merge(
                apex,
                (
                    apex_fit_quality_score=score,
                    apex_fit_quality_log_score=
                        isfinite(score) && score > 0 ? log(score) : NaN,
                    apex_fit_quality_zscore=zscore,
                    apex_fit_quality_zscore_excluded=zexcluded,
                    apex_fit_quality_outlier_z=isnothing(apexfitqualityoutlierz) ?
                        NaN :
                        Float64(apexfitqualityoutlierz),
                    apex_fit_quality_baseline_median=quality.median,
                    apex_fit_quality_baseline_mad=quality.mad,
                    apex_fit_quality_nsteps=quality.nsteps,
                    apex_fit_quality_applied=quality.applied,
                    calibration_excluded=calibrationexcluded,
                    calibration_exclusion_reason=calibrationexcluded ?
                        "apex curve fit residual outlier" :
                        nothing,
                    good_for_calibration=refined && !calibrationexcluded
                )
            )
        )
    end

    updated
end

function alkane_ladder_apex_fit_quality_score(apex)
    hasproperty(apex, :apex_fit_quality_score) &&
        return Float64(apex.apex_fit_quality_score)
    hasproperty(apex, :reduced_normalized_residual) &&
        return Float64(apex.reduced_normalized_residual)

    NaN
end

function alkane_ladder_apex_fit_quality_robust_zscores(
    logscores::AbstractVector{<:Real}
)
    nsteps = length(logscores)
    zscores = fill(NaN, nsteps)
    nsteps >= DEFAULT_ALKANE_APEX_FIT_QUALITY_MIN_STEPS ||
        return (
            zscores=zscores,
            median=NaN,
            mad=NaN,
            nsteps=nsteps,
            applied=false
        )
    values = Float64.(collect(logscores))
    center = median(values)
    deviations = abs.(values .- center)
    mad = median(deviations)
    isfinite(mad) && mad > eps(Float64) ||
        return (
            zscores=zscores,
            median=center,
            mad=mad,
            nsteps=nsteps,
            applied=false
        )

    zscores .= 0.67448975 .* (values .- center) ./ mad
    (
        zscores=zscores,
        median=center,
        mad=mad,
        nsteps=nsteps,
        applied=true
    )
end

"""
    alkaneladderscanorder(msm, variances, pathinfo; mzscanorder=:inferdirection, ...)

Resolve the m/z scan timing model used by alkane ladder apex fitting.

If `mzscanorder` or `mzretentionkwargs.order` is `:ascending`, `:descending`, or
`:simultaneous`, that order is used directly and the returned status is `:provided`.
If it is `:inferdirection`, the function infers only the sequential quadrupole direction by
comparing `:ascending` and `:descending`. Simultaneous acquisition is not part of direction
inference; users with TOF-like/full-spectrum data should pass `mzscanorder=:simultaneous`
or `mzretentionkwargs=(..., order=:simultaneous)`.
"""
function alkaneladderscanorder(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    pathinfo;
    standard=defaultalkanestandard(),
    mzretentionkwargs=nothing,
    mzscanorder::Symbol=:inferdirection,
    scanwindow::Integer=2,
    apexionexcludemzvalues=DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
    apexionminrelativeintensity::Real=DEFAULT_ALKANE_APEX_ION_MIN_RELATIVE_INTENSITY,
    variancefloor::Real=1.0,
    logfloorfraction::Real=1e-3,
    maxapexshiftfromguess::Real=3.0,
    maxpeaks=3,
    minpeaks::Integer=3,
    minapexvarianceratio::Real=1.25,
    shapeioncount::Integer=3,
    shapemzspacing::Integer=14,
    extremeioncount::Integer=5,
    minioncount::Integer=8
)
    validate_alkane_series_variances(msm, variances)
    hasproperty(pathinfo, :path) || throw(ArgumentError("pathinfo must contain path"))
    validate_alkane_ladder_scan_order_settings(
        maxpeaks,
        minpeaks,
        minapexvarianceratio,
        shapeioncount,
        shapemzspacing,
        extremeioncount,
        minioncount
    )

    base_mzkwargs = alkane_ladder_base_mzretention_kwargs(
        msm,
        mzretentionkwargs,
        mzscanorder
    )
    order = alkane_ladder_mzretention_order(base_mzkwargs)
    if order != :inferdirection
        order in (:ascending, :descending, :simultaneous) || throw(ArgumentError(
            "mzretentionkwargs.order must be :ascending, :descending, or :simultaneous"))
        return (
            mzretentionkwargs=base_mzkwargs,
            info=alkane_ladder_provided_scan_order_info(order)
        )
    end

    candidates = alkane_ladder_scan_order_candidates(pathinfo; maxpeaks=maxpeaks)
    initialcount = isnothing(maxpeaks) ?
        length(candidates) :
        min(Int(maxpeaks), length(candidates))
    ascending, descending, info = alkane_ladder_scan_order_expanding_trials(
        msm,
        variances,
        candidates,
        base_mzkwargs;
        initialcount=initialcount,
        minpeaks=minpeaks,
        minapexvarianceratio=minapexvarianceratio,
        standard=standard,
        scanwindow=scanwindow,
        apexionexcludemzvalues=apexionexcludemzvalues,
        apexionminrelativeintensity=apexionminrelativeintensity,
        variancefloor=variancefloor,
        logfloorfraction=logfloorfraction,
        maxapexshiftfromguess=maxapexshiftfromguess,
        shapeioncount=shapeioncount,
        shapemzspacing=shapemzspacing,
        extremeioncount=extremeioncount,
        minioncount=minioncount
    )

    (
        mzretentionkwargs=merge(base_mzkwargs, (order=info.selected_order,)),
        info=info
    )
end

function validate_alkane_ladder_scan_order_settings(
    maxpeaks,
    minpeaks,
    minapexvarianceratio,
    shapeioncount,
    shapemzspacing,
    extremeioncount,
    minioncount
)
    isnothing(maxpeaks) ||
        (maxpeaks isa Integer && maxpeaks >= 1) ||
        throw(ArgumentError(
            "scan-order maxpeaks must be nothing or a positive integer"))
    minpeaks isa Integer && minpeaks >= 1 || throw(ArgumentError(
        "scan-order minpeaks must be a positive integer"))
    isfinite(minapexvarianceratio) && minapexvarianceratio > 0 || throw(ArgumentError(
        "scan-order minapexvarianceratio must be finite and positive"))
    shapeioncount isa Integer && shapeioncount >= 2 || throw(ArgumentError(
        "scan-order shapeioncount must be an integer >= 2"))
    shapemzspacing isa Integer && shapemzspacing >= 1 || throw(ArgumentError(
        "scan-order shapemzspacing must be a positive integer"))
    extremeioncount isa Integer && extremeioncount >= 1 || throw(ArgumentError(
        "scan-order extremeioncount must be a positive integer"))
    minioncount isa Integer && minioncount >= 2 || throw(ArgumentError(
        "scan-order minioncount must be an integer >= 2"))

    nothing
end

function alkane_ladder_base_mzretention_kwargs(
    msm::MassScanMatrix,
    mzretentionkwargs,
    mzscanorder::Symbol
)
    mzscanorder in (:inferdirection, :ascending, :descending, :simultaneous) ||
        throw(ArgumentError(
            "mzscanorder must be :inferdirection, :ascending, :descending, or :simultaneous"))

    if isnothing(mzretentionkwargs)
        return alkane_ladder_default_mzretention_kwargs(msm; order=mzscanorder)
    end

    if hasproperty(mzretentionkwargs, :order)
        mzorder = getproperty(mzretentionkwargs, :order)
        mzorder in (:ascending, :descending, :simultaneous) ||
            throw(ArgumentError(
                "mzretentionkwargs.order must be :ascending, :descending, or :simultaneous"))
        if mzscanorder != :inferdirection && mzscanorder != mzorder
            throw(ArgumentError(
                "mzscanorder=$(repr(mzscanorder)) conflicts with " *
                "mzretentionkwargs.order=$(repr(mzorder))"))
        end

        return mzretentionkwargs
    end

    merge(mzretentionkwargs, (order=mzscanorder,))
end

function alkane_ladder_provided_scan_order_info(order::Symbol)
    (
        selected_order=order,
        requested_order=order,
        status=:provided,
        score_kind=:provided,
        evidence_score=nothing,
        ascending_score=NaN,
        descending_score=NaN,
        ascending_median_apex_variance=NaN,
        descending_median_apex_variance=NaN,
        ascending_median_width=NaN,
        descending_median_width=NaN,
        n_peaks_used=0,
        n_peaks_tried=0,
        initial_n_peaks_tried=0,
        expanded_for_ambiguity=false,
        trials=(ascending=nothing, descending=nothing)
    )
end

function alkane_ladder_no_scan_order_info(reason::Symbol)
    (
        selected_order=:unknown,
        requested_order=:inferdirection,
        status=reason,
        score_kind=:not_applicable,
        evidence_score=nothing,
        ascending_score=nothing,
        descending_score=nothing,
        ascending_median_apex_variance=nothing,
        descending_median_apex_variance=nothing,
        ascending_median_width=nothing,
        descending_median_width=nothing,
        n_peaks_used=0,
        n_peaks_tried=0,
        initial_n_peaks_tried=0,
        expanded_for_ambiguity=false,
        trials=(ascending=nothing, descending=nothing)
    )
end


function alkane_ladder_mzretention_order(mzretentionkwargs)
    hasproperty(mzretentionkwargs, :order) || throw(ArgumentError(
        "mzretentionkwargs must contain an order field"))

    getproperty(mzretentionkwargs, :order)
end

function alkane_ladder_scan_order_candidates(pathinfo; maxpeaks)
    pairs = NamedTuple[]
    for candidate in pathinfo.path
        push!(
            pairs,
            (
                ladderstep=alkane_ladder_candidate_step(candidate),
                scanindex=alkane_ladder_input_scan_index(candidate),
                candidate=candidate
            )
        )
    end
    sort!(pairs; by=pair -> pair.ladderstep)
    isnothing(maxpeaks) && return pairs
    length(pairs) <= maxpeaks && return pairs

    order = alkane_ladder_spread_index_order(length(pairs), Int(maxpeaks))
    pairs[order]
end

function alkane_ladder_spread_index_order(n::Integer, initialcount::Integer)
    n >= 0 || throw(ArgumentError("number of candidates must be nonnegative"))
    n == 0 && return Int[]
    initialcount >= 1 || throw(ArgumentError("initialcount must be positive"))

    selected = Int[]
    seen = Set{Int}()
    firstcount = min(Int(initialcount), Int(n))
    for index in alkane_ladder_spread_indices(Int(n), firstcount)
        push!(selected, index)
        push!(seen, index)
    end

    while length(selected) < Int(n)
        bestindex = 0
        bestdistance = -1
        for index in 1:Int(n)
            index in seen && continue
            distance = minimum(abs(index - selectedindex) for selectedindex in selected)
            if distance > bestdistance
                bestdistance = distance
                bestindex = index
            end
        end
        bestindex == 0 && break
        push!(selected, bestindex)
        push!(seen, bestindex)
    end

    selected
end

function alkane_ladder_spread_indices(n::Integer, count::Integer)
    count = min(Int(count), Int(n))
    count <= 0 && return Int[]
    count == 1 && return [ceil(Int, Int(n) / 2)]
    indices = [
        clamp(round(Int, value), 1, Int(n))
        for value in range(1, Int(n); length=count)
    ]
    seen = Set{Int}()
    output = Int[]
    for index in indices
        index in seen && continue
        push!(output, index)
        push!(seen, index)
    end

    candidate = 1
    while length(output) < count && candidate <= Int(n)
        if !(candidate in seen)
            push!(output, candidate)
            push!(seen, candidate)
        end
        candidate += 1
    end

    output
end

function alkane_ladder_scan_order_expanding_trials(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidates::AbstractVector,
    base_mzretentionkwargs;
    initialcount::Integer,
    minpeaks::Integer,
    minapexvarianceratio::Real,
    standard,
    scanwindow::Integer,
    apexionexcludemzvalues,
    apexionminrelativeintensity::Real,
    variancefloor::Real,
    logfloorfraction::Real,
    maxapexshiftfromguess::Real,
    shapeioncount::Integer,
    shapemzspacing::Integer,
    extremeioncount::Integer,
    minioncount::Integer
)
    tested = NamedTuple[]
    ascending_results = NamedTuple[]
    descending_results = NamedTuple[]
    initialcount = min(Int(initialcount), length(candidates))

    last = nothing
    expanded_for_ambiguity = false
    for candidateinfo in candidates
        push!(tested, candidateinfo)
        push!(
            ascending_results,
            alkane_ladder_scan_order_trial_candidate(
                msm,
                variances,
                candidateinfo,
                :ascending,
                base_mzretentionkwargs;
                standard=standard,
                scanwindow=scanwindow,
                apexionexcludemzvalues=apexionexcludemzvalues,
                apexionminrelativeintensity=apexionminrelativeintensity,
                variancefloor=variancefloor,
                logfloorfraction=logfloorfraction,
                maxapexshiftfromguess=maxapexshiftfromguess,
                shapeioncount=shapeioncount,
                shapemzspacing=shapemzspacing,
                extremeioncount=extremeioncount,
                minioncount=minioncount
            )
        )
        push!(
            descending_results,
            alkane_ladder_scan_order_trial_candidate(
                msm,
                variances,
                candidateinfo,
                :descending,
                base_mzretentionkwargs;
                standard=standard,
                scanwindow=scanwindow,
                apexionexcludemzvalues=apexionexcludemzvalues,
                apexionminrelativeintensity=apexionminrelativeintensity,
                variancefloor=variancefloor,
                logfloorfraction=logfloorfraction,
                maxapexshiftfromguess=maxapexshiftfromguess,
                shapeioncount=shapeioncount,
                shapemzspacing=shapemzspacing,
                extremeioncount=extremeioncount,
                minioncount=minioncount
            )
        )

        length(tested) >= initialcount || continue

        ascending = alkane_ladder_scan_order_trial_summary(
            :ascending,
            tested,
            ascending_results
        )
        descending = alkane_ladder_scan_order_trial_summary(
            :descending,
            tested,
            descending_results
        )
        info = alkane_ladder_choose_scan_order(
            ascending,
            descending;
            minpeaks=minpeaks,
            minapexvarianceratio=minapexvarianceratio
        )
        info = alkane_ladder_scan_order_info_with_expansion(
            info,
            initialcount,
            expanded_for_ambiguity
        )
        last = (ascending, descending, info)

        if expanded_for_ambiguity
            length(tested) == length(candidates) && return last
            continue
        end

        info.status == :accepted && return last
        if info.status == :ambiguous
            length(tested) == length(candidates) && return last
            expanded_for_ambiguity = true
        end
    end

    if isnothing(last)
        ascending = alkane_ladder_scan_order_trial_summary(
            :ascending,
            tested,
            ascending_results
        )
        descending = alkane_ladder_scan_order_trial_summary(
            :descending,
            tested,
            descending_results
        )
        info = alkane_ladder_choose_scan_order(
            ascending,
            descending;
            minpeaks=minpeaks,
            minapexvarianceratio=minapexvarianceratio
        )
        info = alkane_ladder_scan_order_info_with_expansion(
            info,
            initialcount,
            expanded_for_ambiguity
        )
        return ascending, descending, info
    end

    last
end

function alkane_ladder_scan_order_info_with_expansion(
    info,
    initialcount::Integer,
    expanded_for_ambiguity::Bool
)
    ntried = max(
        isnothing(info.trials.ascending) ? 0 : Int(info.trials.ascending.n_tried),
        isnothing(info.trials.descending) ? 0 : Int(info.trials.descending.n_tried)
    )

    merge(
        info,
        (
            n_peaks_tried=ntried,
            initial_n_peaks_tried=Int(initialcount),
            expanded_for_ambiguity=expanded_for_ambiguity
        )
    )
end

function alkane_ladder_scan_order_trial(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidates::AbstractVector,
    order::Symbol,
    base_mzretentionkwargs;
    standard,
    scanwindow::Integer,
    apexionexcludemzvalues,
    apexionminrelativeintensity::Real,
    variancefloor::Real,
    logfloorfraction::Real,
    maxapexshiftfromguess::Real,
    shapeioncount::Integer,
    shapemzspacing::Integer,
    extremeioncount::Integer,
    minioncount::Integer
)
    trial_mzkwargs = merge(base_mzretentionkwargs, (order=order,))
    apex_variances = Float64[]
    shape_widths = Float64[]
    ion_counts = Int[]
    mz_spans = Float64[]
    failures = NamedTuple[]

    for candidateinfo in candidates
        candidate = candidateinfo.candidate
        step = candidateinfo.ladderstep
        try
            shape_selection = alkane_ladder_scan_order_shape_ion_selection(
                msm,
                standard,
                step;
                ioncount=shapeioncount,
                mzspacing=shapemzspacing
            )
            shape = alkane_ladder_ion_apex(
                msm,
                variances,
                candidate;
                scanwindow=scanwindow,
                standard=standard,
                ladderstep=step,
                apexionexcludemzvalues=apexionexcludemzvalues,
                apexionmzvalues=shape_selection.mzvalues,
                apexionminrelativeintensity=apexionminrelativeintensity,
                minioncount=length(shape_selection.mzvalues),
                mzretentionkwargs=trial_mzkwargs,
                variancefloor=variancefloor,
                logfloorfraction=logfloorfraction,
                maxapexshiftfromguess=maxapexshiftfromguess
            )
            if !alkane_ladder_scan_order_shape_is_usable(shape)
                push!(
                    failures,
                    (
                        ladderstep=step,
                        scanindex=candidateinfo.scanindex,
                        reason="shape fit failed: $(alkane_ladder_apex_failure_reason(shape))",
                        ion_count=0,
                        mz_span=NaN,
                        shape_width=alkane_ladder_apex_width(shape),
                        shape_mzvalues=shape_selection.mzvalues
                    )
                )
                continue
            end

            score = alkane_ladder_scan_order_step_apex_variance(
                msm,
                variances,
                candidate,
                shape,
                trial_mzkwargs;
                standard=standard,
                extremeioncount=extremeioncount,
                minioncount=minioncount,
                logfloorfraction=logfloorfraction,
                variancefloor=variancefloor,
                maxapexshiftfromguess=maxapexshiftfromguess
            )
            if score.success
                push!(apex_variances, score.apex_variance)
                push!(shape_widths, score.shape_width)
                push!(ion_counts, score.ion_count)
                push!(mz_spans, score.mz_span)
            else
                push!(
                    failures,
                    (
                        ladderstep=step,
                        scanindex=candidateinfo.scanindex,
                        reason=score.reason,
                        ion_count=score.ion_count,
                        mz_span=score.mz_span,
                        shape_width=score.shape_width
                    )
                )
            end
        catch err
            push!(
                failures,
                (
                    ladderstep=step,
                    scanindex=candidateinfo.scanindex,
                    reason=sprint(showerror, err),
                    ion_count=0,
                    mz_span=NaN,
                    shape_width=NaN
                )
            )
        end
    end

    positive_variances = [max(variance, eps(Float64)) for variance in apex_variances]
    score = isempty(positive_variances) ? Inf : median(log.(positive_variances))
    finite_spans = filter(isfinite, mz_spans)
    (
        order=order,
        score=score,
        score_kind=:fixed_shape_ion_apex_variance,
        median_apex_variance=isempty(apex_variances) ? Inf : median(apex_variances),
        apex_variances=apex_variances,
        median_width=isempty(shape_widths) ? Inf : median(shape_widths),
        shape_widths=shape_widths,
        widths=shape_widths,
        ion_counts=ion_counts,
        mz_spans=mz_spans,
        median_ion_count=isempty(ion_counts) ? 0.0 : median(Float64.(ion_counts)),
        median_mz_span=isempty(finite_spans) ? NaN : median(finite_spans),
        n_success=length(apex_variances),
        n_tried=length(candidates),
        failures=failures
    )
end

function alkane_ladder_scan_order_trial_candidate(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidateinfo,
    order::Symbol,
    base_mzretentionkwargs;
    standard,
    scanwindow::Integer,
    apexionexcludemzvalues,
    apexionminrelativeintensity::Real,
    variancefloor::Real,
    logfloorfraction::Real,
    maxapexshiftfromguess::Real,
    shapeioncount::Integer,
    shapemzspacing::Integer,
    extremeioncount::Integer,
    minioncount::Integer
)
    trial_mzkwargs = merge(base_mzretentionkwargs, (order=order,))
    candidate = candidateinfo.candidate
    step = candidateinfo.ladderstep
    try
        shape_selection = alkane_ladder_scan_order_shape_ion_selection(
            msm,
            standard,
            step;
            ioncount=shapeioncount,
            mzspacing=shapemzspacing
        )
        shape = alkane_ladder_ion_apex(
            msm,
            variances,
            candidate;
            scanwindow=scanwindow,
            standard=standard,
            ladderstep=step,
            apexionexcludemzvalues=apexionexcludemzvalues,
            apexionmzvalues=shape_selection.mzvalues,
            apexionminrelativeintensity=apexionminrelativeintensity,
            minioncount=length(shape_selection.mzvalues),
            mzretentionkwargs=trial_mzkwargs,
            variancefloor=variancefloor,
            logfloorfraction=logfloorfraction,
            maxapexshiftfromguess=maxapexshiftfromguess
        )
        if !alkane_ladder_scan_order_shape_is_usable(shape)
            return (
                success=false,
                failure=(
                    ladderstep=step,
                    scanindex=candidateinfo.scanindex,
                    reason="shape fit failed: $(alkane_ladder_apex_failure_reason(shape))",
                    ion_count=0,
                    mz_span=NaN,
                    shape_width=alkane_ladder_apex_width(shape),
                    shape_mzvalues=shape_selection.mzvalues
                )
            )
        end

        score = alkane_ladder_scan_order_step_apex_variance(
            msm,
            variances,
            candidate,
            shape,
            trial_mzkwargs;
            standard=standard,
            extremeioncount=extremeioncount,
            minioncount=minioncount,
            logfloorfraction=logfloorfraction,
            variancefloor=variancefloor,
            maxapexshiftfromguess=maxapexshiftfromguess
        )
        score.success && return (
            success=true,
            apex_variance=score.apex_variance,
            shape_width=score.shape_width,
            ion_count=score.ion_count,
            mz_span=score.mz_span
        )

        (
            success=false,
            failure=(
                ladderstep=step,
                scanindex=candidateinfo.scanindex,
                reason=score.reason,
                ion_count=score.ion_count,
                mz_span=score.mz_span,
                shape_width=score.shape_width
            )
        )
    catch err
        (
            success=false,
            failure=(
                ladderstep=step,
                scanindex=candidateinfo.scanindex,
                reason=sprint(showerror, err),
                ion_count=0,
                mz_span=NaN,
                shape_width=NaN
            )
        )
    end
end

function alkane_ladder_scan_order_trial_summary(
    order::Symbol,
    candidates::AbstractVector,
    results::AbstractVector
)
    apex_variances = [result.apex_variance for result in results if result.success]
    shape_widths = [result.shape_width for result in results if result.success]
    ion_counts = [result.ion_count for result in results if result.success]
    mz_spans = [result.mz_span for result in results if result.success]
    failures = [result.failure for result in results if !result.success]
    positive_variances = [max(variance, eps(Float64)) for variance in apex_variances]
    score = isempty(positive_variances) ? Inf : median(log.(positive_variances))
    finite_spans = filter(isfinite, mz_spans)

    (
        order=order,
        score=score,
        score_kind=:fixed_shape_ion_apex_variance,
        median_apex_variance=isempty(apex_variances) ? Inf : median(apex_variances),
        apex_variances=apex_variances,
        median_width=isempty(shape_widths) ? Inf : median(shape_widths),
        shape_widths=shape_widths,
        widths=shape_widths,
        ion_counts=ion_counts,
        mz_spans=mz_spans,
        median_ion_count=isempty(ion_counts) ? 0.0 : median(Float64.(ion_counts)),
        median_mz_span=isempty(finite_spans) ? NaN : median(finite_spans),
        n_success=length(apex_variances),
        n_tried=length(candidates),
        failures=failures
    )
end

function alkane_ladder_scan_order_step_apex_variance(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidate,
    shape,
    mzretentionkwargs;
    standard,
    extremeioncount::Integer,
    minioncount::Integer,
    logfloorfraction::Real,
    variancefloor::Real,
    maxapexshiftfromguess::Real
)
    step = alkane_ladder_candidate_step(candidate)
    ion_selection = alkane_ladder_scan_order_extreme_apex_mzvalues(
        msm,
        standard,
        step;
        extremeioncount=extremeioncount
    )
    required_ion_count = min(Int(minioncount), length(ion_selection.mzvalues))
    mzindices = alkane_ladder_apex_mz_indices(
        msm,
        ion_selection.mzvalues;
        minioncount=required_ion_count
    )
    fits = alkane_ladder_fixed_shape_ion_apices(
        msm,
        shape,
        mzindices,
        variances,
        mzretentionkwargs;
        inputscanindex=alkane_ladder_input_scan_index(candidate),
        logfloorfraction=logfloorfraction,
        variancefloor=variancefloor,
        maxapexshiftfromguess=maxapexshiftfromguess
    )
    successful = [fit for fit in fits if fit.success]
    ion_count = length(successful)
    mz_values = [fit.mz_value for fit in successful]
    mz_span = ion_count >= 2 ? maximum(mz_values) - minimum(mz_values) : NaN
    shape_width = alkane_ladder_apex_width(shape)
    if ion_count < required_ion_count
        return (
            success=false,
            reason="fewer than $(required_ion_count) fixed-shape ion apices were inferred",
            ion_count=ion_count,
            mz_span=mz_span,
            shape_width=shape_width,
            apex_variance=NaN
        )
    end

    apex_scan_indices = [fit.apex_scan_index for fit in successful]
    apex_variance = alkane_ladder_robust_variance(apex_scan_indices)
    isfinite(apex_variance) || return (
        success=false,
        reason="fixed-shape ion apex variance is not finite",
        ion_count=ion_count,
        mz_span=mz_span,
        shape_width=shape_width,
        apex_variance=apex_variance
    )

    (
        success=true,
        reason=nothing,
        ion_count=ion_count,
        mz_span=mz_span,
        shape_width=shape_width,
        apex_variance=apex_variance,
        apex_scan_indices=apex_scan_indices,
        ion_selection=ion_selection,
        ion_fits=fits
    )
end

function alkane_ladder_choose_scan_order(
    ascending,
    descending;
    minpeaks::Integer,
    minapexvarianceratio::Real
)
    if !isfinite(ascending.score) && !isfinite(descending.score)
        return (
            selected_order=:unknown,
            requested_order=:inferdirection,
            status=:failed,
            score_kind=:fixed_shape_ion_apex_variance,
            evidence_score=NaN,
            ascending_score=ascending.score,
            descending_score=descending.score,
            ascending_median_apex_variance=ascending.median_apex_variance,
            descending_median_apex_variance=descending.median_apex_variance,
            ascending_median_width=ascending.median_width,
            descending_median_width=descending.median_width,
            n_peaks_used=0,
            n_peaks_tried=max(ascending.n_tried, descending.n_tried),
            initial_n_peaks_tried=0,
            expanded_for_ambiguity=false,
            trials=(ascending=ascending, descending=descending)
        )
    end

    best = descending.score <= ascending.score ? descending : ascending
    other = best.order == :descending ? ascending : descending
    evidence_score = isfinite(other.score) && isfinite(best.score) ?
        exp(other.score - best.score) :
        Inf
    status = best.n_success < minpeaks ?
        :insufficient_peaks :
        evidence_score >= minapexvarianceratio ?
            :accepted :
            :ambiguous

    (
        selected_order=best.order,
        requested_order=:inferdirection,
        status=status,
        score_kind=:fixed_shape_ion_apex_variance,
        evidence_score=evidence_score,
        ascending_score=ascending.score,
        descending_score=descending.score,
        ascending_median_apex_variance=ascending.median_apex_variance,
        descending_median_apex_variance=descending.median_apex_variance,
        ascending_median_width=ascending.median_width,
        descending_median_width=descending.median_width,
        n_peaks_used=best.n_success,
        n_peaks_tried=max(ascending.n_tried, descending.n_tried),
        initial_n_peaks_tried=0,
        expanded_for_ambiguity=false,
        trials=(ascending=ascending, descending=descending)
    )
end

function alkane_ladder_fixed_shape_ion_apices(
    msm::MassScanMatrix,
    shape,
    mzindices,
    variances::AbstractMatrix{<:Real},
    mzretentionkwargs;
    inputscanindex::Integer,
    logfloorfraction::Real,
    variancefloor::Real,
    maxapexshiftfromguess::Real
)
    alkane_ladder_scan_order_shape_is_usable(shape) || throw(ArgumentError(
        "shape apex is not usable for fixed-shape ion apex inference"))
    validate_alkane_series_variances(msm, variances)
    isfinite(logfloorfraction) && logfloorfraction > 0 || throw(ArgumentError(
        "logfloorfraction must be finite and positive"))
    validate_alkane_abundance_variancefloor(variancefloor)
    isfinite(maxapexshiftfromguess) && maxapexshiftfromguess >= 0 ||
        throw(ArgumentError(
            "maxapexshiftfromguess must be finite and nonnegative"))

    raw_scan_retentions = Float64.(rawretentions(msm))
    scan_retentions = retentions(msm)
    Xraw = rawintensities(msm)
    fit_center_retention = alkane_ladder_fit_center_retention(shape)
    fit_center_scan_index = alkane_ladder_fractional_scan_index(
        raw_scan_retentions,
        fit_center_retention
    )
    nscans = length(shape.scan_indices)
    localmargin = max(3, fld(nscans, 2) + 2)
    gamma = Float64(shape.gamma)
    xscale = Float64(shape.x_scale)
    fits = NamedTuple[]

    for mzindex in mzindices
        selected_scans = alkane_ladder_scan_indices_for_mz(
            msm,
            scan_retentions,
            fit_center_retention,
            nscans,
            mzindex,
            mzretentionkwargs,
            fit_center_scan_index,
            localmargin
        )
        intensities = [
            max(Float64(Xraw[scanindex, mzindex]), 0.0)
            for scanindex in selected_scans
        ]
        local_variances = [
            Float64(variances[scanindex, mzindex])
            for scanindex in selected_scans
        ]
        observation_retentions = [
            alkane_ladder_observation_retention(
                msm,
                scan_retentions,
                scanindex,
                mzindex,
                mzretentionkwargs
            )
            for scanindex in selected_scans
        ]
        fit = alkane_ladder_fit_fixed_shape_ion_apex(
            raw_scan_retentions,
            fit_center_retention,
            xscale,
            gamma,
            observation_retentions,
            intensities,
            local_variances;
            inputscanindex=inputscanindex,
            logfloorfraction=logfloorfraction,
            variancefloor=variancefloor,
            maxapexshiftfromguess=maxapexshiftfromguess
        )
        push!(
            fits,
            merge(
                fit,
                (
                    mz_index=Int(mzindex),
                    mz_value=Float64(rawmzvalues(msm)[mzindex]),
                    scan_indices=selected_scans,
                    observation_retentions=observation_retentions,
                    intensities=intensities,
                    variances=local_variances
                )
            )
        )
    end

    fits
end

function alkane_ladder_fit_fixed_shape_ion_apex(
    raw_scan_retentions::AbstractVector{<:Real},
    fit_center_retention::Real,
    xscale::Real,
    gamma::Real,
    observation_retentions,
    intensities,
    variances;
    inputscanindex::Integer,
    logfloorfraction::Real,
    variancefloor::Real,
    maxapexshiftfromguess::Real
)
    gamma_value = Float64(gamma)
    xscale_value = Float64(xscale)
    gamma_value < 0 && isfinite(gamma_value) || return (
        success=false,
        reason="fixed shape is not concave",
        apex_scan_index=NaN,
        apex_retention=NaN,
        apex_x=NaN,
        beta=NaN,
        intercept=NaN,
        reduced_normalized_residual=NaN
    )
    xscale_value > 0 && isfinite(xscale_value) || throw(ArgumentError(
        "fixed shape xscale must be finite and positive"))

    y = Float64.(intensities)
    local_variances = Float64.(variances)
    observation_values = Float64.(observation_retentions)
    nobs = length(y)
    nobs == length(local_variances) == length(observation_values) ||
        throw(DimensionMismatch(
            "fixed-shape ion fit vectors must have equal length"))
    nobs >= 3 || throw(ArgumentError(
        "at least three observations are required for fixed-shape ion apex inference"))
    ymax = maximum(y)
    ymax > 0 || return (
        success=false,
        reason="ion has no positive local signal",
        apex_scan_index=NaN,
        apex_retention=NaN,
        apex_x=NaN,
        beta=NaN,
        intercept=NaN,
        reduced_normalized_residual=NaN
    )

    logfloor = Float64(logfloorfraction) * ymax
    adjusted = y .+ logfloor
    x = (observation_values .- Float64(fit_center_retention)) ./ xscale_value
    logy = log.(adjusted)
    target = logy .- gamma_value .* abs2.(x)
    weights = [
        alkane_ladder_log_weight(adjusted[index], local_variances[index], variancefloor)
        for index in eachindex(adjusted)
    ]
    any(>(0), weights) || return (
        success=false,
        reason="ion fixed-shape fit has no positive weights",
        apex_scan_index=NaN,
        apex_retention=NaN,
        apex_x=NaN,
        beta=NaN,
        intercept=NaN,
        reduced_normalized_residual=NaN
    )

    design = hcat(ones(Float64, nobs), x)
    sqrtweights = sqrt.(weights)
    coefficients = (design .* reshape(sqrtweights, :, 1)) \ (target .* sqrtweights)
    intercept = coefficients[1]
    beta = coefficients[2]
    fitted_logy = intercept .+ beta .* x .+ gamma_value .* abs2.(x)
    apexx = -beta / (2 * gamma_value)
    apex_retention = Float64(fit_center_retention) + apexx * xscale_value
    apex_scan_index = isfinite(apex_retention) ?
        alkane_ladder_fractional_scan_index(raw_scan_retentions, apex_retention) :
        NaN
    apex_in_window = isfinite(apexx) && minimum(x) <= apexx <= maximum(x)
    apex_shift_from_guess_scans = apex_scan_index - Int(inputscanindex)
    apex_within_allowed_shift = isfinite(apex_shift_from_guess_scans) &&
        abs(apex_shift_from_guess_scans) <= Float64(maxapexshiftfromguess) + 1e-9
    success = isfinite(apex_scan_index) &&
        apex_in_window &&
        apex_within_allowed_shift
    ssres = sum(weights .* abs2.(logy .- fitted_logy))
    dof = nobs - 2
    reduced_normalized_residual = dof > 0 ? ssres / dof : NaN
    reason = success ?
        nothing :
        !isfinite(apex_scan_index) ?
            "fixed-shape apex scan index is not finite" :
        !apex_in_window ?
            "fixed-shape apex is outside the local window" :
        !apex_within_allowed_shift ?
            "fixed-shape apex exceeds maxapexshiftfromguess from abundance-trace guess" :
            "fixed-shape ion apex fit failed"

    (
        success=success,
        reason=reason,
        apex_scan_index=apex_scan_index,
        apex_retention=apex_retention,
        apex_x=apexx,
        beta=beta,
        intercept=intercept,
        reduced_normalized_residual=reduced_normalized_residual,
        apex_shift_from_guess_scans=apex_shift_from_guess_scans,
        apex_within_allowed_shift=apex_within_allowed_shift,
        apex_in_window=apex_in_window
    )
end

function alkane_ladder_scan_order_shape_is_usable(shape)
    hasproperty(shape, :success) && shape.success || return false
    width = alkane_ladder_apex_width(shape)

    isfinite(width) &&
        width > 0 &&
        hasproperty(shape, :scan_indices) &&
        length(shape.scan_indices) >= 3 &&
        isfinite(alkane_ladder_fit_center_retention(shape))
end

function alkane_ladder_fit_center_retention(apex)
    hasproperty(apex, :fit_center_retention) &&
        return Float64(apex.fit_center_retention)

    Float64(apex.input_retention)
end

function alkane_ladder_robust_variance(values)
    finite_values = [Float64(value) for value in values if isfinite(value)]
    length(finite_values) >= 2 || return NaN
    center = median(finite_values)
    deviations = abs.(finite_values .- center)
    mad = median(deviations)
    if isfinite(mad) && mad > 0
        return abs2(1.4826 * mad)
    end

    arithmetic_mean = sum(finite_values) / length(finite_values)
    sum(abs2(value - arithmetic_mean) for value in finite_values) /
        (length(finite_values) - 1)
end

function alkane_ladder_apex_width(apex)
    hasproperty(apex, :gamma) && hasproperty(apex, :x_scale) || return NaN
    gamma = Float64(apex.gamma)
    xscale = Float64(apex.x_scale)
    gamma < 0 && isfinite(gamma) && isfinite(xscale) && xscale > 0 || return NaN

    xscale / sqrt(-gamma)
end

function alkane_ladder_apex_failure_reason(apex)
    apex.success && return "success"
    !isfinite(apex.apex_scan_index) && return "nonfinite apex scan index"
    !isfinite(apex.apex_retention) && return "nonfinite apex retention"
    hasproperty(apex, :apex_within_allowed_shift) &&
        !apex.apex_within_allowed_shift &&
        return "apex exceeds maxapexshiftfromguess from abundance-trace guess"
    apex.gamma >= 0 && return "log-quadratic fit is not concave"
    !apex.apex_in_window && return "continuous apex is outside scan window"

    "apex fit failed"
end

function alkane_ladder_discrete_apex(
    candidate,
    retentions::AbstractVector{Float64},
    abundance::AbstractVector{Float64};
    scanindices,
    reason,
    failurereason=String(reason)
)
    scanindex = alkane_ladder_input_scan_index(candidate)
    (
        success=false,
        reason=reason,
        failurereason=failurereason,
        refined=false,
        fallback=true,
        ladderstep=alkane_ladder_candidate_step(candidate),
        scanindex=scanindex,
        retention=retentions[scanindex],
        apex_scan_index=Float64(scanindex),
        apex_retention=retentions[scanindex],
        apex_offset_scans=0.0,
        apex_offset_retention=0.0,
        apexscanindex=Float64(scanindex),
        apexretention=retentions[scanindex],
        output_scan_index=Float64(scanindex),
        output_retention=retentions[scanindex],
        outputscanindex=Float64(scanindex),
        outputretention=retentions[scanindex],
        apexoffsetscans=0.0,
        apexoffsetretention=0.0,
        apexabundance=abundance[scanindex],
        scan_indices=collect(scanindices),
        scanindices=collect(scanindices),
        retentions=retentions[scanindices],
        abundance=abundance[scanindices],
        beta=NaN,
        gamma=NaN,
        apex_x=NaN,
        x_scale=NaN,
        xscale=NaN,
        r2=NaN,
        weighted_log_residual_sum_squares=NaN,
        weightedlogresidualsumsquares=NaN,
        fit_coefficients=Float64[],
        log_floor=NaN,
        logfloor=NaN,
        fit_degrees_of_freedom=0,
        reduced_normalized_residual=NaN,
        source=hasproperty(candidate, :source) ? candidate.source : :molecular_ion_dp,
        gapfilled=hasproperty(candidate, :gapfilled) ? Bool(candidate.gapfilled) : false,
        edgeextended=hasproperty(candidate, :edgeextended) ?
            Bool(candidate.edgeextended) :
            false,
        mass_spectrum_cosine=hasproperty(candidate, :massspectrumcosine) ?
            Float64(candidate.massspectrumcosine) :
            NaN,
        required_cosine=hasproperty(candidate, :requiredcosine) ?
            Float64(candidate.requiredcosine) :
            NaN,
        apex_fit_quality_score=NaN,
        apex_fit_quality_log_score=NaN,
        apex_fit_quality_zscore=NaN,
        apex_fit_quality_zscore_excluded=false,
        apex_fit_quality_outlier_z=NaN,
        apex_fit_quality_baseline_median=NaN,
        apex_fit_quality_baseline_mad=NaN,
        apex_fit_quality_nsteps=0,
        apex_fit_quality_applied=false,
        calibration_excluded=false,
        calibration_exclusion_reason=nothing,
        good_for_calibration=false,
        candidate=candidate
    )
end

function alkane_ladder_candidate_step(candidate)
    hasproperty(candidate, :ladderstep) && return Int(candidate.ladderstep)
    hasproperty(candidate, :step) && return Int(candidate.step)
    throw(ArgumentError("candidate must contain ladderstep"))
end

function alkane_ladder_input_scan_index(candidate)
    if hasproperty(candidate, :scanindex)
        return Int(candidate.scanindex)
    elseif hasproperty(candidate, :apexindex)
        return Int(candidate.apexindex)
    elseif hasproperty(candidate, :scan_index)
        return Int(candidate.scan_index)
    end
    throw(ArgumentError(
        "candidate must contain scanindex, apexindex, or scan_index"))
end

function alkane_ladder_candidate_abundance(abundanceinfo, step, retentions)
    haskey(abundanceinfo.abundances, step) || throw(ArgumentError(
        "abundanceinfo.abundances does not contain C$(step)"))
    abundance = alkane_abundance_values(abundanceinfo.abundances[step], step)
    length(abundance) == length(retentions) || throw(DimensionMismatch(
        "abundance vector for C$(step) must match retention length"))

    abundance
end

function alkane_ladder_fallback_scanindices(
    retentions::AbstractVector{Float64},
    scanindex::Integer,
    scanwindow::Integer,
    candidate
)
    leftindex = max(1, Int(scanindex) - Int(scanwindow))
    rightindex = min(length(retentions), Int(scanindex) + Int(scanwindow))
    if hasproperty(candidate, :window)
        leftindex = max(leftindex, Int(candidate.window.leftindex))
        rightindex = min(rightindex, Int(candidate.window.rightindex))
    end

    collect(leftindex:rightindex)
end

function alkane_ladder_reference_filtered_apex_mzvalues(
    msm::MassScanMatrix,
    standard,
    ladderstep;
    excludemzvalues,
    minrelativeintensity::Real,
    minioncount::Integer
)
    reference_key, spectrum = alkane_ladder_reference_spectrum(standard, ladderstep)
    ref_mzs = Float64.(mzvalues(spectrum))
    ref_intensities = Float64.(intensities(spectrum))
    selected = alkane_ladder_reference_threshold_apex_ion_candidates(
        msm,
        ref_mzs,
        ref_intensities;
        excludemzvalues=excludemzvalues,
        minrelativeintensity=minrelativeintensity
    )
    length(selected) >= minioncount || throw(ArgumentError(
        "fewer than $(minioncount) reference apex ions are available for C$(ladderstep)"))

    (
        selection=:reference_relative_intensity_threshold,
        excluded_mzvalues=alkane_ladder_mzvalue_vector(
            excludemzvalues;
            allowempty=true
        ),
        mzvalues=[candidate.mzvalue for candidate in selected],
        reference_mzvalues=ref_mzs,
        reference_relative_intensities=[
            candidate.relativeintensity for candidate in selected
        ],
        reference_key=reference_key
    )
end

function alkane_ladder_scan_order_extreme_apex_mzvalues(
    msm::MassScanMatrix,
    standard,
    ladderstep;
    extremeioncount::Integer
)
    extremeioncount >= 1 || throw(ArgumentError(
        "scan-order extremeioncount must be a positive integer"))

    reference_key, spectrum = alkane_ladder_reference_spectrum(standard, ladderstep)
    ref_mzs = Float64.(mzvalues(spectrum))
    ref_intensities = Float64.(intensities(spectrum))
    candidates = alkane_ladder_reference_all_present_ion_candidates(
        msm,
        ref_mzs,
        ref_intensities
    )
    selected, excluded_mz = alkane_ladder_scan_order_extreme_ion_candidates(
        candidates,
        ladderstep;
        extremeioncount=extremeioncount
    )
    length(selected) >= 2 || throw(ArgumentError(
        "fewer than two scan-order contrast ions are present on the m/z grid " *
        "for C$(ladderstep)"))

    (
        selection=:reference_alkane_series_extreme_grid_ions,
        extremeioncount=Int(extremeioncount),
        excluded_mzvalues=isfinite(excluded_mz) ? [excluded_mz] : Float64[],
        mzvalues=[candidate.mzvalue for candidate in selected],
        reference_mzvalues=ref_mzs,
        reference_relative_intensities=[
            candidate.relativeintensity for candidate in selected
        ],
        reference_key=reference_key
    )
end

function alkane_ladder_scan_order_extreme_ion_candidates(
    candidates::AbstractVector,
    ladderstep;
    extremeioncount::Integer
)
    step = Int(ladderstep)
    molecular_mz = alkane_ladder_nominal_molecular_ion_mz(step)
    excluded_mz = alkane_ladder_nominal_molecular_minus_14_series_mz(step)

    eligible = NamedTuple[]
    for candidate in candidates
        nominal_mz = round(Int, Float64(candidate.mzvalue))
        abs(Float64(candidate.mzvalue) - nominal_mz) <= 0.5 || continue
        is_fragment = alkane_ladder_is_alkane_fragment_series_mz(nominal_mz)
        is_molecular = nominal_mz == molecular_mz
        (is_fragment || is_molecular) || continue
        nominal_mz == excluded_mz && continue
        push!(eligible, candidate)
    end
    sort!(eligible; by=candidate -> candidate.mzvalue)

    selected = NamedTuple[]
    for pairindex in 1:Int(extremeioncount)
        leftindex = pairindex
        rightindex = length(eligible) - pairindex + 1
        leftindex < rightindex || break

        push!(selected, eligible[leftindex])
        push!(selected, eligible[rightindex])
    end
    sort!(selected; by=candidate -> candidate.mzvalue)

    selected, Float64(excluded_mz)
end

alkane_ladder_nominal_molecular_ion_mz(ladderstep::Integer) = 14 * Int(ladderstep) + 2

alkane_ladder_nominal_molecular_minus_14_series_mz(ladderstep::Integer) =
    1 + 14 * (Int(ladderstep) - 1)

alkane_ladder_is_alkane_fragment_series_mz(mz::Integer) =
    mz >= 1 && mod(mz - 1, 14) == 0

function alkane_ladder_explicit_apex_mzvalue_selection(
    msm::MassScanMatrix,
    mzvalues;
    minioncount::Integer
)
    selected_mzindices = alkane_ladder_apex_mz_indices(
        msm,
        mzvalues;
        minioncount=minioncount
    )
    (
        selection=:explicit_mzvalues,
        excluded_mzvalues=Float64[],
        mzvalues=collect(rawmzvalues(msm)[selected_mzindices]),
        reference_mzvalues=Float64[],
        reference_relative_intensities=Float64[],
        reference_key=nothing
    )
end

function alkane_ladder_reference_threshold_apex_ion_candidates(
    msm::MassScanMatrix,
    ref_mzs::AbstractVector{<:Real},
    ref_intensities::AbstractVector{<:Real};
    excludemzvalues,
    minrelativeintensity::Real
)
    length(ref_mzs) == length(ref_intensities) || throw(DimensionMismatch(
        "reference spectrum m/z and intensity vectors must have the same length"))
    isempty(ref_mzs) && throw(ArgumentError("reference spectrum must not be empty"))
    max_intensity = maximum(abs, Float64.(ref_intensities))
    max_intensity > 0 || throw(ArgumentError(
        "reference spectrum must contain at least one nonzero intensity"))

    grid = Float64.(alkane_mz_bins(msm))
    excluded = alkane_ladder_mzvalue_vector(excludemzvalues; allowempty=true)
    intensity_by_grid_index = Dict{Int, Float64}()
    for (mz, intensity) in zip(ref_mzs, ref_intensities)
        alkane_ladder_mz_is_excluded(mz, excluded) && continue
        grid_index = alkane_ladder_nearest_grid_mz_index(grid, mz)
        isnothing(grid_index) && continue
        intensity_by_grid_index[grid_index] =
            get(intensity_by_grid_index, grid_index, 0.0) +
            abs(Float64(intensity))
    end

    candidates = [
        (
            mzvalue=Float64(grid[grid_index]),
            mzindex=grid_index,
            relativeintensity=intensity / max_intensity
        )
        for (grid_index, intensity) in intensity_by_grid_index
        if intensity / max_intensity >= Float64(minrelativeintensity)
    ]
    sort!(candidates; by=candidate -> (-candidate.relativeintensity, candidate.mzvalue))

    candidates
end

function alkane_ladder_reference_all_present_ion_candidates(
    msm::MassScanMatrix,
    ref_mzs::AbstractVector{<:Real},
    ref_intensities::AbstractVector{<:Real}
)
    length(ref_mzs) == length(ref_intensities) || throw(DimensionMismatch(
        "reference spectrum m/z and intensity vectors must have the same length"))
    isempty(ref_mzs) && throw(ArgumentError("reference spectrum must not be empty"))
    max_intensity = maximum(abs, Float64.(ref_intensities))
    max_intensity > 0 || throw(ArgumentError(
        "reference spectrum must contain at least one nonzero intensity"))

    grid = Float64.(alkane_mz_bins(msm))
    intensity_by_grid_index = Dict{Int, Float64}()
    for (mz, intensity) in zip(ref_mzs, ref_intensities)
        intensity_value = abs(Float64(intensity))
        intensity_value > 0 || continue
        grid_index = alkane_ladder_nearest_grid_mz_index(grid, mz)
        isnothing(grid_index) && continue
        intensity_by_grid_index[grid_index] =
            get(intensity_by_grid_index, grid_index, 0.0) + intensity_value
    end

    candidates = [
        (
            mzvalue=Float64(grid[grid_index]),
            mzindex=grid_index,
            relativeintensity=intensity / max_intensity
        )
        for (grid_index, intensity) in intensity_by_grid_index
    ]
    sort!(candidates; by=candidate -> candidate.mzvalue)

    candidates
end

function alkane_ladder_scan_order_shape_ion_selection(
    msm::MassScanMatrix,
    standard,
    ladderstep;
    ioncount::Integer,
    mzspacing::Integer
)
    ioncount >= 2 || throw(ArgumentError("shape ioncount must be at least 2"))
    mzspacing >= 1 || throw(ArgumentError("shape mzspacing must be at least 1"))
    reference_key, spectrum = alkane_ladder_reference_spectrum(standard, ladderstep)
    ref_mzs = Float64.(mzvalues(spectrum))
    ref_intensities = Float64.(intensities(spectrum))
    candidates = alkane_ladder_reference_all_present_ion_candidates(
        msm,
        ref_mzs,
        ref_intensities
    )
    length(candidates) >= ioncount || throw(ArgumentError(
        "fewer than $(ioncount) nonzero reference ions are present on the m/z grid " *
        "for C$(ladderstep)"))
    best, score = alkane_ladder_best_spaced_shape_ion_candidates(
        candidates;
        ioncount=ioncount,
        mzspacing=mzspacing
    )

    (
        selection=:reference_best_spaced_shape_ions,
        mzspacing=Int(mzspacing),
        mzvalues=[candidate.mzvalue for candidate in best],
        reference_relative_intensities=[
            candidate.relativeintensity for candidate in best
        ],
        score=score,
        reference_key=reference_key,
        reference_mzvalues=ref_mzs
    )
end

function alkane_ladder_best_spaced_shape_ion_candidates(
    candidates;
    ioncount::Integer,
    mzspacing::Integer
)
    by_nominal_mz = Dict{Int, NamedTuple}()
    for candidate in candidates
        nominal_mz = round(Int, Float64(candidate.mzvalue))
        abs(Float64(candidate.mzvalue) - nominal_mz) <= 0.5 || continue
        previous = get(by_nominal_mz, nominal_mz, nothing)
        if isnothing(previous) ||
                candidate.relativeintensity > previous.relativeintensity
            by_nominal_mz[nominal_mz] = candidate
        end
    end

    best = NamedTuple[]
    best_score = -Inf
    for start_mz in sort(collect(keys(by_nominal_mz)))
        selected = NamedTuple[]
        for offset in 0:(ioncount - 1)
            candidate = get(
                by_nominal_mz,
                start_mz + offset * mzspacing,
                nothing
            )
            isnothing(candidate) && break
            push!(selected, candidate)
        end
        length(selected) == ioncount || continue
        score = sum(candidate.relativeintensity for candidate in selected)
        if score > best_score
            best = selected
            best_score = score
        end
    end

    isempty(best) && throw(ArgumentError(
        "no $(ioncount)-ion reference shape set with $(mzspacing) m/z spacing " *
        "is present on the measured m/z grid"))

    best, best_score
end

function alkane_ladder_reference_spectrum(standard, ladderstep)
    target = Int(ladderstep)
    spectra = alkane_standard_spectra(standard)
    for (index, spectrum) in pairs(spectra)
        key = if hasproperty(attrs(spectrum), :order)
            Int(attrs(spectrum).order)
        else
            Int(index)
        end
        key == target && return key, spectrum
    end
    throw(ArgumentError(
        "standard does not contain a reference spectrum for C$(ladderstep)"))
end

function alkane_ladder_mzvalue_vector(values; allowempty::Bool=false)
    targets = Float64.(collect(values))
    if isempty(targets)
        allowempty && return Float64[]
        throw(ArgumentError("m/z value collections must not be empty"))
    end
    all(isfinite, targets) || throw(ArgumentError(
        "m/z value collections must contain only finite values"))

    unique(targets)
end

function alkane_ladder_apex_mz_indices(
    msm::MassScanMatrix,
    mzvalues;
    minioncount::Integer
)
    targets = alkane_ladder_mzvalue_vector(mzvalues)
    grid = Float64.(alkane_mz_bins(msm))
    selected = Int[]
    for target in targets
        index = alkane_ladder_nearest_grid_mz_index(grid, target)
        if !isnothing(index) && !(index in selected)
            push!(selected, index)
        end
    end
    length(selected) >= minioncount || throw(ArgumentError(
        "fewer than $(minioncount) apex ion m/z values are present on the msm grid"))

    selected
end

function alkane_ladder_nearest_grid_mz_index(
    grid::AbstractVector{<:Real},
    target_mz::Real
)
    best_index = firstindex(grid)
    best_distance = Inf
    for index in eachindex(grid)
        distance = abs(Float64(grid[index]) - Float64(target_mz))
        if distance < best_distance
            best_distance = distance
            best_index = index
        end
    end

    best_distance <= 0.5 ? Int(best_index) : nothing
end

function alkane_ladder_mz_is_excluded(
    mz::Real,
    excluded_mzvalues::AbstractVector{<:Real}
)
    any(excluded -> abs(Float64(mz) - Float64(excluded)) <= 0.5, excluded_mzvalues)
end

function alkane_ladder_default_mzretention_kwargs(
    msm::MassScanMatrix;
    order::Symbol=:inferdirection
)
    order in (:inferdirection, :ascending, :descending, :simultaneous) ||
        throw(ArgumentError(
            "order must be :inferdirection, :ascending, :descending, or :simultaneous"))
    retention_values = retentions(msm)
    length(retention_values) >= 2 || throw(ArgumentError(
        "at least two scans are needed to infer mzretention scaninterval"))
    intervals = diff(retention_values)
    scaninterval = sum(intervals) / length(intervals)

    (
        retentionref=:start,
        scaninterval=scaninterval,
        mzcount=mzcount(msm),
        order=order,
        dwellref=:middle,
        dwell=:homogeneous
    )
end

function alkane_ladder_scan_indices_by_mz(
    msm::MassScanMatrix,
    scan_retentions,
    targetretention::Real,
    nscans::Integer,
    mzindices::AbstractVector{<:Integer},
    mzretentionkwargs::NamedTuple,
    centerscanindex::Real,
    localmargin::Integer
)
    selected = Matrix{Int}(undef, nscans, length(mzindices))
    for (mzcol, mzindex) in pairs(mzindices)
        selected[:, mzcol] .= alkane_ladder_scan_indices_for_mz(
            msm,
            scan_retentions,
            targetretention,
            nscans,
            mzindex,
            mzretentionkwargs,
            centerscanindex,
            localmargin
        )
    end

    selected
end

function alkane_ladder_scan_indices_for_mz(
    msm::MassScanMatrix,
    scan_retentions,
    targetretention::Real,
    nscans::Integer,
    mzindex::Integer,
    mzretentionkwargs::NamedTuple,
    centerscanindex::Real,
    localmargin::Integer
)
    scanrange = alkane_ladder_local_scan_search_range(
        scancount(msm),
        centerscanindex,
        nscans,
        localmargin
    )
    candidates = Vector{Tuple{Float64, Float64, Int}}(undef, length(scanrange))
    for (candidateindex, scanindex) in pairs(scanrange)
        observation_retention = alkane_ladder_observation_retention(
            msm,
            scan_retentions,
            scanindex,
            mzindex,
            mzretentionkwargs
        )
        candidates[candidateindex] = (
            abs(observation_retention - targetretention),
            observation_retention,
            scanindex
        )
    end
    sort!(candidates; by=candidate -> candidate[1])
    chosen = candidates[1:nscans]
    sort!(chosen; by=candidate -> candidate[2])

    [candidate[3] for candidate in chosen]
end

function alkane_ladder_local_scan_search_range(
    totalscans::Integer,
    centerscanindex::Real,
    nscans::Integer,
    localmargin::Integer
)
    totalscans >= 1 || throw(ArgumentError("totalscans must be positive"))
    nscans >= 1 || throw(ArgumentError("nscans must be positive"))
    localmargin >= 0 || throw(ArgumentError("localmargin must be nonnegative"))
    isfinite(centerscanindex) || throw(ArgumentError("centerscanindex must be finite"))

    center = clamp(round(Int, centerscanindex), 1, Int(totalscans))
    left = max(1, center - Int(localmargin))
    right = min(Int(totalscans), center + Int(localmargin))

    needed = min(Int(nscans), Int(totalscans))
    while right - left + 1 < needed
        if left > 1
            left -= 1
        elseif right < totalscans
            right += 1
        else
            break
        end
    end

    left:right
end

function alkane_ladder_observation_retention(
    msm::MassScanMatrix,
    scan_retentions,
    scanindex::Integer,
    mzindex::Integer,
    mzretentionkwargs::NamedTuple
)
    kwargs = merge(mzretentionkwargs, (mzindex=Int(mzindex),))
    retention = Core.kwcall(kwargs, mzretention, scan_retentions[scanindex])

    alkane_ladder_raw_retention_value(msm, retention)
end

function alkane_ladder_raw_retention_value(msm::MassScanMatrix, retention)
    if retention isa Real
        isfinite(retention) || throw(ArgumentError("retention must be finite"))
        return Float64(retention)
    end

    unit = retentionunit(msm)
    isnothing(unit) && throw(ArgumentError(
        "unitful retentions require a unitful msm retention axis"))
    try
        return Float64(ustrip(unit, retention))
    catch
        throw(ArgumentError("retention must be compatible with retentionunit(msm)"))
    end
end

function alkane_ladder_log_weight(
    adjusted_intensity::Real,
    variance::Real,
    variancefloor::Real
)
    v = max(Float64(variance), Float64(variancefloor))
    isfinite(v) && v > 0 || return 0.0

    abs2(Float64(adjusted_intensity)) / v
end

function alkane_ladder_max_center_passes(maxapexshiftfromguess::Real)
    max(1, ceil(Int, Float64(maxapexshiftfromguess)) + 1)
end

function alkane_ladder_best_apex_attempt_index(attempts::AbstractVector)
    isempty(attempts) && throw(ArgumentError("at least one apex attempt is required"))
    centered = findfirst(
        apex -> apex.success && alkane_ladder_apex_is_centered(apex),
        attempts
    )
    !isnothing(centered) && return centered
    successful = findlast(apex -> apex.success, attempts)
    !isnothing(successful) && return successful
    finite = findlast(
        apex -> isfinite(apex.apex_scan_index) && isfinite(apex.apex_retention),
        attempts
    )
    !isnothing(finite) && return finite

    length(attempts)
end

function alkane_ladder_has_centered_success(attempts::AbstractVector)
    any(apex -> apex.success && alkane_ladder_apex_is_centered(apex), attempts)
end

function alkane_ladder_center_already_attempted(
    attempts::AbstractVector,
    center_scan_index::Real;
    tolerance::Real=0.25
)
    any(
        apex -> hasproperty(apex, :fit_center_scan_index) &&
            isfinite(apex.fit_center_scan_index) &&
            abs(Float64(apex.fit_center_scan_index) - Float64(center_scan_index)) <=
                Float64(tolerance),
        attempts
    )
end

function alkane_ladder_apex_is_centered(apex; tolerance::Real=0.25)
    isfinite(tolerance) && tolerance >= 0 || throw(ArgumentError(
        "tolerance must be finite and nonnegative"))
    offset = if hasproperty(apex, :apex_offset_from_fit_center_scans)
        Float64(apex.apex_offset_from_fit_center_scans)
    elseif hasproperty(apex, :apex_scan_index) && hasproperty(apex, :fit_center_scan_index)
        Float64(apex.apex_scan_index) - Float64(apex.fit_center_scan_index)
    else
        NaN
    end

    isfinite(offset) && abs(offset) <= Float64(tolerance)
end

function alkane_ladder_next_center_retention(
    apex,
    raw_scan_retentions::AbstractVector{<:Real},
    input_scanindex::Integer,
    maxapexshiftfromguess::Real
)
    hasproperty(apex, :apex_scan_index) || return nothing
    isfinite(apex.apex_scan_index) || return nothing
    lower_scan = max(1.0, Float64(input_scanindex) - Float64(maxapexshiftfromguess))
    upper_scan = min(
        Float64(length(raw_scan_retentions)),
        Float64(input_scanindex) + Float64(maxapexshiftfromguess)
    )
    target_scan = clamp(Float64(apex.apex_scan_index), lower_scan, upper_scan)
    hasproperty(apex, :fit_center_scan_index) &&
        abs(target_scan - Float64(apex.fit_center_scan_index)) < 0.25 &&
        return nothing

    alkane_ladder_retention_at_fractional_scan_index(raw_scan_retentions, target_scan)
end

function alkane_ladder_center_search_retentions(
    raw_scan_retentions::AbstractVector{<:Real},
    input_scanindex::Integer,
    maxapexshiftfromguess::Real
)
    maxshift = Float64(maxapexshiftfromguess)
    maxshift == 0 && return Float64[]
    offsets = Float64[]
    whole_scan_shift = floor(Int, maxshift)
    for shift in 1:whole_scan_shift
        push!(offsets, Float64(shift))
        push!(offsets, -Float64(shift))
    end
    if maxshift - whole_scan_shift > 1e-9
        push!(offsets, maxshift)
        push!(offsets, -maxshift)
    end

    center_scan_indices = Float64[]
    for offset in offsets
        center_scan = clamp(
            Float64(input_scanindex) + offset,
            1.0,
            Float64(length(raw_scan_retentions))
        )
        any(existing -> abs(existing - center_scan) <= 0.25, center_scan_indices) &&
            continue
        push!(center_scan_indices, center_scan)
    end

    [
        alkane_ladder_retention_at_fractional_scan_index(
            raw_scan_retentions,
            center_scan
        )
        for center_scan in center_scan_indices
    ]
end

function alkane_validate_retention_axis(retentions::AbstractVector{Float64})
    isempty(retentions) && throw(ArgumentError("retention axis must not be empty"))
    all(isfinite, retentions) || throw(ArgumentError("retention axis must be finite"))
    for index in 2:length(retentions)
        retentions[index] > retentions[index - 1] || throw(ArgumentError(
            "retention axis must be strictly increasing"))
    end

    nothing
end

function alkane_ladder_fractional_scan_index(
    retentions::AbstractVector{<:Real},
    retention::Real
)
    all(isfinite, retentions) || throw(ArgumentError("scan retentions must be finite"))
    isfinite(retention) || return NaN
    target = Float64(retention)
    length(retentions) == 1 && return 1.0
    increasing = last(retentions) >= first(retentions)
    increasing && return alkane_ladder_fractional_scan_index_increasing(
        retentions,
        target
    )
    reversed_index = alkane_ladder_fractional_scan_index_increasing(
        reverse(retentions),
        target
    )

    length(retentions) - reversed_index + 1
end

function alkane_ladder_fractional_scan_index_increasing(
    retentions::AbstractVector{<:Real},
    retention::Real
)
    retention <= first(retentions) && return 1.0
    retention >= last(retentions) && return Float64(length(retentions))
    right = searchsortedfirst(retentions, retention)
    left = right - 1
    left_retention = Float64(retentions[left])
    right_retention = Float64(retentions[right])
    right_retention == left_retention && return Float64(left)

    left + (Float64(retention) - left_retention) /
        (right_retention - left_retention)
end

alkane_fractional_scan_index(
    retentions::AbstractVector{Float64},
    retention::Real
) = alkane_ladder_fractional_scan_index(retentions, retention)

function alkane_ladder_retention_at_fractional_scan_index(
    retentions::AbstractVector{<:Real},
    scanindex::Real
)
    isempty(retentions) && throw(ArgumentError("scan retentions must not be empty"))
    all(isfinite, retentions) || throw(ArgumentError("scan retentions must be finite"))
    clamped_scanindex = clamp(Float64(scanindex), 1.0, Float64(length(retentions)))
    left = max(1, floor(Int, clamped_scanindex))
    right = min(length(retentions), ceil(Int, clamped_scanindex))
    left == right && return Float64(retentions[left])
    fraction = clamped_scanindex - left

    Float64(retentions[left]) +
        fraction * (Float64(retentions[right]) - Float64(retentions[left]))
end

function alkane_ladder_scan_indices_around_retention(
    retentions::AbstractVector{<:Real},
    targetretention::Real,
    nscans::Integer
)
    nscans >= 1 || throw(ArgumentError("nscans must be at least 1"))
    isempty(retentions) && throw(ArgumentError("scan retentions must not be empty"))
    all(isfinite, retentions) || throw(ArgumentError("scan retentions must be finite"))
    isfinite(targetretention) || throw(ArgumentError(
        "targetretention must be finite"))
    requested = min(Int(nscans), length(retentions))
    candidates = [
        (
            distance=abs(Float64(retentions[index]) - Float64(targetretention)),
            scanindex=Int(index)
        )
        for index in eachindex(retentions)
    ]
    sort!(candidates; by=candidate -> (candidate.distance, candidate.scanindex))
    selected = [candidate.scanindex for candidate in candidates[1:requested]]
    sort!(selected)

    selected
end
