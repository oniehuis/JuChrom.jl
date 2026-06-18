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

struct AlkaneLadderScanOrderTrials{T1, T2}
    ascending::T1
    descending::T2
end

struct AlkaneLadderScanOrderCandidate{T}
    ladderstep::Int
    scanindex::Int
    candidate::T
end

struct AlkaneLadderScanOrderTrialFailure
    ladderstep::Int
    scanindex::Int
    reason::String
    ion_count::Int
    mz_span::Float64
    shape_width::Float64
    shape_mzvalues::Vector{Float64}
end

struct AlkaneLadderScanOrderTrialResult{T}
    success::Bool
    apex_variance::Float64
    shape_width::Float64
    ion_count::Int
    mz_span::Float64
    failure::Union{Nothing,T}
end

struct AlkaneLadderScanOrderStepApexVariance{T1, T2}
    success::Bool
    reason::Union{Nothing,String}
    ion_count::Int
    mz_span::Float64
    shape_width::Float64
    apex_variance::Float64
    apex_scan_indices::Vector{Float64}
    ion_selection::T1
    ion_fits::T2
end

struct AlkaneLadderApexIonCandidate
    mzvalue::Float64
    mzindex::Int
    relativeintensity::Float64
end

struct AlkaneLadderApexIonSelection{
    T1,
    T2<:Real,
    T3<:Real,
    T4<:Real,
    T5<:Real
}
    selection::Symbol
    excluded_mzvalues::Vector{T2}
    mzvalues::Vector{T3}
    reference_mzvalues::Vector{T4}
    reference_relative_intensities::Vector{T5}
    reference_key::T1
    mzspacing::Union{Nothing, Int}
    score::Float64
    extremeioncount::Union{Nothing, Int}
end

struct AlkaneLadderFixedShapeIonFit
    success::Bool
    reason::Union{Nothing, String}
    apex_scan_index::Float64
    apex_retention::Float64
    apex_x::Float64
    beta::Float64
    intercept::Float64
    reduced_normalized_residual::Float64
    apex_shift_from_guess_scans::Float64
    apex_within_allowed_shift::Bool
    apex_in_window::Bool
end

struct AlkaneLadderFixedShapeIonApexFit
    success::Bool
    reason::Union{Nothing, String}
    apex_scan_index::Float64
    apex_retention::Float64
    apex_x::Float64
    beta::Float64
    intercept::Float64
    reduced_normalized_residual::Float64
    apex_shift_from_guess_scans::Float64
    apex_within_allowed_shift::Bool
    apex_in_window::Bool
    mz_index::Int
    mz_value::Float64
    scan_indices::Vector{Int}
    observation_retentions::Vector{Float64}
    intensities::Vector{Float64}
    variances::Vector{Float64}
end

struct AlkaneLadderScanOrderTrialSummary{T<:AbstractVector}
    order::Symbol
    score::Float64
    score_kind::Symbol
    median_apex_variance::Float64
    apex_variances::Vector{Float64}
    median_width::Float64
    shape_widths::Vector{Float64}
    widths::Vector{Float64}
    ion_counts::Vector{Int}
    mz_spans::Vector{Float64}
    median_ion_count::Float64
    median_mz_span::Float64
    n_success::Int
    n_tried::Int
    failures::T
end

struct AlkaneLadderScanOrderInfo{T<:AlkaneLadderScanOrderTrials}
    selected_order::Symbol
    requested_order::Symbol
    status::Symbol
    score_kind::Symbol
    evidence_score::Union{Nothing, Float64}
    ascending_score::Union{Nothing, Float64}
    descending_score::Union{Nothing, Float64}
    ascending_median_apex_variance::Union{Nothing, Float64}
    descending_median_apex_variance::Union{Nothing, Float64}
    ascending_median_width::Union{Nothing, Float64}
    descending_median_width::Union{Nothing, Float64}
    n_peaks_used::Int
    n_peaks_tried::Int
    initial_n_peaks_tried::Int
    expanded_for_ambiguity::Bool
    expanded_ions::Bool
    trials::T
end

struct AlkaneLadderScanOrderResult{
    T1<:NamedTuple,
    T2<:AlkaneLadderScanOrderInfo
}
    mzretentionkwargs::T1
    info::T2
end

struct AlkaneLadderApexSettings{
    T1<:AlkaneStandard,
    T2<:Integer,
    T3<:Real,
    T4<:Real,
    T5<:Union{Nothing, NamedTuple},
    T6<:Symbol,
    T7<:Union{AbstractVector, Tuple},
    T8<:Union{Nothing, AbstractVector, Tuple},
    T9<:Real,
    T10<:Integer,
    T11<:Real,
    T12<:Real,
    T13<:Union{Nothing, Real},
    T14<:Integer,
    T15<:Union{Nothing, Integer},
    T16<:Integer,
    T17<:Real,
    T18<:Integer,
    T19<:Integer,
    T20<:Integer,
    T21<:Integer
}
    standard::T1
    scanwindow::T2
    variancefloor::T3
    logfloorfraction::T4
    mzretentionkwargs::T5
    mzscanorder::T6
    apexionexcludemzvalues::T7
    apexionmzvalues::T8
    apexionminrelativeintensity::T9
    minioncount::T10
    maxapexshiftfromguess::T11
    apexcenteredscantolerance::T12
    apexfitqualityoutlierz::T13
    apexfitqualityminsteps::T14
    mzscanordermaxpeaks::T15
    mzscanorderminpeaks::T16
    mzscanorderminapexvarianceratio::T17
    mzscanordershapeioncount::T18
    mzscanordershapemzspacing::T19
    mzscanorderextremeioncount::T20
    mzscanorderminioncount::T21
end

struct AlkaneLadderApexFitQualityStats
    zscores::Vector{Float64}
    median::Float64
    mad::Float64
    nsteps::Int
    applied::Bool
end

struct AlkaneLadderIonApexAttempt{
    T1<:Union{Nothing, Integer},
    T2<:NamedTuple,
    T3<:Real,
    T4<:Real,
    T5<:Real,
    T6<:Real,
    T7<:Real,
    T8<:Real,
    T9<:Real,
    T10<:Real
}
    success::Bool
    apex_retention::Float64
    apex_time::Float64
    apex_scan_index::Float64
    input_scan_index::Int
    input_retention::T6
    apex_offset_retention::Float64
    apex_offset_scans::Float64
    fit_center_retention::T7
    fit_center_scan_index::T8
    apex_offset_from_fit_center_retention::Float64
    apex_offset_from_fit_center_scans::Float64
    apex_shift_from_guess_scans::Float64
    apex_within_allowed_shift::Bool
    apex_in_window::Bool
    attemptindex::Int
    beta::Float64
    gamma::Float64
    apex_x::Float64
    x_scale::Float64
    r2::Float64
    weighted_log_residual_sum_squares::Float64
    fit_degrees_of_freedom::Int
    reduced_normalized_residual::Float64
    normalized_log_residuals::Matrix{Float64}
    fitted_log_intensities::Matrix{Float64}
    scanwindow::Int
    scan_indices::Vector{Int}
    retentions::Vector{T9}
    peak_model::Vector{Float64}
    observation_peak_model::Matrix{Float64}
    scan_indices_by_mz::Matrix{Int}
    raw_retentions_by_mz::Matrix{Float64}
    observation_retentions::Matrix{Float64}
    fit_intensities::Matrix{Float64}
    fit_variances::Matrix{Float64}
    mz_indices::Vector{Int}
    mz_values::Vector{T3}
    apex_ion_selection::Symbol
    apex_ion_excluded_mzvalues::Vector{T10}
    apex_ion_reference_mzvalues::Vector{T4}
    apex_ion_reference_relative_intensities::Vector{T5}
    apex_ion_min_relative_intensity::Float64
    apex_ion_reference_key::T1
    ion_intercepts::Vector{Float64}
    log_floor::Float64
    logfloorfraction::Float64
    variance_weighted::Bool
    variancefloor::Float64
    mzretentionkwargs::T2
end

struct AlkaneLadderIonApexResult{
    T1<:AlkaneLadderIonApexAttempt,
    T2<:AlkaneLadderIonApexAttempt,
    T3<:AbstractVector{<:AlkaneLadderIonApexAttempt}
}
    fit::T1
    initial_apex::T2
    all_apex_attempts::T3
    recenter_attempts::Int
    recenter_used::Bool
    maxapexshiftfromguess::Float64
end

struct AlkaneLadderApex{
    T1<:Union{AbstractAlkaneLadderCandidate, NamedTuple},
    T2<:Union{Nothing, AlkaneLadderIonApexResult},
    T3<:Real
} <: AbstractAlkaneLadderApex
    success::Bool
    reason::Symbol
    failurereason::Union{Nothing, String}
    refined::Bool
    fallback::Bool
    ladderstep::Int
    scanindex::Int
    retention::T3
    apexscanindex::Float64
    apexretention::Float64
    output_scan_index::Float64
    output_retention::Float64
    outputscanindex::Float64
    outputretention::Float64
    apexoffsetscans::Float64
    apexoffsetretention::Float64
    apexabundance::Float64
    scanindices::Vector{Int}
    abundance::Vector{Float64}
    xscale::Float64
    weightedlogresidualsumsquares::Float64
    fit_coefficients::Vector{Float64}
    logfloor::Float64
    source::Symbol
    gapfilled::Bool
    edgeextended::Bool
    mass_spectrum_cosine::Float64
    required_cosine::Float64
    apex_fit_quality_score::Float64
    apex_fit_quality_log_score::Float64
    apex_fit_quality_zscore::Float64
    apex_fit_quality_zscore_excluded::Bool
    apex_fit_quality_outlier_z::Float64
    apex_fit_quality_baseline_median::Float64
    apex_fit_quality_baseline_mad::Float64
    apex_fit_quality_nsteps::Int
    apex_fit_quality_applied::Bool
    calibration_excluded::Bool
    calibration_exclusion_reason::Union{Nothing, String}
    good_for_calibration::Bool
    candidate::T1
    fit::T2
end

struct AlkaneLadderApexInfo{
    T1<:AbstractVector{<:AlkaneLadderApex},
    T2<:AbstractDict{Int, <:AlkaneLadderApex},
    T3<:AlkaneLadderApexSettings,
    T4<:AlkaneLadderScanOrderInfo
} <: AbstractAlkaneLadderApexInfo
    status::Symbol
    reason::Symbol
    apexes::T1
    bycarbon::T2
    settings::T3
    scanorderinfo::T4
    calibrationexcluded::Vector{Bool}
    goodforcalibration::Vector{Bool}
    apexfitqualityscores::Vector{Float64}
    apexfitqualityzscores::Vector{Float64}
end

"""
    alkaneladderapexes(msm, variances, abundanceinfo, pathinfo, settings)

Refine the apex of each selected alkane ladder path candidate from fitted ion traces.

`settings` is an `AlkaneLadderApexSettings` object built by the ladder finder from the
user-facing apex keywords. The returned `AlkaneLadderApexInfo` stores the resolved
m/z-retention timing settings used for apex fitting.
"""
function alkaneladderapexes(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    pathinfo::AlkaneLadderPathInfo,
    settings::AlkaneLadderApexSettings
)
    validate_alkane_ladder_apex_settings(
        msm,
        variances,
        settings.scanwindow,
        settings.minioncount,
        settings.apexionminrelativeintensity,
        settings.variancefloor,
        settings.logfloorfraction,
        settings.maxapexshiftfromguess,
        settings.apexcenteredscantolerance,
        settings.apexfitqualityoutlierz,
        settings.apexfitqualityminsteps
    )

    if isempty(pathinfo.path)
        base_mzkwargs = alkane_ladder_base_mzretention_kwargs(
            msm,
            settings.mzretentionkwargs,
            settings.mzscanorder
        )
        return alkane_ladder_empty_apexinfo(
            :no_path,
            base_mzkwargs,
            alkane_ladder_no_scan_order_info(:no_path),
            settings
        )
    end

    scanorder = alkaneladderscanorder(msm, variances, pathinfo, settings)
    resolved_mzretentionkwargs = scanorder.mzretentionkwargs
    resolved_settings = alkane_ladder_apex_settings_with_mzretentionkwargs(
        settings,
        resolved_mzretentionkwargs
    )
    if scanorder.info.selected_order ≡ :unknown
        return alkane_ladder_empty_apexinfo(
            :scan_order_inference_failed,
            resolved_mzretentionkwargs,
            scanorder.info,
            resolved_settings
        )
    end

    candidates = pathinfo.path
    apexes = Vector{AlkaneLadderApex}(undef, length(candidates))
    Base.Threads.@threads for index in eachindex(candidates)
        apexes[index] = alkaneladderapex(
            msm,
            variances,
            abundanceinfo,
            candidates[index],
            resolved_settings
        )
    end

    apexes = alkane_ladder_annotate_apex_fit_quality(
        apexes,
        settings.apexfitqualityoutlierz,
        settings.apexfitqualityminsteps
    )
    bycarbon = Dict(apex.ladderstep => apex for apex in apexes)
    AlkaneLadderApexInfo(
        isempty(apexes) ? :failed : :success,
        isempty(apexes) ? :no_path : :success,
        apexes,
        bycarbon,
        resolved_settings,
        scanorder.info,
        [apex.calibration_excluded for apex in apexes],
        [apex.good_for_calibration for apex in apexes],
        [apex.apex_fit_quality_score for apex in apexes],
        [apex.apex_fit_quality_zscore for apex in apexes]
    )
end

function alkane_ladder_empty_apexinfo(
    reason::Symbol,
    mzretentionkwargs,
    scanorderinfo,
    settings::AlkaneLadderApexSettings
)
    AlkaneLadderApexInfo(
        :failed,
        reason,
        AlkaneLadderApex[],
        Dict{Int,AlkaneLadderApex}(),
        alkane_ladder_apex_settings_with_mzretentionkwargs(settings, mzretentionkwargs),
        scanorderinfo,
        Bool[],
        Bool[],
        Float64[],
        Float64[]
    )
end

function alkane_ladder_apex_settings_with_mzretentionkwargs(
    settings::AlkaneLadderApexSettings,
    mzretentionkwargs::NamedTuple
)
    AlkaneLadderApexSettings(
        settings.standard,
        settings.scanwindow,
        settings.variancefloor,
        settings.logfloorfraction,
        mzretentionkwargs,
        settings.mzscanorder,
        settings.apexionexcludemzvalues,
        settings.apexionmzvalues,
        settings.apexionminrelativeintensity,
        settings.minioncount,
        settings.maxapexshiftfromguess,
        settings.apexcenteredscantolerance,
        settings.apexfitqualityoutlierz,
        settings.apexfitqualityminsteps,
        settings.mzscanordermaxpeaks,
        settings.mzscanorderminpeaks,
        settings.mzscanorderminapexvarianceratio,
        settings.mzscanordershapeioncount,
        settings.mzscanordershapemzspacing,
        settings.mzscanorderextremeioncount,
        settings.mzscanorderminioncount
    )
end

function alkane_ladder_apex_settings_with_ion_selection(
    settings::AlkaneLadderApexSettings,
    mzretentionkwargs::NamedTuple,
    apexionmzvalues::Union{Nothing, AbstractVector{<:Real}, Tuple{Vararg{Real}}},
    minioncount::Integer
)
    AlkaneLadderApexSettings(
        settings.standard,
        settings.scanwindow,
        settings.variancefloor,
        settings.logfloorfraction,
        mzretentionkwargs,
        get(mzretentionkwargs, :order, settings.mzscanorder),
        settings.apexionexcludemzvalues,
        apexionmzvalues,
        settings.apexionminrelativeintensity,
        minioncount,
        settings.maxapexshiftfromguess,
        settings.apexcenteredscantolerance,
        settings.apexfitqualityoutlierz,
        settings.apexfitqualityminsteps,
        settings.mzscanordermaxpeaks,
        settings.mzscanorderminpeaks,
        settings.mzscanorderminapexvarianceratio,
        settings.mzscanordershapeioncount,
        settings.mzscanordershapemzspacing,
        settings.mzscanorderextremeioncount,
        settings.mzscanorderminioncount
    )
end

function alkaneladderapexes(
    msm::MassScanMatrix,
    abundanceinfo::AlkaneAbundanceInfo,
    pathinfo::AlkaneLadderPathInfo,
    settings::AlkaneLadderApexSettings
)
    alkaneladderapexes(
        msm,
        ones(size(rawintensities(msm))),
        abundanceinfo,
        pathinfo,
        settings
    )
end

function validate_alkane_ladder_apex_settings(
    scanwindow,
    variancefloor,
    logfloorfraction
)
    scanwindow isa Integer || throw(ArgumentError("scanwindow must be an integer"))
    scanwindow ≥ 1 || throw(ArgumentError("scanwindow must be at least 1"))
    validate_alkane_abundance_variancefloor(variancefloor)
    isfinite(logfloorfraction) && logfloorfraction > 0 || throw(ArgumentError(
        "logfloorfraction must be finite and positive"))

    nothing
end

function validate_alkane_ladder_apex_settings(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    scanwindow::Integer,
    minioncount::Integer,
    apexionminrelativeintensity::Real,
    variancefloor::Real,
    logfloorfraction::Real,
    maxapexshiftfromguess::Real,
    apexcenteredscantolerance::Real,
    apexfitqualityoutlierz::Union{Nothing, Real},
    apexfitqualityminsteps::Integer
)
    validate_alkane_series_variances(msm, variances)
    validate_alkane_ladder_apex_settings(scanwindow, variancefloor, logfloorfraction)
    minioncount ≥ 1 || throw(ArgumentError("minioncount must be at least 1"))
    isfinite(apexionminrelativeintensity) && 0 ≤ apexionminrelativeintensity < 1 ||
        throw(ArgumentError("apexionminrelativeintensity must be finite and in [0, 1)"))
    isfinite(maxapexshiftfromguess) && maxapexshiftfromguess ≥ 0 || 
        throw(ArgumentError("maxapexshiftfromguess must be finite and nonnegative"))
    isfinite(apexcenteredscantolerance) && apexcenteredscantolerance ≥ 0 ||
        throw(ArgumentError(
            "apexcenteredscantolerance must be finite and nonnegative"))
    if !isnothing(apexfitqualityoutlierz)
        isfinite(apexfitqualityoutlierz) && apexfitqualityoutlierz > 0 ||
            throw(ArgumentError(
                "apexfitqualityoutlierz must be finite and positive, or nothing"))
    end
    apexfitqualityminsteps ≥ 2 || throw(ArgumentError(
        "apexfitqualityminsteps must be at least 2"))

    nothing
end

function alkaneladderapex(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    candidate::Union{AbstractAlkaneLadderCandidate, NamedTuple},
    settings::AlkaneLadderApexSettings
)
    validate_alkane_ladder_apex_settings(
        msm,
        variances,
        settings.scanwindow,
        settings.minioncount,
        settings.apexionminrelativeintensity,
        settings.variancefloor,
        settings.logfloorfraction,
        settings.maxapexshiftfromguess,
        settings.apexcenteredscantolerance,
        settings.apexfitqualityoutlierz,
        settings.apexfitqualityminsteps
    )

    retentions = rawretentions(msm)
    alkane_validate_retention_axis(retentions)
    step = alkane_ladder_candidate_step(candidate)
    abundance = alkane_ladder_candidate_abundance(abundanceinfo, step, retentions)
    inputscan = alkane_ladder_input_scan_index(candidate)

    try
        apex = alkane_ladder_ion_apex(msm, variances, candidate, step, settings)
        return alkane_ladder_apex_public_result(apex, candidate, abundance, retentions)
    catch err
        (err isa ArgumentError || err isa DimensionMismatch) || rethrow()
        return alkane_ladder_discrete_apex(
            candidate,
            retentions,
            abundance,
            alkane_ladder_fallback_scanindices(
                retentions,
                inputscan,
                settings.scanwindow,
                candidate
            ),
            :apex_fit_failed,
            sprint(showerror, err)
        )
    end
end

function alkaneladderapex(
    msm::MassScanMatrix,
    abundanceinfo::AlkaneAbundanceInfo,
    candidate::Union{AbstractAlkaneLadderCandidate, NamedTuple},
    settings::AlkaneLadderApexSettings
)
    alkaneladderapex(
        msm,
        ones(size(rawintensities(msm))),
        abundanceinfo,
        candidate,
        settings
    )
end

function alkane_ladder_ion_apex(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidate::Union{AbstractAlkaneLadderCandidate, NamedTuple},
    ladderstep::Integer,
    settings::AlkaneLadderApexSettings
)
    input_scanindex = alkane_ladder_input_scan_index(candidate)
    raw_scan_retentions = rawretentions(msm)
    1 ≤ input_scanindex ≤ length(raw_scan_retentions) || throw(ArgumentError(
        "candidate scan index $(input_scanindex) is outside 1:$(length(raw_scan_retentions))"))

    attempts = AlkaneLadderIonApexAttempt[]
    function fit_at_center!(center_retention)
        apex = alkane_ladder_ion_apex_once(
            msm,
            variances,
            candidate,
            ladderstep,
            center_retention,
            length(attempts) + 1,
            settings
        )
        push!(attempts, apex)

        apex
    end

    center_retention = raw_scan_retentions[input_scanindex]
    for _ in 1:alkane_ladder_max_center_passes(settings.maxapexshiftfromguess)
        apex = fit_at_center!(center_retention)
        apex.success &&
            alkane_ladder_apex_is_centered(
                apex,
                settings.apexcenteredscantolerance
            ) &&
            break

        next_center = alkane_ladder_next_center_retention(
            apex,
            raw_scan_retentions,
            input_scanindex,
            settings.maxapexshiftfromguess,
            settings.apexcenteredscantolerance
        )
        isnothing(next_center) && break
        next_center_scan = alkane_ladder_fractional_scan_index(
            raw_scan_retentions,
            next_center
        )
        abs(next_center_scan - apex.fit_center_scan_index) <
            settings.apexcenteredscantolerance &&
            break
        center_retention = next_center
    end

    if !alkane_ladder_has_centered_success(
            attempts,
            settings.apexcenteredscantolerance
        )
        for fallback_center in alkane_ladder_center_search_retentions(
            raw_scan_retentions,
            input_scanindex,
            settings.maxapexshiftfromguess,
            settings.apexcenteredscantolerance
        )
            fallback_scan = alkane_ladder_fractional_scan_index(
                raw_scan_retentions,
                fallback_center
            )
            alkane_ladder_center_already_attempted(
                attempts,
                fallback_scan,
                settings.apexcenteredscantolerance
            ) && continue

            apex = fit_at_center!(fallback_center)
            apex.success &&
                alkane_ladder_apex_is_centered(
                    apex,
                    settings.apexcenteredscantolerance
                ) &&
                break
        end
    end

    final_index = alkane_ladder_best_apex_attempt_index(
        attempts,
        settings.apexcenteredscantolerance
    )
    AlkaneLadderIonApexResult(
        attempts[final_index],
        first(attempts),
        attempts,
        length(attempts) - 1,
        final_index > 1,
        settings.maxapexshiftfromguess
    )
end

function alkane_ladder_ion_apex_once(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidate::Union{AbstractAlkaneLadderCandidate, NamedTuple},
    ladderstep::Integer,
    fitcenterretention::Real,
    attemptindex::Integer,
    settings::AlkaneLadderApexSettings
)
    input_scanindex = alkane_ladder_input_scan_index(candidate)
    raw_scan_retentions = rawretentions(msm)
    input_retention = raw_scan_retentions[input_scanindex]
    fit_center_retention = fitcenterretention
    isfinite(fit_center_retention) || throw(ArgumentError(
        "fitcenterretention must be finite"))
    fit_center_scan_index = alkane_ladder_fractional_scan_index(
        raw_scan_retentions,
        fit_center_retention
    )
    localmargin = max(3, settings.scanwindow + 2)

    nscans = min(scancount(msm), 2 * settings.scanwindow + 1)
    model_scanindices = alkane_ladder_scan_indices_around_retention(
        raw_scan_retentions,
        fit_center_retention,
        nscans
    )
    length(model_scanindices) ≥ 3 || throw(ArgumentError(
        "at least three scans are needed for a log-quadratic apex fit"))

    ion_selection = isnothing(settings.apexionmzvalues) ?
        alkane_ladder_reference_filtered_apex_mzvalues(
            msm,
            settings.standard,
            ladderstep,
            settings.apexionexcludemzvalues,
            settings.apexionminrelativeintensity,
            settings.minioncount
        ) :
        alkane_ladder_explicit_apex_mzvalue_selection(
            msm,
            settings.apexionmzvalues,
            settings.minioncount
        )
    selected_mzindices = alkane_ladder_apex_mz_indices(
        msm,
        ion_selection.mzvalues,
        settings.minioncount
    )

    mzkwargs = isnothing(settings.mzretentionkwargs) ?
        alkane_ladder_default_mzretention_kwargs(msm, :inferdirection) :
        settings.mzretentionkwargs
    scan_retentions = raw_scan_retentions
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
        intensities_by_mz[scancol, mzcol] = max(Xraw[selected_scan, selected_mz], 0.0)
        variances_by_mz[scancol, mzcol] = variances[selected_scan, selected_mz]
        raw_retentions_by_mz[scancol, mzcol] = raw_scan_retentions[selected_scan]
        observation_retentions[scancol, mzcol] = alkane_ladder_observation_retention(
            scan_retentions,
            selected_scan,
            selected_mz,
            mzkwargs
        )
    end

    ymax = maximum(intensities_by_mz)
    ymax > 0 || throw(ArgumentError("selected ions have no positive local signal"))
    logfloor = settings.logfloorfraction * ymax

    xvalues = observation_retentions .- fit_center_retention
    xscale = maximum(abs, xvalues)
    xscale > 0 || throw(ArgumentError("scan window has zero retention span"))

    nobs = nscans * nmz
    design = zeros(Float64, nobs, nmz + 2)
    logy = Vector{Float64}(undef, nobs)
    weights = Vector{Float64}(undef, nobs)
    row = 0
    for mzcol in 1:nmz, scancol in 1:nscans
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
            settings.variancefloor
        )
    end
    any(>(0), weights) || throw(ArgumentError("peak apex fit has no positive weights"))

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
        minimum(xvalues) / xscale ≤ apexx ≤ maximum(xvalues) / xscale
    apex_shift_from_guess_scans = apex_scan_index - input_scanindex
    apex_within_allowed_shift = isfinite(apex_shift_from_guess_scans) &&
        abs(apex_shift_from_guess_scans) ≤ settings.maxapexshiftfromguess + 1e-9
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
    for mzcol in 1:nmz, scancol in 1:nscans
        row += 1
        normalized_log_residuals[scancol, mzcol] =
            sqrt(weights[row]) * (logy[row] - fitted_logy[row])
        fitted_log_intensities[scancol, mzcol] = fitted_logy[row]
    end

    model_x =
        (raw_scan_retentions[model_scanindices] .- fit_center_retention) ./ xscale
    peak_model = success ?
        exp.(gamma .* abs2.(model_x .- apexx)) :
        fill(Float64(NaN), length(model_x))
    observation_peak_model = success ?
        exp.(gamma .* abs2.(xvalues ./ xscale .- apexx)) :
        fill(Float64(NaN), size(xvalues))

    AlkaneLadderIonApexAttempt(
        success,
        apex_retention,
        apex_retention,
        apex_scan_index,
        input_scanindex,
        input_retention,
        apex_retention - input_retention,
        apex_scan_index - input_scanindex,
        fit_center_retention,
        fit_center_scan_index,
        apex_retention - fit_center_retention,
        apex_scan_index - fit_center_scan_index,
        apex_shift_from_guess_scans,
        apex_within_allowed_shift,
        apex_in_window,
        attemptindex,
        beta,
        gamma,
        apexx,
        xscale,
        r2,
        ssres,
        fit_degrees_of_freedom,
        reduced_normalized_residual,
        normalized_log_residuals,
        fitted_log_intensities,
        settings.scanwindow,
        model_scanindices,
        raw_scan_retentions[model_scanindices],
        peak_model,
        observation_peak_model,
        scanindices_by_mz,
        raw_retentions_by_mz,
        observation_retentions,
        intensities_by_mz,
        variances_by_mz,
        selected_mzindices,
        rawmzvalues(msm)[selected_mzindices],
        ion_selection.selection,
        ion_selection.excluded_mzvalues,
        ion_selection.reference_mzvalues,
        ion_selection.reference_relative_intensities,
        settings.apexionminrelativeintensity,
        ion_selection.reference_key,
        coefficients[1:nmz],
        logfloor,
        settings.logfloorfraction,
        true,
        settings.variancefloor,
        mzkwargs
    )
end

function alkane_ladder_apex_public_result(
    apex::AlkaneLadderIonApexResult,
    candidate::Union{AbstractAlkaneLadderCandidate, NamedTuple},
    abundance::AbstractVector{Float64},
    retentions::AbstractVector{<:Real}
)
    step = alkane_ladder_candidate_step(candidate)
    fit = apex.fit
    scanindex = fit.input_scan_index
    reason = fit.success ? :success : Symbol(
        replace(alkane_ladder_apex_failure_reason(apex), ' ' => '_')
    )
    failurereason = fit.success ? nothing : alkane_ladder_apex_failure_reason(apex)
    qualityscore = fit.reduced_normalized_residual
    AlkaneLadderApex(
        fit.success,
        reason,
        failurereason,
        fit.success,
        !fit.success,
        step,
        scanindex,
        retentions[scanindex],
        fit.apex_scan_index,
        fit.apex_retention,
        fit.apex_scan_index,
        fit.apex_retention,
        fit.apex_scan_index,
        fit.apex_retention,
        fit.apex_offset_scans,
        fit.apex_offset_retention,
        abundance[scanindex],
        fit.scan_indices,
        abundance[fit.scan_indices],
        fit.x_scale,
        fit.weighted_log_residual_sum_squares,
        vcat(fit.ion_intercepts, fit.beta, fit.gamma),
        fit.log_floor,
        alkane_ladder_candidate_source(candidate),
        alkane_ladder_candidate_is_gapfilled(candidate),
        alkane_ladder_candidate_is_edgeextended(candidate),
        alkane_ladder_candidate_mass_spectrum_cosine(candidate),
        alkane_ladder_candidate_required_cosine(candidate),
        qualityscore,
        isfinite(qualityscore) && qualityscore > 0 ? log(qualityscore) : NaN,
        NaN,
        false,
        NaN,
        NaN,
        NaN,
        0,
        false,
        false,
        nothing,
        fit.success,
        candidate,
        apex
    )
end

function alkane_ladder_annotate_apex_fit_quality(
    apexes::AbstractVector{<:AlkaneLadderApex},
    apexfitqualityoutlierz::Union{Nothing, Real},
    apexfitqualityminsteps::Integer
)
    if !isnothing(apexfitqualityoutlierz)
        isfinite(apexfitqualityoutlierz) && apexfitqualityoutlierz > 0 ||
            throw(ArgumentError(
                "apexfitqualityoutlierz must be finite and positive, or nothing"))
    end
    apexfitqualityminsteps ≥ 2 || throw(ArgumentError(
        "apexfitqualityminsteps must be at least 2"))

    qualityindices = Int[]
    logscores = Float64[]
    for (index, apex) in pairs(apexes)
        score = alkane_ladder_apex_fit_quality_score(apex)
        if apex.refined && isfinite(score) && score > 0
            push!(qualityindices, index)
            push!(logscores, log(score))
        end
    end

    quality = alkane_ladder_apex_fit_quality_robust_zscores(
        logscores,
        apexfitqualityminsteps
    )
    zbyindex = Dict{Int, Float64}()
    for (index, zscore) in zip(qualityindices, quality.zscores)
        zbyindex[index] = zscore
    end

    updated = AlkaneLadderApex[]
    for (index, apex) in pairs(apexes)
        score = alkane_ladder_apex_fit_quality_score(apex)
        zscore = get(zbyindex, index, NaN)
        zexcluded = !isnothing(apexfitqualityoutlierz) &&
            apex.refined &&
            quality.applied &&
            isfinite(zscore) &&
            zscore > apexfitqualityoutlierz
        push!(updated, alkane_ladder_apex_with_fit_quality(
                apex,
                score,
                zscore,
                zexcluded,
                quality,
                apexfitqualityoutlierz
            )
        )
    end

    updated
end

function alkane_ladder_apex_with_fit_quality(
    apex::AlkaneLadderApex,
    score::Real,
    zscore::Real,
    zexcluded::Bool,
    quality::AlkaneLadderApexFitQualityStats,
    apexfitqualityoutlierz::Union{Nothing, Real}
)
    AlkaneLadderApex(
        apex.success,
        apex.reason,
        apex.failurereason,
        apex.refined,
        apex.fallback,
        apex.ladderstep,
        apex.scanindex,
        apex.retention,
        apex.apexscanindex,
        apex.apexretention,
        apex.output_scan_index,
        apex.output_retention,
        apex.outputscanindex,
        apex.outputretention,
        apex.apexoffsetscans,
        apex.apexoffsetretention,
        apex.apexabundance,
        apex.scanindices,
        apex.abundance,
        apex.xscale,
        apex.weightedlogresidualsumsquares,
        apex.fit_coefficients,
        apex.logfloor,
        apex.source,
        apex.gapfilled,
        apex.edgeextended,
        apex.mass_spectrum_cosine,
        apex.required_cosine,
        score,
        isfinite(score) && score > 0 ? log(score) : NaN,
        zscore,
        zexcluded,
        isnothing(apexfitqualityoutlierz) ? NaN : apexfitqualityoutlierz,
        quality.median,
        quality.mad,
        quality.nsteps,
        quality.applied,
        zexcluded,
        zexcluded ? "apex curve fit residual outlier" : nothing,
        apex.refined && !zexcluded,
        apex.candidate,
        apex.fit
    )
end

alkane_ladder_apex_fit_quality_score(apex::AlkaneLadderApex) = apex.apex_fit_quality_score

function alkane_ladder_apex_fit_quality_robust_zscores(
    logscores::AbstractVector{Float64},
    minsteps::Integer
)
    nsteps = length(logscores)
    zscores = fill(NaN, nsteps)
    nsteps ≥ minsteps ||
        return AlkaneLadderApexFitQualityStats(zscores, NaN, NaN, nsteps, false)
    values = logscores
    center = median(values)
    deviations = abs.(values .- center)
    mad = median(deviations)
    isfinite(mad) && mad > eps(Float64) ||
        return AlkaneLadderApexFitQualityStats(zscores, center, mad, nsteps, false)

    zscores .= 0.67448975 .* (values .- center) ./ mad
    AlkaneLadderApexFitQualityStats(zscores, center, mad, nsteps, true)
end

"""
    alkaneladderscanorder(msm, variances, pathinfo, settings)

Resolve the m/z scan timing model used by alkane ladder apex fitting.

`settings` is the `AlkaneLadderApexSettings` object built by the ladder finder from the
user-facing apex keywords. If `settings.mzscanorder` or
`settings.mzretentionkwargs.order` is `:ascending`, `:descending`, or `:simultaneous`,
that order is used directly. If it is `:inferdirection`, the function compares only the
sequential quadrupole directions `:ascending` and `:descending`.
"""
function alkaneladderscanorder(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    pathinfo::AlkaneLadderPathInfo,
    settings::AlkaneLadderApexSettings
)
    validate_alkane_series_variances(msm, variances)
    validate_alkane_ladder_scan_order_settings(
        settings.mzscanordermaxpeaks,
        settings.mzscanorderminpeaks,
        settings.mzscanorderminapexvarianceratio,
        settings.mzscanordershapeioncount,
        settings.mzscanordershapemzspacing,
        settings.mzscanorderextremeioncount,
        settings.mzscanorderminioncount
    )

    base_mzkwargs = alkane_ladder_base_mzretention_kwargs(
        msm,
        settings.mzretentionkwargs,
        settings.mzscanorder
    )
    order = alkane_ladder_mzretention_order(base_mzkwargs)
    if order ≢ :inferdirection
        order in (:ascending, :descending, :simultaneous) || throw(ArgumentError(
            "mzretentionkwargs.order must be :ascending, :descending, or :simultaneous"))
        return AlkaneLadderScanOrderResult(
            base_mzkwargs,
            alkane_ladder_provided_scan_order_info(order)
        )
    end

    maxpeaks = settings.mzscanordermaxpeaks
    candidates = alkane_ladder_scan_order_candidates(pathinfo, maxpeaks)
    initialcount = isnothing(maxpeaks) ?
        length(candidates) :
        min(maxpeaks, length(candidates))
    _, _, info = alkane_ladder_scan_order_expanding_trials(
        msm,
        variances,
        candidates,
        base_mzkwargs,
        initialcount,
        settings
    )

    AlkaneLadderScanOrderResult(merge(base_mzkwargs, (order=info.selected_order,)), info)
end

function validate_alkane_ladder_scan_order_settings(
    maxpeaks::Union{Nothing, Integer},
    minpeaks::Integer,
    minapexvarianceratio::Real,
    shapeioncount::Integer,
    shapemzspacing::Integer,
    extremeioncount::Integer,
    minioncount::Integer
)
    isnothing(maxpeaks) || maxpeaks ≥ 1 || throw(ArgumentError(
        "scan-order maxpeaks must be nothing or a positive integer"))
    minpeaks ≥ 1 || throw(ArgumentError(
        "scan-order minpeaks must be a positive integer"))
    isfinite(minapexvarianceratio) && minapexvarianceratio > 0 || throw(ArgumentError(
        "scan-order minapexvarianceratio must be finite and positive"))
    shapeioncount ≥ 2 || throw(ArgumentError(
        "scan-order shapeioncount must be an integer ≥ 2"))
    shapemzspacing ≥ 1 || throw(ArgumentError(
        "scan-order shapemzspacing must be a positive integer"))
    extremeioncount ≥ 1 || throw(ArgumentError(
        "scan-order extremeioncount must be a positive integer"))
    minioncount ≥ 2 || throw(ArgumentError(
        "scan-order minioncount must be an integer ≥ 2"))

    nothing
end

function alkane_ladder_base_mzretention_kwargs(
    msm::MassScanMatrix,
    mzretentionkwargs::Union{Nothing, NamedTuple},
    mzscanorder::Symbol
)
    mzscanorder in (:inferdirection, :ascending, :descending, :simultaneous) ||
        throw(ArgumentError(
            "mzscanorder must be :inferdirection, :ascending, :descending, or :simultaneous"))

    if isnothing(mzretentionkwargs)
        return alkane_ladder_default_mzretention_kwargs(msm, mzscanorder)
    end

    if :order in keys(mzretentionkwargs)
        mzorder = mzretentionkwargs.order
        mzorder in (:ascending, :descending, :simultaneous) ||
            throw(ArgumentError(
                "mzretentionkwargs.order must be :ascending, :descending, or :simultaneous"))
        if mzscanorder ≢ :inferdirection && mzscanorder ≢ mzorder
            throw(ArgumentError(
                "mzscanorder=$(repr(mzscanorder)) conflicts with " *
                "mzretentionkwargs.order=$(repr(mzorder))"))
        end

        return mzretentionkwargs
    end

    merge(mzretentionkwargs, (order=mzscanorder,))
end

function alkane_ladder_provided_scan_order_info(order::Symbol)
    AlkaneLadderScanOrderInfo(
        order,
        order,
        :provided,
        :provided,
        nothing,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        0,
        0,
        0,
        false,
        false,
        AlkaneLadderScanOrderTrials(nothing, nothing)
    )
end

function alkane_ladder_no_scan_order_info(reason::Symbol)
    AlkaneLadderScanOrderInfo(
        :unknown,
        :inferdirection,
        reason,
        :not_applicable,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        0,
        0,
        0,
        false,
        false,
        AlkaneLadderScanOrderTrials(nothing, nothing)
    )
end

function alkane_ladder_mzretention_order(mzretentionkwargs::NamedTuple)
    :order in keys(mzretentionkwargs) || throw(ArgumentError(
        "mzretentionkwargs must contain an order field"))

    mzretentionkwargs.order
end

function alkane_ladder_scan_order_candidates(
    pathinfo::AlkaneLadderPathInfo,
    maxpeaks
)
    pairs = AlkaneLadderScanOrderCandidate[]
    for candidate in pathinfo.path
        push!(pairs, AlkaneLadderScanOrderCandidate(
                alkane_ladder_candidate_step(candidate),
                alkane_ladder_input_scan_index(candidate),
                candidate
            )
        )
    end
    sort!(pairs; by=pair -> pair.ladderstep)
    isnothing(maxpeaks) && return pairs
    length(pairs) ≤ maxpeaks && return pairs

    order = alkane_ladder_spread_index_order(length(pairs), maxpeaks)
    pairs[order]
end

function alkane_ladder_spread_index_order(n::Integer, initialcount::Integer)
    n ≥ 0 || throw(ArgumentError("number of candidates must be nonnegative"))
    n == 0 && return Int[]
    initialcount ≥ 1 || throw(ArgumentError("initialcount must be positive"))

    selected = Int[]
    seen = Set{Int}()
    firstcount = min(initialcount, n)
    for index in alkane_ladder_spread_indices(n, firstcount)
        push!(selected, index)
        push!(seen, index)
    end

    while length(selected) < n
        bestindex = 0
        bestdistance = -1
        for index in 1:n
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
    count = min(count, n)
    count ≤ 0 && return Int[]
    count == 1 && return [ceil(Int, n / 2)]

    unique([clamp(round(Int, value), 1, n) for value in range(1, n; length=count)])
end

function alkane_ladder_scan_order_expanding_trials(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidates::AbstractVector{<:AlkaneLadderScanOrderCandidate},
    base_mzretentionkwargs::NamedTuple,
    initialcount::Integer,
    settings::AlkaneLadderApexSettings
)
    tested = AlkaneLadderScanOrderCandidate[]
    ascending_results = AlkaneLadderScanOrderTrialResult[]
    descending_results = AlkaneLadderScanOrderTrialResult[]
    initialcount = min(initialcount, length(candidates))

    last = nothing
    expanded_for_ambiguity = false
    for candidateinfo in candidates
        push!(tested, candidateinfo)
        push!(ascending_results, alkane_ladder_scan_order_trial_candidate(
                msm,
                variances,
                candidateinfo,
                :ascending,
                base_mzretentionkwargs,
                settings,
                :extreme_reference_ions
            )
        )
        push!(descending_results, alkane_ladder_scan_order_trial_candidate(
                msm,
                variances,
                candidateinfo,
                :descending,
                base_mzretentionkwargs,
                settings,
                :extreme_reference_ions
            )
        )

        length(tested) ≥ initialcount || continue

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
            descending,
            settings.mzscanorderminpeaks,
            settings.mzscanorderminapexvarianceratio
        )
        info = alkane_ladder_scan_order_info_with_expansion(
            info,
            initialcount,
            expanded_for_ambiguity,
            false
        )
        last = (ascending, descending, info)

        if expanded_for_ambiguity
            length(tested) == length(candidates) && break
            continue
        end

        info.status ≡ :accepted && return last
        if info.status ≡ :ambiguous
            length(tested) == length(candidates) && break
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
            descending,
            settings.mzscanorderminpeaks,
            settings.mzscanorderminapexvarianceratio
        )
        info = alkane_ladder_scan_order_info_with_expansion(
            info,
            initialcount,
            expanded_for_ambiguity,
            false
        )
        last = (ascending, descending, info)
    end

    last[3].status ≡ :accepted && return last

    ascending, descending, info = alkane_ladder_scan_order_all_ion_trials(
        msm,
        variances,
        candidates,
        base_mzretentionkwargs,
        initialcount,
        settings
    )

    ascending, descending, info
end

function alkane_ladder_scan_order_info_with_expansion(
    info::AlkaneLadderScanOrderInfo,
    initialcount::Integer,
    expanded_for_ambiguity::Bool,
    expanded_ions::Bool
)
    ntried = max(
        isnothing(info.trials.ascending) ? 0 : info.trials.ascending.n_tried,
        isnothing(info.trials.descending) ? 0 : info.trials.descending.n_tried
    )

    AlkaneLadderScanOrderInfo(
        info.selected_order,
        info.requested_order,
        info.status,
        info.score_kind,
        info.evidence_score,
        info.ascending_score,
        info.descending_score,
        info.ascending_median_apex_variance,
        info.descending_median_apex_variance,
        info.ascending_median_width,
        info.descending_median_width,
        info.n_peaks_used,
        ntried,
        initialcount,
        expanded_for_ambiguity,
        expanded_ions,
        info.trials
    )
end

function alkane_ladder_scan_order_all_ion_trials(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidates::AbstractVector{<:AlkaneLadderScanOrderCandidate},
    base_mzretentionkwargs::NamedTuple,
    initialcount::Integer,
    settings::AlkaneLadderApexSettings
)
    ascending_results = AlkaneLadderScanOrderTrialResult[]
    descending_results = AlkaneLadderScanOrderTrialResult[]

    for candidateinfo in candidates
        push!(ascending_results, alkane_ladder_scan_order_trial_candidate(
                msm,
                variances,
                candidateinfo,
                :ascending,
                base_mzretentionkwargs,
                settings,
                :all_usable_reference_ions
            )
        )
        push!(descending_results, alkane_ladder_scan_order_trial_candidate(
                msm,
                variances,
                candidateinfo,
                :descending,
                base_mzretentionkwargs,
                settings,
                :all_usable_reference_ions
            )
        )
    end

    ascending = alkane_ladder_scan_order_trial_summary(
        :ascending,
        candidates,
        ascending_results
    )
    descending = alkane_ladder_scan_order_trial_summary(
        :descending,
        candidates,
        descending_results
    )
    info = alkane_ladder_choose_scan_order(
        ascending,
        descending,
        settings.mzscanorderminpeaks,
        settings.mzscanorderminapexvarianceratio
    )
    info = alkane_ladder_scan_order_info_with_expansion(
        info,
        initialcount,
        true,
        true
    )

    ascending, descending, info
end

function alkane_ladder_scan_order_trial_candidate(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidateinfo::AlkaneLadderScanOrderCandidate,
    order::Symbol,
    base_mzretentionkwargs::NamedTuple,
    settings::AlkaneLadderApexSettings,
    ionstrategy::Symbol
)
    trial_mzkwargs = merge(base_mzretentionkwargs, (order=order,))
    candidate = candidateinfo.candidate
    step = candidateinfo.ladderstep
    try
        shape_selection = alkane_ladder_scan_order_shape_ion_selection(
            msm,
            settings.standard,
            step,
            settings.mzscanordershapeioncount,
            settings.mzscanordershapemzspacing
        )
        shape_settings = alkane_ladder_apex_settings_with_ion_selection(
            settings,
            trial_mzkwargs,
            shape_selection.mzvalues,
            length(shape_selection.mzvalues)
        )
        shape = alkane_ladder_ion_apex(
            msm,
            variances,
            candidate,
            step,
            shape_settings
        )
        if !alkane_ladder_scan_order_shape_is_usable(shape)
            return AlkaneLadderScanOrderTrialResult{AlkaneLadderScanOrderTrialFailure}(
                false,
                NaN,
                NaN,
                0,
                NaN,
                AlkaneLadderScanOrderTrialFailure(
                    step,
                    candidateinfo.scanindex,
                    "shape fit failed: $(alkane_ladder_apex_failure_reason(shape))",
                    0,
                    NaN,
                    alkane_ladder_apex_width(shape),
                    shape_selection.mzvalues
                )
            )
        end

        score = alkane_ladder_scan_order_step_apex_variance(
            msm,
            variances,
            candidateinfo,
            shape,
            trial_mzkwargs,
            settings,
            ionstrategy
        )
        score.success && return AlkaneLadderScanOrderTrialResult{Nothing}(
            true,
            score.apex_variance,
            score.shape_width,
            score.ion_count,
            score.mz_span,
            nothing
        )

        AlkaneLadderScanOrderTrialResult{AlkaneLadderScanOrderTrialFailure}(
            false,
            score.apex_variance,
            score.shape_width,
            score.ion_count,
            score.mz_span,
            AlkaneLadderScanOrderTrialFailure(
                step,
                candidateinfo.scanindex,
                score.reason,
                score.ion_count,
                score.mz_span,
                score.shape_width,
                Float64[]
            )
        )
    catch err
        AlkaneLadderScanOrderTrialResult{AlkaneLadderScanOrderTrialFailure}(
            false,
            NaN,
            NaN,
            0,
            NaN,
            AlkaneLadderScanOrderTrialFailure(
                step,
                candidateinfo.scanindex,
                sprint(showerror, err),
                0,
                NaN,
                NaN,
                Float64[]
            )
        )
    end
end

function alkane_ladder_scan_order_trial_summary(
    order::Symbol,
    candidates::AbstractVector{<:AlkaneLadderScanOrderCandidate},
    results::AbstractVector{<:AlkaneLadderScanOrderTrialResult}
)
    apex_variances = [result.apex_variance for result in results if result.success]
    shape_widths = [result.shape_width for result in results if result.success]
    ion_counts = [result.ion_count for result in results if result.success]
    mz_spans = [result.mz_span for result in results if result.success]
    failures = [result.failure for result in results if !result.success]
    positive_variances = [max(variance, eps(Float64)) for variance in apex_variances]
    score = isempty(positive_variances) ? Inf : median(log.(positive_variances))
    finite_spans = filter(isfinite, mz_spans)

    AlkaneLadderScanOrderTrialSummary(
        order,
        score,
        :fixed_shape_ion_apex_variance,
        isempty(apex_variances) ? Inf : median(apex_variances),
        apex_variances,
        isempty(shape_widths) ? Inf : median(shape_widths),
        shape_widths,
        shape_widths,
        ion_counts,
        mz_spans,
        isempty(ion_counts) ? 0.0 : median(Float64.(ion_counts)),
        isempty(finite_spans) ? NaN : median(finite_spans),
        length(apex_variances),
        length(candidates),
        failures
    )
end

function alkane_ladder_scan_order_step_apex_variance(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    candidateinfo::AlkaneLadderScanOrderCandidate,
    shape::AlkaneLadderIonApexResult,
    mzretentionkwargs::NamedTuple,
    settings::AlkaneLadderApexSettings,
    ionstrategy::Symbol
)
    candidate = candidateinfo.candidate
    step = alkane_ladder_candidate_step(candidate)
    ion_selection = if ionstrategy ≡ :extreme_reference_ions
        alkane_ladder_scan_order_extreme_apex_mzvalues(
            msm,
            settings.standard,
            step,
            settings.mzscanorderextremeioncount
        )
    elseif ionstrategy ≡ :all_usable_reference_ions
        alkane_ladder_scan_order_all_usable_apex_mzvalues(
            msm,
            settings.standard,
            step,
            settings.apexionexcludemzvalues
        )
    else
        throw(ArgumentError("unknown scan-order ion strategy $(repr(ionstrategy))"))
    end
    required_ion_count = min(
        settings.mzscanorderminioncount,
        length(ion_selection.mzvalues)
    )
    mzindices = alkane_ladder_apex_mz_indices(
        msm,
        ion_selection.mzvalues,
        required_ion_count
    )
    fits = alkane_ladder_fixed_shape_ion_apices(
        msm,
        shape,
        mzindices,
        variances,
        mzretentionkwargs,
        alkane_ladder_input_scan_index(candidate),
        settings
    )
    successful = [fit for fit in fits if fit.success]
    ion_count = length(successful)
    mz_values = [fit.mz_value for fit in successful]
    mz_span = ion_count ≥ 2 ? maximum(mz_values) - minimum(mz_values) : NaN
    shape_width = alkane_ladder_apex_width(shape)
    if ion_count < required_ion_count
        return AlkaneLadderScanOrderStepApexVariance(
            false,
            "fewer than $(required_ion_count) fixed-shape ion apices were inferred",
            ion_count,
            mz_span,
            shape_width,
            NaN,
            Float64[],
            ion_selection,
            fits
        )
    end

    apex_scan_indices = [fit.apex_scan_index for fit in successful]
    apex_variance = alkane_ladder_robust_variance(apex_scan_indices)
    isfinite(apex_variance) || return AlkaneLadderScanOrderStepApexVariance(
        false,
        "fixed-shape ion apex variance is not finite",
        ion_count,
        mz_span,
        shape_width,
        apex_variance,
        apex_scan_indices,
        ion_selection,
        fits
    )

    AlkaneLadderScanOrderStepApexVariance(
        true,
        nothing,
        ion_count,
        mz_span,
        shape_width,
        apex_variance,
        apex_scan_indices,
        ion_selection,
        fits
    )
end

function alkane_ladder_choose_scan_order(
    ascending::AlkaneLadderScanOrderTrialSummary,
    descending::AlkaneLadderScanOrderTrialSummary,
    minpeaks::Integer,
    minapexvarianceratio::Real
)
    if !isfinite(ascending.score) && !isfinite(descending.score)
        return AlkaneLadderScanOrderInfo(
            :unknown,
            :inferdirection,
            :failed,
            :fixed_shape_ion_apex_variance,
            NaN,
            ascending.score,
            descending.score,
            ascending.median_apex_variance,
            descending.median_apex_variance,
            ascending.median_width,
            descending.median_width,
            0,
            max(ascending.n_tried, descending.n_tried),
            0,
            false,
            false,
            AlkaneLadderScanOrderTrials(ascending, descending)
        )
    end

    best = descending.score ≤ ascending.score ? descending : ascending
    other = best.order ≡ :descending ? ascending : descending
    evidence_score = isfinite(other.score) && isfinite(best.score) ?
        exp(other.score - best.score) :
        Inf
    status = best.n_success < minpeaks ?
        :insufficient_peaks :
        evidence_score ≥ minapexvarianceratio ?
            :accepted :
            :ambiguous

    AlkaneLadderScanOrderInfo(
        best.order,
        :inferdirection,
        status,
        :fixed_shape_ion_apex_variance,
        evidence_score,
        ascending.score,
        descending.score,
        ascending.median_apex_variance,
        descending.median_apex_variance,
        ascending.median_width,
        descending.median_width,
        best.n_success,
        max(ascending.n_tried, descending.n_tried),
        0,
        false,
        false,
        AlkaneLadderScanOrderTrials(ascending, descending)
    )
end

function alkane_ladder_fixed_shape_ion_apices(
    msm::MassScanMatrix,
    shape::AlkaneLadderIonApexResult,
    mzindices::AbstractVector{<:Integer},
    variances::AbstractMatrix{<:Real},
    mzretentionkwargs::NamedTuple,
    inputscanindex::Integer,
    settings::AlkaneLadderApexSettings
)
    alkane_ladder_scan_order_shape_is_usable(shape) || throw(ArgumentError(
        "shape apex is not usable for fixed-shape ion apex inference"))
    validate_alkane_series_variances(msm, variances)
    isfinite(settings.logfloorfraction) && settings.logfloorfraction > 0 ||
        throw(ArgumentError("logfloorfraction must be finite and positive"))
    validate_alkane_abundance_variancefloor(settings.variancefloor)
    isfinite(settings.maxapexshiftfromguess) && settings.maxapexshiftfromguess ≥ 0 ||
        throw(ArgumentError("maxapexshiftfromguess must be finite and nonnegative"))

    raw_scan_retentions = rawretentions(msm)
    scan_retentions = raw_scan_retentions
    Xraw = rawintensities(msm)
    shapefit = shape.fit
    fit_center_retention = alkane_ladder_fit_center_retention(shape)
    fit_center_scan_index = alkane_ladder_fractional_scan_index(
        raw_scan_retentions,
        fit_center_retention
    )
    nscans = length(shapefit.scan_indices)
    localmargin = max(3, fld(nscans, 2) + 2)
    gamma = shapefit.gamma
    xscale = shapefit.x_scale
    fits = AlkaneLadderFixedShapeIonApexFit[]

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
        intensities = [max(Xraw[scanindex, mzindex], 0.0) for scanindex in selected_scans]
        local_variances = [variances[scanindex, mzindex] for scanindex in selected_scans]
        observation_retentions = [
            alkane_ladder_observation_retention(
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
            local_variances,
            inputscanindex,
            settings
        )
        push!(fits, AlkaneLadderFixedShapeIonApexFit(
                fit.success,
                fit.reason,
                fit.apex_scan_index,
                fit.apex_retention,
                fit.apex_x,
                fit.beta,
                fit.intercept,
                fit.reduced_normalized_residual,
                fit.apex_shift_from_guess_scans,
                fit.apex_within_allowed_shift,
                fit.apex_in_window,
                mzindex,
                Float64(rawmzvalues(msm)[mzindex]),
                selected_scans,
                observation_retentions,
                intensities,
                local_variances
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
    observation_retentions::AbstractVector{<:Real},
    intensities::AbstractVector{<:Real},
    variances::AbstractVector{<:Real},
    inputscanindex::Integer,
    settings::AlkaneLadderApexSettings
)
    gamma < 0 && isfinite(gamma) || return AlkaneLadderFixedShapeIonFit(
        false,
        "fixed shape is not concave",
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        false,
        false
    )
    xscale > 0 && isfinite(xscale) || throw(ArgumentError(
        "fixed shape xscale must be finite and positive"))

    nobs = length(intensities)
    nobs == length(variances) == length(observation_retentions) || throw(DimensionMismatch(
        "fixed-shape ion fit vectors must have equal length"))
    nobs ≥ 3 || throw(ArgumentError(
        "at least three observations are required for fixed-shape ion apex inference"))
    ymax = maximum(intensities)
    ymax > 0 || return AlkaneLadderFixedShapeIonFit(
        false,
        "ion has no positive local signal",
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        false,
        false
    )

    logfloor = settings.logfloorfraction * ymax
    adjusted = intensities .+ logfloor
    x = (observation_retentions .- fit_center_retention) ./ xscale
    logy = log.(adjusted)
    target = logy .- gamma .* abs2.(x)
    weights = [
        alkane_ladder_log_weight(
            adjusted[index],
            variances[index],
            settings.variancefloor
        )
        for index in eachindex(adjusted)
    ]
    any(>(0), weights) || return AlkaneLadderFixedShapeIonFit(
        false,
        "ion fixed-shape fit has no positive weights",
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        false,
        false
    )

    design = hcat(ones(Float64, nobs), x)
    sqrtweights = sqrt.(weights)
    coefficients = (design .* reshape(sqrtweights, :, 1)) \ (target .* sqrtweights)
    intercept = coefficients[1]
    beta = coefficients[2]
    fitted_logy = intercept .+ beta .* x .+ gamma .* abs2.(x)
    apexx = -beta / (2 * gamma)
    apex_retention = fit_center_retention + apexx * xscale
    apex_scan_index = isfinite(apex_retention) ?
        alkane_ladder_fractional_scan_index(raw_scan_retentions, apex_retention) :
        NaN
    apex_in_window = isfinite(apexx) && minimum(x) ≤ apexx ≤ maximum(x)
    apex_shift_from_guess_scans = apex_scan_index - inputscanindex
    apex_within_allowed_shift = isfinite(apex_shift_from_guess_scans) &&
        abs(apex_shift_from_guess_scans) ≤ settings.maxapexshiftfromguess + 1e-9
    success = isfinite(apex_scan_index) && apex_in_window && apex_within_allowed_shift
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

    AlkaneLadderFixedShapeIonFit(
        success,
        reason,
        apex_scan_index,
        apex_retention,
        apexx,
        beta,
        intercept,
        reduced_normalized_residual,
        apex_shift_from_guess_scans,
        apex_within_allowed_shift,
        apex_in_window
    )
end

function alkane_ladder_scan_order_shape_is_usable(shape::AlkaneLadderIonApexAttempt)
    shape.success || return false
    width = alkane_ladder_apex_width(shape)

    isfinite(width) &&
        width > 0 &&
        length(shape.scan_indices) ≥ 3 &&
        isfinite(alkane_ladder_fit_center_retention(shape))
end

function alkane_ladder_scan_order_shape_is_usable(shape::AlkaneLadderIonApexResult)
    alkane_ladder_scan_order_shape_is_usable(shape.fit)
end

alkane_ladder_fit_center_retention(apex::AlkaneLadderIonApexAttempt) =
    apex.fit_center_retention

function alkane_ladder_fit_center_retention(apex::AlkaneLadderIonApexResult)
    alkane_ladder_fit_center_retention(apex.fit)
end

function alkane_ladder_robust_variance(values::AbstractVector{<:Real})
    finite_values = [value for value in values if isfinite(value)]
    length(finite_values) ≥ 2 || return NaN
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

function alkane_ladder_apex_width(apex::AlkaneLadderIonApexAttempt)
    apex.gamma < 0 &&
        isfinite(apex.gamma) &&
        isfinite(apex.x_scale) &&
        apex.x_scale > 0 ||
        return NaN

    apex.x_scale / sqrt(-apex.gamma)
end

function alkane_ladder_apex_width(apex::AlkaneLadderIonApexResult)
    alkane_ladder_apex_width(apex.fit)
end

function alkane_ladder_apex_failure_reason(apex::AlkaneLadderIonApexAttempt)
    apex.success && return "success"
    !isfinite(apex.apex_scan_index) && return "nonfinite apex scan index"
    !isfinite(apex.apex_retention) && return "nonfinite apex retention"
    !apex.apex_within_allowed_shift &&
        return "apex exceeds maxapexshiftfromguess from abundance-trace guess"
    apex.gamma ≥ 0 && return "log-quadratic fit is not concave"
    !apex.apex_in_window && return "continuous apex is outside scan window"

    "apex fit failed"
end

function alkane_ladder_apex_failure_reason(apex::AlkaneLadderIonApexResult)
    alkane_ladder_apex_failure_reason(apex.fit)
end

function alkane_ladder_discrete_apex(
    candidate::Union{AbstractAlkaneLadderCandidate, NamedTuple},
    retentions::AbstractVector{<:Real},
    abundance::AbstractVector{Float64},
    scanindices::AbstractVector{<:Integer},
    reason::Symbol,
    failurereason::String
)
    scanindex = alkane_ladder_input_scan_index(candidate)
    apexscanindex = Float64(scanindex)
    scanindicesvector = collect(scanindices)
    AlkaneLadderApex(
        false,
        reason,
        failurereason,
        false,
        true,
        alkane_ladder_candidate_step(candidate),
        scanindex,
        retentions[scanindex],
        apexscanindex,
        retentions[scanindex],
        apexscanindex,
        retentions[scanindex],
        apexscanindex,
        retentions[scanindex],
        0.0,
        0.0,
        abundance[scanindex],
        scanindicesvector,
        abundance[scanindicesvector],
        NaN,
        NaN,
        Float64[],
        NaN,
        alkane_ladder_candidate_source(candidate),
        alkane_ladder_candidate_is_gapfilled(candidate),
        alkane_ladder_candidate_is_edgeextended(candidate),
        alkane_ladder_candidate_mass_spectrum_cosine(candidate),
        alkane_ladder_candidate_required_cosine(candidate),
        NaN,
        NaN,
        NaN,
        false,
        NaN,
        NaN,
        NaN,
        0,
        false,
        false,
        nothing,
        false,
        candidate,
        nothing
    )
end

alkane_ladder_candidate_step(candidate::AlkaneLadderPathCandidate) =
    candidate.ladderstep

function alkane_ladder_candidate_step(candidate::NamedTuple)
    :ladderstep in keys(candidate) && return candidate.ladderstep
    :step in keys(candidate) && return candidate.step
    throw(ArgumentError("candidate must contain ladderstep"))
end

alkane_ladder_candidate_source(candidate::AlkaneLadderPathCandidate) = candidate.source

function alkane_ladder_candidate_source(candidate::NamedTuple)
    :source in keys(candidate) && return candidate.source

    :molecular_ion_dp
end

function alkane_ladder_candidate_is_gapfilled(candidate::AlkaneLadderPathCandidate)
    candidate.gapfilled
end

function alkane_ladder_candidate_is_gapfilled(candidate::NamedTuple)
    :gapfilled in keys(candidate) && return candidate.gapfilled
    :source in keys(candidate) && return candidate.source ≡ :gapfilled

    false
end

function alkane_ladder_candidate_is_edgeextended(candidate::AlkaneLadderPathCandidate)
    candidate.source in (:leftextended, :rightextended)
end

function alkane_ladder_candidate_is_edgeextended(candidate::NamedTuple)
    :edgeextended in keys(candidate) && return candidate.edgeextended
    :source in keys(candidate) &&
        return candidate.source in (:leftextended, :rightextended)

    false
end

function alkane_ladder_candidate_mass_spectrum_cosine(candidate::AlkaneLadderPathCandidate)
    candidate.massspectrumcosine
end

function alkane_ladder_candidate_mass_spectrum_cosine(candidate::NamedTuple)
    :massspectrumcosine in keys(candidate) && return candidate.massspectrumcosine

    NaN
end

function alkane_ladder_candidate_required_cosine(candidate::AlkaneLadderPathCandidate)
    candidate.requiredcosine
end

function alkane_ladder_candidate_required_cosine(candidate::NamedTuple)
    :requiredcosine in keys(candidate) && return candidate.requiredcosine

    NaN
end

alkane_ladder_input_scan_index(candidate::AlkaneLadderPathCandidate) = candidate.scanindex

function alkane_ladder_input_scan_index(candidate::NamedTuple)
    if :scanindex in keys(candidate)
        return candidate.scanindex
    elseif :apexindex in keys(candidate)
        return candidate.apexindex
    elseif :scan_index in keys(candidate)
        return candidate.scan_index
    end
    throw(ArgumentError("candidate must contain scanindex, apexindex, or scan_index"))
end

function alkane_ladder_candidate_abundance(
    abundanceinfo::AlkaneAbundanceInfo,
    step::Integer,
    retentions::AbstractVector{<:Real}
)
    haskey(abundanceinfo.abundances, step) || throw(ArgumentError(
        "abundanceinfo.abundances does not contain C$(step)"))
    abundance = alkane_abundance_values(abundanceinfo.abundances[step], step)
    length(abundance) == length(retentions) || throw(DimensionMismatch(
        "abundance vector for C$(step) must match retention length"))

    abundance
end

function alkane_ladder_fallback_scanindices(
    retentions::AbstractVector{<:Real},
    scanindex::Integer,
    scanwindow::Integer,
    candidate::AlkaneLadderPathCandidate
)
    leftindex = max(1, scanindex - scanwindow)
    rightindex = min(length(retentions), scanindex + scanwindow)
    leftindex = max(leftindex, candidate.window.leftindex)
    rightindex = min(rightindex, candidate.window.rightindex)

    collect(leftindex:rightindex)
end

function alkane_ladder_fallback_scanindices(
    retentions::AbstractVector{<:Real},
    scanindex::Integer,
    scanwindow::Integer,
    candidate::NamedTuple
)
    leftindex = max(1, scanindex - scanwindow)
    rightindex = min(length(retentions), scanindex + scanwindow)
    if :window in keys(candidate)
        leftindex = max(leftindex, candidate.window.leftindex)
        rightindex = min(rightindex, candidate.window.rightindex)
    end

    collect(leftindex:rightindex)
end

function alkane_ladder_reference_filtered_apex_mzvalues(
    msm::MassScanMatrix,
    standard::AlkaneStandard,
    ladderstep::Integer,
    excludemzvalues::Union{AbstractVector{<:Real}, Tuple{Vararg{Real}}},
    minrelativeintensity::Real,
    minioncount::Integer
)
    reference_key, spectrum = alkane_ladder_reference_spectrum(standard, ladderstep)
    ref_mzs = mzvalues(spectrum)
    ref_intensities = intensities(spectrum)
    selected = alkane_ladder_reference_threshold_apex_ion_candidates(
        msm,
        ref_mzs,
        ref_intensities,
        excludemzvalues,
        minrelativeintensity
    )
    length(selected) ≥ minioncount || throw(ArgumentError(
        "fewer than $(minioncount) reference apex ions are available for C$(ladderstep)"))

    AlkaneLadderApexIonSelection(
        :reference_relative_intensity_threshold,
        alkane_ladder_mzvalue_vector(
            excludemzvalues,
            true
        ),
        [candidate.mzvalue for candidate in selected],
        ref_mzs,
        [
            candidate.relativeintensity for candidate in selected
        ],
        reference_key,
        nothing,
        NaN,
        nothing
    )
end

function alkane_ladder_scan_order_extreme_apex_mzvalues(
    msm::MassScanMatrix,
    standard::AlkaneStandard,
    ladderstep::Integer,
    extremeioncount::Integer
)
    extremeioncount ≥ 1 || throw(ArgumentError(
        "scan-order extremeioncount must be a positive integer"))

    reference_key, spectrum = alkane_ladder_reference_spectrum(standard, ladderstep)
    ref_mzs = mzvalues(spectrum)
    ref_intensities = intensities(spectrum)
    candidates = alkane_ladder_reference_all_present_ion_candidates(
        msm,
        ref_mzs,
        ref_intensities
    )
    selected, excluded_mz = alkane_ladder_scan_order_extreme_ion_candidates(
        candidates,
        ladderstep,
        extremeioncount
    )
    length(selected) ≥ 2 || throw(ArgumentError(
        "fewer than two scan-order contrast ions are present on the m/z grid " *
        "for C$(ladderstep)"))

    AlkaneLadderApexIonSelection(
        :reference_alkane_series_extreme_grid_ions,
        isfinite(excluded_mz) ? [excluded_mz] : Float64[],
        [candidate.mzvalue for candidate in selected],
        ref_mzs,
        [
            candidate.relativeintensity for candidate in selected
        ],
        reference_key,
        nothing,
        NaN,
        extremeioncount
    )
end

function alkane_ladder_scan_order_all_usable_apex_mzvalues(
    msm::MassScanMatrix,
    standard::AlkaneStandard,
    ladderstep::Integer,
    excludemzvalues::Union{AbstractVector{<:Real}, Tuple{Vararg{Real}}}
)
    reference_key, spectrum = alkane_ladder_reference_spectrum(standard, ladderstep)
    ref_mzs = mzvalues(spectrum)
    ref_intensities = intensities(spectrum)
    candidates = alkane_ladder_reference_threshold_apex_ion_candidates(
        msm,
        ref_mzs,
        ref_intensities,
        excludemzvalues,
        0.0
    )
    length(candidates) ≥ 2 || throw(ArgumentError(
        "fewer than two usable reference ions are present on the m/z grid " *
        "for C$(ladderstep)"))

    AlkaneLadderApexIonSelection(
        :reference_all_usable_grid_ions,
        alkane_ladder_mzvalue_vector(excludemzvalues, true),
        [candidate.mzvalue for candidate in candidates],
        ref_mzs,
        [candidate.relativeintensity for candidate in candidates],
        reference_key,
        nothing,
        sum(candidate.relativeintensity for candidate in candidates),
        nothing
    )
end

function alkane_ladder_scan_order_extreme_ion_candidates(
    candidates::AbstractVector{<:AlkaneLadderApexIonCandidate},
    ladderstep::Integer,
    extremeioncount::Integer
)
    molecular_mz = alkane_ladder_nominal_molecular_ion_mz(ladderstep)
    excluded_mz = alkane_ladder_nominal_molecular_minus_14_series_mz(ladderstep)

    eligible = AlkaneLadderApexIonCandidate[]
    for candidate in candidates
        nominal_mz = round(Int, candidate.mzvalue)
        abs(candidate.mzvalue - nominal_mz) ≤ 0.5 || continue
        is_fragment = alkane_ladder_is_alkane_fragment_series_mz(nominal_mz)
        is_molecular = nominal_mz == molecular_mz
        (is_fragment || is_molecular) || continue
        nominal_mz == excluded_mz && continue
        push!(eligible, candidate)
    end
    sort!(eligible; by=candidate -> candidate.mzvalue)

    selected = AlkaneLadderApexIonCandidate[]
    for pairindex in 1:extremeioncount
        leftindex = pairindex
        rightindex = length(eligible) - pairindex + 1
        leftindex < rightindex || break

        push!(selected, eligible[leftindex])
        push!(selected, eligible[rightindex])
    end
    sort!(selected; by=candidate -> candidate.mzvalue)

    selected, excluded_mz
end

alkane_ladder_nominal_molecular_ion_mz(ladderstep::Integer) = 14 * ladderstep + 2

alkane_ladder_nominal_molecular_minus_14_series_mz(ladderstep::Integer) = 
    1 + 14 * (ladderstep - 1)

alkane_ladder_is_alkane_fragment_series_mz(mz::Integer) = mz ≥ 1 && mod(mz - 1, 14) == 0

function alkane_ladder_explicit_apex_mzvalue_selection(
    msm::MassScanMatrix,
    mzvalues::Union{AbstractVector{<:Real}, Tuple{Vararg{Real}}},
    minioncount::Integer
)
    selected_mzindices = alkane_ladder_apex_mz_indices(msm, mzvalues, minioncount)
    AlkaneLadderApexIonSelection(
        :explicit_mzvalues,
        Float64[],
        rawmzvalues(msm)[selected_mzindices],
        Float64[],
        Float64[],
        nothing,
        nothing,
        NaN,
        nothing
    )
end

function alkane_ladder_reference_threshold_apex_ion_candidates(
    msm::MassScanMatrix,
    ref_mzs::AbstractVector{<:Real},
    ref_intensities::AbstractVector{<:Real},
    excludemzvalues::Union{AbstractVector{<:Real}, Tuple{Vararg{Real}}},
    minrelativeintensity::Real
)
    length(ref_mzs) == length(ref_intensities) || throw(DimensionMismatch(
        "reference spectrum m/z and intensity vectors must have the same length"))
    isempty(ref_mzs) && throw(ArgumentError("reference spectrum must not be empty"))
    max_intensity = maximum(abs, ref_intensities)
    max_intensity > 0 || throw(ArgumentError(
        "reference spectrum must contain at least one nonzero intensity"))

    grid = alkane_mz_bins(msm)
    excluded = alkane_ladder_mzvalue_vector(excludemzvalues, true)
    intensity_by_grid_index = Dict{Int, Float64}()
    for (mz, intensity) in zip(ref_mzs, ref_intensities)
        alkane_ladder_mz_is_excluded(mz, excluded) && continue
        grid_index = alkane_ladder_nearest_grid_mz_index(grid, mz)
        isnothing(grid_index) && continue
        intensity_by_grid_index[grid_index] =
            get(intensity_by_grid_index, grid_index, 0.0) +
            abs(intensity)
    end

    candidates = [
        AlkaneLadderApexIonCandidate(
            grid[grid_index],
            grid_index,
            intensity / max_intensity
        )
        for (grid_index, intensity) in intensity_by_grid_index
        if intensity / max_intensity ≥ minrelativeintensity
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
    max_intensity = maximum(abs, ref_intensities)
    max_intensity > 0 || throw(ArgumentError(
        "reference spectrum must contain at least one nonzero intensity"))

    grid = alkane_mz_bins(msm)
    intensity_by_grid_index = Dict{Int, Float64}()
    for (mz, intensity) in zip(ref_mzs, ref_intensities)
        intensity_value = abs(intensity)
        intensity_value > 0 || continue
        grid_index = alkane_ladder_nearest_grid_mz_index(grid, mz)
        isnothing(grid_index) && continue
        intensity_by_grid_index[grid_index] =
            get(intensity_by_grid_index, grid_index, 0.0) + intensity_value
    end

    candidates = [
        AlkaneLadderApexIonCandidate(
            grid[grid_index],
            grid_index,
            intensity / max_intensity
        )
        for (grid_index, intensity) in intensity_by_grid_index
    ]
    sort!(candidates; by=candidate -> candidate.mzvalue)

    candidates
end

function alkane_ladder_scan_order_shape_ion_selection(
    msm::MassScanMatrix,
    standard::AlkaneStandard,
    ladderstep::Integer,
    ioncount::Integer,
    mzspacing::Integer
)
    ioncount ≥ 2 || throw(ArgumentError("shape ioncount must be at least 2"))
    mzspacing ≥ 1 || throw(ArgumentError("shape mzspacing must be at least 1"))
    reference_key, spectrum = alkane_ladder_reference_spectrum(standard, ladderstep)
    ref_mzs = mzvalues(spectrum)
    ref_intensities = intensities(spectrum)
    candidates = alkane_ladder_reference_all_present_ion_candidates(
        msm,
        ref_mzs,
        ref_intensities
    )
    length(candidates) ≥ ioncount || throw(ArgumentError(
        "fewer than $(ioncount) nonzero reference ions are present on the m/z grid " *
        "for C$(ladderstep)"))
    best, score = alkane_ladder_best_spaced_shape_ion_candidates(
        candidates,
        ioncount,
        mzspacing
    )

    AlkaneLadderApexIonSelection(
        :reference_best_spaced_shape_ions,
        Float64[],
        [candidate.mzvalue for candidate in best],
        ref_mzs,
        [candidate.relativeintensity for candidate in best],
        reference_key,
        mzspacing,
        score,
        nothing
    )
end

function alkane_ladder_best_spaced_shape_ion_candidates(
    candidates::AbstractVector{<:AlkaneLadderApexIonCandidate},
    ioncount::Integer,
    mzspacing::Integer
)
    by_nominal_mz = Dict{Int, AlkaneLadderApexIonCandidate}()
    for candidate in candidates
        nominal_mz = round(Int, candidate.mzvalue)
        abs(candidate.mzvalue - nominal_mz) ≤ 0.5 || continue
        previous = get(by_nominal_mz, nominal_mz, nothing)
        if isnothing(previous) ||
                candidate.relativeintensity > previous.relativeintensity
            by_nominal_mz[nominal_mz] = candidate
        end
    end

    best = AlkaneLadderApexIonCandidate[]
    best_score = -Inf
    for start_mz in sort(collect(keys(by_nominal_mz)))
        selected = AlkaneLadderApexIonCandidate[]
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

function alkane_ladder_reference_spectrum(
    standard::AlkaneStandard,
    ladderstep::Integer
)
    spectra = alkane_standard_spectra(standard)
    for (index, spectrum) in pairs(spectra)
        spectrumattrs = attrs(spectrum)
        key = if :order in keys(spectrumattrs)
            spectrumattrs.order
        else
            index
        end
        key == ladderstep && return key, spectrum
    end
    throw(ArgumentError(
        "standard does not contain a reference spectrum for C$(ladderstep)"))
end

function alkane_ladder_mzvalue_vector(
    values::Union{AbstractVector{<:Real}, Tuple{Vararg{Real}}},
    allowempty::Bool
)
    targets = collect(values)
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
    mzvalues::Union{AbstractVector{<:Real}, Tuple{Vararg{Real}}},
    minioncount::Integer
)
    targets = alkane_ladder_mzvalue_vector(mzvalues, false)
    grid = alkane_mz_bins(msm)
    selected = Int[]
    for target in targets
        index = alkane_ladder_nearest_grid_mz_index(grid, target)
        if !isnothing(index) && !(index in selected)
            push!(selected, index)
        end
    end
    length(selected) ≥ minioncount || throw(ArgumentError(
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
        distance = abs(grid[index] - target_mz)
        if distance < best_distance
            best_distance = distance
            best_index = index
        end
    end

    best_distance ≤ 0.5 ? best_index : nothing
end

function alkane_ladder_mz_is_excluded(
    mz::Real,
    excluded_mzvalues::AbstractVector{<:Real}
)
    any(excluded -> abs(mz - excluded) ≤ 0.5, excluded_mzvalues)
end

function alkane_ladder_default_mzretention_kwargs(
    msm::MassScanMatrix,
    order::Symbol
)
    order in (:inferdirection, :ascending, :descending, :simultaneous) ||
        throw(ArgumentError(
            "order must be :inferdirection, :ascending, :descending, or :simultaneous"))
    retention_values = rawretentions(msm)
    length(retention_values) ≥ 2 || throw(ArgumentError(
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
    scan_retentions::AbstractVector{<:Real},
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
    scan_retentions::AbstractVector{<:Real},
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
    totalscans ≥ 1 || throw(ArgumentError("totalscans must be positive"))
    nscans ≥ 1 || throw(ArgumentError("nscans must be positive"))
    localmargin ≥ 0 || throw(ArgumentError("localmargin must be nonnegative"))
    isfinite(centerscanindex) || throw(ArgumentError("centerscanindex must be finite"))

    center = clamp(round(Int, centerscanindex), 1, totalscans)
    left = max(1, center - localmargin)
    right = min(totalscans, center + localmargin)

    needed = min(nscans, totalscans)
    while right - left + 1 < needed
        if left > 1
            left -= 1
        else
            right += 1
        end
    end

    left:right
end

function alkane_ladder_observation_retention(
    scan_retentions::AbstractVector{<:Real},
    scanindex::Integer,
    mzindex::Integer,
    mzretentionkwargs::NamedTuple
)
    kwargs = merge(mzretentionkwargs, (mzindex=mzindex,))
    retention = Core.kwcall(kwargs, mzretention, scan_retentions[scanindex])
    isfinite(retention) || throw(ArgumentError("retention must be finite"))

    retention
end

function alkane_ladder_log_weight(
    adjusted_intensity::Real,
    variance::Real,
    variancefloor::Real
)
    v = max(variance, variancefloor)
    isfinite(v) && v > 0 || return 0.0

    abs2(adjusted_intensity) / v
end

function alkane_ladder_max_center_passes(maxapexshiftfromguess::Real)
    max(1, ceil(Int, maxapexshiftfromguess) + 1)
end

function alkane_ladder_best_apex_attempt_index(
    attempts::AbstractVector{<:AlkaneLadderIonApexAttempt},
    centeredscantolerance::Real
)
    isempty(attempts) && throw(ArgumentError("at least one apex attempt is required"))
    centered = findfirst(
        apex -> apex.success &&
            alkane_ladder_apex_is_centered(apex, centeredscantolerance),
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

function alkane_ladder_has_centered_success(
    attempts::AbstractVector{<:AlkaneLadderIonApexAttempt},
    centeredscantolerance::Real
)
    any(
        apex -> apex.success &&
            alkane_ladder_apex_is_centered(apex, centeredscantolerance),
        attempts
    )
end

function alkane_ladder_center_already_attempted(
    attempts::AbstractVector{<:AlkaneLadderIonApexAttempt},
    center_scan_index::Real,
    tolerance::Real
)
    any(
        apex -> isfinite(apex.fit_center_scan_index) &&
            abs(apex.fit_center_scan_index - center_scan_index) ≤ tolerance,
        attempts
    )
end

function alkane_ladder_apex_is_centered(
    apex::AlkaneLadderIonApexAttempt,
    tolerance::Real
)
    isfinite(tolerance) && tolerance ≥ 0 || throw(ArgumentError(
        "tolerance must be finite and nonnegative"))
    offset = apex.apex_offset_from_fit_center_scans

    isfinite(offset) && abs(offset) ≤ tolerance
end

function alkane_ladder_next_center_retention(
    apex::AlkaneLadderIonApexAttempt,
    raw_scan_retentions::AbstractVector{<:Real},
    input_scanindex::Integer,
    maxapexshiftfromguess::Real,
    centeredscantolerance::Real
)
    isfinite(apex.apex_scan_index) || return nothing
    lower_scan = max(1, input_scanindex - maxapexshiftfromguess)
    upper_scan = min(
        length(raw_scan_retentions),
        input_scanindex + maxapexshiftfromguess
    )
    target_scan = clamp(apex.apex_scan_index, lower_scan, upper_scan)
    abs(target_scan - apex.fit_center_scan_index) < centeredscantolerance &&
        return nothing

    alkane_ladder_retention_at_fractional_scan_index(raw_scan_retentions, target_scan)
end

function alkane_ladder_center_search_retentions(
    raw_scan_retentions::AbstractVector{<:Real},
    input_scanindex::Integer,
    maxapexshiftfromguess::Real,
    centeredscantolerance::Real
)
    maxshift = maxapexshiftfromguess
    maxshift == 0 && return Float64[]
    offsets = Float64[]
    whole_scan_shift = floor(Int, maxshift)
    for shift in 1:whole_scan_shift
        push!(offsets, shift)
        push!(offsets, -shift)
    end
    if maxshift - whole_scan_shift > 1e-9
        push!(offsets, maxshift)
        push!(offsets, -maxshift)
    end

    center_scan_indices = Float64[]
    for offset in offsets
        center_scan = clamp(
            input_scanindex + offset,
            1,
            length(raw_scan_retentions)
        )
        any(
            existing -> abs(existing - center_scan) ≤ centeredscantolerance,
            center_scan_indices
        ) && continue
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

function alkane_validate_retention_axis(retentions::AbstractVector{<:Real})
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
    length(retentions) == 1 && return 1.0
    increasing = last(retentions) ≥ first(retentions)
    increasing && return alkane_ladder_fractional_scan_index_increasing(
        retentions,
        retention
    )
    reversed_index = alkane_ladder_fractional_scan_index_increasing(
        reverse(retentions),
        retention
    )

    length(retentions) - reversed_index + 1
end

function alkane_ladder_fractional_scan_index_increasing(
    retentions::AbstractVector{<:Real},
    retention::Real
)
    retention ≤ first(retentions) && return 1.0
    retention ≥ last(retentions) && return Float64(length(retentions))
    right = searchsortedfirst(retentions, retention)
    left = right - 1
    left_retention = retentions[left]
    right_retention = retentions[right]
    right_retention == left_retention && return Float64(left)

    left + (retention - left_retention) / (right_retention - left_retention)
end

alkane_fractional_scan_index(
    retentions::AbstractVector{<:Real},
    retention::Real
) = alkane_ladder_fractional_scan_index(retentions, retention)

function alkane_ladder_retention_at_fractional_scan_index(
    retentions::AbstractVector{<:Real},
    scanindex::Real
)
    isempty(retentions) && throw(ArgumentError("scan retentions must not be empty"))
    all(isfinite, retentions) || throw(ArgumentError("scan retentions must be finite"))
    clamped_scanindex = clamp(scanindex, 1, length(retentions))
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
    nscans ≥ 1 || throw(ArgumentError("nscans must be at least 1"))
    isempty(retentions) && throw(ArgumentError("scan retentions must not be empty"))
    all(isfinite, retentions) || throw(ArgumentError("scan retentions must be finite"))
    isfinite(targetretention) || throw(ArgumentError(
        "targetretention must be finite"))
    requested = min(nscans, length(retentions))
    candidates = [
        (abs(retentions[index] - targetretention), index)
        for index in eachindex(retentions)
    ]
    sort!(candidates; by=candidate -> (candidate[1], candidate[2]))
    selected = [candidate[2] for candidate in candidates[1:requested]]
    sort!(selected)

    selected
end
