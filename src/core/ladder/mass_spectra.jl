const ALKANE_SERIES_CHECKSUM_VERSION = 1

function alkane_series_datainfo(
    rawmsm::MassScanMatrix,
    signalmsm::MassScanMatrix,
    variances::AbstractMatrix{<:Real}
)
    validate_alkane_series_variances(signalmsm, variances)
    (
        checksumalgorithm=:sha256,
        checksumversion=ALKANE_SERIES_CHECKSUM_VERSION,
        rawmsmchecksum=alkane_msm_checksum(rawmsm),
        signalmsmchecksum=alkane_msm_checksum(signalmsm),
        variancechecksum=alkane_variance_checksum(variances),
        rawscancount=scancount(rawmsm),
        rawmzcount=mzcount(rawmsm),
        signalscancount=scancount(signalmsm),
        signalmzcount=mzcount(signalmsm)
    )
end

function alkane_msm_checksum(msm::MassScanMatrix)
    ctx = SHA256_CTX()
    alkane_checksum_update_string!(ctx, "JuChrom.AlkaneSeries.MassScanMatrix")
    alkane_checksum_update_string!(ctx, string(ALKANE_SERIES_CHECKSUM_VERSION))
    alkane_checksum_update_string!(ctx, string(retentionunit(msm)))
    alkane_checksum_update_string!(ctx, string(mzunit(msm)))
    alkane_checksum_update_string!(ctx, string(intensityunit(msm)))
    alkane_checksum_update_string!(ctx, string(level(msm)))
    alkane_checksum_update_real_array!(ctx, "retentions", rawretentions(msm))
    alkane_checksum_update_real_array!(ctx, "mzvalues", rawmzvalues(msm))
    alkane_checksum_update_real_array!(ctx, "intensities", rawintensities(msm))

    bytes2hex(digest!(ctx))
end

function alkane_variance_checksum(variances::AbstractMatrix{<:Real})
    ctx = SHA256_CTX()
    alkane_checksum_update_string!(ctx, "JuChrom.AlkaneSeries.Variances")
    alkane_checksum_update_string!(ctx, string(ALKANE_SERIES_CHECKSUM_VERSION))
    alkane_checksum_update_real_array!(ctx, "variances", variances)

    bytes2hex(digest!(ctx))
end

function alkane_checksum_update_string!(ctx, value::AbstractString)
    bytes = collect(codeunits(value))
    alkane_checksum_update_length!(ctx, length(bytes))
    update!(ctx, bytes)
    nothing
end

function alkane_checksum_update_length!(ctx, lengthvalue::Integer)
    update!(ctx, reinterpret(UInt8, [hton(UInt64(lengthvalue))]))
    nothing
end

function alkane_checksum_update_real_array!(ctx, label::AbstractString, values)
    alkane_checksum_update_string!(ctx, label)
    alkane_checksum_update_string!(ctx, string(ndims(values)))
    for dim in size(values)
        alkane_checksum_update_length!(ctx, dim)
    end
    numeric_values = Float64.(vec(values))
    canonical_bits = hton.(reinterpret(UInt64, numeric_values))
    update!(ctx, reinterpret(UInt8, canonical_bits))
    nothing
end

"""
    alkaneladdermassspectra(msm, result; nonnegative=true, acceptedonly=true)

Extract peak-model mass spectra for refined alkane ladder steps.

The supplied `msm` must be the original matrix used for `findalkanes`. When baseline
subtraction was used, the analyzed signal is reconstructed as
`msm - result.baselineinfo.baselines` and verified against the stored SHA-256 checksum.
"""
function alkaneladdermassspectra(
    msm::MassScanMatrix,
    result::AlkaneSeriesResult;
    nonnegative::Bool=true,
    variancefloor::Real=1.0,
    threaded::Bool=true,
    acceptedonly::Bool=true,
    validatechecksum::Bool=true,
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true
)
    signal = alkane_ladder_extraction_signal(
        msm,
        result;
        validatechecksum=validatechecksum
    )
    validate_alkane_abundance_variancefloor(variancefloor)

    spectra = Dict{Int, AbstractMassSpectrum}()
    failures = Dict{Int, String}()
    steps = alkaneladdersteps(
        result;
        molecularion=molecularion,
        gapfilled=gapfilled,
        edgeextended=edgeextended
    )

    for step in steps
        if acceptedonly && !step.goodforcalibration
            failures[step.ladderstep] = "ladder step did not pass apex fit quality gate"
            continue
        end

        try
            spectra[step.ladderstep] = alkane_ladder_step_mass_spectrum(
                signal,
                step,
                result.variances;
                variancefloor=variancefloor,
                nonnegative=nonnegative,
                threaded=threaded
            )
        catch err
            failures[step.ladderstep] = sprint(showerror, err)
        end
    end

    (
        spectra=spectra,
        failures=failures,
        settings=(
            nonnegative=nonnegative,
            variancefloor=Float64(variancefloor),
            threaded=threaded,
            acceptedonly=acceptedonly,
            validatechecksum=validatechecksum
        )
    )
end

function alkaneladdermassspectrum(
    msm::MassScanMatrix,
    result::AlkaneSeriesResult,
    ladderstep::Integer;
    nonnegative::Bool=true,
    variancefloor::Real=1.0,
    threaded::Bool=true,
    acceptedonly::Bool=true,
    validatechecksum::Bool=true
)
    signal = alkane_ladder_extraction_signal(
        msm,
        result;
        validatechecksum=validatechecksum
    )
    validate_alkane_abundance_variancefloor(variancefloor)

    step = alkane_ladder_step_by_carbon(result, ladderstep)
    acceptedonly && !step.goodforcalibration && throw(ArgumentError(
        "ladder step $(Int(ladderstep)) did not pass apex fit quality gate"))

    alkane_ladder_step_mass_spectrum(
        signal,
        step,
        result.variances;
        variancefloor=variancefloor,
        nonnegative=nonnegative,
        threaded=threaded
    )
end

function alkane_ladder_extraction_signal(
    msm::MassScanMatrix,
    result::AlkaneSeriesResult;
    validatechecksum::Bool
)
    validatechecksum && alkane_validate_raw_msm_checksum(msm, result)

    signal = if isnothing(result.baselineinfo)
        msm
    else
        hasproperty(result.baselineinfo, :baselines) || throw(ArgumentError(
            "result.baselineinfo does not contain fitted baselines"))
        validate_alkane_series_baselines(msm, result.baselineinfo.baselines)
        msm - result.baselineinfo.baselines
    end

    validate_alkane_series_variances(signal, result.variances)
    validatechecksum && alkane_validate_signal_and_variance_checksums(signal, result)

    signal
end

function alkane_validate_raw_msm_checksum(msm::MassScanMatrix, result::AlkaneSeriesResult)
    datainfo = alkane_ladder_result_datainfo(result)
    checksum = alkane_msm_checksum(msm)
    checksum == datainfo.rawmsmchecksum || throw(ArgumentError(
        "msm does not match the raw data used for alkane ladder detection"))

    nothing
end

function alkane_validate_signal_and_variance_checksums(
    signal::MassScanMatrix,
    result::AlkaneSeriesResult
)
    datainfo = alkane_ladder_result_datainfo(result)
    signalchecksum = alkane_msm_checksum(signal)
    signalchecksum == datainfo.signalmsmchecksum || throw(ArgumentError(
        "reconstructed alkane ladder signal does not match the analyzed data"))
    variancechecksum = alkane_variance_checksum(result.variances)
    variancechecksum == datainfo.variancechecksum || throw(ArgumentError(
        "result variances do not match the variances used for alkane ladder detection"))

    nothing
end

function alkane_ladder_result_datainfo(result::AlkaneSeriesResult)
    isnothing(result.datainfo) && throw(ArgumentError(
        "alkane series result does not contain SHA-256 data checksums"))
    result.datainfo
end

function alkane_ladder_step_by_carbon(result::AlkaneSeriesResult, ladderstep::Integer)
    matches = [
        step for step in alkaneladdersteps(result)
        if step.ladderstep == Int(ladderstep)
    ]
    isempty(matches) && throw(ArgumentError(
        "alkane ladder result contains no refined step C$(Int(ladderstep))"))
    length(matches) == 1 || throw(ArgumentError(
        "alkane ladder result contains multiple refined steps for C$(Int(ladderstep))"))

    only(matches)
end

function alkane_ladder_step_mass_spectrum(
    msm::MassScanMatrix,
    step::AlkaneLadderStep,
    variances::AbstractMatrix{<:Real};
    variancefloor::Real,
    nonnegative::Bool,
    threaded::Bool
)
    validate_alkane_series_variances(msm, variances)
    apex = step.apex
    alkane_ladder_apex_is_usable_for_spectrum(apex) || throw(ArgumentError(
        "ladder step C$(step.ladderstep) does not contain a usable peak model apex"))

    fit = alkane_ladder_spectrum_apex_fit(apex)
    nscans = length(fit.scan_indices)
    nscans >= 3 || throw(ArgumentError(
        "at least three scans are needed for peak-model spectrum extraction"))

    X = rawintensities(msm)
    scan_retentions = retentions(msm)
    target_retention = alkane_ladder_spectrum_fit_center_retention(apex)
    mzkwargs = fit.mzretentionkwargs
    n_mz = mzcount(msm)

    spectrum_intensities = Vector{Float64}(undef, n_mz)
    standarderrors = fill(NaN, n_mz)
    zscores = fill(NaN, n_mz)
    fitsuccess = falses(n_mz)
    nobservations = zeros(Int, n_mz)

    if threaded
        @threads for mzindex in 1:n_mz
            alkane_ladder_fit_spectrum_ion!(
                spectrum_intensities,
                standarderrors,
                zscores,
                fitsuccess,
                nobservations,
                X,
                variances,
                msm,
                apex,
                scan_retentions,
                target_retention,
                mzkwargs,
                mzindex,
                nscans,
                variancefloor,
                nonnegative
            )
        end
    else
        for mzindex in 1:n_mz
            alkane_ladder_fit_spectrum_ion!(
                spectrum_intensities,
                standarderrors,
                zscores,
                fitsuccess,
                nobservations,
                X,
                variances,
                msm,
                apex,
                scan_retentions,
                target_retention,
                mzkwargs,
                mzindex,
                nscans,
                variancefloor,
                nonnegative
            )
        end
    end

    MassSpectrum(
        collect(rawmzvalues(msm)),
        mzunit(msm),
        spectrum_intensities,
        intensityunit(msm);
        attrs=(
            source=:alkane_ladder_mass_spectrum,
            model=:fixed_peak_model_wls,
            ladderstep=step.ladderstep,
            step_source=step.source,
            apex_retention=apex.apexretention,
            apex_scan_index=apex.apexscanindex,
            input_scan_index=fit.input_scan_index,
            input_retention=fit.input_retention,
            fit_center_scan_index=alkane_ladder_spectrum_fit_center_scan_index(apex),
            fit_center_retention=target_retention,
            scan_indices=collect(fit.scan_indices),
            scanwindow=(nscans - 1) / 2,
            variance_weighted=true,
            variancefloor=Float64(variancefloor),
            nonnegative=nonnegative,
            allownegative=!nonnegative,
            threaded=threaded,
            mzretentionkwargs=mzkwargs,
            peak_model=collect(fit.peak_model),
            standard_errors=standarderrors,
            z_scores=zscores,
            fit_success=fitsuccess,
            n_observations=nobservations,
            mass_spectrum_cosine=step.massspectrumcosine,
            required_cosine=step.requiredcosine,
            good_for_calibration=step.goodforcalibration
        )
    )
end

function alkane_ladder_fit_spectrum_ion!(
    spectrum_intensities::AbstractVector{Float64},
    standarderrors::AbstractVector{Float64},
    zscores::AbstractVector{Float64},
    fitsuccess::AbstractVector{Bool},
    nobservations::AbstractVector{Int},
    X::AbstractMatrix{<:Real},
    variances::AbstractMatrix{<:Real},
    msm::MassScanMatrix,
    apex,
    scan_retentions,
    target_retention::Real,
    mzretentionkwargs::NamedTuple,
    mzindex::Integer,
    nscans::Integer,
    variancefloor::Real,
    nonnegative::Bool
)
    selected_scanindices = alkane_ladder_scan_indices_for_mz(
        msm,
        scan_retentions,
        target_retention,
        nscans,
        mzindex,
        mzretentionkwargs,
        alkane_ladder_spectrum_fit_center_scan_index(apex),
        max(3, fld(Int(nscans), 2) + 2)
    )

    numerator = 0.0
    denominator = 0.0
    observations = 0
    for selected_scanindex in selected_scanindices
        observation_retention = alkane_ladder_observation_retention(
            scan_retentions,
            selected_scanindex,
            mzindex,
            mzretentionkwargs
        )
        model_value = alkane_ladder_spectrum_shape_at_retention(
            apex,
            observation_retention
        )
        y = Float64(X[selected_scanindex, mzindex])
        variance = max(Float64(variances[selected_scanindex, mzindex]), variancefloor)
        isfinite(model_value) && isfinite(y) && isfinite(variance) && variance > 0 ||
            continue

        weight = inv(variance)
        numerator += weight * model_value * y
        denominator += weight * abs2(model_value)
        observations += 1
    end

    raw_height = denominator > 0 ? numerator / denominator : NaN
    height = nonnegative ? max(raw_height, 0.0) : raw_height
    spectrum_intensities[mzindex] = isfinite(height) ? height : 0.0
    fitsuccess[mzindex] = isfinite(raw_height) && denominator > 0
    nobservations[mzindex] = observations
    if denominator > 0
        standarderror = sqrt(inv(denominator))
        standarderrors[mzindex] = standarderror
        zscores[mzindex] = standarderror > 0 ? raw_height / standarderror : NaN
    end

    nothing
end

function alkane_ladder_apex_is_usable_for_spectrum(apex)
    apex.success || return false
    fitresult = apex.fit
    fitresult isa AlkaneLadderIonApexResult || return false
    fit = fitresult.fit

    isfinite(apex.apexretention) &&
        isfinite(apex.apexscanindex) &&
        isfinite(fit.apex_x) &&
        isfinite(fit.x_scale) &&
        fit.x_scale > 0 &&
        isfinite(fit.beta) &&
        isfinite(fit.gamma) &&
        length(fit.scan_indices) == length(fit.peak_model) &&
        all(isfinite, fit.peak_model)
end

function alkane_ladder_spectrum_fit_center_retention(apex)
    alkane_ladder_spectrum_apex_fit(apex).fit_center_retention
end

function alkane_ladder_spectrum_fit_center_scan_index(apex)
    alkane_ladder_spectrum_apex_fit(apex).fit_center_scan_index
end

function alkane_ladder_spectrum_apex_fit(apex)
    fitresult = apex.fit
    fitresult isa AlkaneLadderIonApexResult || throw(ArgumentError(
        "ladder apex does not contain an ion apex fit result"))

    fitresult.fit
end

function alkane_ladder_spectrum_shape_at_retention(apex, observation_retention::Real)
    fit = alkane_ladder_spectrum_apex_fit(apex)
    x = (Float64(observation_retention) - alkane_ladder_spectrum_fit_center_retention(apex)) /
        fit.x_scale

    exp(
        fit.beta * (x - fit.apex_x) +
        fit.gamma * (abs2(x) - abs2(fit.apex_x))
    )
end
