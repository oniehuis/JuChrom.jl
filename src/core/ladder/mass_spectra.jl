const ALKANE_SERIES_CHECKSUM_VERSION = 1

struct AlkaneSeriesDataInfo <: AbstractAlkaneSeriesDataInfo
    checksumalgorithm::Symbol
    checksumversion::Int
    rawmsmchecksum::String
    signalmsmchecksum::String
    variancechecksum::String
    rawscancount::Int
    rawmzcount::Int
    signalscancount::Int
    signalmzcount::Int
end

struct AlkaneLadderMassSpectrumSettings{T<:Real}
    nonnegative::Bool
    variancefloor::T
    threaded::Bool
    acceptedonly::Bool
    validatechecksum::Bool
    molecularion::Bool
    gapfilled::Bool
    edgeextended::Bool
end

struct AlkaneLadderMassSpectrumExtraction{
    T1<:AbstractDict{Int, <:AbstractMassSpectrum},
    T2<:AbstractDict{Int, <:AbstractString},
    T3<:AlkaneLadderMassSpectrumSettings
}
    spectra::T1
    failures::T2
    settings::T3
end

struct AlkaneLadderMassSpectrumStepExtraction
    ladderstep::Int
    spectrum::Union{Nothing, AbstractMassSpectrum}
    failure::Union{Nothing, String}
end

function Base.show(io::IO, extraction::AlkaneLadderMassSpectrumExtraction)
    print(io, "AlkaneLadderMassSpectrumExtraction(")
    print(io, "spectra=", alkane_ladder_mass_spectrum_count_string(extraction.spectra))
    print(io, ", failures=", alkane_ladder_mass_spectrum_count_string(extraction.failures))
    print(io, ", settings=", alkane_ladder_mass_spectrum_settings_string(extraction.settings))
    print(io, ")")
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    extraction::AlkaneLadderMassSpectrumExtraction
)
    println(io, "AlkaneLadderMassSpectrumExtraction")
    println(io, "  spectra: ", alkane_ladder_mass_spectrum_count_string(extraction.spectra))
    println(io, "  failures: ", alkane_ladder_mass_spectrum_count_string(extraction.failures))
    print(io, "  settings: ",
        alkane_ladder_mass_spectrum_settings_string(extraction.settings))
end

function Base.getindex(
    extraction::AlkaneLadderMassSpectrumExtraction,
    carbon::Integer
)
    extraction.spectra[Int(carbon)]
end

function Base.getindex(
    extraction::AlkaneLadderMassSpectrumExtraction,
    carbons::AbstractVector{<:Integer}
)
    [extraction[carbon] for carbon in carbons]
end

function Base.get(
    extraction::AlkaneLadderMassSpectrumExtraction,
    carbon::Integer,
    default
)
    get(extraction.spectra, Int(carbon), default)
end

function Base.haskey(
    extraction::AlkaneLadderMassSpectrumExtraction,
    carbon::Integer
)
    haskey(extraction.spectra, Int(carbon))
end

function Base.keys(extraction::AlkaneLadderMassSpectrumExtraction)
    sort!(collect(keys(extraction.spectra)))
end

function Base.values(extraction::AlkaneLadderMassSpectrumExtraction)
    [extraction.spectra[carbon] for carbon in keys(extraction)]
end

function Base.pairs(extraction::AlkaneLadderMassSpectrumExtraction)
    [carbon => extraction.spectra[carbon] for carbon in keys(extraction)]
end

Base.length(extraction::AlkaneLadderMassSpectrumExtraction) = length(extraction.spectra)
Base.isempty(extraction::AlkaneLadderMassSpectrumExtraction) = isempty(extraction.spectra)

Base.IteratorSize(::Type{<:AlkaneLadderMassSpectrumExtraction}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AlkaneLadderMassSpectrumExtraction}) = Base.HasEltype()
Base.eltype(::Type{<:AlkaneLadderMassSpectrumExtraction{T}}) where {T} = valtype(T)

function Base.iterate(extraction::AlkaneLadderMassSpectrumExtraction)
    carbons = keys(extraction)
    isempty(carbons) && return nothing

    extraction.spectra[first(carbons)], (carbons, 2)
end

function Base.iterate(extraction::AlkaneLadderMassSpectrumExtraction, state)
    carbons, index = state
    index > length(carbons) && return nothing

    extraction.spectra[carbons[index]], (carbons, index + 1)
end

function alkane_ladder_mass_spectrum_count_string(items::AbstractDict)
    carbons = sort!(collect(keys(items)))
    count = length(carbons)
    noun = count == 1 ? "C spectrum" : "C spectra"

    "$count $noun ($(alkane_ladder_mass_spectrum_carbon_string(carbons)))"
end

function alkane_ladder_mass_spectrum_carbon_string(carbons::AbstractVector{<:Integer})
    isempty(carbons) && return "none"

    ranges = UnitRange{Int}[]
    range_start = Int(first(carbons))
    previous = range_start
    for carbon in Iterators.drop(carbons, 1)
        carbon = Int(carbon)
        if carbon == previous + 1
            previous = carbon
        else
            push!(ranges, range_start:previous)
            range_start = carbon
            previous = carbon
        end
    end
    push!(ranges, range_start:previous)

    join(alkane_ladder_mass_spectrum_carbon_range_string.(ranges), ", ")
end

function alkane_ladder_mass_spectrum_carbon_range_string(range::UnitRange{Int})
    first(range) == last(range) && return "C$(first(range))"

    "C$(first(range))-C$(last(range))"
end

function alkane_ladder_mass_spectrum_settings_string(
    settings::AlkaneLadderMassSpectrumSettings
)
    string(
        "nonnegative=", settings.nonnegative,
        ", variancefloor=", settings.variancefloor,
        ", acceptedonly=", settings.acceptedonly,
        ", threaded=", settings.threaded
    )
end

function alkane_series_datainfo(
    rawmsm::AbstractMassScanMatrix,
    signalmsm::AbstractMassScanMatrix,
    variances::AbstractMatrix{<:Real}
)
    validate_alkane_series_variances(signalmsm, variances)
    AlkaneSeriesDataInfo(
        :sha256,
        ALKANE_SERIES_CHECKSUM_VERSION,
        alkane_msm_checksum(rawmsm),
        alkane_msm_checksum(signalmsm),
        alkane_variance_checksum(variances),
        scancount(rawmsm),
        mzcount(rawmsm),
        scancount(signalmsm),
        mzcount(signalmsm)
    )
end

function alkane_msm_checksum(msm::AbstractMassScanMatrix)::String
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

function alkane_variance_checksum(variances::AbstractMatrix{<:Real})::String
    ctx = SHA256_CTX()
    alkane_checksum_update_string!(ctx, "JuChrom.AlkaneSeries.Variances")
    alkane_checksum_update_string!(ctx, string(ALKANE_SERIES_CHECKSUM_VERSION))
    alkane_checksum_update_real_array!(ctx, "variances", variances)

    bytes2hex(digest!(ctx))
end

function alkane_checksum_update_string!(ctx::SHA256_CTX, value::AbstractString)
    bytes = collect(codeunits(value))
    alkane_checksum_update_length!(ctx, length(bytes))
    update!(ctx, bytes)
    nothing
end

function alkane_checksum_update_length!(ctx::SHA256_CTX, lengthvalue::Integer)
    update!(ctx, reinterpret(UInt8, [hton(UInt64(lengthvalue))]))
    nothing
end

function alkane_checksum_update_real_array!(
    ctx::SHA256_CTX,
    label::AbstractString,
    values::AbstractArray{<:Real}
)
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

function alkane_ladder_mass_spectrum_settings(
    nonnegative::Bool,
    variancefloor::Real,
    threaded::Bool,
    acceptedonly::Bool,
    validatechecksum::Bool,
    molecularion::Bool,
    gapfilled::Bool,
    edgeextended::Bool
)
    validate_alkane_abundance_variancefloor(variancefloor)
    AlkaneLadderMassSpectrumSettings(
        nonnegative,
        variancefloor,
        threaded,
        acceptedonly,
        validatechecksum,
        molecularion,
        gapfilled,
        edgeextended
    )
end

"""
    alkaneladdermassspectra(msm, result; nonnegative=true, acceptedonly=true)

Extract peak-model mass spectra for refined alkane ladder steps.

The supplied `msm` must be the original matrix used for `findalkanes`. When baseline
subtraction was used, the analyzed signal is reconstructed as
`subtractbaseline(msm, result.baselineinfo.baselines)` and verified against the stored
SHA-256 checksum.
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
    settings = alkane_ladder_mass_spectrum_settings(
        nonnegative,
        variancefloor,
        threaded,
        acceptedonly,
        validatechecksum,
        molecularion,
        gapfilled,
        edgeextended
    )

    alkane_ladder_mass_spectra(msm, result, settings)
end

function alkane_ladder_mass_spectra(
    msm::MassScanMatrix,
    result::AlkaneSeriesResult,
    settings::AlkaneLadderMassSpectrumSettings
)
    signal = alkane_ladder_extraction_signal(msm, result, settings.validatechecksum)
    steps = alkane_ladder_mass_spectrum_steps(
        result,
        settings.molecularion,
        settings.gapfilled,
        settings.edgeextended
    )

    use_step_threads = settings.threaded && length(steps) > 1
    ionthreaded = settings.threaded && !use_step_threads
    stepresults = Vector{AlkaneLadderMassSpectrumStepExtraction}(undef, length(steps))

    if use_step_threads
        Base.Threads.@threads for index in eachindex(steps)
            stepresults[index] = alkane_ladder_extract_step_mass_spectrum(
                signal,
                steps[index],
                result.variances,
                settings,
                ionthreaded
            )
        end
    else
        for index in eachindex(steps)
            stepresults[index] = alkane_ladder_extract_step_mass_spectrum(
                signal,
                steps[index],
                result.variances,
                settings,
                ionthreaded
            )
        end
    end

    spectra = Dict{Int, AbstractMassSpectrum}()
    failures = Dict{Int, String}()
    for stepresult in stepresults
        if isnothing(stepresult.failure)
            spectra[stepresult.ladderstep] = stepresult.spectrum
        else
            failures[stepresult.ladderstep] = stepresult.failure
        end
    end

    AlkaneLadderMassSpectrumExtraction(
        spectra,
        failures,
        settings
    )
end

function alkane_ladder_mass_spectrum_steps(
    result::AlkaneSeriesResult,
    molecularion::Bool,
    gapfilled::Bool,
    edgeextended::Bool
)
    steps = AlkaneLadderStep[]
    if molecularion
        append!(steps, alkane_ladder_steps_from_apexes(
            result.apexinfo.apexes,
            :molecularion
        ))
    end

    if gapfilled
        append!(steps, alkane_ladder_steps_from_additions(
            result.additioninfo.gapfilled,
            :gapfilled
        ))
    end

    if edgeextended
        append!(steps, alkane_ladder_steps_from_additions(
            result.additioninfo.leftextended,
            :leftextended
        ))
        append!(steps, alkane_ladder_steps_from_additions(
            result.additioninfo.rightextended,
            :rightextended
        ))
    end

    sort!(steps; by=step -> step.ladderstep)

    steps
end

function alkane_ladder_extract_step_mass_spectrum(
    msm::MassScanMatrix,
    step::AlkaneLadderStep,
    variances::AbstractMatrix{<:Real},
    settings::AlkaneLadderMassSpectrumSettings,
    ionthreaded::Bool
)
    if settings.acceptedonly && !step.goodforcalibration
        return AlkaneLadderMassSpectrumStepExtraction(
            step.ladderstep,
            nothing,
            "ladder step did not pass apex fit quality gate"
        )
    end

    try
        spectrum = alkane_ladder_step_mass_spectrum(
            msm,
            step,
            variances,
            settings,
            ionthreaded
        )
        return AlkaneLadderMassSpectrumStepExtraction(
            step.ladderstep,
            spectrum,
            nothing
        )
    catch err
        return AlkaneLadderMassSpectrumStepExtraction(
            step.ladderstep,
            nothing,
            sprint(showerror, err)
        )
    end
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
    settings = alkane_ladder_mass_spectrum_settings(
        nonnegative,
        variancefloor,
        threaded,
        acceptedonly,
        validatechecksum,
        true,
        true,
        true
    )
    signal = alkane_ladder_extraction_signal(msm, result, settings.validatechecksum)

    step = alkane_ladder_step_by_carbon(result, ladderstep)
    settings.acceptedonly && !step.goodforcalibration && throw(ArgumentError(
        "ladder step $(ladderstep) did not pass apex fit quality gate"))

    alkane_ladder_step_mass_spectrum(
        signal,
        step,
        result.variances,
        settings,
        settings.threaded
    )
end

function alkane_ladder_extraction_signal(
    msm::MassScanMatrix,
    result::AlkaneSeriesResult,
    validatechecksum::Bool
)
    validatechecksum && alkane_validate_raw_msm_checksum(msm, result)

    signal = if isnothing(result.baselineinfo)
        msm
    else
        validate_alkane_series_baselines(msm, result.baselineinfo.baselines)
        _subtractbaseline(msm, result.baselineinfo.baselines)
    end

    validate_alkane_series_variances(signal, result.variances)
    validatechecksum && alkane_validate_signal_and_variance_checksums(signal, result)

    signal
end

function alkane_validate_raw_msm_checksum(
    msm::AbstractMassScanMatrix,
    result::AlkaneSeriesResult
)
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

function alkane_ladder_result_datainfo(result::AlkaneSeriesResult)::AlkaneSeriesDataInfo
    isnothing(result.datainfo) && throw(ArgumentError(
        "alkane series result does not contain SHA-256 data checksums"))
    result.datainfo isa AlkaneSeriesDataInfo || throw(ArgumentError(
        "alkane series result datainfo must be an AlkaneSeriesDataInfo"))
    result.datainfo.checksumalgorithm ≡ :sha256 || throw(ArgumentError(
        "alkane series result datainfo checksum algorithm must be :sha256"))
    result.datainfo.checksumversion == ALKANE_SERIES_CHECKSUM_VERSION || throw(
        ArgumentError(
            "alkane series result datainfo checksum version is not supported"))

    result.datainfo
end

function alkane_ladder_step_by_carbon(result::AlkaneSeriesResult, ladderstep::Integer)
    matches = [
        step for step in alkane_ladder_mass_spectrum_steps(result, true, true, true)
        if step.ladderstep == ladderstep
    ]
    isempty(matches) && throw(ArgumentError(
        "alkane ladder result contains no refined step C$(ladderstep)"))
    length(matches) == 1 || throw(ArgumentError(
        "alkane ladder result contains multiple refined steps for C$(ladderstep)"))

    only(matches)
end

function alkane_ladder_step_mass_spectrum(
    msm::MassScanMatrix,
    step::AlkaneLadderStep,
    variances::AbstractMatrix{<:Real},
    settings::AlkaneLadderMassSpectrumSettings,
    ionthreaded::Bool
)
    validate_alkane_series_variances(msm, variances)
    apex = step.apex
    alkane_ladder_apex_is_usable_for_spectrum(apex) || throw(ArgumentError(
        "ladder step C$(step.ladderstep) does not contain a usable peak model apex"))

    fit = alkane_ladder_spectrum_apex_fit(apex)
    nscans = length(fit.scan_indices)
    nscans ≥ 3 || throw(ArgumentError(
        "at least three scans are needed for peak-model spectrum extraction"))

    X = rawintensities(msm)
    scan_retentions = rawretentions(msm)
    target_retention = fit.fit_center_retention
    mzkwargs = fit.mzretentionkwargs
    n_mz = mzcount(msm)

    spectrum_intensities = Vector{Float64}(undef, n_mz)
    standarderrors = fill(NaN, n_mz)
    zscores = fill(NaN, n_mz)
    fitsuccess = falses(n_mz)
    nobservations = zeros(Int, n_mz)

    if ionthreaded
        Base.Threads.@threads for mzindex in 1:n_mz
            alkane_ladder_fit_spectrum_ion!(
                spectrum_intensities,
                standarderrors,
                zscores,
                fitsuccess,
                nobservations,
                X,
                variances,
                msm,
                fit,
                scan_retentions,
                target_retention,
                mzkwargs,
                mzindex,
                nscans,
                settings.variancefloor,
                settings.nonnegative
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
                fit,
                scan_retentions,
                target_retention,
                mzkwargs,
                mzindex,
                nscans,
                settings.variancefloor,
                settings.nonnegative
            )
        end
    end

    MassSpectrum(
        rawmzvalues(msm),
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
            fit_center_scan_index=fit.fit_center_scan_index,
            fit_center_retention=target_retention,
            scan_indices=fit.scan_indices,
            scanwindow=(nscans - 1) / 2,
            variance_weighted=true,
            variancefloor=settings.variancefloor,
            nonnegative=settings.nonnegative,
            allownegative=!settings.nonnegative,
            threaded=settings.threaded,
            ion_threaded=ionthreaded,
            mzretentionkwargs=mzkwargs,
            peak_model=fit.peak_model,
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
    fit::AlkaneLadderIonApexAttempt,
    scan_retentions::AbstractVector{<:Real},
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
        fit.fit_center_scan_index,
        max(3, fld(nscans, 2) + 2)
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
            fit,
            target_retention,
            observation_retention
        )
        y = X[selected_scanindex, mzindex]
        variance = max(variances[selected_scanindex, mzindex], variancefloor)
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

function alkane_ladder_apex_is_usable_for_spectrum(apex::AlkaneLadderApex)::Bool
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

function alkane_ladder_spectrum_apex_fit(
    apex::AlkaneLadderApex
)::AlkaneLadderIonApexAttempt
    fitresult = apex.fit
    fitresult isa AlkaneLadderIonApexResult || throw(ArgumentError(
        "ladder apex does not contain an ion apex fit result"))

    fitresult.fit
end

function alkane_ladder_spectrum_shape_at_retention(
    fit::AlkaneLadderIonApexAttempt,
    fit_center_retention::Real,
    observation_retention::Real
)
    x = (observation_retention - fit_center_retention) / fit.x_scale
    
    exp(fit.beta * (x - fit.apex_x) + fit.gamma * (abs2(x) - abs2(fit.apex_x)))
end
