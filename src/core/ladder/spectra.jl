"""
    ladder_step_mass_spectrum(msm, apexresult, variances; kwargs...) -> MassSpectrum

Extract a full-grid mass spectrum for one alkane ladder step.

Reference ions that were used to infer the apex model are reconstructed directly from the
shared log-quadratic fit. All other m/z channels are fit with the fixed exponential peak
shape from that apex model, estimating only one height per ion.
"""
function ladder_step_mass_spectrum(
    msm::MassScanMatrix,
    apexresult,
    variances::AbstractMatrix{<:Real};
    scanwindow::Integer=2,
    mzretentionkwargs=nothing,
    variancefloor::Real=1.0,
    allownegative::Bool=true,
    threaded::Bool=true,
)
    validate_alkane_series_variances(msm, variances)
    scanwindow >= 1 || throw(ArgumentError("scanwindow must be at least 1"))
    isfinite(variancefloor) && variancefloor > 0 || throw(ArgumentError(
        "variancefloor must be finite and positive"))

    apex = ladder_spectrum_apex(apexresult)
    ladder_spectrum_apex_is_usable(apex) || throw(ArgumentError(
        "apex result must contain a valid continuous peak model"))
    center_scanindex = Int(apex.input_scan_index)
    1 <= center_scanindex <= scancount(msm) || throw(ArgumentError(
        "apex center scan index $(center_scanindex) is outside 1:$(scancount(msm))"))

    mzkwargs = isnothing(mzretentionkwargs) ?
        apex.mzretentionkwargs :
        mzretentionkwargs
    raw_scan_retentions = Float64.(rawretentions(msm))
    scan_retentions = retentions(msm)
    nscans = apex_scan_window_observation_count(scancount(msm), center_scanindex, scanwindow)
    nscans >= 3 || throw(ArgumentError(
        "at least three scans are needed for ladder step spectrum extraction"))

    mzindices = collect(1:mzcount(msm))

    intensities = Vector{Float64}(undef, mzcount(msm))
    standarderrors = fill(NaN, mzcount(msm))
    zscores = fill(NaN, mzcount(msm))
    fitsuccess = falses(mzcount(msm))
    nobservations = zeros(Int, mzcount(msm))
    heightsources = fill(:fixed_shape_fit, mzcount(msm))

    referencecols = zeros(Int, mzcount(msm))
    for (col, mzindex) in pairs(apex.mz_indices)
        referencecols[Int(mzindex)] = col
    end
    for mzindex in mzindices
        referencecol = referencecols[mzindex]
        referencecol == 0 && continue
        height = ladder_reference_ion_height(apex, referencecol)
        height = allownegative ? height : max(height, 0.0)
        intensities[mzindex] = isfinite(height) ? height : 0.0
        fitsuccess[mzindex] = isfinite(intensities[mzindex])
        nobservations[mzindex] = size(apex.fit_intensities, 1)
        heightsources[mzindex] = :apex_model
    end

    if threaded && Threads.nthreads() > 1
        Threads.@threads for mzindex in mzindices
            fit_ladder_step_nonreference_ion!(
                intensities,
                standarderrors,
                zscores,
                fitsuccess,
                nobservations,
                msm,
                variances,
                apex,
                scan_retentions,
                raw_scan_retentions[center_scanindex],
                mzkwargs,
                referencecols,
                mzindex,
                nscans,
                variancefloor,
                allownegative,
            )
        end
    else
        for mzindex in mzindices
            fit_ladder_step_nonreference_ion!(
                intensities,
                standarderrors,
                zscores,
                fitsuccess,
                nobservations,
                msm,
                variances,
                apex,
                scan_retentions,
                raw_scan_retentions[center_scanindex],
                mzkwargs,
                referencecols,
                mzindex,
                nscans,
                variancefloor,
                allownegative,
            )
        end
    end

    MassSpectrum(
        collect(rawmzvalues(msm)),
        mzunit(msm),
        intensities,
        intensityunit(msm);
        attrs=(
            source=:alkane_ladder_step_spectrum,
            model=:fixed_apex_peak_shape,
            carbon=Int(apex.carbon),
            apex_retention=apex.apex_retention,
            apex_scan_index=apex.apex_scan_index,
            center_scan_index=center_scanindex,
            ladder_retention=hasproperty(apex, :ladder_retention) ? apex.ladder_retention : missing,
            ladder_scan_index=hasproperty(apex, :ladder_scan_index) ? apex.ladder_scan_index : missing,
            scanwindow=Int(scanwindow),
            variance_weighted=true,
            variancefloor=Float64(variancefloor),
            allownegative=allownegative,
            threaded=threaded,
            mzretentionkwargs=mzkwargs,
            reference_mz_indices=collect(Int.(apex.mz_indices)),
            reference_mz_values=collect(apex.mz_values),
            reference_intensities=collect(apex.reference_intensities),
            referenceattrs=apex.referenceattrs,
            beta=apex.beta,
            gamma=apex.gamma,
            apex_x=apex.apex_x,
            x_scale=apex.x_scale,
            log_floor=apex.log_floor,
            height_sources=heightsources,
            standard_errors=standarderrors,
            z_scores=zscores,
            fit_success=fitsuccess,
            n_observations=nobservations,
        ),
    )
end

function alkane_series_step_spectra(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    apexes;
    scanwindow::Integer=2,
    carbons=:all,
    mzretentionkwargs=nothing,
    variancefloor::Real=1.0,
    allownegative::Bool=true,
)
    validate_alkane_series_variances(msm, variances)
    selectedcarbons = normalize_alkane_step_spectrum_carbons(carbons)
    settings = (
        scanwindow=Int(scanwindow),
        carbons=isnothing(selectedcarbons) ? :all : sort!(collect(selectedcarbons)),
        variance_weighted=true,
        variancefloor=Float64(variancefloor),
        allownegative=allownegative,
        mzretentionkwargs=mzretentionkwargs,
    )

    apexresults = collect(
        hasproperty(apexes, :results) ? getproperty(apexes, :results) : apexes,
    )
    if isempty(apexresults)
        reason = if hasproperty(apexes, :failurereason) && !isnothing(apexes.failurereason)
            "apex refinement failed: $(apexes.failurereason)"
        else
            "no apex results available"
        end
        return failed_alkane_series_step_spectra_result(reason, settings)
    end
    if !isnothing(selectedcarbons)
        filter!(
            apexresult -> Int(ladder_spectrum_apex(apexresult).carbon) in selectedcarbons,
            apexresults,
        )
        if isempty(apexresults)
            requested = join(("C$(carbon)" for carbon in settings.carbons), ", ")
            return failed_alkane_series_step_spectra_result(
                "no apex results available for requested carbon numbers $(requested)",
                settings,
            )
        end
    end

    results = Vector{NamedTuple}(undef, length(apexresults))
    Threads.@threads for index in eachindex(apexresults)
        apexresult = apexresults[index]
        apex = ladder_spectrum_apex(apexresult)
        spectrum = ladder_step_mass_spectrum(
            msm,
            apexresult,
            variances;
            scanwindow=scanwindow,
            mzretentionkwargs=mzretentionkwargs,
            variancefloor=variancefloor,
            allownegative=allownegative,
            threaded=false,
        )
        results[index] = (
            carbon=Int(apex.carbon),
            apex=apex,
            spectrum=spectrum,
            success=all(attrs(spectrum).fit_success),
        )
    end

    spectra = Dict{Int, AbstractMassSpectrum}()
    for result in results
        spectra[result.carbon] = result.spectrum
    end

    (
        success=all(result -> result.success, results),
        failurereason=nothing,
        carbonnumbers=[result.carbon for result in results],
        spectra=spectra,
        results=results,
        settings=settings,
    )
end

function normalize_alkane_step_spectrum_carbons(carbons)
    carbons === :all && return nothing
    carbons isa Integer &&
        return Set{Int}([validated_alkane_step_spectrum_carbon(carbons)])

    selected = Set{Int}()
    try
        for carbon in carbons
            push!(selected, validated_alkane_step_spectrum_carbon(carbon))
        end
    catch err
        err isa MethodError || rethrow()
        throw(ArgumentError(
            "step spectrum carbons must be :all, an integer carbon number, " *
            "or a collection of integer carbon numbers"))
    end
    isempty(selected) && throw(ArgumentError(
        "step spectrum carbons must not be empty"))

    selected
end

function validated_alkane_step_spectrum_carbon(carbon)
    carbon isa Integer || throw(ArgumentError(
        "step spectrum carbons must be integer carbon numbers"))
    carbon > 0 || throw(ArgumentError(
        "step spectrum carbons must be positive integer carbon numbers"))

    Int(carbon)
end

function failed_alkane_series_step_spectra_result(reason, settings)
    (
        success=false,
        failurereason=reason,
        carbonnumbers=Int[],
        spectra=Dict{Int, AbstractMassSpectrum}(),
        results=NamedTuple[],
        settings=settings,
    )
end

function ladder_spectrum_apex(apexresult)
    hasproperty(apexresult, :apex) && return getproperty(apexresult, :apex)
    apexresult
end

function ladder_spectrum_apex_is_usable(apex)
    required = (
        :carbon,
        :apex_retention,
        :input_scan_index,
        :input_retention,
        :beta,
        :gamma,
        :apex_x,
        :x_scale,
        :log_floor,
        :mz_indices,
        :ion_intercepts,
        :fit_intensities,
        :mzretentionkwargs,
    )
    all(name -> hasproperty(apex, name), required) || return false
    isfinite(apex.apex_retention) && isfinite(apex.apex_x) &&
        isfinite(apex.x_scale) && apex.x_scale > 0 &&
        isfinite(apex.beta) && isfinite(apex.gamma)
end

function fit_ladder_step_nonreference_ion!(
    intensities::AbstractVector{Float64},
    standarderrors::AbstractVector{Float64},
    zscores::AbstractVector{Float64},
    fitsuccess::AbstractVector{Bool},
    nobservations::AbstractVector{Int},
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    apex,
    scan_retentions,
    targetretention::Real,
    mzretentionkwargs::NamedTuple,
    referencecols::AbstractVector{<:Integer},
    mzindex::Integer,
    nscans::Integer,
    variancefloor::Real,
    allownegative::Bool,
)
    referencecols[mzindex] != 0 && return nothing

    selectedscanindices = apex_scan_indices_for_mz(
        msm,
        scan_retentions,
        targetretention,
        nscans,
        mzindex,
        mzretentionkwargs,
    )
    I = rawintensities(msm)
    numerator = 0.0
    denominator = 0.0
    observations = 0
    for selectedscan in selectedscanindices
        obsretention = apex_observation_retention(
            msm,
            scan_retentions,
            selectedscan,
            mzindex,
            mzretentionkwargs,
        )
        shape = ladder_apex_peak_shape(apex, obsretention)
        y = Float64(I[selectedscan, mzindex])
        variance = max(Float64(variances[selectedscan, mzindex]), Float64(variancefloor))
        isfinite(shape) && isfinite(y) && isfinite(variance) && variance > 0 || continue

        weight = 1.0 / variance
        numerator += weight * shape * y
        denominator += weight * abs2(shape)
        observations += 1
    end

    rawheight = denominator > 0 ? numerator / denominator : NaN
    height = allownegative ? rawheight : max(rawheight, 0.0)
    intensities[mzindex] = isfinite(height) ? height : 0.0
    fitsuccess[mzindex] = isfinite(rawheight) && denominator > 0
    nobservations[mzindex] = observations
    if denominator > 0
        standarderror = sqrt(1.0 / denominator)
        standarderrors[mzindex] = standarderror
        zscores[mzindex] = standarderror > 0 ? height / standarderror : NaN
    end

    nothing
end

function ladder_reference_ion_height(apex, referencecol::Integer)
    z = Float64(apex.ion_intercepts[referencecol]) +
        Float64(apex.beta) * Float64(apex.apex_x) +
        Float64(apex.gamma) * abs2(Float64(apex.apex_x))
    exp(z) - Float64(apex.log_floor)
end

function ladder_apex_peak_shape(apex, observationretention::Real)
    x = (Float64(observationretention) - Float64(apex.input_retention)) /
        Float64(apex.x_scale)
    apex_x = Float64(apex.apex_x)
    exp(
        Float64(apex.beta) * (x - apex_x) +
        Float64(apex.gamma) * (abs2(x) - abs2(apex_x)),
    )
end
