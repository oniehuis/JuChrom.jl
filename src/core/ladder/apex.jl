"""
    ladder_peak_apex(msm, peak, variances; channelinfo, kwargs...)

Estimate the chromatographic apex of one alkane ladder peak with a shared
variance-weighted log-quadratic fit.

The fitted ions are the reference-spectrum channels selected for the peak's carbon number
in `channelinfo`. `peak` can be a scan index when `carbon` is supplied, or a path
candidate with `carbon` and `scanindex` fields.
"""
function ladder_peak_apex(
    msm::MassScanMatrix,
    peak,
    variances::AbstractMatrix{<:Real};
    channelinfo,
    carbon=nothing,
    scanindex_override=nothing,
    scanwindow::Integer=2,
    mzretentionkwargs=nothing,
    logfloorfraction::Real=1e-3,
    variancefloor::Real=1.0,
)
    scanindex = isnothing(scanindex_override) ?
        apex_peak_scan_index(peak) :
        Int(scanindex_override)
    1 <= scanindex <= scancount(msm) || throw(ArgumentError(
        "peak scan index $(scanindex) is outside 1:$(scancount(msm))"))
    scanwindow >= 1 || throw(ArgumentError("scanwindow must be at least 1"))
    isfinite(logfloorfraction) && logfloorfraction > 0 || throw(ArgumentError(
        "logfloorfraction must be finite and positive"))
    isfinite(variancefloor) && variancefloor > 0 || throw(ArgumentError(
        "variancefloor must be finite and positive"))
    size(variances) == size(rawintensities(msm)) || throw(DimensionMismatch(
        "variances must have size $(size(rawintensities(msm))), matching rawintensities(msm)"))

    peakcarbon = apex_peak_carbon(peak, carbon)
    reference = alkane_channel_reference(channelinfo, peakcarbon)
    mzindices = Int.(reference.mzindices)
    isempty(mzindices) && throw(ArgumentError(
        "reference spectrum for C$(peakcarbon) has no selected ions on the msm m/z grid"))

    nominal_scanindices = collect(
        max(1, scanindex - scanwindow):min(scancount(msm), scanindex + scanwindow),
    )
    length(nominal_scanindices) >= 3 || throw(ArgumentError(
        "at least three scans are needed for a quadratic apex fit"))

    mzkwargs = isnothing(mzretentionkwargs) ?
        default_alkane_mzretention_kwargs(msm) :
        mzretentionkwargs
    raw_scan_retentions = Float64.(rawretentions(msm))
    scan_retentions = retentions(msm)
    inputretention = raw_scan_retentions[scanindex]
    nscans = length(nominal_scanindices)
    nmz = length(mzindices)

    scanindicesbymz = apex_scan_indices_by_mz(
        msm,
        scan_retentions,
        raw_scan_retentions,
        scanindex,
        nominal_scanindices,
        mzindices,
        mzkwargs,
    )

    I = rawintensities(msm)
    y = Matrix{Float64}(undef, nscans, nmz)
    fitvariances = Matrix{Float64}(undef, nscans, nmz)
    observationretentions = Matrix{Float64}(undef, nscans, nmz)
    retentionsbymz = Matrix{Float64}(undef, nscans, nmz)
    for mzcol in 1:nmz, scancol in 1:nscans
        selectedscan = scanindicesbymz[scancol, mzcol]
        selectedmz = mzindices[mzcol]
        y[scancol, mzcol] = max(Float64(I[selectedscan, selectedmz]), 0.0)
        fitvariances[scancol, mzcol] = Float64(variances[selectedscan, selectedmz])
        retentionsbymz[scancol, mzcol] = raw_scan_retentions[selectedscan]
        observationretentions[scancol, mzcol] = apex_observation_retention(
            msm,
            scan_retentions,
            selectedscan,
            selectedmz,
            mzkwargs,
        )
    end

    ymax = maximum(y)
    ymax > 0 || throw(ArgumentError("selected reference ions have no positive signal"))
    logfloor = Float64(logfloorfraction) * ymax

    xvalues = observationretentions .- inputretention
    scale = maximum(abs, xvalues)
    scale > 0 || throw(ArgumentError("scan window has zero retention span"))

    nobs = nscans * nmz
    X = zeros(Float64, nobs, nmz + 2)
    z = Vector{Float64}(undef, nobs)
    weights = Vector{Float64}(undef, nobs)

    row = 0
    for mzcol in 1:nmz
        for scancol in 1:nscans
            row += 1
            x = xvalues[scancol, mzcol] / scale
            X[row, mzcol] = 1.0
            X[row, nmz + 1] = x
            X[row, nmz + 2] = abs2(x)
            adjusted = y[scancol, mzcol] + logfloor
            z[row] = log(adjusted)
            weights[row] = apex_variance_log_weight(
                adjusted,
                fitvariances[scancol, mzcol],
                variancefloor,
            )
        end
    end
    any(>(0), weights) || throw(ArgumentError("apex fit has no positive weights"))

    sqrtweights = sqrt.(weights)
    coef = (X .* sqrtweights) \ (z .* sqrtweights)
    fitted = X * coef
    beta = coef[nmz + 1]
    gamma = coef[nmz + 2]

    apexx = gamma < 0 ? -beta / (2 * gamma) : NaN
    apexretention = isfinite(apexx) ? inputretention + apexx * scale : NaN
    apexscanindex = isfinite(apexretention) ?
        apex_fractional_scan_index(raw_scan_retentions, apexretention) :
        NaN
    apexinwindow = isfinite(apexx) &&
        minimum(xvalues) / scale <= apexx <= maximum(xvalues) / scale

    zmean = sum(weights .* z) / sum(weights)
    ssres = sum(weights .* abs2.(z .- fitted))
    sstot = sum(weights .* abs2.(z .- zmean))
    r2 = sstot > 0 ? 1 - ssres / sstot : NaN

    (
        carbon=peakcarbon,
        apex_retention=apexretention,
        apex_scan_index=apexscanindex,
        input_scan_index=scanindex,
        input_retention=inputretention,
        apex_offset_retention=apexretention - inputretention,
        apex_offset_scans=apexscanindex - scanindex,
        success=gamma < 0 && apexinwindow,
        apex_in_window=apexinwindow,
        beta=beta,
        gamma=gamma,
        apex_x=apexx,
        x_scale=scale,
        r2=r2,
        scan_indices=nominal_scanindices,
        retentions=raw_scan_retentions[nominal_scanindices],
        scan_indices_by_mz=scanindicesbymz,
        retentions_by_mz=retentionsbymz,
        observation_retentions=observationretentions,
        fit_intensities=y,
        fit_variances=fitvariances,
        mz_values=collect(mzvalues(msm)[mzindices]),
        mz_indices=mzindices,
        reference_intensities=collect(reference.referenceintensities),
        referenceattrs=reference.referenceattrs,
        ion_intercepts=coef[1:nmz],
        log_floor=logfloor,
        logfloorfraction=Float64(logfloorfraction),
        variance_weighted=true,
        variancefloor=Float64(variancefloor),
        mzretentionkwargs=mzkwargs,
    )
end

"""
    ladder_peak_apex_twopass(msm, peak, variances; channelinfo, kwargs...)

Estimate one ladder apex, then shift the center scan by up to two scans when the first
continuous fit places the apex closer to a neighboring scan.
"""
function ladder_peak_apex_twopass(
    msm::MassScanMatrix,
    peak,
    variances::AbstractMatrix{<:Real};
    scanwindow::Integer=2,
    maxapexoffsetscans::Real=1.0,
    kwargs...,
)
    scanwindow >= 1 || throw(ArgumentError("scanwindow must be at least 1"))
    isfinite(maxapexoffsetscans) && maxapexoffsetscans >= 0 || throw(ArgumentError(
        "maxapexoffsetscans must be finite and nonnegative"))

    raw_scan_retentions = Float64.(rawretentions(msm))
    ladder_scanindex = apex_peak_scan_index(peak)
    ladder_retention = apex_peak_retention(msm, peak, raw_scan_retentions)
    ladder_fractional_scanindex =
        apex_fractional_scan_index(raw_scan_retentions, ladder_retention)
    initial_center_scanindex = nearest_valid_apex_scan_index(
        ladder_fractional_scanindex,
        scancount(msm),
        scanwindow,
    )

    initial_apex = ladder_peak_apex(
        msm,
        peak,
        variances;
        scanwindow=scanwindow,
        kwargs...,
        scanindex_override=initial_center_scanindex,
    )
    initial_continuous_valid = continuous_apex_is_valid(
        initial_apex;
        maxapexoffsetscans=maxapexoffsetscans,
    )

    shift_attempted = false
    shift_direction = 0
    shift_count = 0
    shifted_apex = nothing
    shifted_continuous_valid = false
    final_fit = initial_apex
    continuous_refinement = initial_continuous_valid
    best_valid_fit = initial_continuous_valid ? initial_apex : nothing
    best_valid_center_offset = initial_continuous_valid ?
        apex_center_offset_abs(initial_apex) :
        Inf

    if !apex_center_scan_is_closest(initial_apex)
        shift_direction = apex_one_scan_shift_direction(initial_apex)
        current_apex = initial_apex
        current_center_offset = apex_center_offset_abs(current_apex)
        for _ in 1:2
            next_shift_direction = apex_one_scan_shift_direction(current_apex)
            next_shift_direction == shift_direction || break
            shift_direction != 0 || break

            shifted_center_scanindex = current_apex.input_scan_index + shift_direction
            1 <= shifted_center_scanindex <= scancount(msm) || break
            apex_scan_window_observation_count(
                scancount(msm),
                shifted_center_scanindex,
                scanwindow,
            ) >= 3 || break

            shift_attempted = true
            shift_count += 1
            shifted_apex = ladder_peak_apex(
                msm,
                peak,
                variances;
                scanwindow=scanwindow,
                kwargs...,
                scanindex_override=shifted_center_scanindex,
            )
            shifted_fit_valid = continuous_apex_is_valid(
                shifted_apex;
                maxapexoffsetscans=maxapexoffsetscans,
            )
            shifted_continuous_valid |= shifted_fit_valid
            shifted_center_offset = apex_center_offset_abs(shifted_apex)
            if shifted_fit_valid && shifted_center_offset < best_valid_center_offset
                best_valid_fit = shifted_apex
                best_valid_center_offset = shifted_center_offset
            end

            shifted_center_offset < current_center_offset || break
            apex_center_scan_is_closest(shifted_apex) && break
            current_apex = shifted_apex
            current_center_offset = shifted_center_offset
        end
    end

    if best_valid_fit !== nothing
        final_fit = best_valid_fit
        continuous_refinement = true
    end

    final_apex = continuous_refinement ?
        final_fit :
        snapshot_apex_result(
            final_fit;
            apexretention=ladder_retention,
            apexscanindex=ladder_fractional_scanindex,
        )
    center_scan_is_closest = apex_center_scan_is_closest(final_apex)

    merge(
        final_apex,
        (
            two_pass=true,
            ladder_scan_index=ladder_scanindex,
            ladder_retention=ladder_retention,
            ladder_fractional_scan_index=ladder_fractional_scanindex,
            scanwindow=scanwindow,
            initial_apex=initial_apex,
            initial_center_scan_index=initial_center_scanindex,
            initial_center_retention=raw_scan_retentions[initial_center_scanindex],
            initial_continuous_valid=initial_continuous_valid,
            shift_attempted=shift_attempted,
            shift_direction=shift_direction,
            shift_count=shift_count,
            shifted_apex=shifted_apex,
            shifted_continuous_valid=shifted_continuous_valid,
            continuous_refinement=continuous_refinement,
            fallback_to_snapshot=!continuous_refinement,
            maxapexoffsetscans=Float64(maxapexoffsetscans),
            continuous_fit=final_fit,
            continuous_apex_retention=final_fit.apex_retention,
            continuous_apex_scan_index=final_fit.apex_scan_index,
            center_scan_is_closest=center_scan_is_closest,
            unresolved_center_problem=!continuous_refinement,
            fine_center_scan_index=final_apex.input_scan_index,
            fine_center_retention=raw_scan_retentions[final_apex.input_scan_index],
            apex_offset_from_ladder_retention=final_apex.apex_retention - ladder_retention,
            apex_offset_from_ladder_scans=final_apex.apex_scan_index - ladder_scanindex,
        ),
    )
end

function refine_alkane_series_apexes(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    seriespath;
    channelinfo,
    scanwindow::Integer=2,
    maxapexoffsetscans::Real=1.0,
    mzretentionkwargs=nothing,
    logfloorfraction::Real=1e-3,
    variancefloor::Real=1.0,
)
    validate_alkane_series_variances(msm, variances)
    settings = (
        scanwindow=Int(scanwindow),
        maxapexoffsetscans=Float64(maxapexoffsetscans),
        logfloorfraction=Float64(logfloorfraction),
        variancefloor=Float64(variancefloor),
        variance_weighted=true,
        mzretentionkwargs=mzretentionkwargs,
    )

    if hasproperty(seriespath, :success) && !getproperty(seriespath, :success)
        return failed_alkane_series_apex_result(
            "series path failed: $(getproperty(seriespath, :failurereason))",
            settings,
        )
    end

    path = collect(hasproperty(seriespath, :path) ? getproperty(seriespath, :path) : seriespath)
    isempty(path) && return failed_alkane_series_apex_result(
        "series path contains no peaks",
        settings,
    )

    results = Vector{NamedTuple}(undef, length(path))
    Threads.@threads for index in eachindex(path)
        peak = path[index]
        carbon = apex_peak_carbon(peak, nothing)
        apex = ladder_peak_apex_twopass(
            msm,
            peak,
            variances;
            channelinfo=channelinfo,
            scanwindow=scanwindow,
            maxapexoffsetscans=maxapexoffsetscans,
            mzretentionkwargs=mzretentionkwargs,
            logfloorfraction=logfloorfraction,
            variancefloor=variancefloor,
        )
        result = (
            carbon=carbon,
            peak=peak,
            apex=apex,
            apex_retention=apex.apex_retention,
            apex_scan_index=apex.apex_scan_index,
            ladder_retention=apex.ladder_retention,
            ladder_scan_index=apex.ladder_scan_index,
            input_scan_index=apex.input_scan_index,
            input_retention=apex.input_retention,
            continuous_refinement=apex.continuous_refinement,
            fallback_to_snapshot=apex.fallback_to_snapshot,
            success=apex.success,
        )
        results[index] = result
    end

    bycarbon = Dict{Int, NamedTuple}()
    for result in results
        bycarbon[result.carbon] = result
    end

    (
        success=all(result -> result.success, results),
        failurereason=nothing,
        carbonnumbers=[result.carbon for result in results],
        apexretentions=[result.apex_retention for result in results],
        apexscanindices=[result.apex_scan_index for result in results],
        ladderretentions=[result.ladder_retention for result in results],
        ladderscanindices=[result.ladder_scan_index for result in results],
        results=results,
        bycarbon=bycarbon,
        settings=settings,
    )
end

function failed_alkane_series_apex_result(reason, settings)
    (
        success=false,
        failurereason=reason,
        carbonnumbers=Int[],
        apexretentions=Float64[],
        apexscanindices=Float64[],
        ladderretentions=Float64[],
        ladderscanindices=Int[],
        results=NamedTuple[],
        bycarbon=Dict{Int, NamedTuple}(),
        settings=settings,
    )
end

apex_peak_scan_index(peak::Integer) = Int(peak)

function apex_peak_scan_index(peak)
    if hasproperty(peak, :scanindex)
        return Int(getproperty(peak, :scanindex))
    elseif hasproperty(peak, :scan_index)
        return Int(getproperty(peak, :scan_index))
    end
    throw(ArgumentError("peak must be a scan index or have a scanindex field"))
end

function apex_peak_carbon(peak, carbon)
    if !isnothing(carbon)
        return Int(carbon)
    elseif hasproperty(peak, :carbon)
        return Int(getproperty(peak, :carbon))
    elseif hasproperty(peak, :c_count)
        return Int(getproperty(peak, :c_count))
    end
    throw(ArgumentError("peak carbon must be supplied or present as a carbon field"))
end

function apex_peak_retention(
    msm::MassScanMatrix,
    peak,
    raw_scan_retentions::AbstractVector{<:Real},
)
    hasproperty(peak, :retention) &&
        return raw_retention_value(msm, getproperty(peak, :retention))
    raw_scan_retentions[apex_peak_scan_index(peak)]
end

function apex_scan_indices_by_mz(
    msm::MassScanMatrix,
    scan_retentions,
    raw_scan_retentions::AbstractVector{<:Real},
    center_scanindex::Integer,
    nominal_scanindices::AbstractVector{<:Integer},
    mzindices::AbstractVector{<:Integer},
    mzretentionkwargs::NamedTuple,
)
    nscans = length(nominal_scanindices)
    nscans >= 3 || throw(ArgumentError(
        "at least three scans are needed for a quadratic apex fit"))
    targetretention = raw_scan_retentions[center_scanindex]
    selected = Matrix{Int}(undef, nscans, length(mzindices))

    for (mzcol, mzindex) in pairs(mzindices)
        selected[:, mzcol] .= apex_scan_indices_for_mz(
            msm,
            scan_retentions,
            targetretention,
            nscans,
            mzindex,
            mzretentionkwargs,
        )
    end

    selected
end

function apex_scan_indices_for_mz(
    msm::MassScanMatrix,
    scan_retentions,
    targetretention::Real,
    nscans::Integer,
    mzindex::Integer,
    mzretentionkwargs::NamedTuple,
)
    candidates = Vector{Tuple{Float64, Float64, Int}}(undef, scancount(msm))
    for scanindex in 1:scancount(msm)
        obsretention = apex_observation_retention(
            msm,
            scan_retentions,
            scanindex,
            mzindex,
            mzretentionkwargs,
        )
        candidates[scanindex] = (
            abs(obsretention - targetretention),
            obsretention,
            scanindex,
        )
    end
    sort!(candidates; by=candidate -> candidate[1])
    chosen = candidates[1:nscans]
    sort!(chosen; by=candidate -> candidate[2])

    [candidate[3] for candidate in chosen]
end

function apex_observation_retention(
    msm::MassScanMatrix,
    scan_retentions,
    scanindex::Integer,
    mzindex::Integer,
    mzretentionkwargs::NamedTuple,
)
    raw_mzretention(msm, scan_retentions[scanindex], mzindex, mzretentionkwargs)
end

function apex_variance_log_weight(
    adjusted_intensity::Real,
    variance::Real,
    variancefloor::Real,
)
    v = max(Float64(variance), Float64(variancefloor))
    isfinite(v) && v > 0 || return 0.0
    abs2(Float64(adjusted_intensity)) / v
end

function continuous_apex_is_valid(apex; maxapexoffsetscans::Real)
    apex.success || return false
    isfinite(apex.apex_scan_index) || return false
    abs(Float64(apex.apex_scan_index) - Float64(apex.input_scan_index)) <=
        Float64(maxapexoffsetscans)
end

function apex_center_scan_is_closest(apex)
    isfinite(apex.apex_retention) || return false
    center_distance = abs(Float64(apex.apex_retention) - Float64(apex.input_retention))
    for (scanindex, retention) in zip(apex.scan_indices, apex.retentions)
        scanindex == apex.input_scan_index && continue
        if center_distance > abs(Float64(apex.apex_retention) - Float64(retention))
            return false
        end
    end
    true
end

function apex_center_offset_abs(apex)
    isfinite(apex.apex_scan_index) || return Inf
    abs(Float64(apex.apex_scan_index) - Float64(apex.input_scan_index))
end

function apex_one_scan_shift_direction(apex)
    isfinite(apex.apex_scan_index) || return 0
    offset = Float64(apex.apex_scan_index) - Float64(apex.input_scan_index)
    offset > 0 ? 1 : offset < 0 ? -1 : 0
end

function snapshot_apex_result(apex; apexretention::Real, apexscanindex::Real)
    merge(
        apex,
        (
            apex_retention=Float64(apexretention),
            apex_scan_index=Float64(apexscanindex),
            apex_offset_retention=Float64(apexretention) - apex.input_retention,
            apex_offset_scans=Float64(apexscanindex) - apex.input_scan_index,
            success=true,
            apex_in_window=true,
        ),
    )
end

function nearest_valid_apex_scan_index(
    fractional_scanindex::Real,
    nscans::Integer,
    scanwindow::Integer,
)
    nscans >= 3 || throw(ArgumentError("at least three scans are needed"))
    bestindex = 0
    bestdistance = Inf
    for scanindex in 1:nscans
        apex_scan_window_observation_count(nscans, scanindex, scanwindow) >= 3 || continue
        distance = abs(Float64(scanindex) - Float64(fractional_scanindex))
        if distance < bestdistance
            bestindex = scanindex
            bestdistance = distance
        end
    end
    bestindex == 0 && throw(ArgumentError(
        "no scan index has at least three observations for scanwindow=$(scanwindow)"))
    bestindex
end

function apex_scan_window_observation_count(
    nscans::Integer,
    scanindex::Integer,
    scanwindow::Integer,
)
    firstscan = max(1, scanindex - scanwindow)
    lastscan = min(nscans, scanindex + scanwindow)
    lastscan - firstscan + 1
end

function apex_fractional_scan_index(
    retentions::AbstractVector{<:Real},
    retention::Real,
)
    retention <= first(retentions) && return 1.0
    retention >= last(retentions) && return Float64(length(retentions))

    right = searchsortedfirst(retentions, retention)
    left = right - 1
    leftretention = Float64(retentions[left])
    rightretention = Float64(retentions[right])
    rightretention == leftretention && return Float64(left)
    left + (Float64(retention) - leftretention) / (rightretention - leftretention)
end
