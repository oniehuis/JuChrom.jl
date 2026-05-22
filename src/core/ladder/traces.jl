"""
    alkanematchtrace(msm::MassScanMatrix, carbon::Integer;
        standard=defaultalkanestandard(), channelinfo=nothing,
        minrelativeintensity=0.05, spectralpower=0.5,
        positiveonly=true, smoothing=3) -> ChromScanSeries

Infer a cosine-match trace for the n-alkane with carbon number `carbon`.

The trace is computed only from integer m/z channels shared by `msm` and the selected
reference spectrum. Scan intensities and reference intensities are transformed by
`x -> max(x, 0)^spectralpower` by default before cosine similarity is computed.
"""
function alkanematchtrace(
    msm::MassScanMatrix,
    carbon::Integer;
    standard=defaultalkanestandard(),
    channelinfo=nothing,
    minrelativeintensity=0.05,
    spectralpower=0.5,
    positiveonly::Bool=true,
    smoothing::Integer=3
)
    carbon > 0 || throw(ArgumentError("carbon must be positive"))
    validate_alkane_channel_settings(minrelativeintensity)
    isfinite(spectralpower) && spectralpower > 0 || throw(ArgumentError(
        "spectralpower must be finite and positive"))
    smoothing >= 0 || throw(ArgumentError("smoothing must be nonnegative"))
    validate_alkane_mz_unit(mzunit(msm))

    channels = isnothing(channelinfo) ?
        alkane_mz_channels(
            msm;
            standard=standard,
            carbonrange=carbon:carbon,
            minrelativeintensity=minrelativeintensity,
        ) :
        channelinfo
    reference = alkane_channel_reference(channels, carbon)

    mzindices = reference.mzindices
    referencevalues = spectral_power_transform(
        reference.referenceintensities,
        spectralpower,
        positiveonly,
    )
    maximum(referencevalues) > 0 || throw(ArgumentError(
        "reference spectrum for C$(carbon) has no positive intensity on the msm m/z grid"))

    I = rawintensities(msm)
    scanvalues = Vector{Float64}(undef, length(mzindices))
    match = Vector{Float64}(undef, size(I, 1))
    for scanindex in axes(I, 1)
        @inbounds for (k, mzindex) in pairs(mzindices)
            scanvalues[k] = spectral_power_transform_value(
                I[scanindex, mzindex],
                spectralpower,
                positiveonly,
            )
        end
        similarity = cossim(scanvalues, referencevalues, positiveonly)
        match[scanindex] = isnan(similarity) ? 0.0 : similarity
    end
    match = trapezoid_smooth(match, smoothing)

    matched_mz_unit = hasproperty(channels, :mzunit) ? channels.mzunit : mzunit(msm)
    traceattrs = (
        carbon=Int(carbon),
        referenceattrs=reference.referenceattrs,
        matchedmz=reference.mzvalues,
        mzunit=matched_mz_unit,
        minrelativeintensity=Float64(minrelativeintensity),
        spectralpower=Float64(spectralpower),
        positiveonly=positiveonly,
        smoothing=smoothing
    )
    chroms = [
        ChromScan(rt, retentionunit(msm), y, nothing; attrs=traceattrs)
        for (rt, y) in zip(rawretentions(msm), match)
    ]

    ChromScanSeries(
        chroms;
        instrument=instrument(msm),
        acquisition=acquisition(msm),
        user=user(msm),
        sample=sample(msm),
        extras=merge(
            copy(extras(msm)),
            Dict{String, Any}(
                "trace_type" => "alkane_match",
                "carbon" => Int(carbon),
                "reference_attrs" => reference.referenceattrs,
                "matched_mz" => reference.mzvalues,
                "mz_unit" => matched_mz_unit,
                "min_relative_intensity" => Float64(minrelativeintensity),
                "spectral_power" => Float64(spectralpower),
                "positive_only" => positiveonly,
                "smoothing" => smoothing,
            ),
        ),
    )
end

"""
    alkane_molecular_ion_trace(msm::MassScanMatrix, carbon::Integer;
        ion_window=2, step_mass=14, normalize=true, smoothing=3) -> ChromScanSeries

Infer a molecular-ion contrast trace for the n-alkane with carbon number `carbon`.

The nominal molecular ion is `M = 14 * carbon + 2`. For every scan this computes the
signal in the `M` window minus the corresponding windows one alkane mass step below and
above. Negative contrasts are clipped to zero before smoothing and optional normalization.
"""
function alkane_molecular_ion_trace(
    msm::MassScanMatrix,
    carbon::Integer;
    ion_window::Integer=2,
    step_mass::Real=14,
    normalize::Bool=true,
    smoothing::Integer=3
)
    carbon > 0 || throw(ArgumentError("carbon must be positive"))
    ion_window >= 0 || throw(ArgumentError("ion_window must be nonnegative"))
    isfinite(step_mass) && step_mass > 0 || throw(ArgumentError(
        "step_mass must be finite and positive"))
    smoothing >= 0 || throw(ArgumentError("smoothing must be nonnegative"))

    integer_mzs = alkane_mz_bins(msm)
    molecular_ion = 14 * Int(carbon) + 2
    I = rawintensities(msm)

    center_idx = mz_window_indices(integer_mzs, molecular_ion, ion_window)
    lower_idx = mz_window_indices(integer_mzs, molecular_ion - step_mass, ion_window)
    upper_idx = mz_window_indices(integer_mzs, molecular_ion + step_mass, ion_window)

    trace = max.(
        row_sum(I, center_idx) .- row_sum(I, lower_idx) .- row_sum(I, upper_idx),
        0.0,
    )
    trace = trapezoid_smooth(trace, smoothing)

    if normalize
        trace_max = maximum(trace)
        trace_max > 0 && (trace ./= trace_max)
    end

    display_mzs = mzvalues(msm)
    traceattrs = (
        carbon=Int(carbon),
        molecularion=molecular_ion,
        ionwindow=ion_window,
        stepmass=Float64(step_mass),
        mzunit=mzunit(msm),
        centermz=collect(display_mzs[center_idx]),
        lowermz=collect(display_mzs[lower_idx]),
        uppermz=collect(display_mzs[upper_idx]),
        normalized=normalize,
        smoothing=smoothing
    )
    traceunit = normalize ? nothing : intensityunit(msm)
    chroms = [
        ChromScan(rt, retentionunit(msm), y, traceunit; attrs=traceattrs)
        for (rt, y) in zip(rawretentions(msm), trace)
    ]

    ChromScanSeries(
        chroms;
        instrument=instrument(msm),
        acquisition=acquisition(msm),
        user=user(msm),
        sample=sample(msm),
        extras=merge(
            copy(extras(msm)),
            Dict{String, Any}(
                "trace_type" => "alkane_molecular_ion",
                "carbon" => Int(carbon),
                "molecular_ion" => molecular_ion,
                "ion_window" => ion_window,
                "step_mass" => Float64(step_mass),
                "mz_unit" => mzunit(msm),
                "normalized" => normalize,
                "smoothing" => smoothing,
            ),
        ),
    )
end

"""
    alkanemzpeakdistancetrace(msm::MassScanMatrix, carbon::Integer;
        variances, standard=defaultalkanestandard(), channelinfo=nothing,
        minrelativeintensity=0.05, mzretentionkwargs=default, smoothing=3,
        peakzmin=4.0, peakeps=1e-12, peakvariancefloor=1.0)
        -> ChromScanSeries

Infer an m/z-peak-distance evidence trace for the n-alkane with carbon number `carbon`.

For each m/z trace, local maxima with z-scores at least `peakzmin` are converted to
parabolic apex retentions. For every scan, selected reference-spectrum m/z channels are
shifted to their within-scan sampling retention via [`mzretention`](@ref). Their distances
to the nearest m/z peak apex are summarized by the median and converted to unitless evidence with
`1 / (1 + d_scan^2)`. Missing or nonfinite distances contribute zero evidence.
"""
function alkanemzpeakdistancetrace(
    msm::MassScanMatrix,
    carbon::Integer;
    variances,
    standard=defaultalkanestandard(),
    channelinfo=nothing,
    minrelativeintensity=0.05,
    mzretentionkwargs=default_alkane_mzretention_kwargs(msm),
    smoothing::Integer=3,
    peakzmin::Real=4.0,
    peakeps::Real=1e-12,
    peakvariancefloor::Real=1.0
)
    carbon > 0 || throw(ArgumentError("carbon must be positive"))
    validate_alkane_series_variances(msm, variances)
    validate_alkane_channel_settings(minrelativeintensity)
    smoothing >= 0 || throw(ArgumentError("smoothing must be nonnegative"))
    validate_mzpeak_parameters(peakzmin, peakeps, peakvariancefloor)

    channels = isnothing(channelinfo) ?
        alkane_mz_channels(
            msm;
            standard=standard,
            carbonrange=carbon:carbon,
            minrelativeintensity=minrelativeintensity,
        ) :
        channelinfo
    reference = alkane_channel_reference(channels, carbon)

    selected_mz_indices = reference.mzindices
    isempty(selected_mz_indices) && throw(ArgumentError(
        "reference spectrum for C$(carbon) has no selected ions on the msm m/z grid"))

    scan_retention_values = retentions(msm)
    scan_interval = median_scan_interval(msm)
    peak_times_by_mz_raw = alkane_mzpeak_times_by_mz(
        msm,
        variances;
        mzretentionkwargs=mzretentionkwargs,
        peakzmin=peakzmin,
        peakeps=peakeps,
        peakvariancefloor=peakvariancefloor,
    )
    selected_peak_counts = [length(peak_times_by_mz_raw[index]) for index in selected_mz_indices]

    distances = Float64[]
    sizehint!(distances, length(selected_mz_indices))
    evidence = Vector{Float64}(undef, scancount(msm))
    for scan_index in 1:scancount(msm)
        empty!(distances)
        for mz_index in selected_mz_indices
            peak_times = peak_times_by_mz_raw[mz_index]
            isempty(peak_times) && continue

            scan_ion_retention = raw_mzretention(
                msm,
                scan_retention_values[scan_index],
                mz_index,
                mzretentionkwargs,
            )
            push!(distances, nearest_sorted_distance(peak_times, scan_ion_retention))
        end

        median_distance = isempty(distances) ? Inf : median(distances)
        distance_in_scans = median_distance / scan_interval
        evidence[scan_index] = isfinite(distance_in_scans) && distance_in_scans >= 0 ?
            clamp(1 / (1 + abs2(distance_in_scans)), 0.0, 1.0) :
            0.0
    end
    evidence = trapezoid_smooth(evidence, smoothing)

    mz_unit = hasproperty(channels, :mzunit) ? channels.mzunit : mzunit(msm)
    traceattrs = (
        carbon=Int(carbon),
        referenceattrs=reference.referenceattrs,
        selectedmz=reference.mzvalues,
        selectedmzindices=selected_mz_indices,
        referenceintensities=reference.referenceintensities,
        mzunit=mz_unit,
        minrelativeintensity=Float64(minrelativeintensity),
        mzretentionkwargs=mzretentionkwargs,
        scaninterval=scan_interval,
        retentionunit=retentionunit(msm),
        selectedpeakcounts=selected_peak_counts,
        peakzmin=Float64(peakzmin),
        peakeps=Float64(peakeps),
        peakvariancefloor=Float64(peakvariancefloor),
        transform=:inverse_quadratic,
        smoothing=smoothing,
        scaledbetweenzeroandone=true
    )
    chroms = [
        ChromScan(rt, retentionunit(msm), y, nothing; attrs=traceattrs)
        for (rt, y) in zip(rawretentions(msm), evidence)
    ]

    ChromScanSeries(
        chroms;
        instrument=instrument(msm),
        acquisition=acquisition(msm),
        user=user(msm),
        sample=sample(msm),
        extras=merge(
            copy(extras(msm)),
            Dict{String, Any}(
                "trace_type" => "alkane_mzpeak_distance_evidence",
                "carbon" => Int(carbon),
                "reference_attrs" => reference.referenceattrs,
                "selected_mz" => reference.mzvalues,
                "selected_mz_indices" => selected_mz_indices,
                "mz_unit" => mz_unit,
                "min_relative_intensity" => Float64(minrelativeintensity),
                "mzretention_kwargs" => mzretentionkwargs,
                "scan_interval" => scan_interval,
                "retention_unit" => retentionunit(msm),
                "selected_peak_counts" => selected_peak_counts,
                "peak_z_min" => Float64(peakzmin),
                "peak_eps" => Float64(peakeps),
                "peak_variance_floor" => Float64(peakvariancefloor),
                "transform" => :inverse_quadratic,
                "smoothing" => smoothing,
                "scaled_between_zero_and_one" => true,
            ),
        ),
    )
end

"""
    alkaneevidencetrace(moleculariontrace, matchtrace; mzpeakdistancetrace=nothing)
        -> ChromScanSeries

Combine alkane evidence traces by multiplying molecular-ion evidence, spectral-match
evidence, and optional m/z-peak-distance evidence.
"""
function alkaneevidencetrace(
    moleculariontrace::AbstractChromScanSeries,
    matchtrace::AbstractChromScanSeries;
    mzpeakdistancetrace=nothing,
)
    validate_product_trace_alignment(moleculariontrace, matchtrace, "matchtrace")
    if !isnothing(mzpeakdistancetrace)
        validate_product_trace_alignment(
            moleculariontrace,
            mzpeakdistancetrace,
            "mzpeakdistancetrace",
        )
    end

    ion = Float64.(rawintensities(moleculariontrace))
    match = Float64.(rawintensities(matchtrace))
    distance = isnothing(mzpeakdistancetrace) ?
        ones(Float64, length(ion)) :
        Float64.(rawintensities(mzpeakdistancetrace))
    evidence = ion .* match .* distance

    carbon = trace_carbon(moleculariontrace)
    traceattrs = (
        trace_type="alkane_evidence_product",
        carbon=carbon,
        distanceevidenceincluded=!isnothing(mzpeakdistancetrace),
    )
    chroms = [
        ChromScan(rt, retentionunit(moleculariontrace), y, nothing; attrs=traceattrs)
        for (rt, y) in zip(rawretentions(moleculariontrace), evidence)
    ]

    ChromScanSeries(
        chroms;
        instrument=instrument(moleculariontrace),
        acquisition=acquisition(moleculariontrace),
        user=user(moleculariontrace),
        sample=sample(moleculariontrace),
        extras=merge(
            copy(extras(moleculariontrace)),
            Dict{String, Any}(
                "trace_type" => "alkane_evidence_product",
                "carbon" => carbon,
                "distance_evidence_included" => !isnothing(mzpeakdistancetrace),
            ),
        ),
    )
end

function alkane_series_traces(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    carbon_numbers::AbstractVector{<:Integer};
    standard,
    channelinfo,
    minrelativeintensity::Real,
    molecularionwindow::Integer,
    spectralpower::Real,
    smoothing::Integer,
    includedistance::Bool,
    mzpeakzmin::Real,
    mzpeakeps::Real,
    mzpeakvariancefloor::Real,
    mzretentionkwargs,
)
    match_traces = Dict{Int, ChromScanSeries}()
    molecularion_traces = Dict{Int, ChromScanSeries}()
    mzpeakdistance_traces = Dict{Int, ChromScanSeries}()
    evidence_traces = Dict{Int, ChromScanSeries}()

    for carbon in carbon_numbers
        carbon = Int(carbon)
        match_traces[carbon] = alkanematchtrace(
            msm,
            carbon;
            standard=standard,
            channelinfo=channelinfo,
            minrelativeintensity=minrelativeintensity,
            spectralpower=spectralpower,
            smoothing=smoothing,
        )
        molecularion_traces[carbon] = alkane_molecular_ion_trace(
            msm,
            carbon;
            ion_window=molecularionwindow,
            normalize=true,
            smoothing=smoothing,
        )
        if includedistance
            distance_kwargs = isnothing(mzretentionkwargs) ?
                NamedTuple() :
                (mzretentionkwargs=mzretentionkwargs,)
            mzpeakdistance_traces[carbon] = alkanemzpeakdistancetrace(
                msm,
                carbon;
                variances=variances,
                standard=standard,
                channelinfo=channelinfo,
                minrelativeintensity=minrelativeintensity,
                smoothing=smoothing,
                peakzmin=mzpeakzmin,
                peakeps=mzpeakeps,
                peakvariancefloor=mzpeakvariancefloor,
                distance_kwargs...,
            )
        end
        evidence_traces[carbon] = alkaneevidencetrace(
            molecularion_traces[carbon],
            match_traces[carbon];
            mzpeakdistancetrace=get(mzpeakdistance_traces, carbon, nothing),
        )
    end

    (
        match=match_traces,
        molecularion=molecularion_traces,
        mzpeakdistance=mzpeakdistance_traces,
        evidence=evidence_traces,
        carbonrange=Int.(carbon_numbers),
    )
end

function mz_window_indices(
    mzs::AbstractVector{<:Real},
    center::Real,
    halfwidth::Integer
)
    findall(mz -> abs(Float64(mz) - Float64(center)) <= halfwidth, mzs)
end

function validate_product_trace_alignment(
    reference::AbstractChromScanSeries,
    candidate::AbstractChromScanSeries,
    label::AbstractString
)
    length(rawintensities(candidate)) == length(rawintensities(reference)) || throw(
        DimensionMismatch("$(label) must have the same length as moleculariontrace"))
    rawretentions(candidate) == rawretentions(reference) || throw(DimensionMismatch(
        "$(label) must have the same retention values as moleculariontrace"))
    retentionunit(candidate) == retentionunit(reference) || throw(DimensionMismatch(
        "$(label) must have the same retention unit as moleculariontrace"))

    nothing
end

function trace_carbon(trace::AbstractChromScanSeries)
    trace_attrs = attrs(first(trace))
    hasproperty(trace_attrs, :carbon) && return Int(trace_attrs.carbon)

    get(extras(trace), "carbon", nothing) isa Integer && return Int(extras(trace)["carbon"])
    missing
end

function default_alkane_mzretention_kwargs(msm::MassScanMatrix)
    retention_values = retentions(msm)
    length(retention_values) >= 2 || throw(ArgumentError(
        "at least two scans are needed to infer mzretention scaninterval"))

    (
        retentionref=:start,
        scaninterval=mean(diff(retention_values)),
        mzcount=mzcount(msm),
        order=:descending,
        dwellref=:middle,
        dwell=:homogeneous,
    )
end

function median_scan_interval(msm::MassScanMatrix)
    retention_values = Float64.(rawretentions(msm))
    length(retention_values) >= 2 || throw(ArgumentError(
        "at least two scans are needed to convert m/z-peak distances to evidence"))
    intervals = diff(retention_values)
    all(isfinite, intervals) || throw(ArgumentError(
        "scan intervals must be finite"))
    scan_interval = median(intervals)
    scan_interval > 0 || throw(ArgumentError(
        "median scan interval must be positive"))

    scan_interval
end

function row_sum(I::AbstractMatrix{<:Real}, indices::AbstractVector{<:Integer})
    isempty(indices) && return zeros(Float64, size(I, 1))
    Float64.(vec(sum(@view(I[:, indices]); dims=2)))
end

function alkane_mzpeak_times_by_mz(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real};
    mzretentionkwargs::NamedTuple,
    peakzmin::Real,
    peakeps::Real,
    peakvariancefloor::Real
)
    raw_counts = rawintensities(msm)
    retention_values = retentions(msm)
    peak_times_by_mz = Vector{Vector{Float64}}(undef, mzcount(msm))

    for mz_index in 1:mzcount(msm)
        trace = @view raw_counts[:, mz_index]
        variance_trace = floored_mzpeak_variances(
            @view(variances[:, mz_index]),
            peakvariancefloor,
        )
        z_trace = mzpeak_zscores(trace, variance_trace, peakeps)
        local_maxima = localmaxima(trace)

        peak_times = Float64[]
        sizehint!(peak_times, length(local_maxima))
        for scan_index in local_maxima
            z_trace[scan_index] >= peakzmin || continue
            apex_retention, apex_intensity = parabola_apex(
                retention_values,
                trace,
                scan_index,
            )
            corrected_retention = raw_mzretention(
                msm,
                apex_retention,
                mz_index,
                mzretentionkwargs,
            )
            isfinite(apex_intensity) && push!(peak_times, corrected_retention)
        end
        sort!(peak_times)
        peak_times_by_mz[mz_index] = peak_times
    end

    peak_times_by_mz
end

function validate_mzpeak_parameters(
    peakzmin,
    peakeps,
    peakvariancefloor,
)
    isfinite(peakzmin) && peakzmin >= 0 || throw(ArgumentError(
        "peakzmin must be finite and nonnegative"))
    isfinite(peakeps) && peakeps >= 0 || throw(ArgumentError(
        "peakeps must be finite and nonnegative"))
    isfinite(peakvariancefloor) && peakvariancefloor > 0 || throw(ArgumentError(
        "peakvariancefloor must be finite and positive"))

    nothing
end

function floored_mzpeak_variances(
    variances::AbstractVector{<:Real},
    variancefloor::Real
)
    [
        isfinite(variance) && variance >= variancefloor ? Float64(variance) :
            Float64(variancefloor)
        for variance in variances
    ]
end

function mzpeak_zscores(
    trace::AbstractVector{<:Real},
    variances::AbstractVector{<:Real},
    peakeps::Real
)
    [
        Float64(value) / sqrt(Float64(variance) + Float64(peakeps))
        for (value, variance) in zip(trace, variances)
    ]
end

function parabola_apex(
    x::AbstractVector,
    y::AbstractVector{<:Real},
    index::Integer
)
    n = length(y)
    length(x) == n || throw(ArgumentError("x and y must have the same length"))
    1 < index < n || return (x[index], y[index])

    @inbounds begin
        y_minus = Float64(y[index - 1])
        y_zero = Float64(y[index])
        y_plus = Float64(y[index + 1])
        x_index = x[index]
        dx = x[index + 1] - x_index
    end
    denominator = y_minus - 2 * y_zero + y_plus
    denominator == 0 && return (x_index, y_zero)
    delta = (y_minus - y_plus) / (2 * denominator)
    x_apex = x_index + delta * dx
    y_apex = y_zero - (y_minus - y_plus)^2 / (8 * denominator)

    x_apex, y_apex
end

function raw_mzretention(
    msm::MassScanMatrix,
    scan_retention,
    mz_index::Integer,
    mzretentionkwargs::NamedTuple
)
    kwargs = merge(mzretentionkwargs, (mzindex=Int(mz_index),))
    raw_retention_value(msm, Core.kwcall(kwargs, mzretention, scan_retention))
end

function raw_retention_value(msm::MassScanMatrix, retention)
    if retention isa AbstractQuantity
        unit = retentionunit(msm)
        isnothing(unit) && throw(ArgumentError(
            "unitful peak retentions require a unitful msm retention axis"))
        try
            return Float64(ustrip(unit, retention))
        catch
            throw(ArgumentError(
                "peak retentions must be compatible with retentionunit(msm)"))
        end
    end

    isfinite(retention) || throw(ArgumentError("peak retentions must be finite"))
    Float64(retention)
end

function nearest_sorted_distance(sorted_values::AbstractVector{<:Real}, value::Real)
    isempty(sorted_values) && return Inf
    right = searchsortedfirst(sorted_values, value)
    if right == 1
        return abs(Float64(sorted_values[1]) - Float64(value))
    elseif right > length(sorted_values)
        return abs(Float64(value) - Float64(sorted_values[end]))
    end

    min(
        abs(Float64(sorted_values[right]) - Float64(value)),
        abs(Float64(value) - Float64(sorted_values[right - 1])),
    )
end

function alkane_channel_reference(channelinfo, carbon::Integer)
    hasproperty(channelinfo, :references) || throw(ArgumentError(
        "channelinfo must contain reference channel matches"))
    for reference in getproperty(channelinfo, :references)
        if reference.carbon == carbon
            return reference
        end
    end

    throw(ArgumentError("channelinfo does not contain reference channels for C$(carbon)"))
end

function spectral_power_transform(
    values::AbstractVector{<:Real},
    spectralpower::Real,
    positiveonly::Bool
)
    [spectral_power_transform_value(value, spectralpower, positiveonly) for value in values]
end

function spectral_power_transform_value(
    value::Real,
    spectralpower::Real,
    positiveonly::Bool
)
    positiveonly && return max(Float64(value), 0.0) ^ spectralpower
    value = Float64(value)
    sign(value) * abs(value) ^ spectralpower
end

function trapezoid_smooth(values::AbstractVector{<:Real}, smoothing::Integer)
    smoothing >= 0 || throw(ArgumentError("smoothing must be nonnegative"))
    smoothing == 0 && return Float64.(values)

    kernel = trapezoid_kernel(smoothing)
    offsets = -smoothing:smoothing
    smoothed = Vector{Float64}(undef, length(values))

    @inbounds for index in eachindex(values)
        numerator = 0.0
        denominator = 0.0
        for (k, offset) in enumerate(offsets)
            neighbor = index + offset
            1 <= neighbor <= length(values) || continue
            weight = kernel[k]
            numerator += weight * Float64(values[neighbor])
            denominator += weight
        end
        smoothed[index] = numerator / denominator
    end

    smoothed
end

function trapezoid_kernel(smoothing::Integer)
    weights = [
        Float64(min(smoothing, smoothing + 1 - abs(offset)))
        for offset in -smoothing:smoothing
    ]
    weights ./ sum(weights)
end
