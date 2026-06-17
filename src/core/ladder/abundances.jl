"""
    AlkaneAbundanceWindow

Peak window on an alkane abundance track.
"""
struct AlkaneAbundanceWindow
    ladderstep::Int
    leftindex::Int
    apexindex::Int
    rightindex::Int
    leftabundance::Float64
    apexabundance::Float64
    rightabundance::Float64
    threshold::Float64
    leftstop::Symbol
    rightstop::Symbol
end

"""
    AlkaneAbundanceInfo

Result of the alkane abundance stage.
"""
struct AlkaneAbundanceInfo{S <: NamedTuple}
    abundances::Dict{Int, Vector{Float64}}
    abundancevariances::Dict{Int, Vector{Float64}}
    windows::Dict{Int, Vector{AlkaneAbundanceWindow}}
    settings::S
end

Base.keys(::AlkaneAbundanceInfo) = (
    :abundances,
    :abundancevariances,
    :windows,
    :settings,
)

function alkaneabundanceinfo(
    msm::MassScanMatrix,
    variances,
    channelinfo;
    variancefloor::Real=1.0,
    nonnegative::Bool=false,
    thresholdfraction::Real=0.05,
    minrisez::Real=10.0,
)
    abundances = alkaneabundances(
        msm,
        variances,
        channelinfo;
        variancefloor=variancefloor,
        nonnegative=nonnegative,
    )
    abundancevariances = alkaneabundancevariances(
        variances,
        channelinfo;
        variancefloor=variancefloor,
    )
    windows = alkaneabundancewindows(
        abundances,
        abundancevariances;
        thresholdfraction=thresholdfraction,
        minrisez=minrisez,
    )

    AlkaneAbundanceInfo(
        abundances,
        abundancevariances,
        windows,
        (
            variancefloor=Float64(variancefloor),
            nonnegative=nonnegative,
            thresholdfraction=Float64(thresholdfraction),
            minrisez=Float64(minrisez),
        ),
    )
end

function alkaneabundances(
    msm::MassScanMatrix,
    variances,
    channelinfo;
    variancefloor::Real=1.0,
    nonnegative::Bool=false,
)
    validate_alkane_abundance_variancefloor(variancefloor)
    validate_alkane_series_variances(msm, variances)

    X = rawintensities(msm)
    abundances = Dict{Int, Vector{Float64}}()
    for reference in channelinfo.references
        carbon = Int(reference.carbon)
        mzindices = reference.mzindices
        referenceintensities = alkane_reference_abundance_intensities(reference)
        abundances[carbon] = alkaneabundance(
            X,
            variances,
            mzindices,
            referenceintensities;
            variancefloor=Float64(variancefloor),
            nonnegative=nonnegative,
            carbon=carbon,
        )
    end

    abundances
end

function alkaneabundancevariances(
    variances,
    channelinfo;
    variancefloor::Real=1.0,
)
    validate_alkane_abundance_variancefloor(variancefloor)
    validate_alkane_abundance_variances_matrix(variances)

    abundancevariances = Dict{Int, Vector{Float64}}()
    for reference in channelinfo.references
        carbon = Int(reference.carbon)
        mzindices = reference.mzindices
        referenceintensities = alkane_reference_abundance_intensities(reference)
        abundancevariances[carbon] = alkaneabundancevariance(
            variances,
            mzindices,
            referenceintensities;
            variancefloor=Float64(variancefloor),
            carbon=carbon,
        )
    end

    abundancevariances
end

function alkaneabundancewindows(
    abundances,
    abundancevariances=nothing;
    thresholdfraction::Real=0.05,
    minrisez::Real=10.0,
)
    validate_alkane_abundance_window_settings(thresholdfraction, minrisez)
    isnothing(abundancevariances) || abundancevariances isa AbstractDict || throw(
        ArgumentError("abundancevariances must be nothing or a dictionary"))

    windows = Dict{Int, Vector{AlkaneAbundanceWindow}}()
    for (carbon, abundance) in pairs(abundances)
        step = Int(carbon)
        abundancevalues = alkane_abundance_values(abundance, step)
        variancevalues = isnothing(abundancevariances) ? nothing :
            alkane_abundance_window_variances(abundancevariances, step, abundancevalues)
        maxima = localmaxima(abundancevalues)
        minima = Set(localmaxima(-abundancevalues))

        stepwindows = AlkaneAbundanceWindow[]
        for apexindex in sort(collect(maxima); by=index -> abundancevalues[index], rev=true)
            abundancevalues[apexindex] > 0 || continue
            alkane_abundance_index_is_in_window(apexindex, stepwindows) && continue
            push!(
                stepwindows,
                alkane_abundance_peak_window(
                    step,
                    abundancevalues,
                    minima,
                    maxima,
                    apexindex,
                    Float64(thresholdfraction),
                    variancevalues,
                    Float64(minrisez),
                ),
            )
        end

        windows[step] = stepwindows
    end

    windows
end

function alkaneabundance(
    X::AbstractMatrix{<:Real},
    variances::AbstractMatrix{<:Real},
    mzindices,
    referenceintensities::AbstractVector{Float64};
    variancefloor::Float64,
    nonnegative::Bool,
    carbon::Integer,
)
    validate_alkane_reference_abundance_channels(
        mzindices,
        referenceintensities,
        size(X, 2);
        carbon=carbon,
    )
    size(variances) == size(X) || throw(DimensionMismatch(
        "variances must have size $(size(X)), matching rawintensities(msm)"))

    abundance = Vector{Float64}(undef, size(X, 1))
    for scanindex in axes(X, 1)
        numerator = 0.0
        denominator = 0.0
        for (mzindex, referenceintensity) in zip(mzindices, referenceintensities)
            weight = inv(max(Float64(variances[scanindex, mzindex]), variancefloor))
            numerator += weight * Float64(X[scanindex, mzindex]) * referenceintensity
            denominator += weight * referenceintensity^2
        end

        denominator > 0 || throw(ArgumentError(
            "reference spectrum for C$(carbon) has no positive abundance weight"))
        value = numerator / denominator
        abundance[scanindex] = nonnegative ? max(value, 0.0) : value
    end

    abundance
end

function alkaneabundancevariance(
    variances::AbstractMatrix{<:Real},
    mzindices,
    referenceintensities::AbstractVector{Float64};
    variancefloor::Float64,
    carbon::Integer,
)
    validate_alkane_reference_abundance_channels(
        mzindices,
        referenceintensities,
        size(variances, 2);
        carbon=carbon,
    )

    abundancevariance = Vector{Float64}(undef, size(variances, 1))
    for scanindex in axes(variances, 1)
        denominator = 0.0
        for (mzindex, referenceintensity) in zip(mzindices, referenceintensities)
            weight = inv(max(Float64(variances[scanindex, mzindex]), variancefloor))
            denominator += weight * referenceintensity^2
        end

        denominator > 0 || throw(ArgumentError(
            "reference spectrum for C$(carbon) has no positive abundance weight"))
        abundancevariance[scanindex] = inv(denominator)
    end

    abundancevariance
end

function alkane_reference_abundance_intensities(reference)
    referenceintensities = Float64.(reference.referenceintensities)
    all(isfinite, referenceintensities) || throw(ArgumentError(
        "reference intensities for C$(reference.carbon) must be finite"))
    all(≥(0), referenceintensities) || throw(ArgumentError(
        "reference intensities for C$(reference.carbon) must be nonnegative"))

    referenceintensities
end

function validate_alkane_abundance_variancefloor(variancefloor)
    variancefloor isa Real || throw(ArgumentError("variancefloor must be real"))
    isfinite(variancefloor) && variancefloor > 0 || throw(ArgumentError(
        "variancefloor must be finite and positive"))

    nothing
end

function validate_alkane_abundance_variances_matrix(variances)
    variances isa AbstractMatrix{<:Real} || throw(ArgumentError(
        "variances must be a matrix"))
    for variance in variances
        isfinite(variance) || throw(ArgumentError("variances must be finite"))
        variance ≥ 0 || throw(ArgumentError("variances must be nonnegative"))
    end

    nothing
end

function validate_alkane_reference_abundance_channels(
    mzindices,
    referenceintensities::AbstractVector{Float64},
    mzcount::Integer;
    carbon::Integer,
)
    length(mzindices) == length(referenceintensities) || throw(DimensionMismatch(
        "mzindices and referenceintensities for C$(carbon) must have the same length"))
    isempty(mzindices) && throw(ArgumentError(
        "reference spectrum for C$(carbon) has no m/z channels"))
    all(index -> index isa Integer, mzindices) || throw(ArgumentError(
        "mzindices for C$(carbon) must contain integers"))
    all(index -> 1 ≤ index ≤ mzcount, mzindices) || throw(ArgumentError(
        "mzindices for C$(carbon) must be valid matrix columns"))
    any(>(0), referenceintensities) || throw(ArgumentError(
        "reference spectrum for C$(carbon) must contain a positive intensity"))

    nothing
end

function validate_alkane_abundance_window_settings(thresholdfraction, minrisez)
    thresholdfraction isa Real || throw(ArgumentError(
        "thresholdfraction must be real"))
    isfinite(thresholdfraction) && 0 ≤ thresholdfraction < 1 || throw(
        ArgumentError("thresholdfraction must be finite and in [0, 1)"))

    minrisez isa Real || throw(ArgumentError("minrisez must be real"))
    isfinite(minrisez) && minrisez ≥ 0 || throw(ArgumentError(
        "minrisez must be finite and nonnegative"))

    nothing
end

function alkane_abundance_values(abundance, carbon::Integer)
    abundance isa AbstractVector{<:Real} || throw(ArgumentError(
        "abundance vector for C$(carbon) must be a real vector"))
    all(isfinite, abundance) || throw(ArgumentError(
        "abundance vector for C$(carbon) must be finite"))

    Float64.(abundance)
end

function alkane_abundance_window_variances(
    abundancevariances::AbstractDict,
    carbon::Integer,
    abundance::AbstractVector{Float64},
)
    haskey(abundancevariances, carbon) || throw(ArgumentError(
        "abundancevariances does not contain a vector for C$(carbon)"))
    variance = Float64.(abundancevariances[carbon])
    length(variance) == length(abundance) || throw(DimensionMismatch(
        "abundance variance vector for C$(carbon) must match abundance length"))
    all(isfinite, variance) || throw(ArgumentError(
        "abundance variance vector for C$(carbon) must be finite"))
    all(≥(0), variance) || throw(ArgumentError(
        "abundance variance vector for C$(carbon) must be nonnegative"))

    variance
end

function alkane_abundance_index_is_in_window(index::Integer, windows::AbstractVector)
    for window in windows
        window.leftindex ≤ index ≤ window.rightindex && return true
    end

    false
end

function alkane_abundance_peak_window(
    carbon::Integer,
    abundance::AbstractVector{Float64},
    minima::AbstractSet{<:Integer},
    maxima::AbstractVector{<:Integer},
    apexindex::Integer,
    thresholdfraction::Float64,
    abundancevariances::Union{Nothing, AbstractVector{Float64}},
    minrisez::Float64,
)
    apexabundance = abundance[apexindex]
    threshold = thresholdfraction * apexabundance

    leftindex, leftstop = alkane_abundance_peak_window_left(
        abundance,
        abundancevariances,
        minima,
        maxima,
        apexindex,
        threshold,
        minrisez,
    )
    rightindex, rightstop = alkane_abundance_peak_window_right(
        abundance,
        abundancevariances,
        minima,
        maxima,
        apexindex,
        threshold,
        minrisez,
    )

    AlkaneAbundanceWindow(
        carbon,
        leftindex,
        apexindex,
        rightindex,
        abundance[leftindex],
        apexabundance,
        abundance[rightindex],
        threshold,
        leftstop,
        rightstop,
    )
end

function alkane_abundance_peak_window_left(
    abundance::AbstractVector{Float64},
    abundancevariances::Union{Nothing, AbstractVector{Float64}},
    minima::AbstractSet{<:Integer},
    maxima::AbstractVector{<:Integer},
    apexindex::Integer,
    threshold::Float64,
    minrisez::Float64,
)
    index = Int(apexindex)
    while index > firstindex(abundance)
        index -= 1
        if index in minima && alkane_abundance_local_minimum_stops_window(
                abundance,
                abundancevariances,
                maxima,
                apexindex,
                index,
                :left,
                minrisez,
            )
            return index, :localminimum
        end
        abundance[index] ≤ threshold && return index, :threshold
    end

    index, :boundary
end

function alkane_abundance_peak_window_right(
    abundance::AbstractVector{Float64},
    abundancevariances::Union{Nothing, AbstractVector{Float64}},
    minima::AbstractSet{<:Integer},
    maxima::AbstractVector{<:Integer},
    apexindex::Integer,
    threshold::Float64,
    minrisez::Float64,
)
    index = Int(apexindex)
    while index < lastindex(abundance)
        index += 1
        if index in minima && alkane_abundance_local_minimum_stops_window(
                abundance,
                abundancevariances,
                maxima,
                apexindex,
                index,
                :right,
                minrisez,
            )
            return index, :localminimum
        end
        abundance[index] ≤ threshold && return index, :threshold
    end

    index, :boundary
end

function alkane_abundance_local_minimum_stops_window(
    abundance::AbstractVector{Float64},
    abundancevariances::Union{Nothing, AbstractVector{Float64}},
    maxima::AbstractVector{<:Integer},
    apexindex::Integer,
    minimumindex::Integer,
    direction::Symbol,
    minrisez::Float64,
)
    isnothing(abundancevariances) && return true

    neighboringpeakindex = alkane_neighboring_peak_across_minimum(
        maxima,
        apexindex,
        minimumindex,
        direction,
    )
    isnothing(neighboringpeakindex) && return false

    rise = abundance[neighboringpeakindex] - abundance[minimumindex]
    rise > 0 || return false
    riseerror = sqrt(
        abundancevariances[neighboringpeakindex] + abundancevariances[minimumindex],
    )
    riseerror > 0 || return true

    rise / riseerror ≥ minrisez
end

function alkane_neighboring_peak_across_minimum(
    maxima::AbstractVector{<:Integer},
    apexindex::Integer,
    minimumindex::Integer,
    direction::Symbol,
)
    if direction === :right
        for peakindex in maxima
            peakindex > minimumindex && peakindex != apexindex && return Int(peakindex)
        end
        return nothing
    elseif direction === :left
        for peakindex in Iterators.reverse(maxima)
            peakindex < minimumindex && peakindex != apexindex && return Int(peakindex)
        end
        return nothing
    end

    throw(ArgumentError("direction must be :left or :right"))
end
