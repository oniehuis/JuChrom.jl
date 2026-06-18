"""
    CountVarianceEstimate

Count-based variance estimate returned by [`countvariances`](@ref).

`variances` stores the unitfree numeric variance matrix. `varianceunit` stores the unit of
that matrix, or `nothing` for unitless input intensities. The remaining fields record the
empirical noise estimate and settings used to produce the variance matrix.
"""
struct CountVarianceEstimate{
    T1<:AbstractMatrix{<:Real},
    T2<:Union{Nothing, Unitful.Units}
}
    variances::T1
    varianceunit::T2
    noisesigma::Float64
    variancefactor::Float64
    zerowindowthreshold::Float64
    positivecountthreshold::Float64
    normalizeddeviationcount::Int
    zerowindowcount::Int
    windowsize::Int
    mintransitioncount::Int
    zerothresholdquantile::Float64
    positivecountquantile::Float64
    intensityfloor::Float64

    function CountVarianceEstimate(
        variances::T1,
        varianceunit::T2,
        noisesigma::Real,
        variancefactor::Real,
        zerowindowthreshold::Real,
        positivecountthreshold::Real,
        normalizeddeviationcount::Integer,
        zerowindowcount::Integer,
        windowsize::Integer,
        mintransitioncount::Integer,
        zerothresholdquantile::Real,
        positivecountquantile::Real,
        intensityfloor::Real
    ) where {T1<:AbstractMatrix{<:Real}, T2<:Union{Nothing, Unitful.Units}}
        !isempty(variances) || throw(ArgumentError("No variance value(s) provided."))
        all(isfinite, variances) || throw(ArgumentError("All variances must be finite."))
        all(value -> value ≥ zero(value), variances) || throw(ArgumentError(
            "All variances must be nonnegative."))
        isfinite(noisesigma) && noisesigma ≥ 0 || throw(ArgumentError(
            "noisesigma must be finite and nonnegative"))
        isfinite(variancefactor) && variancefactor ≥ 0 || throw(ArgumentError(
            "variancefactor must be finite and nonnegative"))
        isfinite(intensityfloor) && intensityfloor > 0 || throw(ArgumentError(
            "intensityfloor must be finite and positive"))

        new{T1, T2}(
            variances,
            varianceunit,
            Float64(noisesigma),
            Float64(variancefactor),
            Float64(zerowindowthreshold),
            Float64(positivecountthreshold),
            Int(normalizeddeviationcount),
            Int(zerowindowcount),
            Int(windowsize),
            Int(mintransitioncount),
            Float64(zerothresholdquantile),
            Float64(positivecountquantile),
            Float64(intensityfloor)
        )
    end
end

varianceunit(estimate::CountVarianceEstimate) = estimate.varianceunit

@inline function variances(
    estimate::CountVarianceEstimate;
    unit::Union{Nothing, Unitful.Units}=nothing
)
    vunit = varianceunit(estimate)
    isnothing(vunit) && return _handle_unitless(estimate.variances, unit, "variances")
    _handle_unitful_convert(estimate.variances, vunit, unit)
end

@inline function rawvariances(
    estimate::CountVarianceEstimate;
    unit::Union{Nothing, Unitful.Units}=nothing
)
    vunit = varianceunit(estimate)
    isnothing(vunit) && return _handle_unitless(estimate.variances, unit, "variances")
    _handle_unitful_strip(estimate.variances, vunit, unit)
end

"""
    countvariances(msm::MassScanMatrix; intensityfloor=nothing,
        positivecountquantile=0.01, zerothresholdquantile=0.99,
        mintransitioncount=7, windowsize=13)

    countvariances(rawcounts::AbstractMatrix; keywords...)

Estimate count-based variances for a mass-scan matrix from uncorrected raw ion
counts.

The first method extracts raw ion counts from `msm` with [`rawintensities`](@ref). The
second method accepts an already extracted count matrix, with scans along the first
dimension and m/z channels along the second dimension. Both methods return a
[`CountVarianceEstimate`](@ref JuChrom.CountVarianceEstimate) whose `variances` field is a
unitfree numeric matrix and whose `varianceunit` field is `nothing` or the squared
intensity unit. All counts must be finite and nonnegative. Negative values indicate
baseline-corrected or otherwise transformed intensities and cause an `ArgumentError`.

The empirical model assumes that the variance is proportional to the ion count,
scaled by a robust single-sample count-noise estimate:

```text
variance[i, j] = noisesigma^2 * max(rawcounts[i, j], intensityfloor)
```

If `intensityfloor` is `nothing`, the floor is chosen as
`max(1.0, stats.zerowindowthreshold, stats.positivecountthreshold)`, where
`stats` are obtained from the internal count-noise estimate.
`positivecountquantile`, `zerothresholdquantile`, `mintransitioncount`, and
`windowsize` control that estimate.

The function returns `CountVarianceEstimate`, containing `variances`, `varianceunit`,
`noisesigma`, `variancefactor`, `zerowindowthreshold`, `positivecountthreshold`,
`normalizeddeviationcount`, `zerowindowcount`, the settings used for the noise estimate,
and the selected `intensityfloor`. This estimator is intended as a fallback when
replicate-based variance estimates are not available.

# Examples
```julia
countvars = countvariances(msm)
σ² = countvars.variances
countvars.noisesigma
```
"""
function countvariances(
    msm::AbstractMassScanMatrix;
    intensityfloor::Union{Nothing, Real, AbstractQuantity{<:Real}}=nothing,
    positivecountquantile::Real=0.01,
    zerothresholdquantile::Real=0.99,
    mintransitioncount::Integer=7,
    windowsize::Integer=13
)
    count_variance_estimate(
        rawintensities(msm),
        intensityunit(msm);
        intensityfloor=intensityfloor,
        positivecountquantile=positivecountquantile,
        zerothresholdquantile=zerothresholdquantile,
        mintransitioncount=mintransitioncount,
        windowsize=windowsize
    )
end

function countvariances(
    rawcounts::AbstractMatrix{<:Union{Real, AbstractQuantity{<:Real}}};
    intensityfloor::Union{Nothing, Real, AbstractQuantity{<:Real}}=nothing,
    positivecountquantile::Real=0.01,
    zerothresholdquantile::Real=0.99,
    mintransitioncount::Integer=7,
    windowsize::Integer=13
)
    rawcounts_unitfree, inputunit = strip_units_checked(rawcounts, "rawcounts")
    count_variance_estimate(
        rawcounts_unitfree,
        inputunit;
        intensityfloor=intensityfloor,
        positivecountquantile=positivecountquantile,
        zerothresholdquantile=zerothresholdquantile,
        mintransitioncount=mintransitioncount,
        windowsize=windowsize
    )
end

function count_variance_estimate(
    rawcounts::AbstractMatrix{<:Real},
    inputunit::Union{Nothing, Unitful.Units};
    intensityfloor::Union{Nothing, Real, AbstractQuantity{<:Real}},
    positivecountquantile::Real,
    zerothresholdquantile::Real,
    mintransitioncount::Integer,
    windowsize::Integer
)

    validate_count_matrix(rawcounts)
    stats = count_noise_stats(
        rawcounts;
        positivecountquantile=positivecountquantile,
        zerothresholdquantile=zerothresholdquantile,
        mintransitioncount=mintransitioncount,
        windowsize=windowsize
    )

    intensityfloor_value = count_variance_intensity_floor_value(intensityfloor, inputunit)
    selected_floor = if isnothing(intensityfloor_value)
        max(1.0, stats.zerowindowthreshold, stats.positivecountthreshold)
    else
        isfinite(intensityfloor_value) && intensityfloor_value > 0 || throw(ArgumentError(
            "intensityfloor must be nothing or a finite positive number"))
        intensityfloor_value
    end

    variances = [
        stats.variancefactor * max(Float64(rawcounts[i, j]), selected_floor)
        for i in axes(rawcounts, 1), j in axes(rawcounts, 2)
    ]

    CountVarianceEstimate(
        variances,
        isnothing(inputunit) ? nothing : inputunit^2,
        stats.noisesigma,
        stats.variancefactor,
        stats.zerowindowthreshold,
        stats.positivecountthreshold,
        stats.normalizeddeviationcount,
        stats.zerowindowcount,
        stats.windowsize,
        stats.mintransitioncount,
        stats.zerothresholdquantile,
        stats.positivecountquantile,
        selected_floor
    )
end

function count_variance_intensity_floor_value(
    intensityfloor::Nothing,
    inputunit::Union{Nothing, Unitful.Units}
)
    nothing
end

function count_variance_intensity_floor_value(
    intensityfloor::Real,
    inputunit::Union{Nothing, Unitful.Units}
)
    Float64(intensityfloor)
end

function count_variance_intensity_floor_value(
    intensityfloor::AbstractQuantity{<:Real},
    inputunit::Nothing
)
    throw(ArgumentError("unitful intensityfloor requires unitful raw counts"))
end

function count_variance_intensity_floor_value(
    intensityfloor::AbstractQuantity{<:Real},
    inputunit::Unitful.Units
)
    try
        return Float64(ustrip(inputunit, intensityfloor))
    catch
        throw(ArgumentError("intensityfloor must be compatible with raw count units"))
    end
end

"""
    count_noise_sigma(rawcounts; windowsize=13, mintransitioncount=7,
        zerothresholdquantile=0.99)

Estimate a single empirical noise scale from one raw ion-count matrix.
"""
function count_noise_sigma(
    rawcounts::AbstractMatrix{<:Real};
    windowsize::Integer=13,
    mintransitioncount::Integer=7,
    zerothresholdquantile::Real=0.99,
)
    stats = count_noise_stats(
        rawcounts;
        windowsize=windowsize,
        mintransitioncount=mintransitioncount,
        zerothresholdquantile=zerothresholdquantile)

    stats.noisesigma
end

count_noise_sigma(msm::AbstractMassScanMatrix; kwargs...) =
    count_noise_sigma(rawintensities(msm); kwargs...)

count_noise_stats(msm::AbstractMassScanMatrix; kwargs...) =
    count_noise_stats(rawintensities(msm); kwargs...)

function count_noise_stats(
    rawcounts::AbstractMatrix{<:Real};
    windowsize::Integer=13,
    mintransitioncount::Integer=7,
    zerothresholdquantile::Real=0.99,
    positivecountquantile::Real=0.01,
)

    validate_count_matrix(rawcounts)
    windowsize > 1 || throw(ArgumentError("windowsize must be larger than 1"))
    0 ≤ mintransitioncount < windowsize || throw(ArgumentError(
        "mintransitioncount must be in 0:(windowsize - 1)"))
    isfinite(zerothresholdquantile) &&
        0 ≤ zerothresholdquantile ≤ 1 || throw(ArgumentError(
            "zerothresholdquantile must be finite and in [0, 1]"))
    isfinite(positivecountquantile) &&
        0 ≤ positivecountquantile ≤ 1 || throw(ArgumentError(
            "positivecountquantile must be finite and in [0, 1]"))

    residuals = Float64[]
    expected_values = Float64[]
    zero_window_expected = Float64[]
    scan_axis = axes(rawcounts, 1)
    mz_axis = axes(rawcounts, 2)
    scan_indices = collect(scan_axis)
    n_scans = length(scan_axis)
    n_scans ≥ windowsize || throw(ArgumentError(
        "rawcounts must contain at least windowsize scans"))

    for mz_idx in mz_axis
        for first_scan in 1:windowsize:(n_scans - windowsize + 1)
            last_scan = first_scan + windowsize - 1
            window_scans = scan_indices[first_scan:last_scan]
            window = @view rawcounts[window_scans, mz_idx]
            all(isfinite, window) || continue

            expected = Float64(median(window))
            if any(iszero, window)
                isfinite(expected) && push!(zero_window_expected, expected)
                continue
            end

            transition_count = level_crossings(window, expected)
            transition_count < mintransitioncount && continue

            for intensity in window
                push!(residuals, Float64(intensity) - expected)
                push!(expected_values, expected)
            end
        end
    end

    zero_threshold = zero_window_threshold(zero_window_expected, zerothresholdquantile)
    normalized_deviations = Float64[]
    for (residual, expected) in zip(residuals, expected_values)
        expected > zero_threshold || continue
        expected > 0 || continue
        push!(normalized_deviations, residual / sqrt(expected))
    end

    length(normalized_deviations) > 1 || throw(ArgumentError(
        "found no suitable scan segments to estimate the empirical noise scale"))
    noise_sigma = 1.4826 * median(abs.(normalized_deviations))
    positive_threshold = positive_count_threshold(rawcounts, positivecountquantile)

    (
        noisesigma=noise_sigma,
        variancefactor=abs2(noise_sigma),
        zerowindowthreshold=zero_threshold,
        positivecountthreshold=positive_threshold,
        normalizeddeviationcount=length(normalized_deviations),
        zerowindowcount=length(zero_window_expected),
        windowsize=windowsize,
        mintransitioncount=mintransitioncount,
        zerothresholdquantile=zerothresholdquantile,
        positivecountquantile=positivecountquantile,
    )
end

function validate_count_matrix(rawcounts::AbstractMatrix{<:Real})
    for value in rawcounts
        isfinite(value) || throw(ArgumentError(
            "countvariances expects finite raw ion counts"))
        value ≥ 0 || throw(ArgumentError(
            "countvariances expects uncorrected raw ion counts; negative values " * 
            "indicate baseline-corrected or transformed intensities"))
    end
    nothing
end

function level_crossings(values::AbstractVector{<:Real}, level::Real)
    state = 0
    count = 0
    for value in values
        nextstate = Float64(value) > level ? 1 : Float64(value) < level ? -1 : state
        nextstate == 0 && continue
        if state ≠ 0 && nextstate ≠ state
            count += 1
        end
        state = nextstate
    end
    count
end

function zero_window_threshold(
    zero_window_expected::AbstractVector{<:Real},
    zero_threshold_quantile::Real,
)
    isempty(zero_window_expected) && return 0.0
    sorted_expected = sort(Float64.(zero_window_expected))
    threshold_idx = clamp(
        ceil(Int, length(sorted_expected) * Float64(zero_threshold_quantile)),
        1,
        length(sorted_expected))
    sorted_expected[threshold_idx]
end

function positive_count_threshold(rawcounts::AbstractMatrix{<:Real}, q::Real)
    positives = Float64[]
    for value in rawcounts
        isfinite(value) && value > 0 && push!(positives, Float64(value))
    end
    isempty(positives) && return 0.0
    sort!(positives)
    threshold_idx = clamp(
        ceil(Int, length(positives) * Float64(q)),
        1,
        length(positives))
    positives[threshold_idx]
end
