"""
    countvariances(msm::MassScanMatrix; intensityfloor=nothing,
        positivecountquantile=0.01, zerothresholdquantile=0.99,
        mintransitioncount=7, windowsize=13)

    countvariances(rawcounts::AbstractMatrix; keywords...)

Estimate count-based variances for a mass-scan matrix from uncorrected raw ion
counts.

The first method extracts raw ion counts from `msm` with [`rawintensities`](@ref).
The second method is a lower-level method for already extracted count matrices,
with scans along the first dimension and m/z channels along the second dimension.
All counts must be finite and nonnegative. Negative values indicate
baseline-corrected or otherwise transformed intensities and cause an
`ArgumentError`.

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

The function returns a flat named tuple containing `variances`, `noisesigma`,
`variancefactor`, `zerowindowthreshold`, `positivecountthreshold`,
`normalizeddeviationcount`, `zerowindowcount`, the settings used for the noise
estimate, and the selected `intensityfloor`. This estimator is intended as a
fallback when replicate-based variance estimates are not available.

# Examples
```julia
countvars = countvariances(msm)
σ² = countvars.variances
countvars.noisesigma
```
"""
countvariances(msm::AbstractMassScanMatrix; kwargs...) =
    countvariances(rawintensities(msm); kwargs...)

function countvariances(
    rawcounts::AbstractMatrix{<:Real};
    intensityfloor::Union{Nothing, Real}=nothing,
    positivecountquantile::Real=0.01,
    zerothresholdquantile::Real=0.99,
    mintransitioncount::Integer=7,
    windowsize::Integer=13,
)

    validatecountmatrix(rawcounts)
    stats = countnoisestats(
        rawcounts;
        positivecountquantile=positivecountquantile,
        zerothresholdquantile=zerothresholdquantile,
        mintransitioncount=mintransitioncount,
        windowsize=windowsize)

    selectedfloor = if isnothing(intensityfloor)
        max(1.0, stats.zerowindowthreshold, stats.positivecountthreshold)
    else
        isfinite(intensityfloor) && intensityfloor > 0 || throw(ArgumentError(
            "intensityfloor must be nothing or a finite positive number"))
        Float64(intensityfloor)
    end

    variances = [
        stats.variancefactor * max(Float64(rawcounts[i, j]), selectedfloor)
        for i in axes(rawcounts, 1), j in axes(rawcounts, 2)
    ]

    merge((variances=variances,), stats, (intensityfloor=selectedfloor,))
end

"""
    countnoisesigma(rawcounts; windowsize=13, mintransitioncount=7,
        zerothresholdquantile=0.99)

Estimate a single empirical noise scale from one raw ion-count matrix.
"""
function countnoisesigma(
    rawcounts::AbstractMatrix{<:Real};
    windowsize::Integer=13,
    mintransitioncount::Integer=7,
    zerothresholdquantile::Real=0.99,
)
    stats = countnoisestats(
        rawcounts;
        windowsize=windowsize,
        mintransitioncount=mintransitioncount,
        zerothresholdquantile=zerothresholdquantile)

    stats.noisesigma
end

countnoisesigma(msm::AbstractMassScanMatrix; kwargs...) =
    countnoisesigma(rawintensities(msm); kwargs...)

countnoisestats(msm::AbstractMassScanMatrix; kwargs...) =
    countnoisestats(rawintensities(msm); kwargs...)

function countnoisestats(
    rawcounts::AbstractMatrix{<:Real};
    windowsize::Integer=13,
    mintransitioncount::Integer=7,
    zerothresholdquantile::Real=0.99,
    positivecountquantile::Real=0.01,
)

    validatecountmatrix(rawcounts)
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
    expectedvalues = Float64[]
    zerowindowexpected = Float64[]
    scanaxis = axes(rawcounts, 1)
    mzaxis = axes(rawcounts, 2)
    scanindices = collect(scanaxis)
    nscans = length(scanaxis)
    nscans ≥ windowsize || throw(ArgumentError(
        "rawcounts must contain at least windowsize scans"))

    for mzidx in mzaxis
        for firstscan in 1:windowsize:(nscans - windowsize + 1)
            lastscan = firstscan + windowsize - 1
            windowscans = scanindices[firstscan:lastscan]
            window = @view rawcounts[windowscans, mzidx]
            all(isfinite, window) || continue

            expected = Float64(median(window))
            if any(iszero, window)
                isfinite(expected) && push!(zerowindowexpected, expected)
                continue
            end

            transitioncount = levelcrossings(window, expected)
            transitioncount < mintransitioncount && continue

            for intensity in window
                push!(residuals, Float64(intensity) - expected)
                push!(expectedvalues, expected)
            end
        end
    end

    zerothreshold = zerowindowthreshold(zerowindowexpected, zerothresholdquantile)
    normalizeddeviations = Float64[]
    for (residual, expected) in zip(residuals, expectedvalues)
        expected > zerothreshold || continue
        expected > 0 || continue
        push!(normalizeddeviations, residual / sqrt(expected))
    end

    length(normalizeddeviations) > 1 || throw(ArgumentError(
        "found no suitable scan segments to estimate the empirical noise scale"))
    noisesigma = 1.4826 * median(abs.(normalizeddeviations))
    positivethreshold = positivecountthreshold(rawcounts, positivecountquantile)

    (
        noisesigma=noisesigma,
        variancefactor=abs2(noisesigma),
        zerowindowthreshold=zerothreshold,
        positivecountthreshold=positivethreshold,
        normalizeddeviationcount=length(normalizeddeviations),
        zerowindowcount=length(zerowindowexpected),
        windowsize=windowsize,
        mintransitioncount=mintransitioncount,
        zerothresholdquantile=zerothresholdquantile,
        positivecountquantile=positivecountquantile,
    )
end

function validatecountmatrix(rawcounts::AbstractMatrix{<:Real})
    for value in rawcounts
        isfinite(value) || throw(ArgumentError(
            "countvariances expects finite raw ion counts"))
        value ≥ 0 || throw(ArgumentError(
            "countvariances expects uncorrected raw ion counts; negative values " * 
            "indicate baseline-corrected or transformed intensities"))
    end
    nothing
end

function levelcrossings(values::AbstractVector{<:Real}, level::Real)
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

function zerowindowthreshold(
    zerowindowexpected::AbstractVector{<:Real},
    zerothresholdquantile::Real,
)
    isempty(zerowindowexpected) && return 0.0
    sortedexpected = sort(Float64.(zerowindowexpected))
    thresholdidx = clamp(
        ceil(Int, length(sortedexpected) * Float64(zerothresholdquantile)),
        1,
        length(sortedexpected))
    sortedexpected[thresholdidx]
end

function positivecountthreshold(rawcounts::AbstractMatrix{<:Real}, q::Real)
    positives = Float64[]
    for value in rawcounts
        isfinite(value) && value > 0 && push!(positives, Float64(value))
    end
    isempty(positives) && return 0.0
    sort!(positives)
    thresholdidx = clamp(
        ceil(Int, length(positives) * Float64(q)),
        1,
        length(positives))
    positives[thresholdidx]
end
