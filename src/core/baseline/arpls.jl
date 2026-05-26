"""
    arpls(
        intensities::AbstractVector{<:Real};
        λ::Real=1e5,
        ratio::Real=1e-3,
        maxiter::Integer=10_000,
        nonnegative::Bool=false,
        variances=nothing,
        variancefloor::Real=1e-12,
        peakthreshold::Real=4.0,
        peakslope::Real=1.0,
        zerothreshold::Real=1e-8,
        zeroweight::Real=0.01
    ) -> Vector{Float64}

    arpls(msm::MassScanMatrix; keywords...) -> MassScanMatrix

Estimate an arPLS baseline using asymmetrically reweighted penalized least-squares
smoothing, following the approach in Baek et al. (2015).

For a vector, `intensities` is treated as a signal sampled on an equally spaced grid, 
and the fitted baseline is returned. For a `MassScanMatrix`, one independent baseline is 
fitted for each m/z channel; retention coordinates, axes, units, and metadata are copied 
to the returned matrix, but the retention coordinates are not used in the estimate.

Compared with the original implementation described by Baek et al. (2015), this
implementation:

- accepts optional per-point measurement variances through `variances`; the smoothing
  solve is precision-weighted, the residuals used for the logistic update are
  standardized, and precision weights are normalized so scaling all variances by
  the same factor does not change the result;
- can give zero/dropout points selected by `zerothreshold` a reduced fit weight
  using `zeroweight` and exclude them from the negative-residual noise estimate;
- can clip fitted baseline values to be nonnegative after each smoothing solve;
- exposes the logistic peak threshold and transition slope as keywords rather than
  fixing them in the weighting formula;

`λ` specifies the smoothness penalty, with higher values enforcing smoother baselines
and lower values allowing shorter-scale baseline drift.

`ratio` specifies the relative weight-change tolerance for convergence. At each
iteration, the arPLS weights are updated from their previous values; the iteration
stops when the norm of this update, divided by the norm of the previous weights,
falls below `ratio`. Smaller values impose stricter convergence, which can matter
when late weight changes still affect peak areas or baseline shape; larger values
stop the iteration earlier. `maxiter` specifies the maximum number of reweighting
iterations; a warning is emitted if convergence is not reached.

`peakthreshold` and `peakslope` specify the logistic peak weighting. `peakthreshold` is
measured in the noise scale estimated from negative residuals; residuals are
variance-standardized when `variances` are provided. Higher values keep more positive
residuals baseline-like, and lower values reject peaks more aggressively. If the
estimated baseline follows the lower envelope of the noise floor, increasing
`peakthreshold` can make the fit less aggressive. `peakslope` specifies the steepness
of the transition between baseline-like and peak-like residuals.

`nonnegative` specifies whether negative fitted baseline values are clipped upward
to zero after each smoothing solve.

`variances` specifies optional per-point measurement variances for heteroscedastic noise,
and `variancefloor` specifies the lower bound used to avoid infinite precision for zero
or tiny variances.

`zerothreshold` specifies the threshold below which intensities are treated as
zero/dropout points. For points below this threshold, fit weights are multiplied by
`zeroweight`, and the points are excluded from the negative-residual noise estimate.


# Examples
```julia
estimated_baseline = arpls(intensities; λ=1e5)
baseline_matrix = arpls(msm; variances=variance_matrix)
```

# References
Baek SJ, Park A, Ahn YJ, Choo J (2015): Baseline correction using asymmetrically
reweighted penalized least squares smoothing. Analyst 140: 250-257.
"""
function arpls(
    intensities::AbstractVector{<:Real};
    λ::T1=1e5,
    ratio::T2=1e-3,
    maxiter::T3=10_000,
    nonnegative::Bool=false,
    variances=nothing,
    variancefloor::Real=1e-12,
    peakthreshold::Real=4.0,
    peakslope::Real=1.0,
    zerothreshold::Real=1e-8,
    zeroweight::Real=0.01
    ) where {T1<:Real, T2<:Real, T3<:Integer}

    validatearplsargs(
        intensities, λ, ratio, maxiter, peakthreshold, peakslope,
        zerothreshold, zeroweight, variancefloor)
    varianceweights = arplsvarianceweights(variances, length(intensities), variancefloor)
    arplsfit(
        intensities, λ * arplspenalty(length(intensities)), ratio, maxiter,
        nonnegative, varianceweights.residualscales, varianceweights.precisionweights,
        peakthreshold, peakslope, zerothreshold, zeroweight)
end

function arpls(
    msm::MassScanMatrix;
    λ::T1=1e5,
    ratio::T2=1e-3,
    maxiter::T3=10_000,
    nonnegative::Bool=false,
    variances=nothing,
    variancefloor::Real=1e-12,
    peakthreshold::Real=4.0,
    peakslope::Real=1.0,
    zerothreshold::Real=1e-8,
    zeroweight::Real=0.01
    ) where {T1<:Real, T2<:Real, T3<:Integer}

    baselines = zeros(Float64, scancount(msm), mzcount(msm))
    ints = rawintensities(msm)
    scancount(msm) ≥ 3 || throw(ArgumentError("Need at least three data points"))
    all(isfinite, ints) || throw(ArgumentError("All intensities must be finite"))
    validatearplssettings(
        λ, ratio, maxiter, peakthreshold, peakslope, zerothreshold, zeroweight,
        variancefloor)
    if !isnothing(variances)
        variances isa AbstractMatrix{<:Real} || throw(
            ArgumentError("variances must be a matrix matching intensities"))
        size(variances) == size(ints) || throw(
            ArgumentError("variances must have the same size as intensities"))
    end

    penalty = λ * arplspenalty(scancount(msm))
    mzs = rawmzvalues(msm)
    @threads :dynamic for i in axes(ints, 2)
        varianceweights = isnothing(variances) ?
            (residualscales=nothing, precisionweights=nothing) :
            arplsvarianceweights(@view(variances[:, i]), scancount(msm), variancefloor)
        baselines[:, i] .= arplsfit(
            @view(ints[:, i]), penalty, ratio, maxiter, nonnegative,
            varianceweights.residualscales, varianceweights.precisionweights,
            peakthreshold, peakslope, zerothreshold, zeroweight,
            "m/z channel $i (m/z=$(mzs[i]))")
    end

    MassScanMatrix(
       deepcopy(rawretentions(msm)),
       deepcopy(retentionunit(msm)),
       deepcopy(rawmzvalues(msm)),
       deepcopy(mzunit(msm)),
       baselines,
       deepcopy(intensityunit(msm)),
       level=deepcopy(level(msm)),
       instrument=deepcopy(instrument(msm)),
       acquisition=deepcopy(acquisition(msm)),
       user=deepcopy(user(msm)),
       sample=deepcopy(sample(msm)),
       extras=deepcopy(extras(msm))
    )
end

function arplspenalty(n::Integer)
    rows = Int[]
    cols = Int[]
    vals = Float64[]

    for i in 1:(n - 2)
        push!(rows, i); push!(cols, i); push!(vals, 1.0)
        push!(rows, i); push!(cols, i + 1); push!(vals, -2.0)
        push!(rows, i); push!(cols, i + 2); push!(vals, 1.0)
    end

    D = sparse(rows, cols, vals, n - 2, n)
    D' * D
end

function validatearplsargs(
    intensities, λ, ratio, maxiter, peakthreshold, peakslope, zerothreshold, zeroweight,
    variancefloor
    )

    length(intensities) ≥ 3 || throw(ArgumentError("Need at least three data points"))
    all(isfinite, intensities) || throw(ArgumentError("All intensities must be finite"))
    validatearplssettings(
        λ, ratio, maxiter, peakthreshold, peakslope, zerothreshold, zeroweight,
        variancefloor)
    nothing
end

function validatearplssettings(
    λ, ratio, maxiter, peakthreshold, peakslope, zerothreshold, zeroweight,
    variancefloor
    )

    λ > 0 || throw(ArgumentError("Smoothing parameter must be positive"))
    ratio > 0 || throw(ArgumentError("Stopping ratio must be positive"))
    maxiter > 0 || throw(ArgumentError("Maximum iteration count must be positive"))
    isfinite(peakthreshold) && peakthreshold > 0 ||
        throw(ArgumentError("Peak threshold must be finite and positive"))
    isfinite(peakslope) && peakslope > 0 ||
        throw(ArgumentError("Peak slope must be finite and positive"))
    isfinite(zerothreshold) && zerothreshold ≥ 0 ||
        throw(ArgumentError("Zero threshold must be finite and non-negative"))
    isfinite(zeroweight) && zeroweight > 0 ||
        throw(ArgumentError("Zero weight must be finite and positive"))
    isfinite(variancefloor) && variancefloor > 0 ||
        throw(ArgumentError("Variance floor must be finite and positive"))
    nothing
end

function arplsvarianceweights(variances, n::Integer, variancefloor::Real)
    isnothing(variances) && return (residualscales=nothing, precisionweights=nothing)

    variances isa AbstractVector{<:Real} || throw(
        ArgumentError("variances must be a vector with the same length as intensities"))
    length(variances) == n || throw(
        ArgumentError("variances must have the same length as intensities"))
    all(v -> isfinite(v) && v ≥ 0, variances) || throw(
        ArgumentError("All variances must be finite and non-negative"))

    variancefloor_eff = max(Float64(variancefloor), eps(Float64))
    variancevalues = max.(Float64.(variances), variancefloor_eff)
    precisionweights = @. 1 / variancevalues
    precisionweights ./= mean(precisionweights)
    residualscales = @. 1 / sqrt(precisionweights)
    (residualscales=residualscales, precisionweights=precisionweights)
end

function updatearplsweights!(
    nextweights,
    weightedresiduals,
    zeromask,
    peakthreshold::Real,
    peakslope::Real
    )

    negativeresiduals = weightedresiduals[(.!zeromask) .& (weightedresiduals .< 0)]
    length(negativeresiduals) ≥ 2 || return false

    m = mean(negativeresiduals)
    s = std(negativeresiduals)
    isfinite(s) && s > 0 || return false

    @. nextweights = 1 / (
        1 + exp(peakslope * (weightedresiduals - (peakthreshold * s - m)) / s))

    true
end

function arplsfit(
    intensities::AbstractVector{<:Real},
    penalty::AbstractMatrix{<:Real},
    ratio::Real,
    maxiter::Integer,
    nonnegative::Bool,
    residualscales,
    precisionweights,
    peakthreshold::Real,
    peakslope::Real,
    zerothreshold::Real,
    zeroweight::Real,
    warncontext::Union{Nothing, String}=nothing
    )

    y = Float64.(intensities)
    n = length(y)
    weights = ones(Float64, n)
    nextweights = similar(weights)
    baseline = similar(y)
    residuals = similar(y)
    weightedresiduals = isnothing(residualscales) ? residuals : similar(y)
    fitweights = similar(weights)
    zeromask = y .< zerothreshold
    zerofactors = ones(Float64, n)
    zerofactors[zeromask] .= zeroweight
    rhs = similar(y)
    converged = false

    for _ in 1:maxiter
        @. fitweights = weights * zerofactors
        if !isnothing(precisionweights)
            @. fitweights = weights * precisionweights
            @. fitweights *= zerofactors
        end
        systemmatrix = penalty + Diagonal(fitweights)
        @. rhs = fitweights * y
        baseline .= cholesky(systemmatrix) \ rhs
        if nonnegative
            @. baseline = max(baseline, 0.0)
        end
        @. residuals = y - baseline
        if !isnothing(residualscales)
            @. weightedresiduals = residuals / residualscales
        end

        if !updatearplsweights!(
            nextweights, weightedresiduals, zeromask, peakthreshold, peakslope)

            converged = true
            break
        end

        if norm(weights - nextweights) / norm(weights) < ratio
            converged = true
            break
        end

        weights .= nextweights
    end

    if !converged
        context = isnothing(warncontext) ? "" : " for $warncontext"
        iterationnoun = maxiter == 1 ? "iteration" : "iterations"
        @warn "arPLS baseline estimation did not converge in $maxiter $iterationnoun$context."
    end

    baseline
end
