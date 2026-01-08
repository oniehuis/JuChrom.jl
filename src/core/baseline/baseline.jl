# ── airpls ────────────────────────────────────────────────────────────────────────────────

"""
 airpls(retentions::AbstractVector{<:Real}, intensities::AbstractVector{<:Real}; 
        λ::Real=1e7, max_iter::Integer=10^4, ...) -> Vector{Float64}

Estimate a smooth baseline from spectroscopic/chromatographic data using Adaptive
Iteratively Reweighted Penalized Least Squares (airPLS).

Uses asymmetric weighting to suppress peaks while preserving baseline trend.
Enhanced with weight propagation across peak shoulders - when peak regions are
detected, their low weights are extended to neighboring points following the
peak slope, improving baseline estimation for broad or overlapping peaks.

# Arguments
- `retentions::AbstractVector{<:Real}`: Independent variable (e.g., retention times, 
  retention indices)
- `intensities::AbstractVector{<:Real}`: Signal data containing baseline + peaks
- `λ::Real=1e7`: Smoothing parameter (higher -> smoother baseline)
- `variances::AbstractVector{<:Real}=nothing`: Optional measurement uncertainties for each 
  intensity value
- `improvement_threshold::Real=1e-6`: Minimum improvement to continue optimization
- `max_iter::Integer=10^4`: Maximum iterations
- `no_improvement_limit::Integer=10`: Stop after this many non-improving iterations
- `threshold_factor::Real=1.96`: Peak detection sensitivity (higher = less sensitive). 
  A factor of 1.282 corresponds to 90% coverage of the one-sided noise distribution.
- `zero_threshold::Real=1e-8`: Threshold for zero-valued points
- `zero_weight::Real=0.01`: Weight factor for zero-threshold points

# Algorithm
1. Iteratively solves penalized least squares with asymmetric weights
2. Down-weights positive residuals (above baseline) more than negative ones
3. Propagates low weights across peak shoulders using slope analysis
4. Tracks optimization progress and returns best solution found

Returns estimated baseline vector with same length as input intensities.

# References
- Zhang ZM, Chen S, Liang YZ (2010): Baseline correction using adaptive
 iteratively reweighted penalized least squares. – Analyst 135: 1138-1146.

# Example
```julia
# Basic usage
baseline = airpls(wavelengths, spectrum; λ=1e6)

# With measurement uncertainties
baseline = airpls(times, signal; variances=errors.^2, λ=1e5)

# Sensitive peak detection
baseline = airpls(x, y; threshold_factor=2.0, no_improvement_limit=15)
```
"""
function airpls(
    retentions::AbstractVector{<:Real},
    intensities::AbstractVector{<:Real};
    variances::T1=nothing,
    improvement_threshold::T2=1e-6,
    λ::T3=1e7,
    max_iter::T4=10^4,
    no_improvement_limit::T5=10,
    threshold_factor::T6=1.96,
    zero_threshold::T7=1e-8,
    zero_weight::T8=0.01
    ) where {
    T1<:Union{Nothing, AbstractVector{<:Real}},
    T2<:Real,
    T3<:Real,
    T4<:Integer,
    T5<:Integer,
    T6<:Real,
    T7<:Real,
    T8<:Real}
    
    n = length(retentions)
    length(intensities) == n || throw(
        ArgumentError("retentions and intensities must have the same length"))
    n ≥ 3 || throw(ArgumentError("Need at least three data points"))
    all(isfinite, retentions) || throw(ArgumentError("All retentions must be finite"))
    all(isfinite, intensities) || throw(ArgumentError("All intensities must be finite"))
    λ > 0 || throw(ArgumentError("Smoothing parameter must be positive"))
    no_improvement_limit > 0 || throw(
        ArgumentError("No improvement limit must be positive"))
    threshold_factor > 0 || throw(ArgumentError("Threshold factor must be positive"))
    zero_threshold ≥ 0 || throw(ArgumentError("Zero threshold must be non-negative"))
    zero_weight > 0 || throw(ArgumentError("Zero weight must be positive"))

    if !isnothing(variances)
        length(variances) == n || throw(
            ArgumentError("variances must have same length as data"))
        all(v -> v ≥ 0, variances) || throw(
            ArgumentError("All variances must be ≥ 0"))
    end

    # Initialize objective tracker
    fit_tracker = FitTracker(n, no_improvement_limit, improvement_threshold)
    
    D = build_second_derivative_matrix(retentions)
    smoothness_matrix = λ * (D' * D)
    weights = ones(Float64, n)
    zero_mask = intensities .≤ zero_threshold
    rhs = similar(intensities, Float64)
    baseline_old = zeros(Float64, n)
    baseline_new = similar(intensities, Float64)
    residuals = similar(intensities, Float64)
    converged = false
    convergence_reason = "max_iterations"
    
    std_devs = isnothing(variances) ? nothing : sqrt.(max.(variances, zero_threshold))

    for iter in 1:max_iter
        weights[zero_mask] .*= zero_weight
        system_matrix = smoothness_matrix + Diagonal(weights)
        mul!(rhs, Diagonal(weights), intensities)
        baseline_new .= cholesky(system_matrix) \ rhs
        clamp!(baseline_new, zero_threshold, Inf)
        @. residuals = intensities - baseline_new
        
        # Calculate current objective value
        current_fit = calculate_fit(baseline_new, intensities, weights, smoothness_matrix)
        
        # Update best solution tracker
        update!(fit_tracker, current_fit, baseline_new, iter)
        
        # Check if we should stop (no improvement for too long)
        if stop_optimization(fit_tracker)
           converged = true
           convergence_reason = "no_improvement"
           break
        end
        
        # Update weights for next iteration
        weights = compute_weights(residuals, std_devs, threshold_factor, zero_mask, 
                                  zero_threshold)
        baseline_old .= baseline_new
    end
    
    if !converged
        @warn "airPLS did not converge in $max_iter iterations. " *
              "Returning best solution found at iteration $(fit_tracker.iteration_of_best)."
    end
    
    # Return the best baseline found during optimization
    baseline(fit_tracker)
end

"""
    airpls(msm::MassScanMatrix; ...) -> MassScanMatrix

Estimate baseline from mass scan matrix data using airPLS algorithm. The method assumes
the retrieved intensity counts represent raw Poisson shot-noise data with identical
variance. Do not apply this method to data whose counts have been transformed (e.g., by
applying a retention mapper). Use instead the method `airpls(msm::MassScanMatrix,
variances::AbstractMatrix{<:Real}; ...)`.

Convenience method that extracts retention times and intensities of each mz channel from
the MassScanMatrix, applies baseline correction, and returns the identified baselines as a
MassScanMatrix object.

# Arguments
- `msm::MassScanMatrix`: Mass scan matrix object with retention and intensity values

See airpls(retentions, intensities; ...) for algorithm details and full parameter
documentation.

# Example
```julia
corrected_msm = airpls(msm; λ=1e6, threshold_factor=2.0)
```
"""
function airpls(
    msm::MassScanMatrix;
    improvement_threshold::T1=1e-6,
    λ::T2=1e7,
    max_iter::T3=10^4,
    no_improvement_limit::T4=10,
    threshold_factor::T5=1.96,
    zero_threshold::T6=1e-8,
    zero_weight::T7=0.01
    ) where {
        T1<:Real,
        T2<:Real,
        T3<:Integer,
        T4<:Integer,
        T5<:Real,
        T6<:Real,
        T7<:Real
    }

    # Initialize baseline for each scan
    baselines = zeros(Float64, scancount(msm), mzcount(msm))
    ints = rawintensities(msm)
    # Apply airPLS to each ion in RI space
    for i in axes(ints, 2)
        # Run airPLS in RI space with corrected intensities and weights
        mzints = @view ints[:, i]
        baselines[:, i] = airpls(
                        rawretentions(msm), 
                        mzints; 
                        variances=mzints,
                        improvement_threshold=improvement_threshold,
                        λ=λ,
                        max_iter=max_iter,
                        no_improvement_limit=no_improvement_limit,
                        threshold_factor=threshold_factor,
                        zero_threshold=zero_threshold,
                        zero_weight=zero_weight)
    end

    # Return MassScanMatrix
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

"""
    airpls(msm::MassScanMatrix, variances::AbstractMatrix{<:Real}; ...) -> MassScanMatrix

Estimate baseline from mass scan matrix data using airPLS algorithm. The method weights the 
PLS regression by the inverse of the variances.

Convenience method that extracts retention times and intensities of each mz channel from
the MassScanMatrix, applies baseline correction, and returns the identified baselines as a
MassScanMatrix object.

# Arguments
- `msm::MassScanMatrix::MassScanMatrix`: Mass scan matrix object with retention and 
  intensity values
- `variances::AbstractMatrix{<:Real}`: Variance matrix with dimensions (n_scans, n_mzs)
  containing the variances for each intensity value

See airpls(retentions, intensities; ...) for algorithm details and full parameter
documentation.

# Example
```julia
corrected_msm = airpls(msm, vars; λ=1e6, threshold_factor=2.0)
```
"""
function airpls(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real};
    improvement_threshold::T1=1e-6,
    λ::T2=1e7,
    max_iter::T3=10^4,
    no_improvement_limit::T4=10,
    threshold_factor::T5=1.96,
    zero_threshold::T6=1e-8,
    zero_weight::T7=0.01
    ) where {
        T1<:Real,
        T2<:Real,
        T3<:Integer,
        T4<:Integer,
        T5<:Real,
        T6<:Real,
        T7<:Real
    }

    # Validate input dimensions
    n_scans = scancount(msm)
    n_mzs = mzcount(msm)
    size(variances) == (n_scans, n_mzs) || throw(
        ArgumentError("variances must have dimensions ($n_scans, $n_mzs)"))

    # Initialize baseline for each scan
    baselines = zeros(Float64, scancount(msm), mzcount(msm))
    ints = rawintensities(msm)
    # Apply airPLS to each ion in RI space
    for i in axes(ints, 2)
        # Run airPLS in RI space with corrected intensities and weights
        mzints = @view ints[:, i]
        baselines[:, i] = airpls(
                        rawretentions(msm), 
                        mzints; 
                        variances=mzints,
                        improvement_threshold=improvement_threshold,
                        λ=λ,
                        max_iter=max_iter,
                        no_improvement_limit=no_improvement_limit,
                        threshold_factor=threshold_factor,
                        zero_threshold=zero_threshold,
                        zero_weight=zero_weight)
    end

    # Return MassScanMatrix
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

mutable struct FitTracker{T1<:Integer, T2<:Real}
    no_improvement_limit::T1
    improvement_threshold::T2
    best_fit::Float64
    best_baseline::Vector{Float64}
    no_improvement_count::Int
    iteration_of_best::Int
    
    function FitTracker{T1, T2}(
        n::Int,
        no_improvement_limit::T1,
        improvement_threshold::T2
    ) where {T1<:Integer, T2<:Real}

        @assert no_improvement_limit > 0 "no_improvement_limit must be > 0"
        @assert improvement_threshold ≥ 0 "improvement_threshold must be ≥ 0"
        
        new{T1, T2}(
            no_improvement_limit,
            improvement_threshold,
            Inf,               # Start with worst possible objective
            zeros(Float64, n), # Placeholder for best baseline
            0,                 # No iterations without improvement yet
            0                  # No best iteration yet
        )
    end
end

function FitTracker(n::Int, no_improvement_limit::T1, improvement_threshold::T2
    )where {T1<:Integer, T2<:Real}

    FitTracker{T1, T2}(n, no_improvement_limit, improvement_threshold)
end

function calculate_fit(baseline, intensities, weights, smoothness_matrix)::Float64
    # Weighted residual sum of squares
    residuals = intensities - baseline
    weighted_rss = sum(weights .* (residuals .^ 2))
    
    # Smoothness penalty: λ * ||D * baseline||²
    # Since smoothness_matrix = λ * D' * D, we have:
    # smoothness_penalty = baseline' * smoothness_matrix * baseline
    smoothness_penalty = dot(baseline, smoothness_matrix * baseline)

    weighted_rss + smoothness_penalty
end

function update!(tracker::FitTracker, current_fit, current_baseline, iteration)
    
    # Check if this is a significant improvement
    if tracker.best_fit - current_fit > tracker.improvement_threshold
        # New best found!
        tracker.best_fit = current_fit
        tracker.best_baseline .= current_baseline  # Copy baseline values
        tracker.no_improvement_count = 0
        tracker.iteration_of_best = iteration
        return true  # Improvement found
    else
        # No significant improvement
        tracker.no_improvement_count += 1
        return false  # No improvement
    end
end

stop_optimization(tracker::FitTracker) = 
    tracker.no_improvement_count ≥ tracker.no_improvement_limit

baseline(tracker::FitTracker) = tracker.best_baseline

function compute_weights(residuals, std_devs, threshold_factor, zero_mask, zero_threshold)

    n = length(residuals)
    
    normalized_residuals = isnothing(std_devs) ? residuals : residuals ./ std_devs
    
    # Use only negative residuals to estimate noise (standard airPLS approach)
    # Exclude zero points from noise estimation
    valid_for_noise = (.!zero_mask) .& (normalized_residuals .< 0)
    neg_residuals = normalized_residuals[valid_for_noise]
    
    # Protection for extreme cases
    length(neg_residuals) < 3 && return ones(n)
    
    # Robust noise estimation from one-sided distribution
    m = median(neg_residuals)  # Typical noise magnitude 
    σ = max(abs(m) * 1.4826, zero_threshold)  # Convert to Gaussian std dev scale
    
    # Apply airPLS asymmetric weighting logic
    weights = zeros(Float64, n)

    threshold = threshold_factor * σ
    for i in 1:n
        weights[i] = 1 / (1 + exp(clamp(2 * (normalized_residuals[i] - threshold) / σ, 
                                  -30.0, 30.0)))
    end
    expand_low_weights!(weights, residuals)
    weights
end

function expand_low_weights!(w, res)
    n = length(w)

    @assert length(res) == n "weights and residuals must have same length"
    
     # Expand minimal weights across contiguous peak shoulders
    idcs = sortperm(w)
    mark = falses(n)
    
    @inbounds for i in idcs
        mark[i] && continue
        
        # Move left while residuals are decreasing (following peak slope down)
        l = i; while l > 1 && res[l - 1] < res[l]; l -= 1; end
        r = i; while r < n && res[r + 1] < res[r]; r += 1; end
        @inbounds if r-l+1 ≥ 5
            for j in l:r
                if !mark[j]
                    w[j] = w[i]
                    mark[j]=true
                end
            end
        end
    end
    
    w
end

function build_second_derivative_matrix(x::AbstractVector{<:Real})
    n = length(x)
    row = Int[]; col = Int[]; val = Float64[]
    
    # Loop over internal points (ignoring endpoints)
    for i in 2:(n - 1)
        h₁ = x[i]     - x[i - 1]
        h₂ = x[i + 1] - x[i]
        sum_h = h₁ + h₂
        inv_h₁ = inv(h₁)
        inv_h₂ = inv(h₂)
        coeff_left = 2 * inv_h₁ / sum_h
        coeff_right = 2 * inv_h₂ / sum_h
        coeff_center = -(coeff_left + coeff_right)  # enforce ∑coeffs = 0 exactly
        
        # Finite difference coefficients for point i
        push!(row, i); push!(col, i - 1); push!(val, coeff_left)
        push!(row, i); push!(col, i);     push!(val, coeff_center)
        push!(row, i); push!(col, i + 1); push!(val, coeff_right)
    end
    
    sparse(row, col, val, n, n)
end
