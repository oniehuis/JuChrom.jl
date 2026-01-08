# ── clr ───────────────────────────────────────────────────────────────────────────────────

"""
    clr(x::AbstractArray{<:Real}) -> Array{Float64}
    clr(x::AbstractArray{<:Real},
        variances::AbstractArray{<:Real}) -> Tuple{Array{Float64}, Array{Float64}}

Apply the centered log-ratio (CLR) transformation to compositional data, with optional 
variance propagation.

The CLR transformation converts strictly positive data to log-ratio coordinates
by subtracting the geometric mean (in log space) from each log-transformed value:
`clr(x) = log(x) - mean(log(x))`.

# Single Array Usage
When called with one array, returns only the CLR-transformed data.

# Variance Propagation Usage  
When called with two arrays, the second array contains variances of the original data.
Returns both the CLR-transformed data and the propagated variances using exact 
finite-sample correction.

# Arguments
- `x::AbstractArray{<:Real}`: Input array with strictly positive values
- `variances::AbstractArray{<:Real}`: Element-wise variances of original data (same shape 
  as `x`)

# Returns
- Single array: `Array{Float64}` - CLR-transformed data with same shape as input
- Two arrays: `Tuple{Array{Float64}, Array{Float64}}` - (CLR-transformed data, propagated 
  variances)

# Throws
- `DomainError`: if any element of `x` is ≤ 0
- `DimensionMismatch`: if `x` and `variances` have different shapes

# Examples
```jldoctest
julia> x = [0.1, 0.7, 0.2]; σ² = [1, 9, 3];

julia> clr(x) ≈ [-0.8796857765384194, 1.0662243725168936, -0.18653859597847422]
true

julia> x_clr, σ²_clr = clr(x, σ²);

julia> x_clr == clr(x)
true

julia> σ²_clr ≈ [54.81859410430839, 27.60770975056689, 46.485260770975046]
true
```
"""
function clr(x::AbstractArray{<:Real})
    # Ensure all components are strictly positive (required for log transformation)
    all(>(0), x) || throw(DomainError(x, "components must be > 0"))

    # Apply log transformation to convert multiplicative to additive scale
    logx = log.(x)

    # Calculate geometric mean in log space (arithmetic mean of log values)
    μ̄ = mean(logx)

    # Center by subtracting geometric mean to create zero-sum constraint
    logx .- μ̄
end

function clr(x::AbstractArray{<:Real}, variances::AbstractArray{<:Real})
    # Ensure arrays have matching dimensions for element-wise variance propagation
    size(x) == size(variances) || throw(
        DimensionMismatch("x and variances must have the same shape"))

    # Validate all components are strictly positive for log transformation
    all(>(0), x) || throw(DomainError(x, "components must be > 0"))
    
    # Apply CLR transformation to the data
    clr_x = clr(x)
    
    # Convert variances to standard deviations after log transformation using Delta method
    σ_log = sqrt.(variances) ./ x

    # Store array length for finite-sample correction calculations
    N = length(x)

    # Convert standard deviations back to variances in log space
    σ2_log = σ_log .^ 2

    # Sum of all log-space variances (needed for CLR correction term)
    Σσ2 = sum(σ2_log)

    # Apply exact finite-sample correction for CLR variance propagation
    σ2_clr = @. σ2_log * (1 - 2/N) + Σσ2/N^2
    
    # Return both transformed data and propagated variances
    clr_x, σ2_clr
end
