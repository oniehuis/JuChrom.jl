
# ── whiten ────────────────────────────────────────────────────────────────────────────────

"""
    whiten(x::AbstractArray{<:Real}, σ²::AbstractArray{<:Real}, σ²_floor::Real)

Whiten x using its propagated variances.

Performs element-wise whitening by dividing each value by its corresponding
standard deviation (derived from the variances), with a floor applied to prevent division
by very small values.

# Arguments
- x::AbstractArray{<:Real}: Data to be whitened (e.g. CLR-transformed intensities).
- σ²::AbstractArray{<:Real}: Matching variances array (same size as x)
- σ²_floor::Real: Absolute variance floor (must be positive)

# Returns
- Whitened data

# Examples
```jldoctest
julia> x = [0.1, 0.7, 0.2]; σ² = [1, 9, 3];

julia> x_clr, σ²_clr = clr(x, σ²);

julia> x_whitened = whiten(x_clr, σ²_clr, 0.01);

julia> x_whitened ≈ [-0.11881290739752365, 0.20292400048800466, -0.0273596834199544]
true
```
"""
function whiten(x::AbstractArray{<:Real}, σ²::AbstractArray{<:Real}, σ²_floor::Real)
    # Ensure arrays have matching dimensions for element-wise operations
    size(x) == size(σ²) ||
        throw(DimensionMismatch("x and σ² must have identical sizes"))

    # Validate floor is positive to prevent invalid standard deviations
    σ²_floor > 0 || throw(ArgumentError("σ²_floor must be positive."))
    
    # Apply variance floor to prevent division by near-zero values
    σ²_safe = max.(σ², σ²_floor)

    # Convert variances to standard deviations and perform whitening
    x ./ sqrt.(σ²_safe)
end
