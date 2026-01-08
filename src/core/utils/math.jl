# ── cosdis ───────────────────────────────────────────────────────────────────

"""
    cosdis(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, 
           positive_only::Bool=true) -> Float64

Compute the cosine distance between two real-valued vectors `x` and `y`.

Cosine distance is defined as:

    1 - cossim(x, y)

# Arguments
- `x`, `y`: Vectors of real numbers.
- `positive_only`: If `true` (default), clamps result to the range [0, 1], assuming 
  non-negative inputs (e.g., intensity spectra). If `false`, the range is [0, 2] to 
  account for possible negative similarity.

# Returns
- A `Float64` value representing the cosine distance.
- Returns `NaN` if either input vector is a zero vector.

# Examples
```jldoctest
julia> x = [1.0, 2.0, 3.0];

julia> y = [4.0, 5.0, 6.0];

julia> isapprox(cosdis(x, y), 0.025368153802923787; atol=1e-12)
true

julia> isapprox(cosdis(x, y, false), 0.025368153802923787; atol=1e-12)
true
```
"""
function cosdis(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real}, 
    positive_only::Bool = true)::Float64

    sim = cossim(x, y, positive_only)
    dist = 1.0 - sim

    positive_only ? clamp(dist, 0.0, 1.0) : clamp(dist, 0.0, 2.0)
end

# ── cossim ────────────────────────────────────────────────────────────────────────────────

"""
    cossim(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, 
           positive_only::Bool=true) -> Float64

Compute the cosine similarity between two real-valued vectors `x` and `y`.

Cosine similarity is defined as:

    dot(x, y) / (‖x‖ * ‖y‖)

# Arguments
- `x`, `y`: Vectors of real numbers.
- `positive_only`: If `true` (default), clamps similarity to the range [0, 1], assuming 
  non-negative input vectors (e.g., intensity spectra).  If `false`, similarity may range 
  in [-1, 1].

# Returns
- A `Float64` value representing the cosine similarity between `x` and `y`.
- Returns `NaN` if either vector is a zero vector (norm is zero).

# Examples
```jldoctest
julia> x = [1.0, 2.0, 3.0];

julia> y = [4.0, 5.0, 6.0];

julia> isapprox(cossim(x, y), 0.9746318461970762; atol=1e-12)
true

julia> isapprox(cossim(x, y, false), 0.9746318461970762; atol=1e-12)
true
```
"""
function cossim(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    positive_only::Bool=true)::Float64

    norm_x = norm(x, 2)
    norm_y = norm(y, 2)

    # Avoid division by zero
    if iszero(norm_x) || iszero(norm_y)
        return NaN  # max distance if one vector is zero
    end

    cos_sim = dot(x, y) / (norm_x * norm_y)

    positive_only ? clamp(cos_sim, 0.0, 1.0) : clamp(cos_sim, -1.0, 1.0)
end

# ── inverse_minmax_scale ──────────────────────────────────────────────────────────────────

"""
    inverse_minmax_scale(v_norm::Real, v_min::Real, v_max::Real) -> Real

Convert a min-max normalized value back to its original scale.

This function performs the inverse of min-max normalization, mapping a normalized value
`v_norm ∈ [0, 1]` back to its original scale using the provided minimum (`v_min`) and
maximum (`v_max`) values.

# Arguments
- `v_norm`: The normalized value (typically in the range [0, 1]).
- `v_min`: The minimum value of the original scale.
- `v_max`: The maximum value of the original scale.

# Return
- The corresponding value on the original scale as a `Real`.

See also: [`minmax_scale`](@ref)

# Example
```jldoctest
julia> JuChrom.inverse_minmax_scale(0.5, 10, 20)
15.0
julia> JuChrom.inverse_minmax_scale(0.0, 5, 15)
5.0
julia> JuChrom.inverse_minmax_scale(1.0, 5, 15)
15.0
```
"""
function inverse_minmax_scale(v_norm::Real, v_min::Real, v_max::Real)
    v_norm .* (v_max - v_min) .+ v_min
end

# ── localmaxima ───────────────────────────────────────────────────────────────────────────

"""
    localmaxima(values::AbstractVector{<:Real}) -> Vector{Int}

Return the indices of local maxima in a vector or 1D array of real numbers.

A local maximum is defined as an element that is strictly greater than its immediate 
neighbors. Only strictly greater maxima are detected; plateaus or flat maxima are not 
identified.

# Arguments
- `values`: A one-dimensional vector or array of real numbers.

# Return
- A vector of integer indices corresponding to positions in `values` where local maxima 
  occur.

# Notes
- At least three elements are required in `values` to detect maxima.
- Works with arrays having arbitrary indexing (e.g., `OffsetArrays`).
- Only strictly greater maxima are detected; flat or plateau maxima are ignored.

# Example
```jldoctest
julia> vals = [1, 3, 2, 4, 5, 1];

julia> JuChrom.localmaxima(vals) == [2, 5]
true

julia> JuChrom.localmaxima([1, 2, 3]) == []
true

julia> JuChrom.localmaxima([1, 5, 1]) == [2]
true

julia> JuChrom.localmaxima([1, 2, 2, 1]) == []
true
```
"""
function localmaxima(values::AbstractVector{<:Real})
    indices = Int[]
    length(values) < 3 && return indices
    for idx in Iterators.drop(eachindex(values), 2)
        if values[idx-2] < values[idx-1] > values[idx]
            push!(indices, idx-1)
        end
    end
    indices
end

# ── maxdecimals ───────────────────────────────────────────────────────────────────────────

"""
    maxdecimals(vec::AbstractVector{<:Real}) -> Int

Return the maximum number of decimal places present in any element of the input vector.

This function examines each element of `vec`, converts it to a string, and counts the 
number of digits after the decimal point (ignoring trailing zeros). It returns the largest 
such count found among all elements. If the vector is empty or contains only integers, 
returns `0`.

# Arguments
- `vec`: A vector of real numbers to analyze.

# Return
- `Int`: The maximum number of decimal digits present in any element of `vec`.

# Example
```jldoctest
julia> JuChrom.maxdecimals([1.23, 4.567, 8.9])
3
julia> JuChrom.maxdecimals([1.0, 2.0, 3.0])
0
julia> JuChrom.maxdecimals(Int[])
0
```
"""
function maxdecimals(vec::AbstractVector{<:Real})
    max_decimals = 0
    for x in vec
        # Convert to string with full precision
        s = string(x)
        if occursin('.', s)
            # split integer and decimal parts
            parts = split(s, '.')
            decimal_str = parts[2]
            # Remove trailing zeros from decimal part
            decimal_str = replace(decimal_str, r"0+$" => "")
            max_decimals = max(max_decimals, length(decimal_str))
        end
    end
    max_decimals
end

# ── minmax_scale ──────────────────────────────────────────────────────────────────────────

"""
    minmax_scale(v::Real, v_min::Real, v_max::Real) -> Real

Scale a value to the [0, 1] range using min-max normalization.

This function applies min-max normalization to a scalar value `v`, mapping it to the 
range [0, 1] based on the provided minimum (`v_min`) and maximum (`v_max`) values of the 
original scale.

# Arguments
- `v`: The value to normalize.
- `v_min`: The minimum value of the original scale.
- `v_max`: The maximum value of the original scale.

# Return
- The normalized value in the range [0, 1] as a `Real`.

See also: [`inverse_minmax_scale`](@ref)

# Example
```jldoctest
julia> JuChrom.minmax_scale(15, 10, 20)
0.5
julia> JuChrom.minmax_scale(5, 5, 15)
0.0
julia> JuChrom.minmax_scale(15, 5, 15)
1.0
```
"""
function minmax_scale(v::Real, v_min::Real, v_max::Real)
    (v .- v_min) ./ (v_max - v_min)
end
