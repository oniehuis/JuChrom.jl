@enum TraceAction Match SkipCompound SkipPosition

"""
    gapalign(S::Matrix{<:Real}, γ::Real, τ::Real) -> (Dict{Int, Int}, Real)

Find the optimal alignment between compounds and positions using dynamic programming 
with gap penalties and similarity thresholds.

This function solves a constrained assignment problem where compounds can be matched to 
positions based on similarity scores, with penalties for leaving items unmatched and 
minimum thresholds for valid matches.

# Algorithm
Uses dynamic programming to solve:
- **Objective**: Maximize total similarity score minus gap penalties
- **Constraints**: Only allow matches where similarity ≥ τ
- **Complexity**: O(mn) time and space

# Arguments
- `S::Matrix{<:Real}`: An `m × n` similarity matrix where `S[i,j]` represents the 
  similarity between compound `i` and position `j`. Higher values indicate better matches.
- `γ::Real`: Gap penalty (≥ 0) subtracted for each unmatched compound or position. 
  Higher values encourage more matches.
- `τ::Real`: Minimum similarity threshold for valid matches. Matches with `S[i,j] < τ` 
  are forbidden.

# Returns
- `assignment::Dict{Int, Int}`: Optimal compound-to-position mapping (1-indexed). 
  Keys are compound indices, values are position indices. Unmatched compounds are 
  not included.
- `score::Real`: Total alignment score = Σ(matched similarities) - γ × (unmatched items)

# Examples
```jldoctest
julia> S = [0.9 0.1; 0.2 0.8];  # Simple 2×2 case: perfect matches

julia> assignment, score = gapalign(S, 0.1, 0.5);

julia> assignment == Dict(1 => 1, 2 => 2)
true

julia> score ≈ 1.7  # 0.9 + 0.8 - 2*0.1 (no gaps)
true

julia> S = [0.9 0.1; 0.1 0.2];  # Poor second match

julia> assignment, score = gapalign(S, 0.3, 0.5); # With gap penalty: prefer partial matching

julia> assignment == Dict(1 => 1)  # Skip compound 2
true

julia> score ≈ 0.3  # 0.9 - 0.3 (skip compound 2) - 0.3 (skip position 2)
true

julia> S = [0.9 0.3; 0.2 0.8];  # Threshold filtering: block poor matches

julia> assignment, score = gapalign(S, 0.1, 0.4);  # τ=0.4 blocks S[1,2] and S[2,1]

julia> assignment == Dict(1 => 1, 2 => 2)  # Perfect matches
true
```
"""
function gapalign(S::Matrix{<:Real}, γ::Real, τ::Real)
    m, n = size(S)
    dp = fill(-Inf, m + 1, n + 1)
    trace = Matrix{TraceAction}(undef, m + 1, n + 1)

    # Initialization
    dp[1, 1] = 0.0
    for i in 2:m+1
        dp[i, 1] = dp[i-1, 1] - γ
        trace[i, 1] = SkipCompound
    end
    for j in 2:n+1
        dp[1, j] = dp[1, j-1] - γ
        trace[1, j] = SkipPosition
    end

    # Fill DP table with threshold
    @inbounds for i in 2:m+1
        for j in 2:n+1
            match_score = (S[i-1, j-1] ≥ τ) ? (dp[i-1, j-1] + S[i-1, j-1]) : -Inf
            dp[i, j] = match_score
            trace[i, j] = Match

            skip_compound_score = dp[i-1, j] - γ
            if skip_compound_score > dp[i, j]
                dp[i, j] = skip_compound_score
                trace[i, j] = SkipCompound
            end

            skip_position_score = dp[i, j-1] - γ
            if skip_position_score > dp[i, j]
                dp[i, j] = skip_position_score
                trace[i, j] = SkipPosition
            end
        end
    end

    # Backtracking
    i, j = m + 1, n + 1
    assignment = Dict{Int, Int}()

    while i > 1 || j > 1
        action = trace[i, j]
        if action == Match
            assignment[i - 1] = j - 1
            i -= 1
            j -= 1
        elseif action == SkipCompound
            i -= 1
        elseif action == SkipPosition
            j -= 1
        else
            break  # Safety net
        end
    end

    assignment, dp[m + 1, n + 1]
end
