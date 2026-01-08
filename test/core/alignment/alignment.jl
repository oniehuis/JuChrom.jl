module TestAlignment

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core./alignment.jl
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

using Test

# ── Internal Project Imports ─────────────────────────────────────────────────

using JuChrom

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests
# ─────────────────────────────────────────────────────────────────────────────

# ── gapalign ─────────────────────────────────────────────────────────────────

# Define TraceAction enum here (or import if in your module)
@enum TraceAction Match SkipCompound SkipPosition

# Tests for gapalign
@testset "gapalign(S::Matrix{<:Real}, γ::Real, τ::Real)" begin

    # Simple example: identity matrix with threshold 0 and no gap penalty
    S = [1.0 0.0; 0.0 1.0]
    γ = 0.0
    τ = 0.0
    assignment, score = gapalign(S, γ, τ)
    
    @test assignment == Dict(1 => 1, 2 => 2)
    @test score ≈ 2.0
    
    # Test gap penalty effect
    S = [0.9 0.1; 0.1 0.9]
    γ = 0.5
    τ = 0.5
    assignment, score = gapalign(S, γ, τ)
    
    # Expected: only matches ≥ τ allowed; skips penalized by γ
    # Check if assigned keys and values are consistent and score reasonable
    @test all(k in keys(assignment) for k in 1:2)
    @test all(v in 1:2 for v in values(assignment))
    @test -2.0 ≤ score ≤ 2.0  # rough bounds
    
    # Check trace matrix content indirectly by inspecting assignments on small matrix
    S = [0.7 0.4; 0.4 0.7]
    γ = 1.0
    τ = 0.5
    assignment, score = gapalign(S, γ, τ)
    
    # Both diagonals are above threshold; assignment should prefer matches
    @test length(assignment) == 2
    
    # Test behavior when entire matrix is below threshold (forcing skips)
    S = fill(0.1, 3, 3)
    γ = 0.2
    τ = 0.5
    assignment, score = gapalign(S, γ, τ)
    
    @test isempty(assignment) || length(assignment) < 3  # likely all skips or partial
    
    # Edge case: empty matrix
    S = zeros(0, 0)
    γ = 1.0
    τ = 0.0
    assignment, score = gapalign(S, γ, τ)
    
    @test assignment == Dict()
    @test score == 0.0
end

end