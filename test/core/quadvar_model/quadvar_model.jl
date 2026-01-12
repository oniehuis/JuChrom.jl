module TestQuadVarModel

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core./quadvar_model.jl
# ─────────────────────────────────────────────────────────────────────────────

using Test
using JuChrom
using JuChrom: boundedlsβ₂, cholridge, clampvarparams, fitgains!, fitoffsets!, fitquadvar,
               makediff2, makeDTD, pooledfitstats, pooledresidacf, rescalevarparams, 
               residuals, updatesignal!, residautocorr
import Base: broadcastable
using LinearAlgebra
using Random
using SparseArrays
using Statistics
using Unitful

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests
# ─────────────────────────────────────────────────────────────────────────────

# ── boundedlsβ₂ ──────────────────────────────────────────────────────────────

@testset "boundedlsβ₂()" begin
    # Objective helper
    J(A, b, β) = sum(abs2, A * β .- b)

    # Case 1: Unconstrained solution is inside bounds → should match A\b
    A_free = [2.0 0.5; 0.5 1.0; 1.0 1.5]
    b_free = [1.0, 0.0, 2.0]
    β_uncon = A_free \ b_free
    lb_free = (minimum(β_uncon) - 1.0, minimum(β_uncon) - 1.0)
    β_free = boundedlsβ₂(A_free, b_free, lb_free)
    @test isapprox(collect(β_free), β_uncon; rtol=1e-12, atol=1e-12)

    # Case 2: One bound active (β₁ clipped to 0, β₂ free)
    A_one = [1.0 0.0; 0.0 1.0; 1.0 1.0]
    b_one = [-0.5, 0.8, 0.3]
    lb_one = (0.0, 0.0)
    β_one = boundedlsβ₂(A_one, b_one, lb_one)
    @test β_one[1] ≥ lb_one[1] - 1e-12
    @test J(A_one, b_one, collect(β_one)) ≤ J(A_one, b_one, [lb_one[1], β_one[2]]) + 1e-12

    # Case 3: Both bounds active (corner solution at lb)
    A_corner = [1.0 2.0; 2.0 4.0; 3.0 6.0]   # rank-1
    b_corner = [-10.0, -20.0, -30.0]
    lb_corner = (0.1, 0.2)
    β_corner = boundedlsβ₂(A_corner, b_corner, lb_corner)
    @test isapprox(β_corner[1], lb_corner[1]; atol=1e-12)
    @test isapprox(β_corner[2], lb_corner[2]; atol=1e-12)

    # Case 4: Feasibility slack allows slightly negative values
    A_slack = I(2)
    b_slack = [-1e-13, -1e-13]
    β_slack = boundedlsβ₂(Matrix(A_slack), b_slack, (0.0, 0.0); feas_tol=1e-12)
    @test β_slack[1] ≤ 0 && β_slack[1] ≥ -1e-12
    @test β_slack[2] ≤ 0 && β_slack[2] ≥ -1e-12
    @test J(Matrix(A_slack), b_slack, collect(β_slack)) ≤ J(Matrix(A_slack), b_slack, [0.0, 0.0]) + 1e-20

    # Case 5: Approximate optimality vs grid search
    Random.seed!(123)
    A_grid = randn(6,2)
    b_grid = randn(6)
    lb_grid = (-0.2, 0.1)
    β_grid = boundedlsβ₂(A_grid, b_grid, lb_grid)
    grid1 = range(lb_grid[1], stop=2.0, length=121)
    grid2 = range(lb_grid[2], stop=2.0, length=121)
    Jmin = Inf
    for x in grid1, y in grid2
        Jmin = min(Jmin, J(A_grid, b_grid, [x, y]))
    end
    @test J(A_grid, b_grid, collect(β_grid)) ≤ Jmin + 1e-6

    # Case 6: Return type is always Tuple{Float64,Float64}
    A_type = Float32[1 2; 3 4; 5 6]
    b_type = Float32[1, 0, -1]
    β_type = boundedlsβ₂(A_type, b_type, (0f0, 0f0))
    @test β_type isa Tuple{Float64,Float64}

    # Case 7: Rank-deficient but feasible (stability)
    A_rank = [1.0 1.0; 2.0 2.0; 3.0 3.0]
    b_rank = [1.0, 2.1, 3.2]
    lb_rank = (0.0, 0.0)
    β_rank = boundedlsβ₂(A_rank, b_rank, lb_rank)
    @test all(β_rank .≥ lb_rank .- 1e-12)

    # Case 8: Input validation (wrong shapes should throw)
    @test_throws AssertionError boundedlsβ₂(rand(3,3), rand(3), (0.0, 0.0))  # needs 2 columns
    @test_throws DimensionMismatch boundedlsβ₂(rand(3,2), rand(2), (0.0, 0.0))  # b length mismatch

    # Case 9: No feasible active-set pattern forces fallback to clamped A\b
    A_fallback = I(2)
    b_fallback = [-1.0, -2.0]
    β_fallback = boundedlsβ₂(Matrix(A_fallback), b_fallback, (0.0, 0.0); feas_tol=-1e-6)
    @test β_fallback == (0.0, 0.0)

    # Case 10: Force no feasible active-set with -Inf tolerance
    β_fallback_inf = boundedlsβ₂(Matrix(A_fallback), b_fallback, (0.0, 0.0); feas_tol=-Inf)
    @test β_fallback_inf == (0.0, 0.0)
end

# ── broadcastable(::QuadVarFit) ──────────────────────────────────────────────

@testset "broadcastable(::QuadVarFit)" begin
    # Minimal dummy QuadVarFit (small, consistent)
    function _mkfit()
        mz_ref = [100.0, 101.0]
        mzunit = nothing
        mz_idx = [1, 2]
        mz_vals = mz_ref
        batchcount = 1
        n_reps_per_batch = [2]
        n_scans_per_batch = [3]
        params = [QuadVarParams(1.0, 0.1, 0.2), QuadVarParams(1.1, 0.1, 0.2)]
        signal = [[fill(1.0, 3)] for _ in 1:2]
        offsets = [[zeros(3, 2)] for _ in 1:2]
        gains = [[ones(2)] for _ in 1:2]
        scale_c = [1.0, 1.1]
        acf = [0.2, 0.25]
        acf_lag = [1, 1]
        n_pairs = [10, 12]
        qc_z_rms = [1.0, 1.0]
        qc_cov68 = [0.68, 0.68]
        qc_cov95 = [0.95, 0.95]
        qc_s_min = [0.0, 0.0]
        qc_s_max = [10.0, 10.0]
        qc_nz = [5, 5]
        observed = [[zeros(3, 2) for _ in 1:2]]
        intensityunit = nothing

        QuadVarFit(mzunit, mz_ref, mz_idx, mz_vals, batchcount,
                   n_reps_per_batch, n_scans_per_batch,
                   params, intensityunit, signal, offsets, gains,
                   scale_c, acf, acf_lag, n_pairs,
                   qc_z_rms, qc_cov68, qc_cov95, qc_s_min, qc_s_max, qc_nz,
                   observed)
    end

    q = _mkfit()

    # Broadcastable wraps the object into a Ref (so it behaves like a scalar)
    b = broadcastable(q)
    @test b isa Ref
    @test b[] === q

    # Identity broadcasting should return the object itself (scalar-like behavior)
    @test identity.(q) === q

    # Scalar broadcasting with === should yield a scalar Bool
    @test (q .=== q) === true
end

# ── cholridge ────────────────────────────────────────────────────────────────

@testset "cholridge(w, DT_D, λ; jitter)" begin
    n = 8
    D   = makediff2(n)
    DTD = D' * D
    w   = fill(1.0, n)
    λ   = 0.2
    jitter = 1e-8

    # Sparse DT_D → CHOLMOD factor (or Cholesky if dense); accept either
    F = cholridge(w, DTD, λ; jitter=jitter)
    @test (F isa LinearAlgebra.Cholesky) || (F isa SparseArrays.CHOLMOD.Factor)

    # Reconstruct A (as dense) and verify A * (F \ b) ≈ b
    j = jitter * max(mean(w), eps(Float64))
    A = Symmetric(Diagonal(w .+ j) + λ * DTD, :L) |> Matrix
    b = randn(n)
    x = F \ b
    @test isapprox(A * x, b; rtol=1e-8, atol=1e-10)

    # λ = 0 → diagonal + jitter; solve is still consistent
    F0 = cholridge(w, DTD, 0.0; jitter=jitter)
    A0 = Symmetric(Diagonal(w .+ jitter * max(mean(w), eps(Float64))), :L) |> Matrix
    b0 = randn(n)
    x0 = F0 \ b0
    @test isapprox(A0 * x0, b0; rtol=1e-8, atol=1e-10)

    # All-zero weights: symmetric positive definite via jitter alone
    w0 = zeros(n)
    Fz = cholridge(w0, DTD, 0.0; jitter=1e-6)
    Az = Symmetric(Diagonal(fill(1e-6 * max(mean(w0), eps(Float64)), n)), :L) |> Matrix
    bz = randn(n)
    xz = Fz \ bz
    @test isapprox(Az * xz, bz; rtol=1e-10, atol=1e-12)

    # Sparse vs dense DT_D: compare solutions, not factor matrices
    F_sparse = cholridge(w, DTD, λ)
    F_dense  = cholridge(w, Matrix(DTD), λ)
    bd = randn(n)
    xs = F_sparse \ bd
    xd = F_dense \ bd
    @test isapprox(xs, xd; rtol=1e-10, atol=1e-12)

    # Size mismatch should error (length(w) ≠ size(DTD,1))
    @test_throws DimensionMismatch cholridge(ones(n-1), DTD, λ)
end

# ── _strip_mz_value ──────────────────────────────────────────────────────────

@testset "_strip_mz_value" begin
    @test JuChrom._strip_mz_value(100.0, nothing) === 100.0
    @test JuChrom._strip_mz_value(100.0, u"Th") === 100.0
    @test JuChrom._strip_mz_value(100.0u"Th", u"Th") == 100.0
    @test_throws ArgumentError JuChrom._strip_mz_value(100.0u"Th", nothing)
end

# ── clampvarparams ───────────────────────────────────────────────────────────

@testset "clampvarparams" begin
    # No clamping needed (all above mins) — result equals numerically
    p = QuadVarParams(1.0, 0.2, 0.3)
    q = @inferred clampvarparams(p; σ₀²_min=0.0, ϕ_min=0.0, κ_min=0.0)
    @test q == QuadVarParams(1.0, 0.2, 0.3)

    # Clamp each field individually
    p1 = QuadVarParams(-1.0, 0.2, 0.3)
    q1 = clampvarparams(p1; σ₀²_min=0.0)
    @test q1 == QuadVarParams(0.0, 0.2, 0.3)

    p2 = QuadVarParams(1.0, -0.2, 0.3)
    q2 = clampvarparams(p2; ϕ_min=0.0)
    @test q2 == QuadVarParams(1.0, 0.0, 0.3)

    p3 = QuadVarParams(1.0, 0.2, -0.3)
    q3 = clampvarparams(p3; κ_min=0.0)
    @test q3 == QuadVarParams(1.0, 0.2, 0.0)

    # Clamp multiple fields at once
    p4 = QuadVarParams(-1.0, -0.2, -0.3)
    q4 = clampvarparams(p4; σ₀²_min=0.1, ϕ_min=0.0, κ_min=0.25)
    @test q4 == QuadVarParams(0.1, 0.0, 0.25)

    # Boundary: exactly at the min → unchanged numerically
    p5 = QuadVarParams(0.1, 0.0, 0.25)
    q5 = clampvarparams(p5; σ₀²_min=0.1, ϕ_min=0.0, κ_min=0.25)
    @test q5 == p5

    # Type behavior: Float32 fields upcast to Float64 (mins are Float64)
    p32 = QuadVarParams(Float32(1), Float32(0.1), Float32(0.3))
    q32 = clampvarparams(p32; σ₀²_min=0.0, ϕ_min=0.2, κ_min=0.0)
    @test typeof(q32) ≡ QuadVarParams{Float64,Float64,Float64}
    @test isapprox(q32.σ₀², 1.0; rtol=0, atol=0)  # exact here
    @test isapprox(q32.ϕ, 0.2; rtol=0, atol=0)  # exact here
    @test isapprox(q32.κ, 0.3; rtol=0, atol=eps(Float32))  # tolerate tiny rounding

    # BigFloat fields remain BigFloat (promotion with Float64 → BigFloat)
    pbig = QuadVarParams(big"1.0", big"0.1", big"0.3")
    qbig = clampvarparams(pbig; σ₀²_min=0.5, ϕ_min=0.0, κ_min=0.25)
    @test typeof(qbig) ≡ QuadVarParams{BigFloat,BigFloat,BigFloat}
    @test isapprox(qbig.σ₀², big"1.0"; rtol=0, atol=eps(BigFloat))
    @test isapprox(qbig.ϕ, big"0.1"; rtol=0, atol=eps(BigFloat))
    @test isapprox(qbig.κ, big"0.3"; rtol=0, atol=eps(BigFloat))

    # Mixed field types are allowed and clamp independently
    pmix = QuadVarParams(1//2, 0.0, 0.1)  # Rational, Float64, Float64
    qmix = clampvarparams(pmix; σ₀²_min=0.6, ϕ_min=0.2, κ_min=0.05)
    @test qmix == QuadVarParams(0.6, 0.2, 0.1)

    # NaN handling: max(NaN, x) = NaN → NaN passes through (current behavior)
    pnan = QuadVarParams(NaN, 0.1, 0.2)
    qnan = clampvarparams(pnan; σ₀²_min=0.0)
    @test isnan(qnan.σ₀²) && qnan.ϕ == 0.1 && qnan.κ == 0.2

    # Unitful parameters: bounds without units inherit the parameter unit
    pu = QuadVarParams(1.0u"pA^2", 0.2u"pA", 0.3)
    qu = clampvarparams(pu; σ₀²_min=2.0, ϕ_min=0.1, κ_min=0.1)
    @test qu.σ₀² == 2.0u"pA^2"
    @test qu.ϕ == 0.2u"pA"   # already above floor
    @test qu.κ == 0.3        # already above floor

    # Mismatched units on bounds raise a DimensionError
    @test_throws Unitful.DimensionError clampvarparams(pu; σ₀²_min=1.0u"s")
end

@testset "QuadVarParams unit validation" begin
    @test QuadVarParams(0.1, 0.2, 0.0) isa QuadVarParams
    @test_throws Unitful.DimensionError QuadVarParams(0.1u"pA^2", 0.0, 0.0)
    @test_throws Unitful.DimensionError QuadVarParams(0.1, 0.0u"pA", 0.0)
    @test_throws Unitful.DimensionError QuadVarParams(0.1u"pA^2", 0.1u"s", 0.0)
    @test_throws Unitful.DimensionError QuadVarParams(0.1, 0.2, 1.0u"s")
end

# ── fitgains! ────────────────────────────────────────────────────────────────

@testset "fitgains!()" begin
    # Helper to synthesize perfectly consistent data:
    # Y = o + s * g_true'  (so the LS step can recover g_true exactly)
    function make_data(s::Vector{Float64}, g_true::Vector{Float64};
                       o::Matrix{Float64}=zeros(length(s), length(g_true)),
                       w::Vector{Float64}=ones(length(s)))
        n, r = length(s), length(g_true)
        @assert size(o) == (n, r)
        y = o .+ s * g_true'
        (y, o, copy(s), copy(w))
    end

    # Case 1: exact recovery (uniform weights, zero offsets) + mean-1 normalization
    s0 = collect(1.0:5.0)  # length N=5
    g_true = [0.8, 1.2, 1.5]  # length R=3
    Y, o, s, w = make_data(s0, g_true)
    g = zeros(length(g_true))
    fitgains!(g, Y, o, s, w)  # mutates g and s
    gm_true = mean(g_true)
    @test isapprox(g, g_true ./ gm_true; rtol=1e-12, atol=1e-12)
    @test isapprox(s, s0 .* gm_true; rtol=1e-12, atol=1e-12)

    # Case 2: nonzero offsets — they cancel correctly in the numerator
    o2 = fill(0.3, length(s0), length(g_true))
    y2, o2, s2, w2 = make_data(s0, g_true; o=o2)
    g2 = similar(g_true, Float64); fill!(g2, 0.0)
    fitgains!(g2, y2, o2, s2, w2)
    @test isapprox(g2, g_true ./ mean(g_true); rtol=1e-12, atol=1e-12)
    @test isapprox(s2, s0 .* mean(g_true); rtol=1e-12, atol=1e-12)

    # Case 3: non-uniform weights — still exact on noiseless model
    w3 = [1.0, 2.0, 0.5, 3.0, 4.0]
    y3, o3, s3, w3 = make_data(s0, g_true; w=w3)
    g3 = zeros(length(g_true))
    fitgains!(g3, y3, o3, s3, w3)
    @test isapprox(g3, g_true ./ mean(g_true); rtol=1e-12, atol=1e-12)
    @test isapprox(s3, s0 .* mean(g_true); rtol=1e-12, atol=1e-12)

    # Case 4: ridge_to_1 = 1 forces all gains to 1, s unchanged by normalization
    y4, o4, s4, w4 = make_data(s0, g_true)
    g4 = zeros(length(g_true))
    fitgains!(g4, y4, o4, s4, w4; ridge_to_1=1.0)
    @test all(isapprox.(g4, 1.0; atol=1e-12))
    @test isapprox(s4, s0; rtol=1e-12, atol=1e-12)  # gm == 1 → no rescaling

    # Case 5: partial ridge shrinks dispersion of gains (variance decreases)
    y5, o5, s5a, w5 = make_data(s0, g_true)
    g5a = zeros(length(g_true))
    fitgains!(g5a, y5, o5, s5a, w5; ridge_to_1=0.0)
    var_no_ridge = var(g5a)

    y5b, o5b, s5b, w5b = make_data(s0, g_true)
    g5b = zeros(length(g_true))
    fitgains!(g5b, y5b, o5b, s5b, w5b; ridge_to_1=0.5)
    var_with_ridge = var(g5b)

    @test var_with_ridge ≤ var_no_ridge + 1e-14
    @test isapprox(mean(g5b), 1.0; atol=1e-12)  # always renormalized

    # Case 6: degenerate denominator (s ≡ 0) → gains set to 1 and return
    s6 = zeros(length(s0))
    y6 = s6 * g_true'  # all zeros
    o6 = zeros(size(y6))
    w6 = ones(length(s6))
    g6 = fill(NaN, length(g_true))
    fitgains!(g6, y6, o6, s6, w6)
    @test all(isapprox.(g6, 1.0; atol=0))
    @test s6 == zeros(length(s0))  # s unchanged on early return

    # Case 7: degenerate denominator (w ≡ 0) → gains set to 1 and return
    s7 = copy(s0)
    w7 = zeros(length(s0))
    y7 = s7 * g_true'
    o7 = zeros(size(y7))
    g7 = fill(NaN, length(g_true))
    fitgains!(g7, y7, o7, s7, w7)
    @test all(isapprox.(g7, 1.0; atol=0))
    @test s7 == s0  # unchanged

    # Case 8: return value and in-place mutation
    y8, o8, s8, w8 = make_data(s0, g_true)
    g8 = zeros(length(g_true))
    ret = fitgains!(g8, y8, o8, s8, w8)
    @test ret ≡ nothing
    @test !all(iszero, g8)  # g8 mutated
    @test s8 != s0  # s8 rescaled (because mean(g_true) ≠ 1)
end

# ── fitoffsets! ──────────────────────────────────────────────────────────────

@testset "fitoffsets!()" begin
    # Synthesize consistent data:
    # Y = o_true + s * g_true'
    # Choose o_true with/without row-wise mean = 0 depending on the case.
    function make_data(s::Vector{Float64}, g_true::Vector{Float64};
                       o_true::Matrix{Float64}, w::Vector{Float64})
        n, r = length(s), length(g_true)
        @assert size(o_true) == (n, r)
        @assert length(w) == n
        Y = o_true .+ s * g_true'
        (Y, copy(o_true), copy(s), copy(w))
    end

    # A simple "solver" that matches the equation used in fitoffsets!:
    # (Diag(w)) * o_r = w .* (Y[:,r] - g[r]*s)  →  o_r = Y[:,r] - g[r]*s
    # Using F = Diagonal(w) makes the solution exact (before effect-coding).
    function factor_diagonal(w)
        Diagonal(w)  # supports "\" via LinearAlgebra
    end

    # A regularized variant that mimics A = Diag(w) + λ I (for predictable shrinkage)
    function factor_diag_plus_lambda(w, λ)
        Diagonal(@. w + λ)
    end

    # Exact recovery when offsets already have zero row-wise mean
    n, r = 6, 4
    s = collect(range(0.5, 1.5; length=n))
    g_true = [0.8, 1.2, 1.5, 0.9]
    o0 = [sin(0.3*i + 0.7*j) for i in 1:n, j in 1:r]
    o0 .-= mean(o0, dims=2)                     # enforce row-wise mean zero
    w = ones(n)
    y, o_true, s1, w1 = make_data(s, g_true; o_true=o0, w=w)
    f = factor_diagonal(w1)
    o = zeros(n, r)
    ret = fitoffsets!(o, y, g_true, s1, f, w1)
    @test ret ≡ nothing
    @test isapprox(o, o_true; rtol=1e-12, atol=1e-12)  # exact (effect-coding preserves o_true)

    # Offsets with nonzero row means → effect-coding should remove row means
    o_raw = [0.3 + 0.2 * cos(0.2 * i) + 0.1 * j for i in 1:n, j in 1:r]
    y2, o_raw_copy, s2, w2 = make_data(s, g_true; o_true=o_raw, w=w)
    f2 = factor_diagonal(w2)
    o2 = similar(o_raw)
    fitoffsets!(o2, y2, g_true, s2, f2, w2)
    o2_expected = o_raw .- mean(o_raw, dims=2)
    @test isapprox(o2, o2_expected; rtol=1e-12, atol=1e-12)

    # Non-uniform weights (still exact with F = Diag(w))
    w3 = [1.0, 0.5, 2.0, 3.0, 0.7, 4.0]
    o3_true = [0.1*sin(i) + 0.05*j for i in 1:n, j in 1:r]
    o3_true .-= mean(o3_true, dims=2)
    y3, o3_true_copy, s3, w3c = make_data(s, g_true; o_true=o3_true, w=w3)
    f3 = factor_diagonal(w3c)
    o3 = zeros(n, r)
    fitoffsets!(o3, y3, g_true, s3, f3, w3c)
    @test isapprox(o3, o3_true; rtol=1e-12, atol=1e-12)

    # Regularization (A = Diag(w) + λ I) induces shrinkage, then effect-coding
    λ = 0.5
    w4 = fill(1.0, n)
    o4_true = [0.25*sin(0.4*i) + 0.05*(j-2) for i in 1:n, j in 1:r]
    y4, _, s4, w4c = make_data(s, g_true; o_true=o4_true, w=w4)
    f4 = factor_diag_plus_lambda(w4c, λ)
    o4 = zeros(n, r)
    fitoffsets!(o4, y4, g_true, s4, f4, w4c)

    # Closed-form expected solution BEFORE effect-coding: (w .* (Y - g*s')) ./ (w + λ)
    # Then subtract row-wise mean to match the function's final step.
    o4_unscaled = similar(o4_true)
    @inbounds for r in 1:r
        rhs = @. w4c * (y4[:, r] - g_true[r] * s4)
        o4_unscaled[:, r] = rhs ./ (@. w4c + λ)
    end
    o4_expected = o4_unscaled .- mean(o4_unscaled, dims=2)
    @test isapprox(o4, o4_expected; rtol=1e-12, atol=1e-12)

    # Shapes and in-place behavior (no allocation semantics checked, just mutation)
    o5 = zeros(n, r)
    y5, _, s5, w5 = make_data(s, g_true; o_true=o0, w=ones(n))
    f5 = factor_diagonal(w5)
    fitoffsets!(o5, y5, g_true, s5, f5, w5)
    @test size(o5) == (n, r)
    @test !(o5 ≡ zeros(n, r))  # mutated

    # Robustness to column count (R=1)
    r1 = 1
    g1 = [1.3]
    o1_true = [0.2*sin(0.5*i) for i in 1:n, _ in 1:r1]
    y1, _, s1c, w1c = make_data(s, g1; o_true=o1_true, w=ones(n))
    f1 = factor_diagonal(w1c)
    o1 = zeros(n, r1)
    fitoffsets!(o1, y1, g1, s1c, f1, w1c)
    o1_expected = o1_true .- mean(o1_true, dims=2)
    @test isapprox(o1, o1_expected; rtol=1e-12, atol=1e-12)
end

# ── fitquadvar ───────────────────────────────────────────────────────────────

@testset "fitquadvar()" begin
    # Reproducible synthetic setup
    Random.seed!(42)

    # True quadratic variance parameters used to generate residuals
    ptrue = QuadVarParams(0.30, 0.25, 0.20)

    # Build two batches with different lengths and replicate counts
    t1, r1 = 120, 3
    t2, r2 = 80, 4

    # Smooth, positive signals with different shapes per batch
    s1 = 2.0 .+ 1.5 .* sin.(LinRange(0, 2π, t1))
    s2 = 1.0 .+ 0.8 .* cos.(LinRange(0, 1.5π, t2))

    # Helper: draw residual matrix with row-wise std = sqrt(varpred(s, ptrue))
    function sim_residuals(s::Vector{Float64}, r::Int; p::QuadVarParams)
        t = length(s)
        sigma = sqrt.(varpred(s, p; varfloor=0.0))
        z = randn(t, r)
        @inbounds for i in 1:t, j in 1:r
            z[i, j] *= sigma[i]
        end
        z
    end

    r1mat = sim_residuals(s1, r1; p=ptrue)
    r2mat = sim_residuals(s2, r2; p=ptrue)

    # Add a few outliers to stress robustness (heavy tails)
    r1mat[10, 1] *= 10
    r2mat[25, 3] *= 8

    s_list = [s1, s2]
    r_list = [r1mat, r2mat]

    # Convenience: pooled vectors for curve-comparison diagnostics
    sall  = vcat(s_list...)
    vtrue = varpred(sall, ptrue; varfloor=0.0)

    # Baseline fit (variance statistic = :var, no binning, no ridge)
    phat, chat = fitquadvar(s_list, r_list; var_stat=:var, bins=0)
    @test phat isa QuadVarParams
    @test isfinite(chat) && chat > 0

    vhat = varpred(sall, phat; varfloor=0.0)

    # Curve agreement: relative RMS error and correlation on the observed s-grid
    rel_rms = sqrt(mean(((vhat .- vtrue) ./ max.(vtrue, eps(Float64))).^2))
    @test rel_rms ≤ 0.20  # allow some slack due to robustness/outliers
    rho = cor(vhat, vtrue)
    @test rho ≥ 0.98

    # Using MAD-based per-row variance should give a very similar curve
    pmad, cmad = fitquadvar(s_list, r_list; var_stat=:mad2, bins=0)
    vmad = varpred(sall, pmad; varfloor=0.0)
    rel_rms_mad_vs_var = sqrt(mean(((vmad .- vhat) ./ max.(vhat, eps(Float64))).^2))
    @test rel_rms_mad_vs_var ≤ 0.45
    @test isfinite(cmad) && cmad > 0

    # Optional equal-frequency binning should preserve the curve up to small error
    pbin, cbin = fitquadvar(s_list, r_list; var_stat=:var, bins=10)
    vbin = varpred(sall, pbin; varfloor=0.0)
    rel_rms_bin_vs_var = sqrt(mean(((vbin .- vhat) ./ max.(vhat, eps(Float64))).^2))
    @test rel_rms_bin_vs_var ≤ 0.15
    @test isfinite(cbin) && cbin > 0

    # Lower bounds are enforced
    pbounded, _ = fitquadvar(s_list, r_list; σ₀²_min=0.05, ϕ_min=0.05, κ_min=0.02, var_stat=:var)
    @test pbounded.σ₀² ≥ 0.05 - 1e-12
    @test pbounded.ϕ ≥ 0.05 - 1e-12
    @test pbounded.κ ≥ 0.02 - 1e-12

    # Strong ridge shrinks ϕ and κ toward 0 in the scaled stage
    pnoridge, _ = fitquadvar(s_list, r_list; ϕ_ridge=0.0, κ_ridge=0.0)
    pridge, _ = fitquadvar(s_list, r_list; ϕ_ridge=1e2, κ_ridge=1e2)
    @test pridge.ϕ ≤ pnoridge.ϕ + 1e-12
    @test pridge.κ ≤ pnoridge.κ + 1e-12

    # Winsorization of v tampers extreme rows but keeps finite parameters
    pclip, cclip = fitquadvar(s_list, r_list; v_clip_q=0.80)
    @test all(isfinite.((pclip.σ₀², pclip.ϕ, pclip.κ, cclip)))

    # Rank-deficient design stress: collapse replicates to make collinearity in variance 
    # summaries still expect finite parameters and a reasonable curve
    s_rd = copy(s1)
    r_rd = [r1mat[:,1] r1mat[:,1] r1mat[:,1]]  # identical columns
    prd, _ = fitquadvar([s_rd], [r_rd])
    vrd = varpred(s_rd, prd; varfloor=0.0)
    @test all(isfinite, vrd)

    # Small problem sizes should still run
    s_small = 0.5 .+ 0.4 .* sin.(LinRange(0, π, 12))
    r_small = sim_residuals(s_small, 2; p=ptrue)
    psmall, csmall = fitquadvar([s_small], [r_small])
    @test isfinite(csmall) && csmall > 0
    @test all(isfinite, (psmall.σ₀², psmall.ϕ, psmall.κ))

    # Empty variance summary returns default params
    s_empty = [1.0, 2.0, 3.0]
    r_empty = zeros(0, 1)
    pempty, cempty = fitquadvar([s_empty], [r_empty])
    expected_c = max(quantile(s_empty, 0.97), 1e-12)
    @test pempty == QuadVarParams(1.0, 0.0, 0.0)
    @test isapprox(cempty, expected_c; rtol=0, atol=1e-12)

    s_empty2 = [0.5, 1.5]
    r_empty2 = zeros(0, 2)
    pempty2, cempty2 = fitquadvar([s_empty2], [r_empty2])
    expected_c2 = max(quantile(s_empty2, 0.97), 1e-12)
    @test pempty2 == QuadVarParams(1.0, 0.0, 0.0)
    @test isapprox(cempty2, expected_c2; rtol=0, atol=1e-12)

    # The fitted curve should be nonnegative on the observed grid
    @test all(varpred(sall, phat; varfloor=0.0) .≥ 0)

    # Increasing z_max should not break the fit and should keep curves close
    ploz, _ = fitquadvar(s_list, r_list; z_max=2.0)
    phiz, _ = fitquadvar(s_list, r_list; z_max=10.0)
    vloz = varpred(sall, ploz; varfloor=0.0)
    vhiz = varpred(sall, phiz; varfloor=0.0)
    rel_rms_z = sqrt(mean(((vloz .- vhiz) ./ max.(vhiz, eps(Float64))).^2))
    @test rel_rms_z ≤ 0.20
end

# ── fitquadvarmodel(Ms::Vector{<:AbstractMatrix}) ─────────────────────────────

@testset "fitquadvarmodel(Ms::Vector{<:AbstractMatrix})" begin

    # Helper: simulate heteroscedastic residuals with std = sqrt(varpred(s, p))
    function sim_residuals(s::Vector{Float64}, r::Int; p::QuadVarParams)
        t = length(s)
        σ = sqrt.(varpred(s, p; varfloor=0.0))
        Z = randn(t, r)
        @inbounds for i in 1:t, j in 1:r
            Z[i, j] *= σ[i]
        end
        Z
    end

    # Helper: generate a batch matrix Y = o + s * g' + ε
    function make_batch(s::Vector{Float64}, g::Vector{Float64};
                        o::Matrix{Float64}, p::QuadVarParams)
        t, r = length(s), length(g)
        @assert size(o) == (t, r)
        ε = sim_residuals(s, r; p=p)
        Float64.(o .+ s * g' .+ ε)
    end

    # Helper: signal roughness via 2nd-difference penalty
    roughness(s::Vector{Float64}) = begin
        D = makediff2(length(s))
        sum(abs2, D * s)
    end

    Random.seed!(2025)

    # Input validation
    Ms_bad_scans = [randn(2, 2)]  # <3 scans
    @test_throws ErrorException fitquadvarmodel(Ms_bad_scans)

    Ms_bad_reps = [randn(3, 1)]   # <2 replicates
    @test_throws ErrorException fitquadvarmodel(Ms_bad_reps)

    # Synthetic setup
    ptrue = QuadVarParams(0.30, 0.25, 0.20)
    t1, r1 = 120, 3
    t2, r2 =  90, 4

    s1 = 2.0 .+ 1.5 .* sin.(LinRange(0, 2π, t1))
    s2 = 1.0 .+ 0.8 .* cos.(LinRange(0, 1.7π, t2))
    g1 = [0.8, 1.2, 1.5]
    g2 = [0.9, 1.0, 1.1, 1.2]

    o1 = [0.15*sin(0.15*i) + 0.05*(j-2) for i in 1:t1, j in 1:r1]
    o1 .-= mean(o1, dims=2)
    o2 = [0.12*cos(0.12*i) - 0.03*(j-2) for i in 1:t2, j in 1:r2]
    o2 .-= mean(o2, dims=2)

    Y1 = make_batch(s1, g1; o=o1, p=ptrue)
    Y2 = make_batch(s2, g2; o=o2, p=ptrue)
    Ms = [Y1, Y2]

    # Return contract and types
    r = fitquadvarmodel(Ms)
    @test haskey(r, :params) && r.params isa QuadVarParams
    @test haskey(r, :s)      && r.s isa Vector{Vector{Float64}} && length(r.s) == 2
    @test haskey(r, :o)      && r.o isa Vector{Matrix{Float64}} && length(r.o) == 2
    @test haskey(r, :g)      && r.g isa Vector{Vector{Float64}} && length(r.g) == 2
    @test haskey(r, :c)      && r.c isa Float64 && isfinite(r.c) && r.c > 0
    @test haskey(r, :stats)  && r.stats isa NamedTuple
    @test haskey(r, :acf)    && r.acf isa Float64
    @test haskey(r, :acf_lag) && r.acf_lag isa Int
    @test haskey(r, :n_pairs) && r.n_pairs isa Int

    # At least one variance update (even if main loop exits early)
    r0 = fitquadvarmodel(Ms; maxiter=0)
    @test all(isfinite, (r0.params.σ₀², r0.params.ϕ, r0.params.κ, r0.c))

    # Basic statistical sanity
    st = r.stats
    @test all(haskey(st, k) for k in (:z_rms, :cov68, :cov95, :s_min, :s_max))
    @test isfinite(st.z_rms) && isfinite(st.cov68) && isfinite(st.cov95)
    @test st.s_min ≤ st.s_max

    # λ_signal smoothing check
    r_lo = fitquadvarmodel(Ms; λ_signal=0.0)
    r_hi = fitquadvarmodel(Ms; λ_signal=1.0)
    rough_lo = sum(roughness, r_lo.s)
    rough_hi = sum(roughness, r_hi.s)
    @test rough_hi ≤ rough_lo + 1e-8

    # gain_ridge effect
    r_nr = fitquadvarmodel(Ms; gain_ridge=0.0)
    r_rg = fitquadvarmodel(Ms; gain_ridge=0.5)
    gains_nr = vcat(r_nr.g...)
    gains_rg = vcat(r_rg.g...)
    @test var(gains_rg) ≤ var(gains_nr) + 1e-12
    @test isapprox(mean(gains_rg), 1.0; atol=1e-10)

    # Offsets are effect-coded (row means ≈ 0)
    for (Yo, oo) in zip(Ms, r.o)
        m = vec(mean(oo, dims=2))
        @test maximum(abs, m) ≤ 1e-8
        @test size(oo) == size(Yo)
    end

    # Bounds and ridge for variance parameters
    rb = fitquadvarmodel(Ms; σ₀²_min=0.05, ϕ_min=0.05, κ_min=0.02, ϕ_ridge=1e-2, 
                             κ_ridge=1e-2)
    @test rb.params.σ₀² ≥ 0.05 - 1e-12
    @test rb.params.ϕ   ≥ 0.05 - 1e-12
    @test rb.params.κ   ≥ 0.02 - 1e-12
    rnr = fitquadvarmodel(Ms; ϕ_ridge=0.0, κ_ridge=0.0)
    rrg = fitquadvarmodel(Ms; ϕ_ridge=1e2, κ_ridge=1e2)
    @test rrg.params.ϕ ≤ rnr.params.ϕ + 1e-12
    @test rrg.params.κ ≤ rnr.params.κ + 1e-12

    # Floors and robustness: force finite ACF by disabling ACF restriction
    rr = fitquadvarmodel(Ms; τ_floor=1e-2, qc_varfloor=1e-10, acf_restrict=:none)
    @test all(isfinite, (rr.params.σ₀², rr.params.ϕ, rr.params.κ, rr.c, rr.acf))
    @test all(isfinite, (rr.stats.z_rms, rr.stats.cov68, rr.stats.cov95, rr.stats.s_min, 
                         rr.stats.s_max))
    @test rr.n_pairs > 0

    # ACF wiring
    ra = fitquadvarmodel(Ms; ac_lag=2, acf_restrict=:none)
    @test ra.acf_lag == 2
    @test isfinite(ra.acf)
    @test ra.n_pairs > 0

    # Invariance to zero-mean replicate shifts (effect-coding)
    Ms_shift = deepcopy(Ms)

    # Per-batch replicate shifts with zero mean
    δ1 = [0.20, -0.20, 0.00]                    # mean == 0
    δ2 = [0.10, -0.10, 0.05, -0.05]             # mean == 0
    Ms_shift[1] .+= ones(size(Ms[1], 1)) * δ1'  # add column-wise constants
    Ms_shift[2] .+= ones(size(Ms[2], 1)) * δ2'
    rs = fitquadvarmodel(Ms_shift)

    # Signal and gains should remain essentially unchanged
    @test sum(norm, (rs.s .- r.s)) ≤ 1e-3
    @test norm(vcat(rs.g...) .- vcat(r.g...)) ≤ 1e-3

    # Small-problem path
    Ms_small = [randn(3,2), randn(4,2)]
    r_small = fitquadvarmodel(Ms_small)
    @test all(isfinite, (r_small.params.σ₀², r_small.params.ϕ, r_small.params.κ, r_small.c))
    @test length(r_small.s) == 2 && length(r_small.o) == 2 && length(r_small.g) == 2
end

# ── fitquadvarmodel(series_batches::Vector{<:Vector}) ─────────────────────────────────────

@testset "fitquadvarmodel(series_batches::Vector{<:Vector})" begin

    # Minimal mock to satisfy the interface used by fitquadvarmodel(series_batches)
    struct MockSeries{T<:Number}
        mz::Vector{T}        # reference m/z grid (unitless or unitful)
        M::Matrix{Float64}  # n_scans × n_mz intensities
    end

    JuChrom.uniquemzvalues(x::MockSeries) = x.mz
    JuChrom.mscanmatrix(x::MockSeries)    = x  # passthrough; rawintensities unwraps
    JuChrom.rawintensities(x::MockSeries) = x.M

    Random.seed!(2025)

    # Reference m/z grid and helpers
    mz_ref = collect(100.0:1.0:105.0)  # length 6
    n_mz   = length(mz_ref)

    # Build a batch with R replicates and T scans from a shared clean signal
    function _make_batch(T::Int, R::Int; σ=0.02)
        s = @. 1.0 + 0.2 * sin((1:T) / T * 2π)
        g = 0.9 .+ 0.2 .* rand(R)
        o = [0.05 .* sin.(0.1 .* (1:T)) for _ in 1:R] |> x -> hcat(x...)
        o .-= mean(o, dims=2)  # effect-coded offsets (row mean 0)
        Y = s * g' .+ o .+ σ .* randn(T, R)  # n_scans × R (per selected m/z later)
        
        # Expand to n_scans × n_mz by repeating columns with small jitter
        M = hcat([Y[:, mod1(j, R)] .+ 1e-4 .* randn(T) for j in 1:n_mz]...)
        reps = [MockSeries(mz_ref, M) for _ in 1:R]  # each "replicate series"
        (reps, T, R)
    end

    # Two batches with different sizes/replicate counts
    reps1, T1, R1 = _make_batch(80, 3)
    reps2, T2, R2 = _make_batch(60, 4)
    series_batches = [reps1, reps2]

    # Basic success path (default mzsel = all)
    q = fitquadvarmodel(series_batches; show_progress=false)
    @test q isa QuadVarFit
    @test q.mzunit === nothing
    @test length(q.mz_ref) == n_mz
    @test q.mzvalues == mz_ref
    @test q.batchcount == 2
    @test q.n_reps_per_batch == [R1, R2]
    @test q.n_scans_per_batch == [T1, T2]
    @test length(q.params) == n_mz
    @test length(q.signal) == n_mz && length(q.offsets) == n_mz && length(q.gains) == n_mz
    @test length(q.scale_c) == n_mz && length(q.acf) == n_mz && length(q.qc_z_rms) == n_mz

    # Observed contains raw data restricted to selected m/z (default = all)
    @test length(q.observed) == 2
    @test length.(q.observed) == [R1, R2]
    @test size(q.observed[1][1]) == (T1, n_mz)

    # mzsel by indices
    qidx = fitquadvarmodel(series_batches; mzsel=[2, 4, 6], show_progress=false)
    @test qidx.mz_idx == [2, 4, 6]
    @test qidx.mzvalues == mz_ref[[2, 4, 6]]
    @test length(qidx.params) == 3
    @test size(qidx.observed[2][1]) == (T2, 3)

    # show_progress logs first-selected m/z info
    @test_logs (:info, r"m/z") begin
        fitquadvarmodel(series_batches; mzsel=[1], show_progress=true)
    end

    @test_logs (:info, r"m/z") (:info, r"m/z") begin
        fitquadvarmodel(series_batches; mzsel=[1, 2], show_progress=true, progress_every=2)
    end

    # Observed really matches the selected columns from the raw data
    # Compare batch 1, replicate 2, column 2 vs 4 vs 6
    M12 = rawintensities(mscanmatrix(series_batches[1][2]))
    @test qidx.observed[1][2][:, 1] ≈ M12[:, 2]
    @test qidx.observed[1][2][:, 2] ≈ M12[:, 4]
    @test qidx.observed[1][2][:, 3] ≈ M12[:, 6]

    # mzsel by values (exact matches in grid)
    qval = fitquadvarmodel(series_batches; mzsel=[mz_ref[1], mz_ref[5]], 
        show_progress=false)
    @test qval.mz_idx == [1, 5]
    @test qval.mzvalues == [mz_ref[1], mz_ref[5]]
    @test_logs (:info, r"m/z") begin
        fitquadvarmodel(series_batches; mzsel=[mz_ref[1]], show_progress=true)
    end

    # Unitful m/z grids are preserved
    mz_ref_unit = mz_ref .* u"Th"
    reps1_unit = [MockSeries(mz_ref_unit, rep.M) for rep in reps1]
    reps2_unit = [MockSeries(mz_ref_unit, rep.M) for rep in reps2]
    series_batches_unit = [reps1_unit, reps2_unit]

    qunit = fitquadvarmodel(series_batches_unit; show_progress=false)
    @test qunit.mzunit == u"Th"
    @test qunit.mz_ref == mz_ref
    @test all(x -> !JuChrom.isunitful(x), qunit.mz_ref)
    @test qunit.mzvalues == mz_ref

    @test_logs (:info, r"m/z") begin
        fitquadvarmodel(series_batches_unit; mzsel=[mz_ref_unit[1]], show_progress=true)
    end

    @test_logs (:info, r"m/z") (:info, r"m/z") begin
        fitquadvarmodel(series_batches_unit; mzsel=[mz_ref_unit[1], mz_ref_unit[3]],
            show_progress=true, progress_every=2)
    end

    qunit_sel = fitquadvarmodel(series_batches_unit; mzsel=[mz_ref_unit[2], mz_ref_unit[4]], show_progress=false)
    @test qunit_sel.mzunit == u"Th"
    @test qunit_sel.mzvalues == [mz_ref[2], mz_ref[4]]

    # Forwarding of kwargs to inner solver: change λ_signal to smoother pooled signal
    function _roughness_sum(qv::QuadVarFit)
        sum(begin
            # D for each batch length
            D = JuChrom.makediff2(length(s)); sum(abs2, D * s)
        end for s in vcat(qv.signal...))
    end
    qa = fitquadvarmodel(series_batches; λ_signal=0.0, show_progress=false)
    qb = fitquadvarmodel(series_batches; λ_signal=1.0, show_progress=false)
    @test _roughness_sum(qb) ≤ _roughness_sum(qa) + 1e-6

    # Error: no batches provided
    @test_throws ArgumentError fitquadvarmodel(Vector{Vector{MockSeries}}(); 
        show_progress=false)

    # Error: progress_every must be ≥ 1
    @test_throws ArgumentError fitquadvarmodel(series_batches; progress_every=0, 
        show_progress=false)

    # Error: batch with < 2 replicates
    bad_one_rep = [[MockSeries(mz_ref, randn(T1, n_mz))]]
    @test_throws ArgumentError fitquadvarmodel(bad_one_rep; show_progress=false)

    # Error: m/z grid mismatch across replicates
    mz_bad = copy(mz_ref); mz_bad[end] += 0.5
    reps_bad = deepcopy(reps1); reps_bad[2] = MockSeries(mz_bad, reps_bad[2].M)
    @test_throws ArgumentError fitquadvarmodel([reps_bad]; show_progress=false)

    # Error: scan count mismatch across replicates
    reps_mixed = deepcopy(reps2)
    reps_mixed[1] = MockSeries(mz_ref, randn(T2 + 1, n_mz))
    @test_throws ArgumentError fitquadvarmodel([reps_mixed]; show_progress=false)

    # Error: mzsel out of bounds (indices)
    @test_throws ArgumentError fitquadvarmodel(series_batches; mzsel=[0, 1], 
        show_progress=false)
    @test_throws ArgumentError fitquadvarmodel(series_batches; mzsel=[n_mz + 1], 
        show_progress=false)

    # Error: mzsel with values not in grid
    @test_throws ArgumentError fitquadvarmodel(series_batches; mzsel=[999.0], 
        show_progress=false)

    # Error: unsupported mzsel type
    @test_throws ArgumentError fitquadvarmodel(series_batches; mzsel="bad", 
        show_progress=false)
    @test_throws ArgumentError fitquadvarmodel(series_batches; mzsel=Set([1, 2]), 
        show_progress=false)
end

# ── makediff2 ────────────────────────────────────────────────────────────────

@testset "makediff2" begin
    n = 8
    D = makediff2(n)

    # Basic structure checks
    @test D isa SparseMatrixCSC{Float64, Int}  # should return sparse Float64
    @test size(D) == (n - 2, n)  # (n-2) × n size
    @test nnz(D) == 3 * (n - 2)  # 3 nonzeros per row

    # Row pattern: [1, -2, 1] stencil
    for i in 1:(n-2)
        cols = findall(!iszero, D[i, :])  # nonzero column indices
        @test length(cols) == 3  # exactly 3 per row
        @test cols == [i, i + 1, i + 2]  # located at i, i+1, i+2
        @test D[i, i]   == 1.0
        @test D[i, i + 1] == -2.0
        @test D[i, i + 2] == 1.0
    end

    # Invalid n should throw
    @test_throws ArgumentError makediff2(2)
    @test_throws ArgumentError makediff2(0)

    # -Other element types
    D32 = makediff2(7, Float32)
    @test D32 isa SparseMatrixCSC{Float32, Int}
    allowed32 = Set(Float32[1, -2])
    @test all(x -> x in allowed32, nonzeros(D32))  # all nonzeros are 1 or -2

    Db = makediff2(6, BigFloat)
    @test Db isa SparseMatrixCSC{BigFloat, Int}
    @test eltype(nonzeros(Db)) ≡ BigFloat  # coefficients promoted

    # Functional behavior checks
    c = fill(3.14, n)  # constant vector
    @test D * c ≈ zeros(n - 2)  # 2nd diff of constant = 0

    a, b = 2.5, -1.1
    lin = @. a * (1:n) + b  # linear function
    @test D * lin ≈ zeros(n - 2)  # 2nd diff of linear = 0

    aq, bq, cq = 0.75, -0.3, 10.0
    quad = @. aq * (1:n)^2 + bq * (1:n) + cq  # quadratic function
    @test D * quad ≈ fill(2 * aq, n - 2)  # 2nd diff of quad = 2a

    # Laplacian check: D'D should be positive semi-definite
    L = Symmetric(D' * D, :L)
    vals = eigvals(Matrix(L))
    @test minimum(vals) ≥ -1e-12  # eigenvalues nonnegative (up to tolerance)
end

# ── makeDTD ──────────────────────────────────────────────────────────────────

@testset "makeDTD" begin
    # Basic structure and types (default Float64)
    n = 10
    DTD = makeDTD(n)  # = (makediff2(n))' * makediff2(n)
    @test DTD isa SparseMatrixCSC{Float64, Int}
    @test size(DTD) == (n, n)
    @test issymmetric(Matrix(DTD))

    # Positive semi-definite property (nonnegative eigenvalues up to numerical noise)
    vals = eigvals(Matrix(Symmetric(DTD, :L)))
    @test minimum(vals) ≥ -1e-12

    # Nullspace contains constants and linear trends (D annihilates them)
    # Because D is a second-difference, D' D is PSD and has nullspace span{1, i}
    v_const = ones(n)
    v_lin = collect(1:n)
    @test norm(DTD * v_const) ≤ 1e-12
    @test norm(DTD * v_lin) ≤ 1e-12

    # Interior stencil check: for 3 ≤ k ≤ n-2, the row has 5-point pattern [1, -4, 6, -4, 1]
    # D is [1, -2, 1] on columns i..i+2, so D' D is the autocorrelation: [-2..2] => [1, -4, 6, -4, 1]
    for k in 3:(n - 2)
        @test DTD[k, k] == 6.0
        @test DTD[k, k - 1] == -4.0
        @test DTD[k, k + 1] == -4.0
        @test DTD[k, k - 2] == 1.0
        @test DTD[k, k + 2] == 1.0
        # zero outside the 5-band (cheap spot check)
        if k - 3 ≥ 1
            @test DTD[k, k - 3] == 0.0
        end
        if k + 3 ≤ n
            @test DTD[k, k + 3] == 0.0
        end
    end

    # Edge rows have truncated patterns (don’t assert exact numbers there, just bandedness)
    # First two and last two rows should have nonzeros only within distance ≤ 2
    for k in (1, 2, n - 1, n)
        nzcols = findall(!iszero, DTD[k, :])
        @test all(abs.(nzcols .- k) .≤ 2)
    end

    # Float32 element type path
    n32 = 9
    DTD32 = makeDTD(n32, Float32)
    @test DTD32 isa SparseMatrixCSC{Float32, Int}
    @test size(DTD32) == (n32, n32)
    @test DTD32[5, 5] == Float32(6)
    @test DTD32[5, 4] == Float32(-4)
    @test DTD32[5, 6] == Float32(-4)
    @test DTD32[5, 3] == Float32(1)
    @test DTD32[5, 7] == Float32(1)

    # BigFloat element type path
    nb = 8
    DTDb = makeDTD(nb, BigFloat)
    @test DTDb isa SparseMatrixCSC{BigFloat, Int}
    @test DTDb[4, 4] == BigFloat(6)
    @test DTDb[4, 3] == BigFloat(-4)
    @test DTDb[4, 5] == BigFloat(-4)
    @test DTDb[4, 2] == BigFloat(1)
    @test DTDb[4, 6] == BigFloat(1)

    # Error propagation for small n (makediff2 throws; makeDTD should too)
    @test_throws ArgumentError makeDTD(2)
    @test_throws ArgumentError makeDTD(0)

    # Adding tiny diagonal jitter should make it symmetric positive definite (strictly positive eigs)
    ϵ = 1e-8
    A_spd = Matrix(DTD) + ϵ * I
    vals_spd = eigvals(Symmetric(A_spd, :L))
    @test minimum(vals_spd) > 0
end

# ── pooledfitstats ───────────────────────────────────────────────────────────

@testset "pooledfitstats()" begin

    # Reproducible synthetic setup
    Random.seed!(2025)

    # True variance parameters and helper
    p_true = QuadVarParams(0.10, 0.20, 0.15)
    varσ(s) = sqrt.(varpred(s, p_true; varfloor=0.0))

    # Two batches with different sizes/replicates
    t1, r1 = 150, 3
    t2, r2 = 90, 4

    # Smooth, positive signals per batch
    s1 = 1.0 .+ 0.7 .* sin.(LinRange(0, 2π, t1))
    s2 = 0.8 .+ 0.5 .* cos.(LinRange(0, 1.5π, t2))

    # Per-replicate gains and simple offsets (vary across replicates)
    g1 = [0.9, 1.0, 1.1]
    g2 = [0.95, 1.0, 1.05, 1.10]
    o1 = [fill(-0.05, t1) fill(0.0, t1) fill(0.03, t1)]
    o2 = [fill(0.02, t2) fill(0.0, t2) fill(-0.01, t2) fill(0.04, t2)]

    # Generate residuals with row-wise std = sqrt(σ²(s))
    r1 = randn(t1, r1) .* varσ(s1)
    r2 = randn(t2, r2) .* varσ(s2)

    # Build observations: Y = o + g .* s + residual
    y1 = o1 .+ s1 * g1' .+ r1
    y2 = o2 .+ s2 * g2' .+ r2

    Ms = [y1, y2]
    o = [o1, o2]
    g = [g1, g2]
    s = [s1, s2]

    # Basic sanity and type checks
    stats = pooledfitstats(Ms, o, g, s, p_true; qc_varfloor=1e-12)
    @test stats isa NamedTuple
    @test keys(stats) == (:z_rms, :cov68, :cov95, :s_min, :s_max, :nz)
    @test stats.nz > 0
    @test isfinite(stats.s_min) && isfinite(stats.s_max)
    @test stats.s_min ≈ minimum(vcat(s...)) atol=0 rtol=0
    @test stats.s_max ≈ maximum(vcat(s...)) atol=0 rtol=0

    # Distributional checks on pooled z-scores
    # With correct variance, pooled RMS ≈ 1 and 68/95% coverage near nominal.
    @test isapprox(stats.z_rms, 1.0; atol=0.08)  # allow some Monte Carlo spread
    @test isapprox(stats.cov68, 0.6827; atol=0.05)
    @test isapprox(stats.cov95, 0.95;   atol=0.03)

    # Robustness to non-finite entries: inject some NaNs/Infs and ensure they are ignored
    y1_bad = copy(y1)
    y1_bad[5, 2]  = NaN
    y1_bad[20, 1] = Inf
    Ms_bad = [y1_bad, y2]
    stats_bad = pooledfitstats(Ms_bad, o, g, s, p_true; qc_varfloor=1e-12)
    @test stats_bad.nz == stats.nz - 2  # two corrupted entries removed
    @test isfinite(stats_bad.z_rms)
    @test 0.0 ≤ stats_bad.cov68 ≤ 1.0 && 0.0 ≤ stats_bad.cov95 ≤ 1.0

    # Effect of qc_varfloor: forcing a larger floor should not produce non-finite values
    stats_floor = pooledfitstats(Ms, o, g, s, p_true; qc_varfloor=1e-4)
    @test all(isfinite, (stats_floor.z_rms, stats_floor.cov68, stats_floor.cov95))

    # Edge case: no finite z-scores → nz = 0, z_rms = NaN, covariances = 0
    t_small, r_small = 10, 2
    y_allnan = fill(NaN, t_small, r_small)
    ms_empty = [y_allnan]
    o_empty = [zeros(t_small, r_small)]
    g_empty = [ones(r_small)]
    s_empty = [ones(t_small)]
    stats_empty = pooledfitstats(ms_empty, o_empty, g_empty, s_empty, p_true)
    @test stats_empty.nz == 0
    @test isnan(stats_empty.z_rms)
    @test stats_empty.cov68 == 0.0
    @test stats_empty.cov95 == 0.0
    @test isfinite(stats_empty.s_min) && isfinite(stats_empty.s_max)
end

# ── pooledresidacf ───────────────────────────────────────────────────────────

@testset "pooledresidacf()" begin

    # Reproducible synthetic setup
    Random.seed!(7)

    # Variance model used to standardize residuals inside pooledresidacf
    p_true = QuadVarParams(0.05, 0.15, 0.10)
    varσ(s) = sqrt.(varpred(s, p_true; varfloor=0.0))

    # Helper: generate R (T×R) with AR(1) structure per column and unit marginal var
    function ar1_errors(t::Int, r::Int, ρ::Float64)
        e = zeros(t, r)
        ϵ = randn(t, r)
        α = sqrt(1 - ρ^2)
        @inbounds for j in 1:r
            e[1, j] = randn()
            for i in 2:t
                e[i, j] = ρ * e[i-1, j] + α * ϵ[i, j]
            end
        end
        e
    end

    # Two batches with different sizes and replicate counts
    t1, r1 = 200, 3
    t2, r2 = 160, 4

    # Smooth positive signals
    s1 = 1.0 .+ 0.8 .* sin.(LinRange(0, 2π, t1))
    s2 = 0.7 .+ 0.6 .* cos.(LinRange(0, 1.5π, t2))

    # Per-replicate gains and offsets
    g1 = [0.95, 1.0, 1.05]
    g2 = [0.9, 1.0, 1.05, 1.1]
    o1 = [fill(-0.03, t1) fill(0.0, t1) fill(0.02, t1)]
    o2 = [fill(0.01, t2) fill(0.0, t2) fill(-0.02, t2) fill(0.03, t2)]

    # Target AR(1) autocorrelation
    ρ_true = 0.40

    # Build residuals with heteroscedastic scaling by sqrt(σ²(s))
    r1 = ar1_errors(t1, r1, ρ_true) .* varσ(s1)
    r2 = ar1_errors(t2, r2, ρ_true) .* varσ(s2)

    # Observations Y = o + s*g' + residual
    y1 = o1 .+ s1 * g1' .+ r1
    y2 = o2 .+ s2 * g2' .+ r2

    y = [y1, y2]
    o = [o1, o2]
    g = [g1, g2]
    s = [s1, s2]

    # Baseline: lag=1, no restriction → acf close to ρ_true; positive #pairs
    acf_base = pooledresidacf(y, o, g, s, p_true; lag=1, restrict=:none)
    @test isfinite(acf_base.acf)
    @test acf_base.n_pairs > 0
    @test isapprox(acf_base.acf, ρ_true; atol=0.05)

    # Lag=2 on AR(1) → approximately ρ_true^2
    acf_lag2 = pooledresidacf(y, o, g, s, p_true; lag=2, restrict=:none)
    @test isfinite(acf_lag2.acf)
    @test acf_lag2.n_pairs > 0
    @test isapprox(acf_lag2.acf, ρ_true^2; atol=0.06)

    # Low-intensity restriction: fewer pairs, still finite ACF
    acf_low = pooledresidacf(y, o, g, s, p_true; restrict=:low, z_quant=0.90, z_max_keep=0.30)
    @test isfinite(acf_low.acf)
    @test 0 < acf_low.n_pairs < acf_base.n_pairs

    # Flat restriction: fewer pairs, still finite ACF
    acf_flat = pooledresidacf(y, o, g, s, p_true; restrict=:flat, flat_q=0.20)
    @test isfinite(acf_flat.acf)
    @test 0 < acf_flat.n_pairs < acf_base.n_pairs

    # Combined low+flat: typically the fewest pairs
    acf_lowflat = pooledresidacf(y, o, g, s, p_true; restrict=:lowflat, z_quant=0.90, 
                                 z_max_keep=0.30, flat_q=0.20)
    @test isfinite(acf_lowflat.acf)
    @test 0 < acf_lowflat.n_pairs ≤ min(acf_low.n_pairs, acf_flat.n_pairs)

    # Trimming: should not explode and usually changes ACF only slightly
    acf_trim0 = pooledresidacf(y, o, g, s, p_true; lag=1, restrict=:none, trimfrac=0.0)
    acf_trim10 = pooledresidacf(y, o, g, s, p_true; lag=1, restrict=:none, trimfrac=0.10)
    @test isfinite(acf_trim10.acf)
    @test 0 ≤ acf_trim10.n_pairs ≤ acf_trim0.n_pairs
    @test abs(acf_trim10.acf - acf_trim0.acf) ≤ 0.15

    # Robustness to NaNs/Infs in Y → fewer pairs, finite ACF
    y1_bad = copy(y1)
    y1_bad[5, 2]  = NaN
    y1_bad[30, 1] = Inf
    acf_bad = pooledresidacf([y1_bad, y2], o, g, s, p_true; lag=1)
    @test isfinite(acf_bad.acf)
    @test acf_bad.n_pairs == acf_base.n_pairs - 2 * 1  # each bad entry kills one pair at lag=1

    # Degenerate: too-large lag leaves no pairs → acf = NaN, n_pairs = 0
    t_small, r_small = 5, 2
    s_small = 1 .+ 0.1 .* sin.(LinRange(0, π, t_small))
    y_small = randn(t_small, r_small) .* varσ(s_small)
    acf_none = pooledresidacf([y_small], [zeros(t_small, r_small)],
                              [ones(r_small)], [s_small], p_true; lag=10)
    @test isnan(acf_none.acf)
    @test acf_none.n_pairs == 0
end

# ── rescalevarparams ─────────────────────────────────────────────────────────

@testset "rescalevarparams()" begin

    # Reproducible synthetic setup
    Random.seed!(7)

    # True params only for simulating residuals (not used by rescale itself)
    p_true = QuadVarParams(0.30, 0.25, 0.20)

    # Two batches with different sizes/replicates
    t1, r1 = 120, 3
    t2, r2 = 90, 4
    s1 = 1.2 .+ 0.9 .* sin.(LinRange(0, 2π, t1))
    s2 = 2.0 .+ 1.1 .* cos.(LinRange(0, 1.7π, t2))

    # Draw residuals with row std = sqrt(varpred(s, p_true))
    function sim_residuals(s::Vector{Float64}, r::Int; p::QuadVarParams)
        t = length(s)
        σ = sqrt.(varpred(s, p; varfloor=0.0))
        z = randn(t, r)
        @inbounds for i in 1:t, j in 1:r
            z[i, j] *= σ[i]
        end
        z
    end
    r1mat = sim_residuals(s1, r1; p=p_true)
    r2mat = sim_residuals(s2, r2; p=p_true)

    s_list = [s1, s2]
    r_list = [r1mat, r2mat]

    # helper: compute expected β = clamp(1 / median(zr)^2, clip_lo, clip_hi)
    # using the same definition of per-row z_rms as the implementation
    function expected_beta(p::QuadVarParams, s_list, r_list;
                           varfloor=1e-12, trim_lo=0.20, trim_hi=0.80,
                           clip_lo=0.7, clip_hi=1.4)
        zr = Float64[]
        @inbounds for (s, r) in zip(s_list, r_list)
            v  = varpred(s, p; varfloor=varfloor)  # length T
            r2 = vec(sum(abs2, r; dims=2))  # length T
            z_row = sqrt.(r2 ./ size(r, 2)) ./ sqrt.(v)  # per-row z_rms
            append!(zr, z_row)
        end
        isempty(zr) && return 1.0
        sort!(zr)
        n = length(zr)
        i1 = clamp(floor(Int, trim_lo * n), 1, n)
        i2 = clamp(ceil(Int,  trim_hi * n), 1, n)
        i1, i2 = min(i1, i2), max(i1, i2)
        m = median(@view zr[i1:i2])
        (isfinite(m) && m > 0) || return 1.0
        clamp(1.0 / (m * m), clip_lo, clip_hi)
    end

    # Baseline: default trimming and clips
    p_start1 = QuadVarParams(0.40, 0.28, 0.21)
    p1 = rescalevarparams(p_start1, s_list, r_list)
    β1_exp = expected_beta(p_start1, s_list, r_list)
    @test isapprox([p1.σ₀², p1.ϕ, p1.κ], β1_exp .* [p_start1.σ₀², p_start1.ϕ, p_start1.κ];
                   rtol=1e-10, atol=1e-12)

    # Ratios between parameters are preserved (shape invariance)
    @test isapprox(p1.σ₀² / p1.ϕ, p_start1.σ₀² / p_start1.ϕ; rtol=1e-12, atol=0)
    @test isapprox(p1.σ₀² / p1.κ, p_start1.σ₀² / p_start1.κ; rtol=1e-12, atol=0)

    # Clip high: choose large starting variances so z is small → β gets clipped to clip_hi
    p_start_big = QuadVarParams(10.0, 8.0, 5.0)
    p_big = rescalevarparams(p_start_big, s_list, r_list; clip_lo=0.7, clip_hi=1.4)
    β_big_exp = expected_beta(p_start_big, s_list, r_list; clip_lo=0.7, clip_hi=1.4)
    @test β_big_exp == 1.4
    @test isapprox([p_big.σ₀², p_big.ϕ, p_big.κ], 1.4 .* [p_start_big.σ₀², p_start_big.ϕ, 
                    p_start_big.κ]; rtol=1e-12, atol=1e-12)

    # Clip low: choose tiny starting variances so z is large → β gets clipped to clip_lo
    p_start_small = QuadVarParams(1e-6, 5e-7, 1e-7)
    p_small = rescalevarparams(p_start_small, s_list, r_list; clip_lo=0.7, clip_hi=1.4)
    β_small_exp = expected_beta(p_start_small, s_list, r_list; clip_lo=0.7, clip_hi=1.4)
    @test β_small_exp == 0.7
    @test isapprox([p_small.σ₀², p_small.ϕ, p_small.κ],
                   0.7 .* [p_start_small.σ₀², p_start_small.ϕ, p_start_small.κ];
                   rtol=1e-12, atol=1e-12)

    # Different trimming window should produce consistent scaling vs. its own β
    p_start2 = QuadVarParams(0.35, 0.22, 0.18)
    p_trim = rescalevarparams(p_start2, s_list, r_list; trim_lo=0.25, trim_hi=0.75, 
                              clip_lo=0.7, clip_hi=1.4)
    β_trim_exp = expected_beta(p_start2, s_list, r_list; trim_lo=0.25, trim_hi=0.75, 
                               clip_lo=0.7, clip_hi=1.4)
    @test isapprox([p_trim.σ₀², p_trim.ϕ, p_trim.κ], β_trim_exp .* [p_start2.σ₀², 
                   p_start2.ϕ, p_start2.κ]; rtol=1e-10, atol=1e-12)

    # 5) empty data: returns the input parameters unchanged
    p0 = QuadVarParams(0.9, 0.3, 0.1)
    p_empty = rescalevarparams(p0, Vector{Vector{Float64}}(), Vector{Matrix{Float64}}())
    @test p_empty == p0

    # 6) finite outputs for mild floors and tiny clip range
    p_start3 = QuadVarParams(0.2, 0.15, 0.12)
    p_tight = rescalevarparams(p_start3, s_list, r_list; varfloor=1e-10, clip_lo=0.95, 
                               clip_hi=1.05)
    @test all(isfinite, (p_tight.σ₀², p_tight.ϕ, p_tight.κ))
    @test 0.95 ≤ p_tight.σ₀² / p_start3.σ₀² ≤ 1.05
    @test 0.95 ≤ p_tight.ϕ / p_start3.ϕ ≤ 1.05
    @test 0.95 ≤ p_tight.κ / p_start3.κ ≤ 1.05
end

# ── updatesignal! ────────────────────────────────────────────────────────────

@testset "updatesignal!()" begin

    # Helper to build a clean synthetic case: Y = s_true * g' + o
    function make_clean_case(n::Int, r::Int)
        s_true = @. 0.5 + 0.3 * sin((1:n) / n * 2π)
        g = [0.8, 1.2, 1.0][1:r]
        o = zeros(n, r)
        y = s_true * g'
        (s_true, g, o, y)
    end

    # Noiseless recovery with damping=1 → s matches s_true exactly
    n, r = 32, 3
    s_true, g, o, y = make_clean_case(n, r)
    w = ones(n)
    s = zeros(n)
    updatesignal!(s, y, o, g, w; damping=1.0)
    @test isapprox(s, s_true; rtol=1e-12, atol=1e-12)

    # Small noise → still close for damping=1
    y_noisy = copy(y)
    Random.seed!(123)
    y_noisy .+= 1e-4 .* randn(size(y_noisy))
    s2 = zeros(n)
    updatesignal!(s2, y_noisy, o, g, w; damping=1.0)
    @test norm(s2 .- s_true) / max(norm(s_true), eps()) ≤ 1e-3

    # Damping=0 leaves s unchanged
    s0 = copy(s_true)
    s_hold = copy(s0)
    updatesignal!(s_hold, y, o, g, w; damping=0.0)
    @test s_hold == s0

    # Zero weights → numerator = 0, denominator ~ 0 → s shrinks by (1 - damping)
    w0 = zeros(n)
    s_start = fill(0.7, n)
    s_shrunk = copy(s_start)
    updatesignal!(s_shrunk, y, o, g, w0; damping=0.4)
    @test isapprox(s_shrunk, (1.0 - 0.4) .* s_start; rtol=0, atol=1e-14)

    # Nonnegativity projection: LS update negative → clamped to 0
    # make o * g dominate Y * g so num < 0 everywhere
    omin = 5.0
    o_big = fill(omin, n, r)
    y_zero = zeros(n, r)
    s_neg = fill(0.1, n)
    updatesignal!(s_neg, y_zero, o_big, g, w; damping=1.0)
    @test all(s_neg .== 0.0)

    # Only s is mutated; Y, o, g, w remain unchanged
    y_keep = copy(y)
    o_keep = copy(o)
    g_keep = copy(g)
    w_keep = copy(w)
    s_mut = zeros(n)
    updatesignal!(s_mut, y_keep, o_keep, g_keep, w_keep; damping=0.7)
    @test y_keep == y
    @test o_keep == o
    @test g_keep == g
    @test w_keep == w
    @test s_mut !== s  # different object
end

# ── residautocorr ────────────────────────────────────────────────────────────

@testset "residautocorr()" begin

    # Deterministic AR(1) generator (unit-variance innovations)
    function ar1(n::Int, ϕ::Float64)
        x = zeros(n)
        ϵ = randn(n)
        @inbounds for t in 2:n
            x[t] = ϕ * x[t-1] + ϵ[t]
        end
        x
    end

    # Simple variance model: constant variance = 1
    var_one(_) = 1.0
    # Degenerate variance model (forces fallback to varfloor)
    var_zero(_) = 0.0

    Random.seed!(2024)

    # Basic setup: three replicates, long series for tight CI
    n, r = 4000, 3
    ϕ_true = 0.6
    y = [ar1(n, ϕ_true) ar1(n, ϕ_true) ar1(n, ϕ_true)]
    o = zeros(n, r)
    g = zeros(r)  # μ = o[:,j] + g[j]*s → here μ ≡ 0
    s = zeros(n)  # unused by var_one/var_zero

    # lag-1 ACF close to ϕ_true for each replicate and pooled
    ρ_all, ρ_series, n_pairs = residautocorr(y, o, g, s, var_one; lag=1)
    @test all(n_pairs .== n - 1)
    @test isfinite(ρ_all) && all(isfinite.(ρ_series))
    @test abs(ρ_all - ϕ_true) ≤ 0.03
    @test all(abs.(ρ_series .- ϕ_true) .≤ 0.05)

    # lag-2 ≈ ϕ_true^2
    ρ_all2, ρ_series2, n_pairs2 = residautocorr(y, o, g, s, var_one; lag=2)
    @test all(n_pairs2 .== n - 2)
    @test abs(ρ_all2 - ϕ_true^2) ≤ 0.04
    @test all(abs.(ρ_series2 .- ϕ_true^2) .≤ 0.06)

    # Mask: keep only latter 75% → same correlation, fewer pairs
    keep = trues(n); keep[1:floor(Int, 0.25n)] .= false
    ρ_all_mask, ρ_series_mask, n_pairs_mask = residautocorr(y, o, g, s, var_one; 
        lag=1, mask=keep)
    @test all(n_pairs_mask .== count(keep) - 1)
    @test abs(ρ_all_mask - ϕ_true) ≤ 0.04
    @test all(abs.(ρ_series_mask .- ϕ_true) .≤ 0.06)

    # Trimming protects against outliers
    y_spiky = copy(y)
    y_spiky[500, 1]  = 1e6
    y_spiky[1500, 2] = -1e6
    y_spiky[2500, 3] = 1e6
    ρ_all_untrim, _, _ = residautocorr(y_spiky, o, g, s, var_one; lag=1, trimfrac=0.0)
    ρ_all_trim,  _, _ = residautocorr(y_spiky, o, g, s, var_one; lag=1, trimfrac=0.05)
    @test abs(ρ_all_trim - ϕ_true) ≤ abs(ρ_all_untrim - ϕ_true)

    # Variance floor: using var_zero with a positive varfloor should match var_one
    ρ_all_floor, ρ_series_floor, n_pairs_floor = residautocorr(y, o, g, s, var_zero; 
        lag=1, varfloor=1e-6)
    @test ρ_all_floor ≈ ρ_all atol=1e-12 rtol=1e-12
    @test ρ_series_floor ≈ ρ_series atol=1e-12 rtol=1e-12
    @test n_pairs_floor == n_pairs

    # Edge case: lag too large → no pairs, NaNs and zeros
    ρ_all_none, ρ_series_none, n_pairs_none = residautocorr(y, o, g, s, var_one; lag=n)
    @test isnan(ρ_all_none)
    @test all(isnan, ρ_series_none)
    @test all(==(0), n_pairs_none)

    # Short series with finite mask ensures shape checks and minimal pairs path
    y_short = y[1:5, :]
    o_short = o[1:5, :]
    s_short = s[1:5]
    ρ_all_s, ρ_series_s, n_pairs_s = residautocorr(y_short, o_short, g, s_short, var_one; 
        lag=1, mask=[true, true, false, true, true])
    @test length(ρ_series_s) == r
    @test length(n_pairs_s) == r
    @test sum(n_pairs_s) ≥ 0  # just exercises the masked branch

    # Input validation
    @test_throws ArgumentError residautocorr(y, o, g, s, var_one; lag=0)
    @test_throws ArgumentError residautocorr(y, o[:,1:2], g, s, var_one)
    @test_throws ArgumentError residautocorr(y, o, g[1:2], s, var_one)
    @test_throws ArgumentError residautocorr(y, o, g, s[1:end-1], var_one)
    @test_throws ArgumentError residautocorr(y, o, g, s, var_one; trimfrac=0.6)
end

# ── QuadVarFit ───────────────────────────────────────────────────────────────

@testset "QuadVarFit()" begin
    # Small synthetic fixture
    mz_ref = [100.0, 101.0, 102.0]
    mzunit = nothing
    mz_idx = [1, 3]
    mzvalues = mz_ref[mz_idx]
    n_sel = length(mz_idx)

    batchcount = 2
    n_reps_per_batch = [2, 1]
    n_scans_per_batch = [4, 3]

    # Concrete Float64 params per selected m/z
    params = QuadVarParams{Float64,Float64,Float64}[
        QuadVarParams(1.0, 0.1, 0.2),
        QuadVarParams(0.5, 0.0, 0.3),
    ]

    # signal[j][b] :: Vector{Float64} (len = scans per batch b)
    signal = [
        [fill(10.0, n_scans_per_batch[1]), fill(11.0, n_scans_per_batch[2])],
        [fill(20.0, n_scans_per_batch[1]), fill(21.0, n_scans_per_batch[2])],
    ]

    # offsets[j][b] :: Matrix{Float64} (scans_b × reps_b)
    offsets = [
        [fill(0.1, n_scans_per_batch[1], n_reps_per_batch[1]),
         fill(0.2, n_scans_per_batch[2], n_reps_per_batch[2])],
        [fill(0.3, n_scans_per_batch[1], n_reps_per_batch[1]),
         fill(0.4, n_scans_per_batch[2], n_reps_per_batch[2])],
    ]

    # gains[j][b] :: Vector{Float64} (len = reps_b)
    gains = [
        [fill(1.1, n_reps_per_batch[1]), fill(1.2, n_reps_per_batch[2])],
        [fill(2.1, n_reps_per_batch[1]), fill(2.2, n_reps_per_batch[2])],
    ]

    # per-selected-m/z scalars/vectors
    scale_c = [1.0, 2.0]
    acf = [0.1, 0.2]
    acf_lag = [1, 1]
    n_pairs = [100, 90]
    qc_z_rms = [0.9, 1.1]
    qc_cov68 = [0.68, 0.70]
    qc_cov95 = [0.95, 0.93]
    qc_s_min = [0.0, 0.0]
    qc_s_max = [25.0, 30.0]
    qc_nz = [400, 380]

    # observed[b][r] :: scans_b × n_sel
    observed = Vector{Vector{Matrix{Float64}}}(undef, batchcount)
    for b in 1:batchcount
        observed[b] = [fill(0.0, n_scans_per_batch[b], n_sel) for _ in 1:n_reps_per_batch[b]]
    end
    intensityunit = nothing

    f = QuadVarFit(
        mzunit, mz_ref, mz_idx, mzvalues,
        batchcount,
        n_reps_per_batch, n_scans_per_batch,
        params, intensityunit, signal, offsets, gains,
        scale_c, acf, acf_lag, n_pairs,
        qc_z_rms, qc_cov68, qc_cov95, qc_s_min, qc_s_max, qc_nz,
        observed
    )

    # Type/param checks
    @test f isa QuadVarFit
    @test f isa QuadVarFit{<:AbstractVector{<:Real}, QuadVarParams{Float64,Float64,Float64}, Nothing, Nothing}

    # Basic shape/consistency
    @test f.batchcount == batchcount
    @test f.mzunit === mzunit
    @test f.mzvalues == f.mz_ref[f.mz_idx]
    @test length(f.mz_idx) == length(f.mzvalues) == n_sel
    @test allunique(f.mz_idx)
    @test f.intensityunit === intensityunit

    # params element type is concrete subtype
    @test f.params isa Vector{<:QuadVarParams}
    @test eltype(f.params) ≡ QuadVarParams{Float64,Float64,Float64}
    @test length(f.params) == n_sel

    # signal/offsets/gains structure
    @test length(f.signal) == n_sel
    @test length(f.offsets) == n_sel
    @test length(f.gains) == n_sel
    for j in 1:n_sel
        @test length(f.signal[j]) == batchcount
        @test length(f.offsets[j]) == batchcount
        @test length(f.gains[j]) == batchcount
        for b in 1:batchcount
            @test f.signal[j][b] isa Vector{Float64}
            @test length(f.signal[j][b]) == n_scans_per_batch[b]
            @test f.offsets[j][b] isa Matrix{Float64}
            @test size(f.offsets[j][b]) == (n_scans_per_batch[b], n_reps_per_batch[b])
            @test f.gains[j][b] isa Vector{Float64}
            @test length(f.gains[j][b]) == n_reps_per_batch[b]
        end
    end

    # qc / summary fields lengths
    @test length.((f.scale_c, f.acf, f.acf_lag, f.n_pairs, f.qc_z_rms, f.qc_cov68, 
                   f.qc_cov95, f.qc_s_min, f.qc_s_max, f.qc_nz)) == ntuple(_ -> n_sel, 10)

    # Observed shape
    @test length(f.observed) == batchcount
    for b in 1:batchcount
        @test length(f.observed[b]) == n_reps_per_batch[b]
        for r in 1:n_reps_per_batch[b]
            @test size(f.observed[b][r]) == (n_scans_per_batch[b], n_sel)
        end
    end

    # Broadcastable behavior (treat as scalar)
    @test Base.broadcastable(f) isa Ref
    bmask = f .== f
    @test size(bmask) == ()
    @test only(bmask) ≡ true

    arr = [1, 2, 3]
    res = ((_, q) -> q).(arr, f)
    @test length(res) == length(arr)
    @test all(x -> x ≡ f, res)

    # Immutability and concreteness
    @test isimmutable(f)
    @test !ismutabletype(typeof(f))
    @test isconcretetype(typeof(f))

    # Field type sanity
    @test f.mz_ref isa AbstractVector{<:Number}
    @test f.mz_idx isa Vector{Int}
    @test f.params isa Vector{QuadVarParams{Float64,Float64,Float64}}
    @test f.signal isa Vector{Vector{Vector{Float64}}}
    @test f.offsets isa Vector{Vector{Matrix{Float64}}}
    @test f.gains  isa Vector{Vector{Vector{Float64}}}
end

# ── QuadVarParams ────────────────────────────────────────────────────────────

@testset "QuadVarParams(σ₀²::Real, ϕ::Real, κ::Real)" begin
    p = QuadVarParams(1.0, Float32(0.2), 3)
    @test typeof(p) ≡ QuadVarParams{Float64, Float32, Int64}
    @test p.σ₀² ≡ 1.0
    @test p.ϕ ≡ Float32(0.2)
    @test p.κ ≡ 3

    # Explicit parametric constructor
    p32 = QuadVarParams{Float32,Float32,Float32}(1, 2, 3)
    @test typeof(p32) ≡ QuadVarParams{Float32,Float32,Float32}
    @test p32.σ₀² ≡ Float32(1)
    @test p32.ϕ ≡ Float32(2)
    @test p32.κ ≡ Float32(3)

    # BigFloat path
    pbig = QuadVarParams(big"1.0", big"0.2", big"0.5")
    @test typeof(pbig) ≡ QuadVarParams{BigFloat,BigFloat,BigFloat}

    # Immutability & concreteness
    @test isimmutable(p)  # instance is immutable
    @test !ismutabletype(typeof(p))  # type is an immutable struct
    @test isconcretetype(typeof(p))
    @test isbitstype(typeof(QuadVarParams(1.0, 0.0, 0.0)))

    # Property names
    @test (:σ₀² in propertynames(p)) &&
          (:ϕ   in propertynames(p)) &&
          (:κ   in propertynames(p))

    # broadcastable: should wrap in Ref so p is treated as a scalar
    @test Base.broadcastable(p) isa Ref
    b = p .== p
    @test size(b) == ()  # 0-dim result (scalar broadcast)
    @test only(b) ≡ true

    # Interaction with broadcast over an array (p should act as scalar)
    xs = [0.0, 1.0, 2.0]
    @test ((x, q)-> q.κ * x^2 + q.ϕ * x + q.σ₀²).(xs, p) ≈
          [p.κ * x^2 + p.ϕ * x + p.σ₀² for x in xs]

    # Mixed-type fields still permitted
    pmix = QuadVarParams(1//2, 0.1, 2)
    @test typeof(pmix) ≡ QuadVarParams{Rational{Int64}, Float64, Int64}
    @test pmix.σ₀² ≡ 1//2

    # Edge coefficients
    pzero = QuadVarParams(0.0, 0.0, 0.0)
    @test pzero.σ₀² == 0.0 && pzero.ϕ == 0.0 && pzero.κ == 0.0

    # Float32 parameters should produce a Float32-typed QuadVarParams
    p32 = QuadVarParams(1f0, 2f0, 3f0)
    @test p32 isa QuadVarParams{Float32, Float32, Float32}

    # BigFloat parameters should produce a BigFloat-typed QuadVarParams
    pBig = QuadVarParams(big(1.0), big(2.0), big(3.0))
    @test pBig isa QuadVarParams{BigFloat, BigFloat, BigFloat}

    # Mixed real types should preserve the individual parameter types
    pMix = QuadVarParams(1f0, 2.0, 3)  # Float32, Float64, Int
    @test pMix isa QuadVarParams{Float32, Float64, Int}

    # Sanity: varpred should honor/promote eltypes correctly with non-Float64 params
    y32 = Float32[0, 1, 2]
    @test eltype(varpred(y32, p32)) == Float32

    y = [0.0, 1.0, 2.0]
    @test eltype(varpred(y, pBig)) == BigFloat
end

# ── residuals ────────────────────────────────────────────────────────────────

@testset "residuals(Y, o, s, g)" begin
    # Basic correctness & shape
    n, r = 5, 3
    Y = [10.0 11 12; 13 14 15; 16 17 18; 19 20 21; 22 23 24]
    o = fill(1.0, n, r)
    s = collect(1.0:n)  # length n
    g = [1.0, 2.0, 3.0]  # length r
    G = s * g'  # n × r rank - 1 matrix
    Rexp = Y .- o .- G
    Rgot = residuals(Y, o, s, g)
    @test size(Rgot) == (n, r)
    @test Rgot ≈ Rexp

    # Zero offsets or unit gains reduce to simpler forms
    @test residuals(Y, zeros(n, r), s, g) ≈ Y .- (s * g')
    @test residuals(Y, o, s, ones(r)) ≈ Y .- o .- s * ones(r)'

    # Single column / single row corner cases
    Y1 = reshape(1.0:5.0, n, 1); o1 = zeros(n, 1); g1 = [2.0]
    @test residuals(Y1, o1, s, g1) ≈ Y1 .- s * g1'
    Yr = reshape(1.0:3.0, 1, 3); or = zeros(1, 3); sr = [4.0]; gr = [1.0, 2.0, 3.0]
    @test residuals(Yr, or, sr, gr) ≈ Yr .- sr * gr'

    # Type behavior (promotion along the arithmetic path)
    Y32 = Float32.(Y); o32 = fill(Float32(0.5), n, r); s32 = Float32.(s); g32 = Float32.(g)
    R32 = residuals(Y32, o32, s32, g32)
    @test eltype(R32) ≡ Float32
    # Mixed Float32 + Float64 → Float64
    Rmix = residuals(Y, o32, s32, g)
    @test eltype(Rmix) ≡ Float64

    # BigFloat path
    Yb = big.(Y); ob = big.(o); sb = big.(s); gb = big.(g)
    Rb = residuals(Yb, ob, sb, gb)
    @test eltype(Rb) ≡ BigFloat
    @test Rb ≈ big.(Rexp)

    # NaN/Inf propagate as standard arithmetic does
    Ynan = copy(Y); Ynan[2,3] = NaN
    @test isnan(residuals(Ynan, o, s, g)[2,3])
    Yinf = copy(Y); Yinf[3,1] = Inf
    @test isinf(residuals(Yinf, o, s, g)[3,1])

    # Dimension checks (should throw)
    @test_throws DimensionMismatch residuals(Y, zeros(n, r + 1), s, g)   # o wrong size
    @test_throws DimensionMismatch residuals(Y, o, s[1:end-1], g)      # s length mismatch
    @test_throws DimensionMismatch residuals(Y, o, s, g[1:end-1])      # g length mismatch
    @test_throws DimensionMismatch residuals(Y[1:end-1, :], o, s, g)   # Y wrong size

    # Immutability of inputs (light sanity: run and compare)
    Yc, oc, sc, gc = copy(Y), copy(o), copy(s), copy(g)
    _ = residuals(Yc, oc, sc, gc)
    @test Yc == Y && oc == o && sc == s && gc == g
end

# ── show(::QuadVarFit) and show(::MIME"text/plain", ::QuadVarFit) ────────────

@testset "show(::QuadVarFit) and show(::MIME\"text/plain\", ::QuadVarFit)" begin

    # Minimal, internally-consistent builder
    function _mkfit(; n_mz=8, idx=collect(2:7), reps=(3,4), scans=(10,12))
        mz_ref = collect(range(100.0, length=n_mz, step=0.5))
        mzunit = nothing
        mz_idx = idx
        mz_vals = mz_ref[mz_idx]
        batchcount = length(reps)
        n_reps_per_batch = collect(reps)
        n_scans_per_batch = collect(scans)
        n_sel = length(mz_idx)

        params = [QuadVarParams(1.0 + 0.01*j, 0.1, 0.2) for j in 1:n_sel]
        signal = [[fill(1.0 + 0.1*b, n_scans_per_batch[b]) for b in 1:batchcount] 
                  for _ in 1:n_sel]
        offsets = [[zeros(n_scans_per_batch[b], n_reps_per_batch[b]) for b in 1:batchcount] 
                   for _ in 1:n_sel]
        gains   = [[ones(n_reps_per_batch[b]) for b in 1:batchcount] for _ in 1:n_sel]
        scale_c = fill(1.0, n_sel)
        acf     = fill(0.2, n_sel)
        acf_lag = fill(1, n_sel)
        n_pairs = fill(100, n_sel)
        qc_z_rms = fill(1.0, n_sel)
        qc_cov68 = fill(0.68, n_sel)
        qc_cov95 = fill(0.95, n_sel)
        qc_s_min = fill(0.0, n_sel)
        qc_s_max = fill(10.0, n_sel)
        qc_nz    = fill(50, n_sel)
        observed = [[zeros(n_scans_per_batch[b], n_sel) for _ in 1:n_reps_per_batch[b]] 
                    for b in 1:batchcount]
        intensityunit = nothing

        QuadVarFit(mzunit, mz_ref, mz_idx, mz_vals, batchcount,
                   n_reps_per_batch, n_scans_per_batch,
                   params, intensityunit, signal, offsets, gains,
                   scale_c, acf, acf_lag, n_pairs,
                   qc_z_rms, qc_cov68, qc_cov95, qc_s_min, qc_s_max, qc_nz,
                   observed)
    end

    # Compact show(io, q) prints exactly summary(q)
    q = _mkfit(; n_mz=8, idx=collect(2:7), reps=(3,4), scans=(10,12))
    s_compact = sprint(io -> show(io, q))
    @test s_compact == summary(q)

    # Pretty show(io, "text/plain", q) includes detailed sections
    s_pretty = sprint(io -> show(io, MIME"text/plain"(), q))
    @test startswith(s_pretty, summary(q) * "\n")
    @test occursin("  batches          : 2", s_pretty)
    @test occursin("    reps per batch : [3, 4]", s_pretty)
    @test occursin("    scans per batch: [10, 12]", s_pretty)
    @test occursin("  ion grid         : 8 total; selected: 6", s_pretty)
    @test occursin("    m/z values     : ", s_pretty)  # ≤ 6 selected → full list printed

    # > 6 selected m/z → head … tail format
    q_many = _mkfit(; n_mz=12, idx=collect(1:10), reps=(2,2), scans=(5,6))
    s_pretty_many = sprint(io -> show(io, MIME"text/plain"(), q_many))
    @test occursin("  ion grid         : 12 total; selected: 10", s_pretty_many)
    @test occursin("    m/z values     : [", s_pretty_many)
    @test occursin("…", s_pretty_many)  # ellipsis present

    # 0 selected m/z → no "m/z values" line
    q_empty = _mkfit(; n_mz=5, idx=Int[], reps=(2,2), scans=(7,7))
    s_pretty_empty = sprint(io -> show(io, MIME"text/plain"(), q_empty))
    @test occursin("QuadVarFit (batches=2, ions=0)", s_pretty_empty)
    @test occursin("  ion grid         : 5 total; selected: 0", s_pretty_empty)
    @test !occursin("m/z values", s_pretty_empty)
end

# ── summary(::QuadVarFit) ────────────────────────────────────────────────────

@testset "summary(::QuadVarFit)" begin

    # Minimal, internally-consistent builder
    function _mkfit(; n_mz=8, idx=collect(2:7), reps=(3,4), scans=(10,12))
        mz_ref = collect(range(100.0, length=n_mz, step=0.5))
        mzunit = nothing
        mz_idx = idx
        mz_vals = mz_ref[mz_idx]
        batchcount = length(reps)
        n_reps_per_batch = collect(reps)
        n_scans_per_batch = collect(scans)
        n_sel = length(mz_idx)

        params  = [QuadVarParams(1.0, 0.1, 0.2) for _ in 1:n_sel]
        signal  = [[fill(1.0, n_scans_per_batch[b]) for b in 1:batchcount] for _ in 1:n_sel]
        offsets  = [[zeros(n_scans_per_batch[b], n_reps_per_batch[b]) for b in 1:batchcount] 
                    for _ in 1:n_sel]
        gains    = [[ones(n_reps_per_batch[b]) for b in 1:batchcount] for _ in 1:n_sel]
        scale_c  = fill(1.0, n_sel)
        acf      = fill(0.0, n_sel)
        acf_lag  = fill(1, n_sel)
        n_pairs  = fill(0, n_sel)
        qc_z_rms = fill(1.0, n_sel)
        qc_cov68 = fill(0.68, n_sel)
        qc_cov95 = fill(0.95, n_sel)
        qc_s_min = fill(0.0, n_sel)
        qc_s_max = fill(1.0, n_sel)
        qc_nz    = fill(1, n_sel)
        observed = [[zeros(n_scans_per_batch[b], n_sel) for _ in 1:n_reps_per_batch[b]] 
                    for b in 1:batchcount]
        intensityunit = nothing

        QuadVarFit(mzunit, mz_ref, mz_idx, mz_vals, batchcount,
                   n_reps_per_batch, n_scans_per_batch,
                   params, intensityunit, signal, offsets, gains,
                   scale_c, acf, acf_lag, n_pairs,
                   qc_z_rms, qc_cov68, qc_cov95, qc_s_min, qc_s_max, qc_nz,
                   observed)
    end

    # Basic cases
    q = _mkfit(; n_mz=8, idx=collect(2:7), reps=(3,4), scans=(10,12))
    @test summary(q) == "QuadVarFit (batches=2, ions=6)"

    q_empty = _mkfit(; n_mz=5, idx=Int[], reps=(2,2), scans=(7,7))
    @test summary(q_empty) == "QuadVarFit (batches=2, ions=0)"

    q_many = _mkfit(; n_mz=12, idx=collect(1:10), reps=(1,3), scans=(9,11))
    @test summary(q_many) == "QuadVarFit (batches=2, ions=10)"
end

# ── varpred ──────────────────────────────────────────────────────────────────

@testset "varpred(y, p::QuadVarParams; varfloor)" begin
    # Helpers (assumes QuadVarParams(σ₀², ϕ, κ))
    p1  = QuadVarParams(1.0, 0.2, 0.5)
    p2  = QuadVarParams(0.0, 0.0, 0.0)  # κ=0 → linear/constant
    p3  = QuadVarParams(-1e-3, 0.0, 0.0)  # negative base to trigger floor
    p32 = QuadVarParams(Float32(1.0), Float32(0.1), Float32(0.3))
    pbig = QuadVarParams(big"1.0", big"0.2", big"0.5")

    # Scalar basics
    y = 3.0
    @test isapprox(varpred(y, p1), p1.σ₀² + p1.ϕ * y + p1.κ * y^2; rtol=1e-12)

    # κ=0 reduces to σ₀² + ϕ y
    @test varpred(2.0, p2) == 0.0

    # Vector broadcasting & values
    ys = [0.0, 1.0, 2.0]
    exp_vec = [p1.σ₀² + p1.ϕ * u + p1.κ * u^2 for u in ys]
    @test size(varpred(ys, p1)) == size(ys)
    @test isapprox(varpred(ys, p1), exp_vec; rtol=1e-12)

    # Matrix broadcasting & values
    Y = reshape(1.0:4.0, 2, 2)
    exp_mat = [p1.σ₀² + p1.ϕ * u + p1.κ * u^2 for u in Y] |> x -> reshape(x, size(Y))
    @test size(varpred(Y, p1)) == (2, 2)
    @test isapprox(varpred(Y, p1), exp_mat; rtol=1e-12)

    # varfloor: lower bound & replacement of non-finite
    @test varpred(0.0, p3; varfloor=0.0) == 0.0  # negative → floored to 0
    @test varpred(0.0, p3; varfloor=1e-6) == 1e-6  # negative → floored to varfloor

    # Replace non-finite y with varfloor
    @test varpred(Inf, p1; varfloor=1e-9) == 1e-9
    @test varpred(-Inf, p1; varfloor=1e-9) == 1e-9
    @test isnan(p1.σ₀² + p1.ϕ * NaN + p1.κ * NaN^2)  # sanity of raw math
    @test varpred(NaN, p1; varfloor=1e-9) == 1e-9

    # Mixed finite and non-finite arrays
    ymix = [0.0, NaN, 2.0, Inf]
    outmix = varpred(ymix, p1; varfloor=1e-8)
    @test outmix[1] ≥ 1e-8 && isfinite(outmix[1])
    @test outmix[2] == 1e-8
    @test isfinite(outmix[3]) && outmix[3] ≥ 1e-8
    @test outmix[4] == 1e-8

    # Monotonicity (elementwise quadratic)
    ymono = 3.0
    p_lowϕ = QuadVarParams(1.0, 0.1, 0.5)
    p_hiϕ  = QuadVarParams(1.0, 0.3, 0.5)
    @test varpred(ymono, p_lowϕ) < varpred(ymono, p_hiϕ)

    p_lowκ = QuadVarParams(1.0, 0.2, 0.1)
    p_hiκ  = QuadVarParams(1.0, 0.2, 0.7)
    @test varpred(ymono, p_lowκ) ≤ varpred(ymono, p_hiκ)

    # Type promotion and element types
    # The result element type is promote_type(eltype(y), typeof(p.σ₀²), typeof(p.ϕ), typeof(p.κ), typeof(varfloor))
    @test typeof(varpred(1, p32)) ≡ Float32
    @test typeof(varpred(1, p32; varfloor=1.0)) ≡ Float64
    @test typeof(varpred(1, p32; varfloor=big"0.0")) ≡ BigFloat
    @test typeof(varpred(Int16(2), QuadVarParams(1.0, 0.0, 0.0))) ≡ Float64
    @test typeof(varpred(Float32[1,2,3], p32)) ≡ Vector{Float32}
    @test typeof(varpred(reshape(Float32.(1:4), 2, 2), p32)) ≡ Matrix{Float32}
    @test typeof(varpred(Float32[1,2,3], p32; varfloor=1.0)) ≡ Vector{Float64}
    @test typeof(varpred(1, pbig)) ≡ BigFloat

    # Float32 stability (avoid exact decimal equality)
    y32 = Float32(5)
    v32 = @inferred varpred(y32, p32; varfloor=Float32(0))  # pass Float32 floor for type-stable inference
    exp32 = p32.σ₀² + p32.ϕ * y32 + p32.κ * y32 * y32
    @test isapprox(v32, exp32; rtol=1e-6, atol=eps(Float32))

    # Shape preservation
    A = randn(3, 4)
    outA = varpred(A, p1)
    @test size(outA) == size(A)
    @test all(isfinite, outA)

    # Unitful inputs propagate squared intensity units and auto-coerce floors
    p_units = QuadVarParams(1.0u"pA^2", 0.25u"pA", 0.1)
    y_unit = 3.0u"pA"
    res_unit = varpred(y_unit, p_units)
    @test res_unit isa Unitful.AbstractQuantity
    @test Unitful.dimension(res_unit) == Unitful.dimension(1.0u"pA^2")
    p_units_floor = QuadVarParams(-1.0u"pA^2", 0.0u"pA", 0.0)
    @test varpred(0.0u"pA", p_units_floor; varfloor=1e-9) == 1e-9u"pA^2"
    @test_throws Unitful.DimensionError varpred(y_unit, p_units; varfloor=1.0u"s")
end

# ── varpred ──────────────────────────────────────────────────────────────────

@testset "varpred(y, σ₀²::Number, ϕ::Number, κ::Number; varfloor::Number=0)" begin
    # Scalar basics
    @test varpred(0.0, 1.0, 0.2, 0.5) == 1.0
    y = 3.0
    @test isapprox(varpred(y, 1.0, 0.2, 0.5), 1.0 + 0.2 * y + 0.5 * y^2; rtol=1e-12)

    # With κ=0 the formula reduces to σ₀² + ϕ y
    @test varpred(2.0, 0.0, 0.0, 0.0) == 0.0

    # Vector broadcasting and values
    ys = [0.0, 1.0, 2.0]
    exp_vec = [(1.0 + 0.2 * u + 0.5 * u^2) for u in ys]
    @test size(varpred(ys, 1.0, 0.2, 0.5)) == size(ys)
    @test isapprox(varpred(ys, 1.0, 0.2, 0.5), exp_vec; rtol=1e-12)

    # Matrix broadcasting & values
    Y = reshape(1.0:4.0, 2, 2)
    exp_mat = [(1.0 + 0.2 * u + 0.5 * u^2) for u in Y]
    exp_mat = reshape(exp_mat, size(Y))
    @test size(varpred(Y, 1.0, 0.2, 0.5)) == (2, 2)
    @test isapprox(varpred(Y, 1.0, 0.2, 0.5), exp_mat; rtol=1e-12)

    # varfloor: lower bound and replacement of non-finite
    @test varpred(0.0, -1e-3, 0.0, 0.0; varfloor=0.0) == 0.0
    @test varpred(0.0, -1e-3, 0.0, 0.0; varfloor=1e-6) == 1e-6

    # Replace non-finite results (Inf/NaN) with varfloor
    @test varpred(Inf, 1.0, 0.2, 0.5; varfloor=1e-9) == 1e-9
    @test varpred(-Inf, 1.0, 0.2, 0.5; varfloor=1e-9) == 1e-9
    @test isnan(1.0 + 0.2 * NaN + 0.5 * NaN^2)
    @test varpred(NaN, 1.0, 0.2, 0.5; varfloor=1e-9) == 1e-9

    # Mixed finite and non-finite in arrays
    ymix = [0.0, NaN, 2.0, Inf]
    outmix = varpred(ymix, 1.0, 0.2, 0.5; varfloor=1e-8)
    @test outmix[1] ≥ 1e-8 && isfinite(outmix[1])
    @test outmix[2] == 1e-8
    @test isfinite(outmix[3]) && outmix[3] ≥ 1e-8
    @test outmix[4] == 1e-8

    # Monotonicity (simple quadratic, elementwise)
    ymono = 3.0
    @test varpred(ymono, 1.0, 0.1, 0.5) <  varpred(ymono, 1.0, 0.3, 0.5)
    @test varpred(ymono, 1.0, 0.2, 0.1) ≤ varpred(ymono, 1.0, 0.2, 0.7)

    # Type promotion and element types
    @test typeof(varpred(1, Float32(1), Float32(0.1), Float32(0.3))) ≡
          promote_type(Int, Float32, Float32, Float32, Int)

    @test typeof(varpred(1, Float32(1), Float32(0.1), Float32(0.3); varfloor=1.0)) ≡
          promote_type(Int, Float32, Float32, Float32, Float64)

    @test typeof(varpred(Float32[1,2,3], Float32(1), Float32(0.1), Float32(0.3))) ≡
          Vector{promote_type(Float32, Float32, Float32, Float32, Int)}

    @test typeof(varpred(Float32[1,2,3], Float32(1), Float32(0.1), Float32(0.3); varfloor=1.0)) ≡
          Vector{promote_type(Float32, Float32, Float32, Float32, Float64)}

    @test typeof(varpred(1, big"1.0", big"0.2", big"0.5")) ≡
          promote_type(Int, BigFloat, BigFloat, BigFloat, Int)

    # Float32 block — ensure type-stable inference by passing a Float32 floor
    y32 = Float32(5)
    σ0²32, ϕ32, κ32 = Float32(1), Float32(0.1), Float32(0.3)
    v32 = @inferred varpred(y32, σ0²32, ϕ32, κ32; varfloor=Float32(0))
    exp32 = (σ0²32 + ϕ32 * y32 + κ32 * y32 * y32)
    @test isapprox(v32, exp32; rtol=1e-6, atol=eps(Float32))

    # Shape preservation
    A = randn(3, 4)
    outA = varpred(A, 1.0, 0.2, 0.5)
    @test size(outA) == size(A)
    @test all(isfinite, outA)

    # Unitful scalars inherit squared intensity units and auto-coerce floors
    res_unit = varpred(2.0u"pA", 1.0u"pA^2", 0.2u"pA", 0.1)
    @test res_unit isa Unitful.AbstractQuantity
    @test Unitful.dimension(res_unit) == Unitful.dimension(1.0u"pA^2")
    @test varpred(0.0u"pA", -1.0u"pA^2", 0.0u"pA", 0.0; varfloor=1e-9) == 1e-9u"pA^2"
end

# ── varpredbias ──────────────────────────────────────────────────────────────

@testset "varpredbias(y, p::QuadVarParams; varfloor::Number=0)" begin
    # Helpers to build params (assumes positional ctor: QuadVarParams(σ₀², ϕ, κ))
    p1 = QuadVarParams(1.0, 0.2, 0.5)  # den = 1 + κ = 1.5
    p2 = QuadVarParams(0.0, 0.0, 0.0)  # identity quadratic, den=1
    p3 = QuadVarParams(-1e-3, 0.0, 0.0)  # negative base to trigger varfloor
    p32 = QuadVarParams(Float32(1.0), Float32(0.1), Float32(0.3))
    pbig = QuadVarParams(big"1.0", big"0.2", big"0.5")

    # Scalar basics
    y = 3.0
    @test isapprox(varpredbias(y, p1), (1.0 + 0.2 * y + 0.5 * y^2) / (1 + 0.5); rtol=1e-12)

    # With κ=0 the formula reduces to σ₀² + ϕ y + κ y^2
    @test varpredbias(2.0, p2) == 0.0

    # Vector broadcasting and values
    ys = [0.0, 1.0, 2.0]
    exp_vec = [(p1.σ₀² + p1.ϕ * u + p1.κ * u^2) / (1 + p1.κ) for u in ys]
    @test size(varpredbias(ys, p1)) == size(ys)
    @test isapprox(varpredbias(ys, p1), exp_vec; rtol=1e-12)

    # Matrix broadcasting & values
    Y = reshape(1.0:4.0, 2, 2)
    exp_mat = [(p1.σ₀² + p1.ϕ * u + p1.κ * u^2) / (1 + p1.κ) for u in Y]
    exp_mat = reshape(exp_mat, size(Y))
    @test size(varpredbias(Y, p1)) == (2, 2)
    @test isapprox(varpredbias(Y, p1), exp_mat; rtol=1e-12)

    # varfloor: lower bound and replacement of non-finite
    @test varpredbias(0.0, p3; varfloor=0.0) == 0.0  # negative → floored to 0
    @test varpredbias(0.0, p3; varfloor=1e-6) == 1e-6  # negative → floored to varfloor

    # Replace non-finite results (Inf/NaN) with varfloor
    @test varpredbias(Inf, p1; varfloor=1e-9) == 1e-9
    @test varpredbias(-Inf, p1; varfloor=1e-9) == 1e-9
    @test isnan((p1.σ₀² + p1.ϕ * NaN + p1.κ * NaN^2) / (1 + p1.κ))  # sanity check of raw math
    @test varpredbias(NaN, p1; varfloor=1e-9) == 1e-9

    # Mixed finite and non-finite in arrays
    ymix = [0.0, NaN, 2.0, Inf]
    outmix = varpredbias(ymix, p1; varfloor=1e-8)
    @test outmix[1] ≥ 1e-8 && isfinite(outmix[1])
    @test outmix[2] == 1e-8
    @test isfinite(outmix[3]) && outmix[3] ≥ 1e-8
    @test outmix[4] == 1e-8

    # Monotonicity when 1 + κ > 0
    # For fixed y > 0 and κ ≥ 0, result increases with ϕ and κ (denominator effect bounded)
    ymono = 3.0
    p_lowκ = QuadVarParams(1.0, 0.2, 0.1)
    p_hiκ = QuadVarParams(1.0, 0.2, 0.7)
    @test varpredbias(ymono, p_lowκ) ≤ varpredbias(ymono, p_hiκ)
    p_lowϕ = QuadVarParams(1.0, 0.1, 0.5)
    p_hiϕ = QuadVarParams(1.0, 0.3, 0.5)
    @test varpredbias(ymono, p_lowϕ) < varpredbias(ymono, p_hiϕ)

    # Type promotion and element types
    # The result element type is promote_type(eltype(y), typeof(p.σ₀²), typeof(p.ϕ), typeof(p.κ), typeof(varfloor))
    @test typeof(varpredbias(1, p32)) ≡ Float32
    @test typeof(varpredbias(1, p32; varfloor=1.0)) ≡ Float64
    @test typeof(varpredbias(1, p32; varfloor=big"0.0")) ≡ BigFloat
    @test typeof(varpredbias(Int16(2), QuadVarParams(1.0, 0.0, 0.0))) ≡ Float64
    @test typeof(varpredbias(Float32[1,2,3], p32)) ≡ Vector{Float32}
    @test typeof(varpredbias(reshape(Float32.(1:4), 2, 2), p32)) ≡ Matrix{Float32}
    @test typeof(varpredbias(1, pbig)) ≡ BigFloat

    # Finite-numeric stability (32-bit tolerances)
    # Don’t assert exact Float32 decimals; use a tolerance
    y32 = Float32(5)
    val32 = varpredbias(y32, p32)
    exp32 = (p32.σ₀² + p32.ϕ * y32 + p32.κ * y32 * y32) / (one(p32.κ) + p32.κ)
    @test isapprox(val32, exp32; rtol=1e-6, atol=eps(Float32))

    # Denominator positivity assumption: 1 + κ > 0
    # Pick κ near -1 but valid; should still compute finite results
    p_neg_close = QuadVarParams(1.0, 0.0, -0.999)
    @test isfinite(varpredbias(2.0, p_neg_close))
    @test varpredbias(2.0, p_neg_close; varfloor=0.0) ≥ 0.0  # floor still applies

    # Shape preservation
    A = randn(3, 4)
    outA = varpredbias(A, p1)
    @test size(outA) == size(A)
    @test all(isfinite, outA)

    # Unitful inputs return squared intensity units
    p_units = QuadVarParams(1.0u"pA^2", 0.25u"pA", 0.1)
    y_unit = 2.0u"pA"
    res_units = varpredbias(y_unit, p_units)
    @test res_units isa Unitful.AbstractQuantity
    @test Unitful.dimension(res_units) == Unitful.dimension(1.0u"pA^2")
    p_units_floor = QuadVarParams(-0.5u"pA^2", 0.0u"pA", 0.0)
    @test varpredbias(0.0u"pA", p_units_floor; varfloor=1e-9) == 1e-9u"pA^2"
end

# ── vif ──────────────────────────────────────────────────────────────────────

@testset "vif()" begin
    # Basic numerics and defaults
    @test vif(0.0, 1) == 1.0
    @test vif(0.5, 5) == 1 + 2 * 0.5 * (5 - 1) / 5  # 1.8
    @test vif(0.8, 10) == 1 + 2 * 0.8 * 9 / 10      # 2.44
    @test vif(-0.3, 7) == 1.0                       # negatives clipped to 0
    @test vif(0.99, 50) == 1 + 2 * 0.8 * 49 / 50    # capped to ρmax=0.8
    @test vif(0.2, 1) == 1.0                        # (n-1)/n = 0 when n=1

    # Keyword: ρmax capping
    @test vif(0.99, 10; ρmax=0.6) == 1 + 2 * 0.6 * 9 / 10
    @test vif(0.1, 10; ρmax=0.0) == 1.0  # everything capped to 0

    # Keyword: nonnegative=false (symmetric cap)
    @test vif(-0.05, 20; nonnegative=false) == 1.0  # clamp to ≥ 1
    @test vif(-0.9, 20; nonnegative=false) == 1.0
    @test vif(0.3, 10; nonnegative=false) == 1 + 2 * 0.3 * 9 / 10

    # Keyword: nmin (effective sample size)
    @test vif(0.5, 3; nmin=10) == 1 + 2 * 0.5 * (10 - 1) / 10
    @test vif(0.5, 12; nmin=10) == 1 + 2 * 0.5 * 11 / 12

    # Input sanitization for ρ
    @test vif(NaN, 8) == 1.0
    @test vif(Inf, 8) == 1.0
    @test vif(-Inf, 8) == 1.0

    # Invariants and errors
    for ρ in (-1.0, -0.5, 0.0, 0.3, 0.8, 1.0), n in (1, 2, 5, 50)
        @test vif(ρ, n) ≥ 1.0
    end
    @test_throws ArgumentError vif(0.1, 10; ρmax=-0.1)
    @test_throws ArgumentError vif(0.1, 10; ρmax=1.1)
    @test_throws ArgumentError vif(0.1, 10; nmin=0)

    # Monotonicity (where applicable)
    @test vif(0.2, 20) < vif(0.4, 20) ≤ vif(0.8, 20)  # increasing in ρ up to cap
    @test vif(0.3, 2) < vif(0.3, 5) < vif(0.3, 50)  # increases with n for ρ>0
    @test all(vif(-0.2, n) == 1.0 for n in (1, 2, 10, 100))  # clipped when nonnegative=true

    # Type behavior / stability
    @test typeof(vif(1, 5)) ≡ Float64
    @test typeof(vif(Float32(0.2), 5)) ≡ Float32
    @test typeof(vif(Float64(0.2), 5)) ≡ Float64
    @test typeof(vif(BigFloat(0.2), 5)) ≡ BigFloat

    # Inference + Float32 value check (tolerant to 32-bit rounding)
    val32 = @inferred vif(Float32(0.2), 5)
    @test isapprox(val32, Float32(1 + 2 * Float32(0.2) * (5 - 1) / 5); atol=eps(Float32))

    # Large-n behavior: check against the finite-n formula (not the limit)
    nbig = 1_000_000
    @test isapprox(vif(0.7,  nbig), 1 + 2 * 0.7 * (nbig - 1) / nbig; rtol=1e-12)
    @test isapprox(vif(0.95, nbig), 1 + 2 * 0.8 * (nbig - 1) / nbig; rtol=1e-12)
    @test vif(-0.95, nbig; nonnegative=false) == 1.0
end

end  # module
