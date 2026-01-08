module BaselineTests

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core/baseline.jl
# ─────────────────────────────────────────────────────────────────────────────

using Test
using LinearAlgebra
using Random
using SparseArrays
using Statistics
using JuChrom
using JuChrom: airpls, FitTracker, calculate_fit, update!, stop_optimization,
               baseline, compute_weights, expand_low_weights!,
               build_second_derivative_matrix

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests
# ─────────────────────────────────────────────────────────────────────────────

# ── build_second_derivative_matrix ───────────────────────────────────────────

@testset "build_second_derivative_matrix" begin
    # Uniform grid → standard [1, -2, 1] stencil on rows 2..n-1
    n = 12
    x = collect(1.0:n)
    D = build_second_derivative_matrix(x)
    @test D isa SparseMatrixCSC{Float64, Int}
    @test size(D) == (n, n)
    @test nnz(D) == 3 * (n - 2)

    # Row pattern spot-check on interior rows
    for i in 2:(n-1)
        cols = findall(!iszero, D[i, :])
        @test cols == [i-1, i, i+1]
        @test D[i, i-1] == 1.0
        @test D[i, i]   == -2.0
        @test D[i, i+1] == 1.0
    end

    # Functional checks: constant and linear → second diff ≈ 0; quadratic → 2a
    c = fill(3.14, n)
    lin = collect(0.5:0.5:0.5n)
    a = 0.75
    quad = @. a * (1:n)^2 + 0.3*(1:n) + 2.0

    @test D * c   ≈ zeros(n)
    @test D * lin ≈ zeros(n)
    y = D * quad
    @test all(isapprox.(y[2:end-1], 2a; atol=1e-12))  # interior rows are 2a
    @test y[1] == 0 && y[end] == 0                    # endpoints unused

    # Nonuniform grid should still annihilate linear functions
    xn = cumsum!(zeros(n), [0; randexp(n-1)]) .+ 1.0
    Dn = build_second_derivative_matrix(xn)
    @test isapprox(Dn * (0.2 .* xn .+ 1.0), zeros(n); atol=1e-12, rtol=0)
end

# ── FitTracker / calculate_fit / update! / stop_optimization / baseline ──────

@testset "FitTracker + helpers" begin
    n = 10
    tracker = FitTracker(n, 3, 1e-6)

    # Initial state
    @test tracker.best_fit === Inf
    @test tracker.iteration_of_best == 0
    @test tracker.no_improvement_count == 0
    @test length(baseline(tracker)) == n
    @test all(iszero, baseline(tracker))

    # calculate_fit matches weighted RSS + smoothness quadratic form
    w = ones(n)
    D = build_second_derivative_matrix(collect(1.0:n))
    λ = 0.5
    S = λ * (D' * D)
    y = @. 1.0 + 0.1 * sin(pi*(1:n)/n)
    b = zeros(n)
    rss = sum(w .* (y .- b).^2)
    sm  = dot(b, S * b)
    @test isapprox(calculate_fit(b, y, w, S), rss + sm; rtol=0, atol=1e-12)

    # update! improves when objective drops by > threshold; tracks iteration
    b1 = fill(0.5, n)
    f1 = calculate_fit(b1, y, w, S)
    @test update!(tracker, f1, b1, 1) === true
    @test tracker.best_fit == f1
    @test tracker.iteration_of_best == 1
    @test baseline(tracker) == b1

    # Smaller improvement than threshold → no improvement; count increases
    @test update!(tracker, f1 - 5e-7, b1, 2) === false
    @test tracker.no_improvement_count == 1

    # Strictly better again → reset count and update baseline
    b2 = fill(0.4, n)
    f2 = f1 - 1e-3
    @test update!(tracker, f2, b2, 3) === true
    @test tracker.no_improvement_count == 0
    @test tracker.iteration_of_best == 3
    @test baseline(tracker) == b2

    # Stop after hitting no-improvement limit
    for k in 1:3
        update!(tracker, f2 - 1e-7, b2, 3 + k)  # not enough improvement
    end
    @test stop_optimization(tracker) === true

    # Bad construction arguments (assertions inside inner ctor)
    @test_throws AssertionError FitTracker(n, 0, 1e-6)
    @test_throws AssertionError FitTracker(n, 2, -1e-9)
end

# ── compute_weights / expand_low_weights! ────────────────────────────────────

@testset "compute_weights / expand_low_weights!" begin
    Random.seed!(1234)
    n = 60
    # Residuals with a single positive peak on otherwise small-noise baseline
    res = 0.02 .* randn(n)
    res[30] += 5.0   # sharp positive residual (strong peak)
    zero_mask = falses(n)

    # Case 1: std_devs = nothing; positive peak gets a low weight
    w = compute_weights(res, nothing, 1.96, zero_mask, 1e-12)
    @test length(w) == n
    @test all(0 .< w .≤ 1)
    @test w[30] == minimum(w)
    @test mean(w) > 0.5  # only a few points are penalized

    # Case 2: providing std_devs rescales residuals; result still in (0,1]
    σ = abs.(0.2 .+ 0.1 .* sin.(range(0, 2π; length=n)))
    w2 = compute_weights(res, σ, 1.5, zero_mask, 1e-12)
    @test all(0 .< w2 .≤ 1)
    @test w2[30] ≤ w[30] + 1e-12  # still among the smallest

    # Case 3: not enough negative residuals → returns ones (early guard)
    res_pos = abs.(randn(n)) .+ 0.1
    wpos = compute_weights(res_pos, nothing, 2.0, zero_mask, 1e-12)
    @test wpos == ones(n)

    # expand_low_weights! should propagate minimal weights over descending shoulders
    wmin = copy(w)
    # Force a very small core weight at index 25 and make residuals decrease to left/right
    wmin .= 1.0
    wmin[25] = 1e-6
    res2 = zeros(n)
    # Create a local "hill" (residual descends away from 25 over ~6 points)
    for j in 1:n
        res2[j] = -abs(j - 25)  # decreases as we move away (peak at center)
    end
    expand_low_weights!(wmin, res2)
    # Expect at least a contiguous block ≥5 around 25 to be set to that minimum
    block = findall(x -> x == minimum(wmin), wmin)
    @test !isempty(block)
    @test maximum(diff(block)) == 1  # contiguous
    @test length(block) ≥ 5
end

# ── airpls(retentions, intensities; ...) ─────────────────────────────────────

@testset "airpls(::AbstractVector, ::AbstractVector; ...)" begin
    # Synthetic baseline + positive peaks
    Random.seed!(7)
    n = 400
    x = range(0.0, 10.0; length=n) |> collect
    # True smooth baseline (quadratic drift + gentle sinusoid)
    btrue = @. 0.5 + 0.02 * x^2 + 0.05 * sin(0.6 * x)
    y = copy(btrue)

    # Add a few positive peaks (Gaussian bumps) + small noise; enforce nonnegativity
    function add_peak!(y, x0, amp, width)
        @. y += amp * exp(-((x - x0)^2) / (2 * width^2))
        nothing
    end
    add_peak!(y, 2.0,  2.0, 0.10)
    add_peak!(y, 5.0,  3.0, 0.15)
    add_peak!(y, 7.5,  1.5, 0.20)
    y .+= 0.01 .* randn(n)
    y .= max.(y, 0.0)

    # Run airPLS with reasonably strong smoothing and default asymmetry
    b̂ = airpls(x, y; λ=1e6, threshold_factor=1.96, max_iter=200, no_improvement_limit=10)

    # Basic contract
    @test length(b̂) == n
    @test all(isfinite, b̂)
    @test all(b̂ .≥ 0)

    # Accuracy on the smooth background (allowing moderate error due to strong peaks)
    rel_rms = sqrt(mean(((b̂ .- btrue) ./ max.(btrue, eps(Float64))).^2))
    @test rel_rms ≤ 0.20

    # Validation errors
    @test_throws ArgumentError airpls(x[1:10], y)                         # length mismatch
    @test_throws ArgumentError airpls(x[1:2], y[1:2])                     # n < 3
    @test_throws ArgumentError airpls(x, y; λ=0.0)                        # λ must be > 0
    @test_throws ArgumentError airpls(x, y; threshold_factor=0.0)         # > 0
    @test_throws ArgumentError airpls(x, y; zero_threshold=-1e-9)         # ≥ 0
    @test_throws ArgumentError airpls(x, y; zero_weight=0.0)              # > 0
    @test_throws ArgumentError airpls(x, y; variances=fill(-1.0, n))      # variances ≥ 0
    @test_throws ArgumentError airpls(x, y; variances=ones(n-1))          # length match

    # With measurement variances supplied (heteroscedastic), should still behave
    vars = @. 1e-4 + 1e-4 * (1 + sin(0.5 * x))^2
    b̂v = airpls(x, y; variances=vars, λ=5e5, max_iter=200)
    @test length(b̂v) == n && all(isfinite, b̂v)
end

end # module
