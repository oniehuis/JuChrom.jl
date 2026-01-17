
# =============================================================================
# OSQP-based nonnegative least squares with OPTIONAL quadratic prior on a
# =============================================================================
"""
    nnls_osqp_prior(F, y; ridge=1e-12, μ=nothing, λ=nothing)

Solve a small nonnegative least squares problem (NNLS) using OSQP:

    minimize_a  ||F*a - y||^2  +  sum_k λ[k]*(a[k] - μ[k])^2
    subject to  a >= 0

This is used to update the per-channel component weights `A[i, :]` given the
current component shapes evaluated at that channel's sample times.

Arguments
- `F::Matrix{Float64}`: (n_scans × Kpeaks) design matrix. Column k is component
  k evaluated at this channel's scan times.
- `y::Vector{Float64}`: (n_scans) observed intensities for one channel.

Keyword arguments
- `ridge`: tiny diagonal added for numerical stability.
- `μ`: optional prior mean for `a` (length Kpeaks). Use when you have an expected
  spectrum for some peaks and want a SOFT pull toward it.
- `λ`: optional nonnegative weights (length Kpeaks). Larger `λ[k]` means stronger
  pull of `a[k]` toward `μ[k]`. Set `λ[k]=0` to disable prior for peak k.

Typical settings
- No prior: leave `μ=nothing, λ=nothing`.
- Soft prior for a subset of peaks: pass `μ` and `λ` with zeros for peaks you
  don't want to regularize.
"""
function nnls_osqp_prior(
    F::Matrix{Float64},
    y::Vector{Float64};
    ridge::Float64 = 1e-12,
    μ::Union{Nothing, Vector{Float64}} = nothing,
    λ::Union{Nothing, Vector{Float64}} = nothing
)
    K = size(F, 2)
    P = F' * F
    q = -(F' * y)

    # numerical stability
    P .+= ridge .* I(K)

    # soft prior: sum_k λ_k (a_k - μ_k)^2
    # contributes: +Λ to P and -Λ*μ to q
    if μ !== nothing && λ !== nothing
        @assert length(μ) == K
        @assert length(λ) == K
        Λ = Diagonal(λ)
        P .+= Λ
        q .-= Λ * μ
    end

    A = sparse(I(K))
    l = zeros(K)
    u = fill(Inf, K)

    model = OSQP.Model()
    OSQP.setup!(model;
        P = sparse(Symmetric(P)),
        q = q,
        A = A,
        l = l,
        u = u,
        verbose = false,
        polish = true
    )
    res = OSQP.solve!(model)
    return res.x
end

# =============================================================================
# Multi-peak deconvolution: shared unimodal B-spline shapes + nonnegative spectra
# =============================================================================
"""
    fit_shared_unimodal_bspline_shifted(Y, t_actual, t0_times;
        A_init=nothing,
        lock_Azeros=nothing,
        A_prior=nothing,
        lambda_peaks=nothing,
        K=30,
        tgrid_n=500,
        iters=30,
        ridge=1e-8,
        nnls_ridge=1e-12)

Fit Kpeaks chromatographic peaks (components) to GC/MS intensity data where each
m/z channel can be measured at a slightly shifted time within each scan.

Model
For channel i and scan j:

    Y[i,j] ≈ sum_{k=1..Kpeaks}  A[i,k] * f_k(t_actual[i,j])  + noise

- `A[i,k] >= 0` is the contribution ("height" / spectrum weight) of peak k in
  channel i.
- `f_k(t) >= 0` is the shared peak shape of component k over continuous time.
- Each f_k is parameterized as a cubic B-spline:
      f_k(t) = B(t) * c_k

Unimodality constraint (per component k)
- Using `t0_times[k]` as the apex switch time:
    f'_k(t) >= 0 for t <= t0_times[k]
    f'_k(t) <= 0 for t >= t0_times[k]
  plus f_k(t) >= 0 on a dense time grid.

Algorithm (alternating optimization)
1) Given shapes {f_k}, update each channel's weights A[i,:] by solving a small
   nonnegative LS (OSQP) possibly with a SOFT prior.
2) Given A, update spline coefficients {c_k} by solving one larger constrained
   QP (OSQP) enforcing positivity + unimodality per component.
3) Normalize each component to max(f_k)=1 and push scale into A[:,k].
Repeat for `iters` iterations.

Arguments
- `Y::Matrix{Float64}`: (n_mz × n_scans) intensity matrix.
- `t_actual::Matrix{Float64}`: (n_mz × n_scans) actual acquisition times per
  channel per scan (accounts for within-scan offsets).
- `t0_times::Vector{Float64}`: length Kpeaks. Approximate apex time for each peak.

Keyword arguments
- `A_init`: optional initial A (n_mz × Kpeaks). Good initialization matters when
  peaks overlap. If `nothing`, starts near-zero.
- `lock_Azeros`: optional BitMatrix (n_mz × Kpeaks). If `true` at (i,k), enforces
  A[i,k]=0 at every iteration. Use for "diagnostic ion exclusivity".
- `A_prior`: optional prior mean for A (n_mz × Kpeaks). This is a SOFT preference,
  not a hard constraint. Use when you have a known/expected spectrum for some peaks.
- `lambda_peaks`: optional Vector length Kpeaks with nonnegative strengths λ[k].
  If provided with `A_prior`, the A-update includes penalty:
      sum_k λ[k]*(A[i,k] - A_prior[i,k])^2
  Typical:
    - Start small: λ ~ 0.1..1.0
    - Strong prior: λ ~ 10..100
    - For peaks without a prior: set λ[k]=0
- `K`: number of knots in the B-spline basis. Larger => more flexibility, but can
  overfit and destabilize deconvolution when peaks overlap. Start 20..40.
- `tgrid_n`: number of grid points used to enforce positivity/unimodality in time.
  Larger => tighter constraints but slower. Start 200..600.
- `iters`: alternating iterations. 20..60 usually enough.
- `ridge`: small ridge added to the spline QP for stability.
- `nnls_ridge`: small ridge added to the A-update QP for stability.

Returns
- `A_hat`: (n_mz × Kpeaks) estimated spectra/weights.
- `basis`: BSplineBasis used for shapes.
- `C_hat`: (p × Kpeaks) spline coefficients.
- `(tgrid, Fhat_grid)`: time grid and fitted shapes on that grid (tgrid_n × Kpeaks).
- `Yfit`: fitted Y at observed times (n_mz × n_scans).
"""
function fit_shared_unimodal_bspline_shifted(
    Y::Matrix{Float64},
    t_actual::Matrix{Float64},
    t0_times::Vector{Float64};
    A_init::Union{Nothing, Matrix{Float64}} = nothing,
    lock_Azeros::Union{Nothing, BitMatrix} = nothing,
    A_prior::Union{Nothing, Matrix{Float64}} = nothing,
    lambda_peaks::Union{Nothing, Vector{Float64}} = nothing,
    K::Int = 30,
    tgrid_n::Int = 500,
    iters::Int = 30,
    ridge::Float64 = 1e-8,
    nnls_ridge::Float64 = 1e-12
)
    n_mz, n_scans = size(Y)
    @assert size(t_actual) == (n_mz, n_scans)

    Kpeaks = length(t0_times)
    @assert Kpeaks >= 1

    if A_init !== nothing
        @assert size(A_init) == (n_mz, Kpeaks)
    end
    if lock_Azeros !== nothing
        @assert size(lock_Azeros) == (n_mz, Kpeaks)
    end
    if A_prior !== nothing
        @assert size(A_prior) == (n_mz, Kpeaks)
    end
    if (A_prior === nothing) != (lambda_peaks === nothing)
        error("Provide both A_prior and lambda_peaks, or neither.")
    end
    if lambda_peaks !== nothing
        @assert length(lambda_peaks) == Kpeaks
        @assert all(lambda_peaks .>= 0.0)
    end

    # Flatten times in the same order as vec(Y) (column-major)
    tvec = vec(t_actual)
    tmin, tmax = minimum(tvec), maximum(tvec)
    Nobs = length(tvec)

    # ---- Spline basis ----
    knots = range(tmin, tmax; length=K)
    basis = BSplineBasis(BSplineOrder(4), knots)  # cubic
    Bobs = collocation_matrix(basis, tvec, Derivative(0), Matrix{Float64})  # (Nobs x p)
    p = size(Bobs, 2)

    # ---- Constraint grid ----
    tgrid = range(tmin, tmax; length=tgrid_n)
    Bgrid  = collocation_matrix(basis, tgrid, Derivative(0), Matrix{Float64})
    Bprime = collocation_matrix(basis, tgrid, Derivative(1), Matrix{Float64})

    # ---- Build block-diagonal constraints for all peaks ----
    Acon_blocks = Vector{SparseMatrixCSC{Float64,Int}}(undef, Kpeaks)
    l_blocks = Vector{Vector{Float64}}(undef, Kpeaks)
    u_blocks = Vector{Vector{Float64}}(undef, Kpeaks)

    for k in 1:Kpeaks
        t0 = t0_times[k]
        left  = tgrid .<= t0
        right = tgrid .>= t0

        Acon_k = vcat(
            sparse(Bgrid),                 # f(t) >= 0
            sparse(Bprime[left, :]),       # f'(t) >= 0 on left
            sparse(-Bprime[right, :])      # f'(t) <= 0 on right
        )
        Acon_blocks[k] = Acon_k
        l_blocks[k] = zeros(size(Acon_k, 1))
        u_blocks[k] = fill(Inf, size(Acon_k, 1))
    end

    function blockdiag(mats::Vector{SparseMatrixCSC{Float64,Int}})
        rows = sum(size(M,1) for M in mats)
        cols = sum(size(M,2) for M in mats)
        out = spzeros(rows, cols)
        r0 = 1
        c0 = 1
        for M in mats
            rr, cc = size(M)
            out[r0:(r0+rr-1), c0:(c0+cc-1)] = M
            r0 += rr
            c0 += cc
        end
        return out
    end

    Acon_all = blockdiag(Acon_blocks)
    lcon = vcat(l_blocks...)
    ucon = vcat(u_blocks...)

    # ---- Initialize A ----
    A_hat = if A_init === nothing
        fill(1e-3, n_mz, Kpeaks)
    else
        max.(A_init, 0.0)
    end
    if lock_Azeros !== nothing
        A_hat[lock_Azeros] .= 0.0
    end

    # ---- Initialize C (p x Kpeaks) ----
    σ0 = (tmax - tmin) / 12
    C = zeros(p, Kpeaks)
    tg = collect(tgrid)
    for k in 1:Kpeaks
        bump = exp.(-0.5 .* ((tg .- t0_times[k]) ./ σ0).^2)
        bump ./= maximum(bump) > 0 ? maximum(bump) : 1.0
        C[:, k] = Bgrid \ bump
        fk = Bgrid * C[:, k]
        fmax = maximum(fk)
        if fmax > 0
            C[:, k] ./= fmax
        end
    end

    yvec = vec(Y)

    for _ in 1:iters
        # ---- Evaluate shapes at observations ----
        Fobs_vec = Bobs * C  # (Nobs x Kpeaks)

        # ---- Update A via per-channel OSQP NNLS (+ optional soft prior) ----
        for i in 1:n_mz
            Fi = zeros(n_scans, Kpeaks)
            for k in 1:Kpeaks
                @inbounds for j in 1:n_scans
                    Fi[j, k] = Fobs_vec[i + (j-1)*n_mz, k]
                end
            end

            yi = vec(Y[i, :])

            if A_prior === nothing
                A_hat[i, :] .= nnls_osqp_prior(Fi, yi; ridge=nnls_ridge)
            else
                μ_i = vec(A_prior[i, :])
                λ_i = vec(lambda_peaks)
                A_hat[i, :] .= nnls_osqp_prior(Fi, yi; ridge=nnls_ridge, μ=μ_i, λ=λ_i)
            end
        end

        if lock_Azeros !== nothing
            A_hat[lock_Azeros] .= 0.0
        end

        # ---- Update C via one big QP in x = vec(C) ----
        BW_all = zeros(Nobs, p * Kpeaks)
        for k in 1:Kpeaks
            w = repeat(view(A_hat, :, k), n_scans)     # length Nobs
            BWk = Bobs .* w                             # (Nobs x p)
            BW_all[:, ((k-1)*p+1):(k*p)] .= BWk
        end

        P = BW_all' * BW_all
        P .+= ridge .* I(size(P, 1))
        q = -(BW_all' * yvec)

        model = OSQP.Model()
        OSQP.setup!(model;
            P = sparse(Symmetric(P)),
            q = q,
            A = sparse(Acon_all),
            l = lcon,
            u = ucon,
            verbose = false,
            polish = true
        )
        res = OSQP.solve!(model)
        x = res.x

        for k in 1:Kpeaks
            C[:, k] = view(x, ((k-1)*p+1):(k*p))
        end

        # ---- Normalize each component (max on grid = 1) ----
        for k in 1:Kpeaks
            fk_grid = Bgrid * C[:, k]
            fmax = maximum(fk_grid)
            if fmax > 0
                C[:, k] ./= fmax
                A_hat[:, k] .*= fmax
            end
        end
    end

    # ---- Final outputs ----
    Fhat_grid = Bgrid * C  # (tgrid_n x Kpeaks)

    Fobs_vec = Bobs * C
    Yfit_vec = zeros(Nobs)
    for k in 1:Kpeaks
        w = repeat(view(A_hat, :, k), n_scans)
        Yfit_vec .+= w .* view(Fobs_vec, :, k)
    end
    Yfit = reshape(Yfit_vec, n_mz, n_scans)

    return A_hat, basis, C, (tgrid, Fhat_grid), Yfit
end

"""
    fit_with_t0_search(
        Y, t_actual;
        t0_guess,
        half_width=1.5,
        ngrid=21,
        min_sep=0.2,
        max_sep=nothing,
        strategy=:single,
        stages=nothing,
        adaptive_window=false,
        shrink=0.35,
        min_half_width=0.05,
        coord_sweeps=2,
        enforce_sorted=true,
        max_iter=50,
        tol_sse=1e-6,
        verbose=true,
        kwargs...
    )

Estimate chromatographic apex times `t0[k]` for an arbitrary number of overlapping
peaks using coordinate-wise grid search combined with constrained B-spline
deconvolution.

This function wraps `fit_shared_unimodal_bspline_shifted` and searches for apex
times that minimize the total squared reconstruction error (SSE) of the model.

The method is designed to scale to any number of peaks K without requiring an
exponential K-dimensional grid.

──────────────────────────────────────────────────────────────────────────────
MODEL
──────────────────────────────────────────────────────────────────────────────

For channel i and scan j:

    Y[i,j] ≈ Σₖ A[i,k] · fₖ(t_actual[i,j])

Each peak k has:
  • a shared unimodal shape fₖ(t) modeled by a constrained cubic B-spline
  • a nonnegative mass spectrum A[:,k]

The spline shapes enforce:
  • fₖ(t) ≥ 0
  • fₖ is unimodal around t0[k]

──────────────────────────────────────────────────────────────────────────────
SEARCH STRATEGY
──────────────────────────────────────────────────────────────────────────────

The search uses **coordinate-wise grid optimization**:

1. Start from an initial guess t0_guess (length K).
2. For each peak k:
     - Search a 1D grid around current t0[k]
     - Keep all other t0 fixed
     - Keep the best candidate
3. Repeat for all k (one sweep).
4. Repeat sweeps until convergence.

Optional coarse→fine refinement is supported.

This avoids the exponential cost of a full K-dimensional grid search.

──────────────────────────────────────────────────────────────────────────────
ARGUMENTS
──────────────────────────────────────────────────────────────────────────────

Required:
- `Y::Matrix{Float64}`  
    Intensity matrix (n_mz × n_scans)

- `t_actual::Matrix{Float64}`  
    Actual acquisition time per ion per scan (n_mz × n_scans)

Keyword arguments:

- `t0_guess::Vector{Float64}`  
    Initial guess of apex times (length K). Must be roughly correct.

- `half_width::Float64`  
    Half-width of the initial search window around each t0.

- `ngrid::Int`  
    Number of grid points per coordinate search.

- `min_sep::Float64`  
    Minimum allowed separation between adjacent peaks.

- `max_sep::Union{Nothing,Float64}`  
    Optional maximum separation constraint.

- `strategy::Symbol`  
    `:single`   → one grid search  
    `:iterative` → multi-stage coarse→fine search

- `stages::Union{Nothing,Vector}`  
    Optional explicit refinement stages:
      [(half_width=2.0, ngrid=11),
       (half_width=0.6, ngrid=17),
       (half_width=0.2, ngrid=21)]

- `adaptive_window::Bool`  
    Automatically shrink half-width after each stage.

- `shrink::Float64`  
    Shrink factor for adaptive window (0 < shrink < 1).

- `min_half_width::Float64`  
    Smallest allowed window size when shrinking.

- `coord_sweeps::Int`  
    Number of coordinate sweeps per stage.

- `enforce_sorted::Bool`  
    Enforce t0 ordering (t0[1] < t0[2] < ...). Recommended.

- `max_iter::Int`  
    Global maximum number of coordinate iterations.

- `tol_sse::Float64`  
    Relative SSE improvement tolerance for early stopping.

- `verbose::Bool`  
    Print progress.

- `kwargs...`  
    Passed directly to `fit_shared_unimodal_bspline_shifted`  
    (e.g. A_init, A_prior, lambda_peaks, lock_Azeros, K, tgrid_n, iters, ...)

──────────────────────────────────────────────────────────────────────────────
RETURNS
──────────────────────────────────────────────────────────────────────────────

Returns:

    best_t0s, best_sse, A_hat, basis_hat, C_hat, (tgrid, Fhat_grid), Y_fit

Where:
- `best_t0s`     → Vector of estimated apex times
- `best_sse`     → final reconstruction SSE
- `A_hat`        → estimated mass spectra (n_mz × K)
- `basis_hat`    → B-spline basis
- `C_hat`        → spline coefficients (p × K)
- `tgrid`        → time grid for plotting
- `Fhat_grid`    → fitted shapes on grid (tgrid_n × K)
- `Y_fit`        → reconstructed signal (n_mz × n_scans)

──────────────────────────────────────────────────────────────────────────────
TYPICAL USAGE
──────────────────────────────────────────────────────────────────────────────

Two peaks, coarse→fine search:

    best_t0s, best_sse, A_hat, basis_hat, C_hat, (tgrid, Fhat_grid), Y_fit =
        fit_with_t0_search(Y, t_actual;
            t0_guess=[15.0, 16.0],
            half_width=2.0,
            ngrid=21,
            min_sep=0.2,
            strategy=:iterative,
            adaptive_window=true,
            coord_sweeps=2,
            A_init=A_init,
            A_prior=A_prior,
            lambda_peaks=lambda_peaks,
            K=35,
            tgrid_n=300,
            iters=40
        )

──────────────────────────────────────────────────────────────────────────────
NOTES
──────────────────────────────────────────────────────────────────────────────

• This is a deterministic optimizer (no randomness).
• Complexity scales linearly with K (number of peaks).
• Works well even for heavily overlapping peaks when spectra differ.
• Soft priors and hard diagnostic-ion constraints integrate naturally.

This function is intended as the primary t0 optimizer for multi-peak
GC/MS deconvolution.
"""
function fit_with_t0_search(
    Y, t_actual;
    t0_guess::Vector{Float64},
    half_width::Float64 = 1.5,
    ngrid::Int = 21,
    min_sep::Float64 = 0.2,
    max_sep::Union{Nothing,Float64} = nothing,
    strategy::Symbol = :single,
    stages = nothing,
    adaptive_window::Bool = false,
    shrink::Float64 = 0.35,
    min_half_width::Float64 = 0.05,
    coord_sweeps::Int = 2,
    enforce_sorted::Bool = true,
    max_iter::Int = 50,           # NEW: total coordinate iterations
    tol_sse::Float64 = 1e-6,      # NEW: relative SSE improvement tolerance
    verbose::Bool = true,
    kwargs...
)
    Kpeaks = length(t0_guess)
    @assert Kpeaks >= 1
    @assert half_width > 0
    @assert ngrid >= 3
    @assert coord_sweeps >= 1
    @assert max_iter >= 1
    @assert tol_sse > 0

    # stage plan
    plan = if strategy == :single
        [(half_width=half_width, ngrid=ngrid)]
    elseif strategy == :iterative
        if stages === nothing
            hw1 = half_width
            hw2 = max(min_half_width, 0.35 * hw1)
            hw3 = max(min_half_width, 0.12 * hw1)
            [(half_width=hw1, ngrid=max(9,  min(ngrid, 13))),
             (half_width=hw2, ngrid=max(13, min(ngrid, 19))),
             (half_width=hw3, ngrid=max(17, ngrid))]
        else
            collect(stages)
        end
    else
        error("strategy must be :single or :iterative")
    end

    # initial t0
    t0 = Float64.(t0_guess)
    if enforce_sorted
        sort!(t0)
    end

    # constraint check
    function feasible(t0s::Vector{Float64})
        if enforce_sorted
            for k in 1:(Kpeaks-1)
                dt = t0s[k+1] - t0s[k]
                if dt < min_sep
                    return false
                end
                if max_sep !== nothing && dt > max_sep
                    return false
                end
            end
        end
        return true
    end

    # model evaluation
    function eval_fit(t0s::Vector{Float64})
        A_hat, basis_hat, C_hat, (tgrid, Fhat_grid), Y_fit =
            fit_shared_unimodal_bspline_shifted(Y, t_actual, t0s; kwargs...)
        sse = sum((Y .- Y_fit).^2)
        return sse, (A_hat, basis_hat, C_hat, (tgrid, Fhat_grid), Y_fit)
    end

    @assert feasible(t0)

    best_sse, best_pack = eval_fit(copy(t0))
    best_t0 = copy(t0)

    if verbose
        println("Initial t0 = ", best_t0, "  SSE = ", best_sse)
    end

    prev_hw = nothing
    iter = 0

    for (si, st) in enumerate(plan)
        hw = hasproperty(st, :half_width) ? Float64(getproperty(st, :half_width)) : half_width
        ng = hasproperty(st, :ngrid) ? Int(getproperty(st, :ngrid)) : ngrid

        if adaptive_window && si > 1 && prev_hw !== nothing && stages === nothing
            hw = max(min_half_width, shrink * prev_hw)
        end
        prev_hw = hw

        if verbose
            println("\nStage $si  (half_width=$(round(hw, digits=4)), ngrid=$ng)")
        end

        for sweep in 1:coord_sweeps
            iter += 1
            if iter > max_iter
                if verbose
                    println("Stopping: reached max_iter = $max_iter")
                end
                return best_t0, best_sse, best_pack...
            end

            improved = false
            prev_sse = best_sse

            for k in 1:Kpeaks
                center = best_t0[k]
                grid = range(center - hw, center + hw; length=ng)

                local_best_sse = best_sse
                local_best_t0 = copy(best_t0)
                local_best_pack = best_pack

                for cand in grid
                    t0_try = copy(best_t0)
                    t0_try[k] = Float64(cand)

                    if enforce_sorted
                        sort!(t0_try)
                    end
                    if !feasible(t0_try)
                        continue
                    end

                    sse, pack = eval_fit(t0_try)

                    if sse < local_best_sse
                        local_best_sse = sse
                        local_best_t0 = t0_try
                        local_best_pack = pack
                    end
                end

                if local_best_sse < best_sse
                    best_sse = local_best_sse
                    best_t0 = local_best_t0
                    best_pack = local_best_pack
                    improved = true
                end
            end

            rel_improve = abs(prev_sse - best_sse) / max(prev_sse, 1e-12)

            if verbose
                println("  sweep $sweep → t0=$(best_t0)  SSE=$(best_sse)  Δrel=$(round(rel_improve, sigdigits=3))")
            end

            if rel_improve < tol_sse
                if verbose
                    println("Stopping: SSE improvement < tol_sse = $tol_sse")
                end
                return best_t0, best_sse, best_pack...
            end

            if !improved
                break
            end
        end
    end

    return best_t0, best_sse, best_pack...
end
