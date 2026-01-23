"""
    nnlspenalized(
        F::AbstractMatrix{<:Real},
        y::AbstractVector{<:Real};
        ridge::Real=1e-12,
        μ=nothing,
        λ=nothing,
        w=nothing
    ) -> Vector{Float64}

Return the optimal nonnegative coefficient vector `a` (length `Kpeaks`) from solving
a small nonnegative least squares problem (NNLS) using OSQP, a quadratic programming solver:

    minimize_a  ||F*a - y||^2  +  sum_k λ[k]*(a[k] - μ[k])^2
    subject to  a >= 0.

The optional quadratic prior is controlled by `μ` and `λ`. When provided, it softly
pulls the solution toward `μ` with per-component strength `λ`. Set `λ[k]=0` to
disable the prior for component `k`.

Optional heteroscedastic weights `w` apply per observation. They are implemented
by row-scaling `F` and `y` by `sqrt.(w)`, which converts the weighted problem to
an equivalent unweighted least squares system.

`F` is the (`n_scans` × `Kpeaks`) design matrix where column `k` is component `k` evaluated
at the retentions, and `y` is the length-`n_scans` observed intensities for one channel.
`ridge` adds a small diagonal regularization for numerical stability.

See also
[`unimodalfit`](@ref),
[`unimodalfit_apexsearch`](@ref).
"""
function nnlspenalized(
    F::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real};
    ridge::T1=1e-12,
    μ::T2=nothing,
    λ::T3=nothing,
    w::T4=nothing
    ) where {
        T1<:Real,
        T2<:Union{Nothing, AbstractVector{<:Real}},
        T3<:Union{Nothing, AbstractVector{<:Real}},
        T4<:Union{Nothing, AbstractVector{<:Real}}
    }

    # Number of components (columns of F).
    K = size(F, 2)

    # Optional row weighting via sqrt(w): converts weighted LS to unweighted LS.
    Fw = F
    yw = y
    if !isnothing(w)
        if length(w) ≠ length(y)
            throw(DimensionMismatch(
                "w must have the same length as y (length(w)=$(length(w)) " *
                "vs length(y)=$(length(y)))."))
        end
        if any(w .< 0.0)
            throw(ArgumentError("w must be nonnegative."))
        end
        sw = sqrt.(w)
        Fw = F .* sw
        yw = y .* sw
    end

    # Quadratic form for ||Fw*a - yw||^2 is (1/2) a'Pa + q'a (OSQP uses this form).
    P = Fw' * Fw
    q = -(Fw' * yw)

    # Numerical stability: small diagonal ridge to keep P PSD in finite precision.
    P .+= ridge .* I(K)

    # Soft prior: sum_k λ_k (a_k - μ_k)^2
    # Contributes: +Λ to P and -Λ*μ to q
    if !isnothing(μ) && !isnothing(λ)
        if length(μ) ≠ K
            throw(DimensionMismatch("μ must have length K=$(K)."))
        end
        if length(λ) ≠ K
            throw(DimensionMismatch("λ must have length K=$(K)."))
        end
        Λ = Diagonal(λ)
        P .+= Λ
        q .-= Λ * μ
    end

    # Nonnegativity constraints: a ≥ 0.
    A = sparse(I(K))
    l = zeros(K)
    u = fill(Inf, K)

    # Solve the NNLS quadratic program with OSQP.
    model = OSQP.Model()
    OSQP.setup!(model;
        P=sparse(Symmetric(P)),
        q=q,
        A=A,
        l=l,
        u=u,
        verbose=false,
        polish=true
    )
    res = OSQP.solve!(model)

    # Return the optimal coefficients a.
    res.x
end

function sigma2weights(
    σ::AbstractVector{<:Real};
    σ_floor_rel::T1=1e-3,
    σ_floor_abs::T2=1e-12,
    w_cap::T3=1e6
    ) where {
        T1<:Real,
        T2<:Real,
        T3<:Real
    }

    # Validate σ vector is finite.
    if any(.!isfinite.(σ))
        throw(ArgumentError("σ must contain only finite values."))
    end
    # Validate σ vector is positive.
    if any(σ .≤ 0.0)
        throw(ArgumentError("σ must be strictly positive."))
    end

    # Largest σ value for relative flooring.
    smax = maximum(σ)
    # Compute effective floor for σ values.
    σ_floor = max(σ_floor_abs, σ_floor_rel * smax)
    # Apply the floor to σ values.
    σ_eff = max.(σ, σ_floor)

    # Convert sigma to weights (inverse variance).
    w = 1.0 ./ (σ_eff .^ 2)

    # Normalize typical magnitude to ~1 (does not change minimizer).
    wm = median(w)
    if !isfinite(wm) || wm ≤ 0.0
        wm = mean(w)
    end
    if !isfinite(wm) || wm ≤ 0.0
        return nothing
    end
    w ./= wm

    # Cap extreme weights for numerical stability.
    w = min.(w, w_cap)
    if !all(isfinite.(w))
        return nothing
    end

    w
end

"""
    unimodalfit(
        Iobs::AbstractMatrix{<:Real},
        R::AbstractMatrix{<:Real},
        peakretentions::AbstractVector{<:Real};
        σ::Union{Nothing, AbstractMatrix{<:Real}}=nothing,
        spectra_init::Union{Nothing, AbstractMatrix{<:Real}}=nothing,
        spectra_zero_mask::Union{Nothing, BitMatrix}=nothing,
        spectra_prior::Union{Nothing, AbstractMatrix{<:Real}}=nothing,
        spectra_prior_weights::Union{Nothing, AbstractVector{<:Real}}=nothing,
        knots_n::Integer=30,
        rgrid_n::Integer=500,
        iters::Integer=30,
        shape_ridge::Real=1e-8,
        spectra_ridge::Real=1e-12,
        shape_couple::Real=0.0,
        shape_couple_mode::Symbol=:d2,
        shape_couple_graph::Symbol=:neighbors,
        shape_couple_window::Real=0.0,
        shape_couple_tau_halfwidth::Real=1.0,
        shape_couple_tau_n::Integer=101,
        apex_localize::Real=0.0,
        apex_localize_scale::Real=1.0
    )

Fit `K` chromatographic peaks (components) to GC/MS intensity data `Iobs` where each m/z 
channel can be measured at a slightly shifted retention within each scan. The model for 
channel `i` and scan `j` is

    Iobs[i,j] ≈ sum_{k=1..K} A[i,k] * f_k(R[i,j]) + noise

where `A[i,k] ≥ 0` is the contribution of peak `k` in channel `i` and `f_k(r) ≥ 0` is the
shared peak shape of component `k` over continuous retention `r`. Each `f_k` is 
parameterized as a cubic B-spline `f_k(r) = B(r) * c_k`.

Unimodality is enforced using `peakretentions` as the apex switch retention. The
constraints  require `f'_k(r) ≥ 0` for `r ≤ peakretentions[k]`, `f'_k(r) ≤ 0` for 
`r ≥ peakretentions[k]`, and `f_k(r) ≥ 0` on a dense retention grid.

The optimization alternates between updating `A` and updating the spline coefficients `C`. 
Given shapes `f_k`, each channel's weights `A[i,:]` are updated by solving a nonnegative 
least squares problem (optionally with a soft prior). Given `A`, the spline coefficients 
are updated by solving a constrained quadratic program that enforces positivity and 
unimodality. Each component is normalized so `max(f_k)=1` and the scale is absorbed into 
`A[:,k]`. This alternation repeats for `iters` steps.

Heteroscedasticity is optional. If `σ` is provided with the same size as `Iobs`, it is used 
in the `A`-update via weights `w = 1 ./ σ.^2`. The shape update remains unweighted to avoid
bias in the peak height when the noise model is approximate.

When `shape_couple > 0`, a convex quadratic penalty is added in the shape-update (the QP
solving for `C`) to discourage nearby peaks from taking very different shapes. 
`shape_couple_mode = :d1` couples first derivatives and `shape_couple_mode = :d2` couples 
second derivatives. `shape_couple_graph` controls which peak pairs are coupled: `:neighbors` 
couples only adjacent peaks in index order (k with k±1). `:window` couples all pairs whose 
apex retentions differ by at most `shape_couple_window` (requires 
`shape_couple_window > 0`). For derivative coupling (`shape_couple_mode = :d1` or `:d2`), 
comparisons are done in an apex-centered coordinate, i.e. it compares 
`f_k^{(d)}(peakretentions[k] + τ)` across peaks for the same offsets `τ`. The offset grid 
is `τ = range(-shape_couple_tau_halfwidth, shape_couple_tau_halfwidth; 
length=shape_couple_tau_n)`; offsets that push `peakretentions[k] + τ` outside the observed 
`[tmin, tmax]` range are skipped pairwise. When `apex_localize > 0`, an additional convex 
quadratic penalty is added in the shape-update to discourage each peak shape from placing 
mass far away from its own apex:

    apex_localize * sum_k  ∑_g  w_k[g] * f_k(rgrid[g])^2

with

    w_k[g] = ((rgrid[g] - peakretentions[k]) / apex_localize_scale)^2.

This tends to reduce broad/flat solutions and makes it harder for one peak to explain 
signal in the neighborhood of another peak. `apex_localize_scale` should be in the same 
retention units as `R`/`peakretentions` (e.g. minutes); it sets how quickly the penalty 
grows with distance.

Returns `A_hat` (estimated spectra), `basis` (B-spline basis), `C_hat` (spline
coefficients), `(rgrid, Fhat_grid)` for plotting shapes, and `Ifit` (reconstructed
signal at observed retentions).

See also [`unimodalfit_apexsearch`](@ref).
"""
function unimodalfit(
    Iobs::AbstractMatrix{<:Real},
    R::AbstractMatrix{<:Real},
    peakretentions::AbstractVector{<:Real};
    σ::T1=nothing,
    spectra_init::T2=nothing,
    spectra_zero_mask::T3=nothing,
    spectra_prior::T4=nothing,
    spectra_prior_weights::T5=nothing,
    knots_n::T6=30,
    rgrid_n::T7=500,
    iters::T8=30,
    shape_ridge::T9=1e-8,
    spectra_ridge::T10=1e-12,
    shape_couple::T11=0.0,
    shape_couple_mode::Symbol=:d2,
    shape_couple_graph::Symbol=:neighbors,
    shape_couple_window::T12=0.0,
    shape_couple_tau_halfwidth::T13=1.0,
    shape_couple_tau_n::T14=101,
    apex_localize::T15=0.0,
    apex_localize_scale::T16=1.0
    ) where {
        T1<:Union{Nothing, AbstractMatrix{<:Real}},
        T2<:Union{Nothing, AbstractMatrix{<:Real}},
        T3<:Union{Nothing, BitMatrix},
        T4<:Union{Nothing, AbstractMatrix{<:Real}},
        T5<:Union{Nothing, AbstractVector{<:Real}},
        T6<:Integer,
        T7<:Integer,
        T8<:Integer,
        T9<:Real,
        T10<:Real,
        T11<:Real,
        T12<:Real,
        T13<:Real,
        T14<:Integer,
        T15<:Real,
        T16<:Real
}
    n_mz, n_scans = size(Iobs)
    if size(R) ≠ (n_mz, n_scans)
        throw(DimensionMismatch("R must have size ($(n_mz), $(n_scans))."))
    end

    if !isnothing(σ)
        if size(σ) ≠ (n_mz, n_scans)
            throw(DimensionMismatch("σ must have size ($(n_mz), $(n_scans))."))
        end
        if any(.!isfinite.(σ))
            throw(ArgumentError("σ must contain only finite values."))
        end
        if any(σ .≤ 0.0)
            throw(ArgumentError("σ must be strictly positive."))
        end
    end

    if shape_couple < 0
        throw(ArgumentError("shape_couple must be nonnegative."))
    end
    if !(shape_couple_mode in (:d1, :d2))
        throw(ArgumentError("shape_couple_mode must be :d1 or :d2."))
    end
    if !(shape_couple_graph in (:neighbors, :window))
        throw(ArgumentError("shape_couple_graph must be :neighbors or :window."))
    end
    if shape_couple_graph == :window && shape_couple > 0 && shape_couple_window ≤ 0
        throw(ArgumentError("shape_couple_window must be > 0 when shape_couple_graph == :window and shape_couple > 0."))
    end
    if shape_couple > 0 && shape_couple_mode in (:d1, :d2)
        if shape_couple_tau_halfwidth ≤ 0
            throw(ArgumentError("shape_couple_tau_halfwidth must be > 0 when using apex-centered coupling."))
        end
        if shape_couple_tau_n < 3
            throw(ArgumentError("shape_couple_tau_n must be at least 3 when using apex-centered coupling."))
        end
    end

    if apex_localize < 0
        throw(ArgumentError("apex_localize must be nonnegative."))
    end
    if apex_localize > 0 && apex_localize_scale ≤ 0
        throw(ArgumentError("apex_localize_scale must be > 0 when apex_localize > 0."))
    end

    Kpeaks = length(peakretentions)
    if Kpeaks < 1
        throw(ArgumentError("peakretentions must contain at least one peak retention."))
    end

    if  !isnothing(spectra_init)
        if size(spectra_init) ≠ (n_mz, Kpeaks)
            throw(DimensionMismatch("spectra_init must have size ($(n_mz), $(Kpeaks))."))
        end
    end
    if !isnothing(spectra_zero_mask)
        if size(spectra_zero_mask) ≠ (n_mz, Kpeaks)
            throw(DimensionMismatch("spectra_zero_mask must have size ($(n_mz), $(Kpeaks))."))
        end
    end
    if !isnothing(spectra_prior)
        if size(spectra_prior) ≠(n_mz, Kpeaks)
            throw(DimensionMismatch("spectra_prior must have size ($(n_mz), $(Kpeaks))."))
        end
    end

    if isnothing(spectra_prior) ≠ isnothing(spectra_prior_weights)
        throw(ArgumentError("Provide both spectra_prior and spectra_prior_weights, or neither."))
    end
    if  !isnothing(spectra_prior_weights)
        if length(spectra_prior_weights) ≠ Kpeaks
            throw(DimensionMismatch("spectra_prior_weights must have length Kpeaks=$(Kpeaks)."))
        end
        if any(spectra_prior_weights .< 0.0)
            throw(ArgumentError("spectra_prior_weights must be nonnegative."))
        end
    end

    # Flatten retentions in the same order as vec(Iobs) (column-major).
    tvec = vec(R)
    tmin, tmax = minimum(tvec), maximum(tvec)
    Nobs = length(tvec)

    # Spline basis (BSplineKit)
    knots = range(tmin, tmax; length=knots_n)
    basis = BSplineBasis(BSplineOrder(4), knots)  # cubic
    Bobs = collocation_matrix(basis, tvec, Derivative(0), Matrix{Float64})  # (Nobs x p)
    p = size(Bobs, 2)

    # Constraint grid
    rgrid = range(tmin, tmax; length=rgrid_n)
    Bgrid  = collocation_matrix(basis, rgrid, Derivative(0), Matrix{Float64})
    Bprime = collocation_matrix(basis, rgrid, Derivative(1), Matrix{Float64})
    B2     = collocation_matrix(basis, rgrid, Derivative(2), Matrix{Float64})
    rg = Float64.(collect(rgrid))  # concrete vector for weights

    # Build block-diagonal constraints for all peaks
    Acon_blocks = Vector{SparseMatrixCSC{Float64, Int}}(undef, Kpeaks)
    l_blocks = Vector{Vector{Float64}}(undef, Kpeaks)
    u_blocks = Vector{Vector{Float64}}(undef, Kpeaks)

    for k in 1:Kpeaks
        t0 = peakretentions[k]
        left  = rgrid .≤ t0
        right = rgrid .≥ t0

        Acon_k = vcat(
            sparse(Bgrid),                 # f(t) >= 0
            sparse(Bprime[left, :]),       # f'(t) >= 0 on left
            sparse(-Bprime[right, :])      # f'(t) <= 0 on right
        )
        Acon_blocks[k] = Acon_k
        l_blocks[k] = zeros(size(Acon_k, 1))
        u_blocks[k] = fill(Inf, size(Acon_k, 1))
    end

    function blockdiag(mats::Vector{SparseMatrixCSC{Float64, Int}})
        rows = sum(size(M,1) for M in mats)
        cols = sum(size(M,2) for M in mats)
        out = spzeros(rows, cols)
        r0 = 1
        c0 = 1
        for M in mats
            rr, cc = size(M)
            out[r0:(r0 + rr - 1), c0:(c0 + cc - 1)] = M
            r0 += rr
            c0 += cc
        end
        out
    end

    Acon_all = blockdiag(Acon_blocks)
    lcon = vcat(l_blocks...)
    ucon = vcat(u_blocks...)

    # --- helper: build list of coupled peak pairs (edges) ---
    function peak_edges(
        peakretentions::AbstractVector{<:Real},
        graph::Symbol,
        window::Real
    )
        K = length(peakretentions)
        edges = Tuple{Int,Int}[]
        if K ≤ 1
            return edges
        end
        if graph == :neighbors
            for k in 1:(K-1)
                push!(edges, (k, k+1))
            end
            return edges
        elseif graph == :window
            t0 = Float64.(peakretentions)
            for k in 1:K
                for ℓ in (k+1):K
                    if abs(t0[k] - t0[ℓ]) ≤ window
                        push!(edges, (k, ℓ))
                    end
                end
            end
            return edges
        else
            error("shape_couple_graph must be :neighbors or :window")
        end
    end

    edges = (shape_couple > 0) ? peak_edges(peakretentions, shape_couple_graph, shape_couple_window) : Tuple{Int,Int}[]

    # --- precompute coupling Hessian term for C-update (independent of A) ---
    P_couple = zeros(Float64, p*Kpeaks, p*Kpeaks)

    if shape_couple > 0 && !isempty(edges)
        d = (shape_couple_mode == :d1) ? 1 : 2

        τ = collect(range(-shape_couple_tau_halfwidth, shape_couple_tau_halfwidth; length=shape_couple_tau_n))

        rshift = Vector{Vector{Float64}}(undef, Kpeaks)
        inrng  = Vector{BitVector}(undef, Kpeaks)
        for k in 1:Kpeaks
            rk = Float64.(peakretentions[k] .+ τ)
            rshift[k] = rk
            inrng[k] = (rk .>= tmin) .& (rk .<= tmax)
        end

        for (k, ℓ) in edges
            mask = inrng[k] .& inrng[ℓ]
            m = count(mask)
            if m < 3
                continue
            end
            rk = rshift[k][mask]
            rℓ = rshift[ℓ][mask]

            Bk = collocation_matrix(basis, rk, Derivative(d), Matrix{Float64})
            Bℓ = collocation_matrix(basis, rℓ, Derivative(d), Matrix{Float64})

            Gkk = Bk' * Bk
            Gℓℓ = Bℓ' * Bℓ
            Gkℓ = Bk' * Bℓ

            rk0 = (k-1)*p + 1
            ck0 = k*p
            rℓ0 = (ℓ-1)*p + 1
            cℓ0 = ℓ*p
            @views begin
                P_couple[rk0:ck0, rk0:ck0] .+= Gkk
                P_couple[rℓ0:cℓ0, rℓ0:cℓ0] .+= Gℓℓ
                P_couple[rk0:ck0, rℓ0:cℓ0] .-= Gkℓ
                P_couple[rℓ0:cℓ0, rk0:ck0] .-= Gkℓ'
            end
        end
    end

    # --- precompute apex localization Hessian term for C-update (independent of A) ---
    P_apex = zeros(Float64, p*Kpeaks, p*Kpeaks)
    if apex_localize > 0
        s = Float64(apex_localize_scale)
        for k in 1:Kpeaks
            t0 = Float64(peakretentions[k])
            wk = ((rg .- t0) ./ s) .^ 2                       # (rgrid_n)
            # Gk = Bgrid' * Diagonal(wk) * Bgrid
            Gk = Bgrid' * (Bgrid .* wk)                       # (p x p)
            rk0 = (k-1)*p + 1
            ck0 = k*p
            @views P_apex[rk0:ck0, rk0:ck0] .+= Gk
        end
    end

    # Initialize A
    A_hat = if isnothing(spectra_init)
        fill(1e-3, n_mz, Kpeaks)
    else
        max.(spectra_init, 0.0)
    end
    if !isnothing(spectra_zero_mask)
        A_hat[spectra_zero_mask] .= 0.0
    end

    # Initialize C (p x Kpeaks)
    σ0 = (tmax - tmin) / 12
    C = zeros(p, Kpeaks)
    tg = collect(rgrid)
    for k in 1:Kpeaks
        bump = exp.(-0.5 .* ((tg .- peakretentions[k]) ./ σ0).^2)
        bump ./= maximum(bump) > 0 ? maximum(bump) : 1.0
        C[:, k] = Bgrid \ bump
        fk = Bgrid * C[:, k]
        fmax = maximum(fk)
        if fmax > 0
            C[:, k] ./= fmax
        end
    end

    # Vectorize Iobs for the global QP (convert to Float64 even if Iobs is Int).
    yvec = Float64.(vec(Iobs))

    # Alternate between A and C updates.
    for _ in 1:iters
        Fobs_vec = Bobs * C  # (Nobs x Kpeaks)

        # Update A via per-channel OSQP NNLS (+ optional soft prior)
        for i in 1:n_mz
            Fi = zeros(n_scans, Kpeaks)
            for k in 1:Kpeaks
                @inbounds for j in 1:n_scans
                    Fi[j, k] = Fobs_vec[i + (j-1)*n_mz, k]
                end
            end

            yi = Float64.(vec(Iobs[i, :]))

            wi = if isnothing(σ)
                nothing
            else
                sigma2weights(vec(σ[i, :]))
            end

            if isnothing(spectra_prior)
                A_hat[i, :] .= nnlspenalized(Fi, yi; ridge=spectra_ridge, w=wi)
            else
                μ_i = vec(spectra_prior[i, :])
                λ_i = vec(spectra_prior_weights)
                A_hat[i, :] .= nnlspenalized(Fi, yi; ridge=spectra_ridge, μ=μ_i, λ=λ_i, w=wi)
            end

            if any(.!isfinite.(A_hat[i, :]))
                throw(ArgumentError(
                    "Non-finite A_hat produced in channel $i (check σ / weights)."))
            end
        end

        if !isnothing(spectra_zero_mask)
            A_hat[spectra_zero_mask] .= 0.0
        end

        # Update C via one big QP in x = vec(C)
        BW_all = zeros(Nobs, p * Kpeaks)
        for k in 1:Kpeaks
            wA = repeat(view(A_hat, :, k), n_scans)      # length Nobs
            BWk = Bobs .* wA                              # (Nobs x p)
            BW_all[:, ((k - 1) * p + 1):(k * p)] .= BWk
        end

        P = BW_all' * BW_all
        P .+= shape_ridge .* I(size(P, 1))
        q = -(BW_all' * yvec)

        # Add the (precomputed) shape coupling penalty.
        if shape_couple > 0 && !isempty(edges)
            P .+= shape_couple .* P_couple
        end

        # Add the (precomputed) apex localization penalty.
        if apex_localize > 0
            P .+= apex_localize .* P_apex
        end

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
            C[:, k] = view(x, ((k - 1) * p + 1):(k * p))
        end

        # Normalize each component (max on grid = 1) and push scale into A.
        for k in 1:Kpeaks
            fk_grid = Bgrid * C[:, k]
            fmax = maximum(fk_grid)
            if fmax > 0
                C[:, k] ./= fmax
                A_hat[:, k] .*= fmax
            end
        end
    end

    # Final outputs
    Fhat_grid = Bgrid * C  # (rgrid_n x Kpeaks)

    Fobs_vec = Bobs * C
    Ifit_vec = zeros(Float64, Nobs)
    for k in 1:Kpeaks
        wA = repeat(view(A_hat, :, k), n_scans)
        Ifit_vec .+= wA .* view(Fobs_vec, :, k)
    end
    Ifit = reshape(Ifit_vec, n_mz, n_scans)

    A_hat, basis, C, (rgrid, Fhat_grid), Ifit
end

"""
    unimodalfit_apexsearch(Iobs, R; peakretentions_guess, kwargs...)

Wrapper around `unimodalfit` that searches apex retentions using coordinate-wise grid
search to minimize the total sum of squared errors (SSE). All keyword arguments accepted
by `unimodalfit` are forwarded through `kwargs` and apply in the same way as in a direct
call, including `σ`, `spectra_init`, `spectra_zero_mask`, `spectra_prior`,
`spectra_prior_weights`, `knots_n`, `rgrid_n`, `iters`, `shape_ridge`, `spectra_ridge`,
shape coupling keywords `shape_couple`, `shape_couple_mode`, `shape_couple_graph`,
`shape_couple_window`, `shape_couple_tau_halfwidth`, `shape_couple_tau_n`, and apex
localization keywords `apex_localize` and `apex_localize_scale`.

The most important keyword to tune here is `peakretentions_guess`, which sets the initial
apex locations and anchors the unimodality constraints. `half_width` and `ngrid` define
the one-dimensional search grid around each apex; increasing `half_width` expands the
search region, and increasing `ngrid` improves resolution at additional cost. `min_sep`
and `max_sep` enforce separation between adjacent peaks when `enforce_sorted=true`, so
adjust them to reflect expected peak spacing. `strategy=:single` performs one grid search
at the specified resolution, while `strategy=:iterative` runs a coarse-to-fine sequence;
in iterative mode `stages` can override the per-stage `half_width` and `ngrid`, falling
back to the top-level defaults for any missing fields. When using iterative search,
`adaptive_window`, `shrink`, and `min_half_width` control how aggressively the window
narrows between stages and are reasonable to tune if you want faster convergence toward a
known apex. `coord_sweeps` sets the number of coordinate passes per stage; increasing it
can help if the search stalls, but the default is typically sufficient. `enforce_sorted`
should usually remain `true` to keep peak order consistent, and `max_iter` and `tol_sse`
are termination controls that generally do not need adjustment unless you want to cap
runtime or force tighter convergence. `verbose` only affects logging.

For the forwarded `unimodalfit` keywords, it is usually safe to keep the defaults unless
you have a specific modeling reason. `σ` controls heteroscedastic weighting in the `A`
update only; `spectra_init`, `spectra_zero_mask`, `spectra_prior`, and
`spectra_prior_weights` can be tuned when you have strong spectral constraints. The spline
and regularization parameters (`knots_n`, `rgrid_n`, `iters`, `shape_ridge`,
`spectra_ridge`, and the shape coupling or apex localization controls) are advanced knobs
that can change smoothness and overlap behavior, so they are best left untouched unless
you are diagnosing a particular fit artifact.

If `σ` is provided through `kwargs`, it is forwarded to `unimodalfit` and applied only in
the per-channel `A` update. The peak-shape update remains unweighted.

Returns the best apex retention vector, the best SSE, and the same outputs as `unimodalfit`
for that retention vector: `A_hat`, `basis`, `C_hat`, `(rgrid, Fhat_grid)`, and `Ifit`.

See also [`unimodalfit`](@ref).
"""
function unimodalfit_apexsearch(
    Iobs::AbstractMatrix{<:Real},
    R::AbstractMatrix{<:Real};
    peakretentions_guess::Vector{Float64},
    half_width::T2=1.5,
    ngrid::T3=21,
    min_sep::T4=0.2,
    max_sep::T5=nothing,
    strategy::Symbol=:single,
    stages::T6=nothing,
    adaptive_window::Bool=false,
    shrink::T7=0.35,
    min_half_width::T8=0.05,
    coord_sweeps::T9=2,
    enforce_sorted::Bool=true,
    max_iter::T10=50,
    tol_sse::T11=1e-6,
    verbose::Bool=true,
    kwargs...
    ) where {
        T2<:Real,
        T3<:Integer,
        T4<:Real,
        T5<:Union{Nothing, Real},
        T6<:Union{Nothing, Any},
        T7<:Real,
        T8<:Real,
        T9<:Integer,
        T10<:Integer,
        T11<:Real
    }

    Kpeaks = length(peakretentions_guess)
    if Kpeaks < 1
        throw(ArgumentError("peakretentions_guess must contain at least one peak retention."))
    end
    if half_width ≤ 0
        throw(ArgumentError("half_width must be positive."))
    end
    if ngrid < 3
        throw(ArgumentError("ngrid must be at least 3."))
    end
    if coord_sweeps < 1
        throw(ArgumentError("coord_sweeps must be at least 1."))
    end
    if max_iter < 1
        throw(ArgumentError("max_iter must be at least 1."))
    end
    if tol_sse ≤ 0
        throw(ArgumentError("tol_sse must be positive."))
    end

    plan = if strategy == :single
        [(half_width=half_width, ngrid=ngrid)]
    elseif strategy == :iterative
        if isnothing(stages)
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

    t0 = Float64.(peakretentions_guess)
    if enforce_sorted
        sort!(t0)
    end

    function feasible(t0s::Vector{Float64})
        if enforce_sorted
            for k in 1:(Kpeaks-1)
                dt = t0s[k + 1] - t0s[k]
                if dt < min_sep
                    return false
                end
                if !isnothing(max_sep) && dt > max_sep
                    return false
                end
            end
        end
        true
    end

    function eval_fit(t0s::Vector{Float64})
        A_hat, basis_hat, C_hat, (rgrid, Fhat_grid), I_fit =
            unimodalfit(Iobs, R, t0s; kwargs...)
        sse = sum((Float64.(Iobs) .- I_fit).^2)
        sse, (A_hat, basis_hat, C_hat, (rgrid, Fhat_grid), I_fit)
    end

    if !feasible(t0)
        throw(ArgumentError("peakretentions_guess violates min_sep/max_sep constraints."))
    end

    best_sse, best_pack = eval_fit(copy(t0))
    best_t0 = copy(t0)

    if verbose
        println("Initial peakretentions = ", best_t0, "  SSE = ", best_sse)
    end

    prev_hw = nothing
    iter = 0

    for (si, st) in enumerate(plan)
        hw = hasproperty(st, :half_width) ? Float64(getproperty(st, :half_width)) : half_width
        ng = hasproperty(st, :ngrid) ? Int(getproperty(st, :ngrid)) : ngrid

        if adaptive_window && si > 1 && !isnothing(prev_hw) && isnothing(stages)
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
                println("  sweep $sweep → peakretentions=$(best_t0)  SSE=$(best_sse)  ",
                    "Δrel=$(round(rel_improve, sigdigits=3))")
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

    best_t0, best_sse, best_pack...
end
