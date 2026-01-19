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
        # Fallback to mean if median is not usable.
        wm = mean(w)
    end
    if !isfinite(wm) || wm ≤ 0.0
        # Abort if no valid scale can be determined.
        return nothing
    end
    # Scale weights by their typical magnitude.
    w ./= wm

    # Cap extreme weights for numerical stability.
    w = min.(w, w_cap)
    if !all(isfinite.(w))
        # Abort if weights are still non-finite after capping.
        return nothing
    end
    
    # Return the stabilized weights.
    w
end

"""
    unimodalfit(
        I::AbstractMatrix{<:Real},
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
        spectra_ridge::Real=1e-12
    )

Fit `Kpeaks` chromatographic peaks (components) to GC/MS intensity data `I` where each
m/z channel can be measured at a slightly shifted retention within each scan. The model
for channel `i` and scan `j` is

    I[i,j] ≈ sum_{k=1..Kpeaks} A[i,k] * f_k(R[i,j]) + noise

where `A[i,k] ≥ 0` is the contribution of peak `k` in channel `i` and `f_k(r) ≥ 0` is the
shared peak shape of component `k` over continuous retention `r`. Each `f_k` is
parameterized as a cubic B-spline `f_k(r) = B(r) * c_k`.

Unimodality is enforced using `peakretentions` as the apex switch retention. The 
constraints require `f'_k(r) ≥ 0` for `r ≤ peakretentions[k]`, `f'_k(r) ≤ 0` for 
`r ≥ peakretentions[k]`, and `f_k(r) ≥ 0` on a dense retention grid.

The optimization alternates between updating `A` and updating the spline coefficients
`C`. Given shapes `f_k`, each channel's weights `A[i,:]` are updated by solving a 
nonnegative least squares problem (optionally with a soft prior). Given `A`, the
spline coefficients are updated by solving a constrained quadratic program that
enforces positivity and unimodality. Each component is normalized so `max(f_k)=1`
and the scale is absorbed into `A[:,k]`. This alternation repeats for `iters` steps.

Heteroscedasticity is optional. If `σ` is provided with the same size as `I`, it is used in 
the `A`-update via weights `w = 1 ./ σ.^2`. The shape update remains unweighted to avoid 
bias in the peak height when the noise model is approximate. If `σ` is `nothing`, all 
observations are equally weighted.

`I` is the (`n_mz` × `n_scans`) intensity matrix, `R` is the matching matrix of acquisition 
retentions, and `peakretentions` gives the approximate apex retention for each peak. `knots_n`
sets the number of spline knots, `rgrid_n` sets the grid used for constraints, `iters` 
controls the number of alternating updates, `shape_ridge` stabilizes the spline QP, and 
`spectra_ridge` stabilizes the per-channel `A` updates. `spectra_init`, `spectra_zero_mask`, `spectra_prior`,
and `spectra_prior_weights` control initialization and optional spectral priors.

Returns `A_hat` (estimated spectra), `basis` (B-spline basis), `C_hat` (spline
coefficients), `(rgrid, Fhat_grid)` for plotting shapes, and `Ifit` (reconstructed
signal at observed retentions).

See also [`unimodalfit_apexsearch`](@ref).
"""
function unimodalfit(
    I::AbstractMatrix{<:Real},
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
    spectra_ridge::T10=1e-12
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
        T10<:Real
}
    # Determine the matrix dimensions of the input signal.
    n_mz, n_scans = size(I)
    # Ensure R aligns with I for all observations.
    if size(R) ≠ (n_mz, n_scans)
        throw(DimensionMismatch("R must have size ($(n_mz), $(n_scans))."))
    end

    # Validate sigma if provided.
    if !isnothing(σ)
        # Sigma must match I in size.
        if size(σ) ≠ (n_mz, n_scans)
            throw(DimensionMismatch("σ must have size ($(n_mz), $(n_scans))."))
        end
        # Sigma must be finite.
        if any(.!isfinite.(σ))
            throw(ArgumentError("σ must contain only finite values."))
        end
        # Sigma must be strictly positive.
        if any(σ .≤ 0.0)
            throw(ArgumentError("σ must be strictly positive."))
        end
    end

    # Number of peaks to fit.
    Kpeaks = length(peakretentions)
    # Require at least one peak.
    if Kpeaks < 1
        throw(ArgumentError("peakretentions must contain at least one peak retention."))
    end

    # Validate optional A initialization.
    if  !isnothing(spectra_init)
        if size(spectra_init) ≠ (n_mz, Kpeaks)
            throw(DimensionMismatch("spectra_init must have size ($(n_mz), $(Kpeaks))."))
        end
    end
    # Validate optional A zero-mask.
    if !isnothing(spectra_zero_mask)
        if size(spectra_zero_mask) ≠ (n_mz, Kpeaks)
            throw(DimensionMismatch("spectra_zero_mask must have size ($(n_mz), $(Kpeaks))."))
        end
    end
    # Validate optional A prior.
    if !isnothing(spectra_prior)
        if size(spectra_prior) ≠(n_mz, Kpeaks)
            throw(DimensionMismatch("spectra_prior must have size ($(n_mz), $(Kpeaks))."))
        end
    end

    # spectra_prior and spectra_prior_weights must be provided together.
    if isnothing(spectra_prior) ≠ isnothing(spectra_prior_weights)
        throw(ArgumentError("Provide both spectra_prior and spectra_prior_weights, or neither."))
    end
    # Validate prior weights if provided.
    if  !isnothing(spectra_prior_weights)
        if length(spectra_prior_weights) ≠ Kpeaks
            throw(DimensionMismatch("spectra_prior_weights must have length Kpeaks=$(Kpeaks)."))
        end
        if any(spectra_prior_weights .< 0.0)
            throw(ArgumentError("spectra_prior_weights must be nonnegative."))
        end
    end

    # Flatten retentions in the same order as vec(I) (column-major).
    tvec = vec(R)
    # Compute retention bounds for spline construction.
    tmin, tmax = minimum(tvec), maximum(tvec)
    # Total number of observations across channels and scans.
    Nobs = length(tvec)

    # Spline basis
    # Choose knot locations across the observed retention range.
    knots = range(tmin, tmax; length=knots_n)
    # Construct a cubic B-spline basis.
    basis = BSplineBasis(BSplineOrder(4), knots)  # cubic
    # Collocation matrix at observation retentions.
    Bobs = collocation_matrix(basis, tvec, Derivative(0), Matrix{Float64})  # (Nobs x p)
    # Number of spline coefficients per component.
    p = size(Bobs, 2)

    # Constraint grid
    # Define a dense grid for enforcing constraints.
    rgrid = range(tmin, tmax; length=rgrid_n)
    # Evaluate basis values on the grid.
    Bgrid  = collocation_matrix(basis, rgrid, Derivative(0), Matrix{Float64})
    # Evaluate basis derivatives on the grid.
    Bprime = collocation_matrix(basis, rgrid, Derivative(1), Matrix{Float64})

    # Build block-diagonal constraints for all peaks
    # Allocate constraint matrices per peak.
    Acon_blocks = Vector{SparseMatrixCSC{Float64, Int}}(undef, Kpeaks)
    # Allocate lower bound vectors per peak.
    l_blocks = Vector{Vector{Float64}}(undef, Kpeaks)
    # Allocate upper bound vectors per peak.
    u_blocks = Vector{Vector{Float64}}(undef, Kpeaks)

    # Build constraint blocks for each peak.
    for k in 1:Kpeaks
        # Apex retention for this peak.
        t0 = peakretentions[k]
        # Left side of the apex for non-decreasing constraint.
        left  = rgrid .≤ t0
        # Right side of the apex for non-increasing constraint.
        right = rgrid .≥ t0

        # Assemble positivity and derivative constraints for this peak.
        Acon_k = vcat(
            sparse(Bgrid),                 # f(t) >= 0
            sparse(Bprime[left, :]),       # f'(t) >= 0 on left
            sparse(-Bprime[right, :])      # f'(t) <= 0 on right
        )
        # Store the constraint matrix.
        Acon_blocks[k] = Acon_k
        # Lower bounds for all constraints.
        l_blocks[k] = zeros(size(Acon_k, 1))
        # Upper bounds for all constraints.
        u_blocks[k] = fill(Inf, size(Acon_k, 1))
    end

    function blockdiag(mats::Vector{SparseMatrixCSC{Float64, Int}})
        # Total number of rows in the block-diagonal matrix.
        rows = sum(size(M,1) for M in mats)
        # Total number of columns in the block-diagonal matrix.
        cols = sum(size(M,2) for M in mats)
        # Allocate the output sparse matrix.
        out = spzeros(rows, cols)
        # Row offset for the next block.
        r0 = 1
        # Column offset for the next block.
        c0 = 1
        # Insert each block into the output matrix.
        for M in mats
            # Size of the current block.
            rr, cc = size(M)
            # Place the block at the current offsets.
            out[r0:(r0 + rr - 1), c0:(c0 + cc - 1)] = M
            # Advance row offset.
            r0 += rr
            # Advance column offset.
            c0 += cc
        end
        # Return the assembled block-diagonal matrix.
        out
    end

    # Concatenate all constraint blocks into a single matrix.
    Acon_all = blockdiag(Acon_blocks)
    # Concatenate all lower bounds.
    lcon = vcat(l_blocks...)
    # Concatenate all upper bounds.
    ucon = vcat(u_blocks...)

    # Initialize A
    # Initialize A from user input or a small positive default.
    A_hat = if isnothing(spectra_init)
        fill(1e-3, n_mz, Kpeaks)
    else
        max.(spectra_init, 0.0)
    end
    # Enforce any fixed zeros in A.
    if !isnothing(spectra_zero_mask)
        A_hat[spectra_zero_mask] .= 0.0
    end

    # Initialize C (p x Kpeaks)
    # Initial width for Gaussian seed shapes.
    σ0 = (tmax - tmin) / 12
    # Allocate spline coefficients.
    C = zeros(p, Kpeaks)
    # Cache grid as a concrete vector.
    tg = collect(rgrid)
    # Seed each component with a unimodal bump around its apex.
    for k in 1:Kpeaks
        # Build a Gaussian bump centered at peakretentions[k].
        bump = exp.(-0.5 .* ((tg .- peakretentions[k]) ./ σ0).^2)
        # Normalize the bump to unit height.
        bump ./= maximum(bump) > 0 ? maximum(bump) : 1.0
        # Fit the bump with the spline basis.
        C[:, k] = Bgrid \ bump
        # Normalize the spline shape to unit max on the grid.
        fk = Bgrid * C[:, k]
        fmax = maximum(fk)
        if fmax > 0
            # Rescale coefficients so the grid maximum is 1.
            C[:, k] ./= fmax
        end
    end

    # Vectorize I for the global QP.
    yvec = vec(I)

    # Alternate between A and C updates.
    for _ in 1:iters
        # Evaluate shapes at observations
        # Compute all component values at observation retentions.
        Fobs_vec = Bobs * C  # (Nobs x Kpeaks)

        # Update A via per-channel OSQP NNLS (+ optional soft prior)
        # Iterate over each m/z channel.
        for i in 1:n_mz
            # Allocate the design matrix for this channel.
            Fi = zeros(n_scans, Kpeaks)
            for k in 1:Kpeaks
                @inbounds for j in 1:n_scans
                    # Fill the design matrix entry for (scan j, component k).
                    Fi[j, k] = Fobs_vec[i + (j-1)*n_mz, k]
                end
            end

            # Observations for this channel.
            yi = vec(I[i, :])

            # Optional per-channel weights derived from sigma.
            wi = if isnothing(σ)
                nothing
            else
                sigma2weights(vec(σ[i, :]))
            end

            # Solve for A with or without a soft prior.
            if isnothing(spectra_prior)
                A_hat[i, :] .= nnlspenalized(Fi, yi; ridge=spectra_ridge, w=wi)
            else
                μ_i = vec(spectra_prior[i, :])
                λ_i = vec(spectra_prior_weights)
                A_hat[i, :] .= nnlspenalized(Fi, yi; ridge=spectra_ridge, μ=μ_i, λ=λ_i, w=wi)
            end

            # Ensure numerical validity of the solution.
            if any(.!isfinite.(A_hat[i, :]))
                throw(ArgumentError(
                    "Non-finite A_hat produced in channel $i (check σ / weights)."))
            end
        end

        # Re-apply any hard-zero constraints on A.
        if !isnothing(spectra_zero_mask)
            A_hat[spectra_zero_mask] .= 0.0
        end

        # Update C via one big QP in x = vec(C)
        # NOTE: UNWEIGHTED shape update (ignores σ by design)
        # Allocate the global design matrix for the shape update.
        BW_all = zeros(Nobs, p * Kpeaks)
        for k in 1:Kpeaks
            # Repeat A for this component across scans.
            wA = repeat(view(A_hat, :, k), n_scans)     # length Nobs
            # Apply component weights to the basis matrix.
            BWk = Bobs .* wA                             # (Nobs x p)
            # Insert this component block into the global design matrix.
            BW_all[:, ((k - 1) * p + 1):(k * p)] .= BWk
        end

        # Build the quadratic objective for the global shape update.
        P = BW_all' * BW_all
        # Add ridge for numerical stability.
        P .+= shape_ridge .* JuChrom.I(size(P, 1))
        # Linear term for the QP.
        q = -(BW_all' * yvec)

        # Solve the constrained QP for spline coefficients.
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

        # Unpack the solution into per-peak coefficient vectors.
        for k in 1:Kpeaks
            C[:, k] = view(x, ((k - 1) * p + 1):(k * p))
        end

        # Normalize each component (max on grid = 1)
        # Normalize shapes and push scale into A.
        for k in 1:Kpeaks
            fk_grid = Bgrid * C[:, k]
            fmax = maximum(fk_grid)
            if fmax > 0
                # Normalize shape and adjust A to preserve scale.
                C[:, k] ./= fmax
                # Push scaling into the spectrum for this component.
                A_hat[:, k] .*= fmax
            end
        end
    end

    # Final outputs
    # Evaluate fitted shapes on the constraint grid.
    Fhat_grid = Bgrid * C  # (rgrid_n x Kpeaks)

    # Evaluate fitted shapes at observation retentions.
    Fobs_vec = Bobs * C
    # Accumulate the fitted signal in vectorized form.
    Ifit_vec = zeros(Nobs)
    for k in 1:Kpeaks
        # Repeat A for this component across scans.
        wA = repeat(view(A_hat, :, k), n_scans)
        # Add this component's contribution to the fit.
        Ifit_vec .+= wA .* view(Fobs_vec, :, k)
    end
    # Reshape the vectorized fit back to n_mz × n_scans.
    Ifit = reshape(Ifit_vec, n_mz, n_scans)

    # Return spectra, basis, coefficients, shapes, and fitted signal.
    A_hat, basis, C, (rgrid, Fhat_grid), Ifit
end

"""
    unimodalfit_apexsearch(I, R; peakretentions_guess, kwargs...)

Wrapper around `unimodalfit` that searches apex retentions using coordinate-wise grid
search to minimize the total sum of squared errors (SSE). All keyword arguments accepted by
`unimodalfit` are forwarded through `kwargs` and apply in the same way as in a direct call,
including `σ`, `spectra_init`, `spectra_zero_mask`, `spectra_prior`, `spectra_prior_weights`,
`knots_n`, `rgrid_n`, `iters`, `shape_ridge`, and `spectra_ridge`.

The search parameters are specific to this wrapper. `peakretentions_guess` provides the
initial apex retentions. `half_width` and `ngrid` define the one-dimensional grid around
each apex. `min_sep` and `max_sep` optionally constrain adjacent apex retentions. `strategy`
selects the search plan: `:single` runs one grid search at the given half-width and grid
size, while `:iterative` runs a coarse-to-fine sequence and optionally accepts explicit
`stages`. When provided, `stages` should be an iterable of stage objects (e.g. a
vector/tuple of `NamedTuple`s) with `half_width` and/or `ngrid` fields; missing fields fall
back to the top-level defaults. `adaptive_window`, `shrink`, and `min_half_width` control
automatic narrowing of the search window between stages. `coord_sweeps` sets the number of
coordinate passes per stage. `enforce_sorted` keeps apex retentions ordered. `max_iter` and
`tol_sse` control termination, and `verbose` controls logging.

If `σ` is provided through `kwargs`, it is forwarded to `unimodalfit` and applied only in
the per-channel `A` update. The peak-shape update remains unweighted.

Returns the best apex retention vector, the best SSE, and the same outputs as `unimodalfit`
for that retention vector: `A_hat`, `basis`, `C_hat`, `(rgrid, Fhat_grid)`, and `Ifit`.

See also [`unimodalfit`](@ref).
"""
function unimodalfit_apexsearch(
    I::AbstractMatrix{<:Real}, 
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

    # Count the number of peaks implied by the initial guesses.
    Kpeaks = length(peakretentions_guess)
    # Require at least one peak.
    if Kpeaks < 1
        throw(ArgumentError("peakretentions_guess must contain at least one peak retention."))
    end
    # Require a positive search half-width.
    if half_width ≤ 0
        throw(ArgumentError("half_width must be positive."))
    end
    # Require a minimum grid resolution.
    if ngrid < 3
        throw(ArgumentError("ngrid must be at least 3."))
    end
    # Require at least one coordinate sweep per stage.
    if coord_sweeps < 1
        throw(ArgumentError("coord_sweeps must be at least 1."))
    end
    # Require at least one total iteration.
    if max_iter < 1
        throw(ArgumentError("max_iter must be at least 1."))
    end
    # Require a positive convergence tolerance.
    if tol_sse ≤ 0
        throw(ArgumentError("tol_sse must be positive."))
    end

    # Build the stage plan for the grid search.
    # Select the stage plan based on the requested strategy.
    plan = if strategy == :single
        # Single-stage plan with the provided parameters.
        [(half_width=half_width, ngrid=ngrid)]
    elseif strategy == :iterative
        if isnothing(stages)
            # Default coarse-to-fine schedule.
            # Start with the user-provided half-width.
            hw1 = half_width
            # Use a medium refinement based on the initial width.
            hw2 = max(min_half_width, 0.35 * hw1)
            # Use a fine refinement based on the initial width.
            hw3 = max(min_half_width, 0.12 * hw1)
            # Assemble the three-stage plan.
            [(half_width=hw1, ngrid=max(9,  min(ngrid, 13))),
             (half_width=hw2, ngrid=max(13, min(ngrid, 19))),
             (half_width=hw3, ngrid=max(17, ngrid))]
        else
            # Use caller-provided stage definitions.
            collect(stages)
        end
    else
        # Reject unsupported strategies.
        error("strategy must be :single or :iterative")
    end

    # Initialize and optionally sort t0.
    # Convert the initial guesses to Float64.
    t0 = Float64.(peakretentions_guess)
    # Enforce sorted order if requested.
    if enforce_sorted
        # Sort the apex retentions in-place.
        sort!(t0)
    end

    function feasible(t0s::Vector{Float64})
        # Check min/max separation constraints.
        if enforce_sorted
            # Scan adjacent apex retentions for violations.
            for k in 1:(Kpeaks-1)
                # Compute separation between neighbors.
                dt = t0s[k + 1] - t0s[k]
                # Enforce minimum separation.
                if dt < min_sep
                    return false
                end
                # Enforce maximum separation if requested.
                if !isnothing(max_sep) && dt > max_sep
                    return false
                end
            end
        end
        # Feasible if no constraint was violated.
        true
    end

    function eval_fit(t0s::Vector{Float64})
        # Fit shapes and spectra at the proposed t0.
        A_hat, basis_hat, C_hat, (rgrid, Fhat_grid), I_fit =
            unimodalfit(I, R, t0s; kwargs...)
        # Compute the total sum of squared errors.
        sse = sum((I .- I_fit).^2)
        # Return the SSE and the fitted outputs.
        sse, (A_hat, basis_hat, C_hat, (rgrid, Fhat_grid), I_fit)
    end

    # Ensure the initial guess is feasible.
    if !feasible(t0)
        throw(ArgumentError("peakretentions_guess violates min_sep/max_sep constraints."))
    end

    # Evaluate the initial solution.
    best_sse, best_pack = eval_fit(copy(t0))
    # Keep the best t0 vector.
    best_t0 = copy(t0)

    # Report the initial state if requested.
    if verbose
        # Print initial apex retentions and SSE.
        println("Initial peakretentions = ", best_t0, "  SSE = ", best_sse)
    end

    # Track adaptive window between stages.
    # Track the previous half-width for adaptive windowing.
    prev_hw = nothing
    # Track the total number of sweeps performed.
    iter = 0

    # Iterate over stages in the search plan.
    for (si, st) in enumerate(plan)
        # Determine search window and grid size for this stage.
        hw = hasproperty(st, :half_width) ? Float64(getproperty(st, :half_width)) : half_width
        # Determine grid resolution for this stage.
        ng = hasproperty(st, :ngrid) ? Int(getproperty(st, :ngrid)) : ngrid

        # Optionally apply adaptive shrinking for default stages.
        if adaptive_window && si > 1 && !isnothing(prev_hw) && isnothing(stages)
            # Shrink the half-width while respecting the minimum.
            hw = max(min_half_width, shrink * prev_hw)
        end
        # Store the current half-width for the next stage.
        prev_hw = hw

        # Report stage settings if requested.
        if verbose
            # Print stage configuration.
            println("\nStage $si  (half_width=$(round(hw, digits=4)), ngrid=$ng)")
        end

        # Coordinate descent sweeps over apex retentions.
        for sweep in 1:coord_sweeps
            # Increment total sweep counter.
            iter += 1
            # Enforce the global iteration cap.
            if iter > max_iter
                if verbose
                    # Report early termination.
                    println("Stopping: reached max_iter = $max_iter")
                end
                # Return the best solution found so far.
                return best_t0, best_sse, best_pack...
            end

            # Track whether this sweep improved the solution.
            improved = false
            # Save the SSE at the start of the sweep.
            prev_sse = best_sse

            # Optimize each t0 coordinate on its local grid.
            for k in 1:Kpeaks
                # Center the search grid at the current apex.
                center = best_t0[k]
                # Build a one-dimensional grid around the center.
                grid = range(center - hw, center + hw; length=ng)

                # Track the best result for this coordinate.
                local_best_sse = best_sse
                # Initialize the local best t0.
                local_best_t0 = copy(best_t0)
                # Initialize the local best fit package.
                local_best_pack = best_pack

                # Try candidate apex positions on the grid.
                for cand in grid
                    # Propose a new t0 vector.
                    t0_try = copy(best_t0)
                    # Update the k-th apex position.
                    t0_try[k] = Float64(cand)

                    # Enforce ordering if requested.
                    if enforce_sorted
                        # Sort the proposed t0 in-place.
                        sort!(t0_try)
                    end
                    # Skip infeasible candidates.
                    if !feasible(t0_try)
                        continue
                    end

                    # Fit the model at this candidate.
                    sse, pack = eval_fit(t0_try)

                    # Update the local best if SSE improves.
                    if sse < local_best_sse
                        # Save the improved SSE.
                        local_best_sse = sse
                        # Save the improved t0 vector.
                        local_best_t0 = t0_try
                        # Save the improved fit package.
                        local_best_pack = pack
                    end
                end

                # Update global best if this coordinate improved.
                if local_best_sse < best_sse
                    # Record the improved SSE.
                    best_sse = local_best_sse
                    # Record the improved t0 vector.
                    best_t0 = local_best_t0
                    # Record the improved fit package.
                    best_pack = local_best_pack
                    # Mark that an improvement occurred.
                    improved = true
                end
            end

            # Track relative improvement for termination checks.
            rel_improve = abs(prev_sse - best_sse) / max(prev_sse, 1e-12)

            # Report progress for this sweep if requested.
            if verbose
                # Print current apex retentions, SSE, and relative change.
                println("  sweep $sweep → peakretentions=$(best_t0)  SSE=$(best_sse)  ",
                    "Δrel=$(round(rel_improve, sigdigits=3))")
            end

            # Stop if improvement falls below tolerance.
            if rel_improve < tol_sse
                if verbose
                    # Report convergence based on tolerance.
                    println("Stopping: SSE improvement < tol_sse = $tol_sse")
                end
                # Return the best solution found so far.
                return best_t0, best_sse, best_pack...
            end

            # Stop this stage if no coordinate improved.
            if !improved
                # Exit the current stage early.
                break
            end
        end
    end

    # Return the best solution after all stages.
    best_t0, best_sse, best_pack...
end
