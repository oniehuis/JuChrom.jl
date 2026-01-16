@inline function _strip_mz_value(val::Number, mzunit::Union{Nothing, Unitful.Units})
    if isnothing(mzunit)
        isunitful(val) && throw(ArgumentError("m/z grid is unitless but the value $val carries units."))
        return val
    elseif isunitful(val)
        return Unitful.ustrip(Unitful.uconvert(mzunit, val))
    else
        return val
    end
end

@inline _strip_mz_values(vals, mzunit::Union{Nothing, Unitful.Units}) =
    map(v -> _strip_mz_value(v, mzunit), vals)

# ── fitquadvarmodel (series → QuadVarFit) ───────────────────────────────────

"""
    fitquadvarmodel(series_batches::AbstractVector{<:AbstractVector};
                    mzsel=nothing, show_progress::Bool=true, progress_every::Int=1, kwargs...)

Fit the signal/offset/gain + quadratic-variance model across a selection of m/z.
Returns a **self-contained** `QuadVarFit` that includes the observed matrices
restricted to the selected m/z.

# Arguments
- `series_batches`: vector of batches; each batch is a vector of replicate series.
  Within a batch, replicates must share the same m/z grid and the same scan count.

# Keyword arguments
- `mzsel`: m/z selection:
    - `nothing` → all m/z in the reference grid,
    - indices (range or `Vector{Int}`) into `mz_ref`,
    - `Vector{<:Number}` of exact m/z values (unitless numbers interpreted in `mzunit`;
      Unitful quantities are converted to `mzunit`) found in `mz_ref`.
- `show_progress` (default `true`), `progress_every` (default `1`).
- `kwargs...` are forwarded to `fitquadvarmodel(Ms::Vector{<:AbstractMatrix})`.

# Returns
A `QuadVarFit` struct containing per-m/z fits and the observed data for the selected m/z.
"""
function fitquadvarmodel(
    series_batches::AbstractVector{<:AbstractVector};
    mzsel=nothing, show_progress::Bool=true, progress_every::Int=1, kwargs...
)

    batchcount = length(series_batches)
    batchcount ≥ 1 || throw(ArgumentError("No batches provided."))
    progress_every ≥ 1 || throw(ArgumentError("progress_every must be ≥ 1"))

    # Track and validate the intensity unit across all inputs (fallback to nothing if absent)
    first_ms = mscanmatrix(series_batches[1][1])
    intensity_unit = applicable(intensityunit, first_ms) ? intensityunit(first_ms) : nothing
    mz_ref_raw = uniquemzvalues(series_batches[1][1])
    mzunit = (!isempty(mz_ref_raw) && isunitful(first(mz_ref_raw))) ? Unitful.unit(first(mz_ref_raw)) : nothing

    # Convert series to n_scans × n_mz matrices (strict checks)
    matrices_batches = Vector{Vector{Matrix{Float64}}}(undef, batchcount)
    mz_ref = collect(_strip_mz_values(mz_ref_raw, mzunit))
    n_mz = length(mz_ref)

    for b in 1:batchcount
        reps = series_batches[b]
        repcount = length(reps)
        repcount ≥ 2 || throw(ArgumentError("Batch $b needs ≥ 2 replicates (got $repcount)."))

        mats = Vector{Matrix{Float64}}(undef, repcount)
        n_scans_ref = nothing

        for (r, obj) in enumerate(reps)
            mz_r = _strip_mz_values(uniquemzvalues(obj), mzunit)
            (length(mz_r) == n_mz && mz_r == mz_ref) ||
                throw(ArgumentError("m/z grid mismatch at batch $b, replicate $r."))

            ms_obj = mscanmatrix(obj)
            M = rawintensities(ms_obj)  # n_scans × n_mz
            size(M, 2) == n_mz ||
                throw(ArgumentError("Column count mismatch at batch $b, rep $r."))

            iu = applicable(intensityunit, ms_obj) ? intensityunit(ms_obj) : nothing
            iu == intensity_unit || throw(ArgumentError("Intensity unit mismatch at batch $b, rep $r."))

            n_scans = size(M, 1)
            if n_scans_ref === nothing
                n_scans_ref = n_scans
            elseif n_scans != n_scans_ref
                throw(ArgumentError("Batch $b has inconsistent scan counts: rep $r has $n_scans, expected $n_scans_ref."))
            end

            mats[r] = Float64.(M)
        end
        matrices_batches[b] = mats
    end

    # Determine selected m/z indices
    mz_idx = Int[]
    sel_mz_vals = nothing

    if mzsel === nothing
        mz_idx = collect(1:n_mz)
    elseif mzsel isa AbstractRange{<:Integer} || mzsel isa AbstractVector{<:Integer}
        mz_idx = collect(mzsel)
        all(1 .≤ mz_idx .≤ n_mz) ||
            throw(ArgumentError("m/z indices out of bounds; valid range is 1:$n_mz."))
    elseif mzsel isa AbstractVector{<:Number}
        vals_raw = collect(mzsel)
        vals = collect(_strip_mz_values(vals_raw, mzunit))
        idxmap = Dict(v => i for (i, v) in enumerate(mz_ref))   # value -> index
        mz_idx = [get(idxmap, v, 0) for v in vals]
        all(mz_idx .> 0) || throw(ArgumentError("At least one requested m/z value not found in grid."))
        sel_mz_vals = vals
    else
        throw(ArgumentError("Unsupported `mzsel` type. Use indices (range/vector) or m/z values vector."))
    end

    n_mz_select = length(mz_idx)
    n_mz_select ≥ 1 || throw(ArgumentError("Empty m/z selection."))

    # Reusable per-batch buffers Ms[b] :: n_scans_b × n_reps_b
    Ms = Vector{Matrix{Float64}}(undef, batchcount)
    for b in 1:batchcount
        n_reps_b  = length(matrices_batches[b])
        n_scans_b = size(matrices_batches[b][1], 1)
        Ms[b] = Matrix{Float64}(undef, n_scans_b, n_reps_b)
    end

    # Fit first selected m/z to establish result type
    i1 = mz_idx[1]
    for b in 1:batchcount
        n_reps_b = size(Ms[b], 2)
        @inbounds for r in 1:n_reps_b
            @views Ms[b][:, r] .= matrices_batches[b][r][:, i1]
        end
    end
    first_result = fitquadvarmodel(Ms; kwargs...)
    results = Vector{typeof(first_result)}(undef, n_mz_select)
    results[1] = first_result

    if show_progress
        r  = results[1]; p = r.params; st = r.stats
        mz_val = isnothing(sel_mz_vals) ? mz_ref[i1] : sel_mz_vals[1]
        mz_fmt = isnothing(mzunit) ? @sprintf("%6.1f", mz_val) :
            @sprintf("%6.1f %s", mz_val, string(mzunit))
        @info @sprintf("m/z %s  z_rms = %5.3f  cov95 = %5.3f  σ₀² = %12.3f  ϕ = %7.3f  κ = %8.4f  acf = %6.3f  [%4d/%4d]",
                       mz_fmt, st.z_rms, st.cov95, p.σ₀², p.ϕ, p.κ, r.acf, 1, n_mz_select)
    end

    # Remaining m/z
    for j in 2:n_mz_select
        mi = mz_idx[j]
        for b in 1:batchcount
            n_reps_b = size(Ms[b], 2)
            @inbounds for r in 1:n_reps_b
                @views Ms[b][:, r] .= matrices_batches[b][r][:, mi]
            end
        end
        results[j] = fitquadvarmodel(Ms; kwargs...)

        if show_progress && (j % progress_every == 0 || j == n_mz_select)
            r  = results[j]; p = r.params; st = r.stats
            mz_val = isnothing(sel_mz_vals) ? mz_ref[mi] : sel_mz_vals[j]
            mz_fmt = isnothing(mzunit) ? @sprintf("%6.1f", mz_val) :
                @sprintf("%6.1f %s", mz_val, string(mzunit))
            @info @sprintf("m/z %s  z_rms = %5.3f  cov95 = %5.3f  σ₀² = %12.3f  ϕ = %7.3f  κ = %8.4f  acf = %6.3f  [%4d/%4d]",
                           mz_fmt, st.z_rms, st.cov95, p.σ₀², p.ϕ, p.κ, r.acf, j, n_mz_select)
        end
    end

    n_scans_per_batch = [size(matrices_batches[b][1], 1) for b in 1:batchcount]
    n_reps_per_batch  = [length(matrices_batches[b])     for b in 1:batchcount]

    # Pack results into QuadVarFit fields
    params   = Vector{QuadVarParams}(undef, n_mz_select)
    signal   = Vector{Vector{Vector{Float64}}}(undef, n_mz_select)
    offsets  = Vector{Vector{Matrix{Float64}}}(undef, n_mz_select)
    gains    = Vector{Vector{Vector{Float64}}}(undef, n_mz_select)
    scale_c  = Vector{Float64}(undef, n_mz_select)
    acf      = Vector{Float64}(undef, n_mz_select)
    acf_lag  = Vector{Int}(undef, n_mz_select)
    n_pairs  = Vector{Int}(undef, n_mz_select)
    qc_z_rms = Vector{Float64}(undef, n_mz_select)
    qc_cov68 = Vector{Float64}(undef, n_mz_select)
    qc_cov95 = Vector{Float64}(undef, n_mz_select)
    qc_s_min = Vector{Float64}(undef, n_mz_select)
    qc_s_max = Vector{Float64}(undef, n_mz_select)
    qc_nz    = Vector{Int}(undef, n_mz_select)

    @inbounds for j in 1:n_mz_select
        r = results[j]
        params[j]   = r.params
        signal[j]   = r.s
        offsets[j]  = r.o
        gains[j]    = r.g
        scale_c[j]  = r.c
        acf[j]      = r.acf
        acf_lag[j]  = r.acf_lag
        n_pairs[j]  = r.n_pairs

        st = r.stats
        qc_z_rms[j] = st.z_rms
        qc_cov68[j] = st.cov68
        qc_cov95[j] = st.cov95
        qc_s_min[j] = st.s_min
        qc_s_max[j] = st.s_max
        qc_nz[j]    = st.nz
    end

    # Store observed data restricted to selected m/z (self-contained)
    observed = Vector{Vector{Matrix{Float64}}}(undef, batchcount)
    for b in 1:batchcount
        n_reps_b = n_reps_per_batch[b]
        obs_b = Vector{Matrix{Float64}}(undef, n_reps_b)
        for r in 1:n_reps_b
            obs_b[r] = copy(matrices_batches[b][r][:, mz_idx])  # n_scans_b × n_mz_select
        end
        observed[b] = obs_b
    end

    mzvalues = isnothing(sel_mz_vals) ? mz_ref[mz_idx] : sel_mz_vals

    QuadVarFit(
        mzunit,
        mz_ref,
        mz_idx,
        mzvalues,
        batchcount,
        n_reps_per_batch,
        n_scans_per_batch,
        params,
        intensity_unit,
        signal,
        offsets,
        gains,
        scale_c,
        acf,
        acf_lag,
        n_pairs,
        qc_z_rms,
        qc_cov68,
        qc_cov95,
        qc_s_min,
        qc_s_max,
        qc_nz,
        observed,
    )
end

"""
    fitquadvarmodel(Ms::Vector{<:AbstractMatrix{<:Real}}; kwargs...)

End-to-end fitting of the **signal/offset/gain + quadratic variance** model via
alternating updates with smoothing, ridge regularization, pooled variance
re-estimation, optional global scale calibration, and diagnostics.

# Arguments
- `Ms`: list of batches; each `Ms[g]` is `N_g × R_g` (N scans, R replicates).

# Keyword arguments
**Signal / offset / gain**
- `λ_signal`: smoothing strength on the shared signal.
- `λ_offset`: smoothing (ridge) on offsets via `D' D`.
- `gain_ridge`: ridge shrinkage of gains toward `1`.
- `damping`: step size for the signal update in `[0, 1]`.

**Variance re-fit**
- `quantlevel`, `bins`, `var_stat`, `σ₀²_min`, `ϕ_min`, `κ_min`,
  `ϕ_ridge`, `κ_ridge`, `z_max`, `c_eps`, `v_clip_q`, `z_low_q`, `feas_tol`.

**Floors and jitter**
- `qc_varfloor`: variance floor for QC/diagnostics.
- `τ_floor`: fraction to set the adaptive floor from per-batch medians.
- `varpred_floor`: hard floor used inside `varpred` when forming weights.
- `chol_jitter`: relative diagonal jitter in the offset solver.

**Iteration / stopping**
- `maxiter`, `var_update_every`, `obj_reltol`.

**Calibration**
- `calibrate_scale`, `cal_trim_lo`, `cal_trim_hi`, `cal_clip_lo`, `cal_clip_hi`.

**Diagnostics**
- `ac_lag::Int=1`: residual autocorrelation lag.
- `acf_restrict::Symbol=:flat`: observation mask used for ACF.
    - `:none`     → use all rows
    - `:low`      → keep rows with low normalized signal `z = s/c′`
    - `:flat`     → keep rows with small local slope (|Δs| below a quantile)
    - `:lowflat`  → apply both low and flat masks
- `acf_z_quant::Float64=quantlevel`: quantile for pooled scale `c′` (used only for `:low`/`:lowflat`).
- `acf_z_max_keep::Float64=0.30`: keep rows with `s/c′ ≤ acf_z_max_keep` when low mask is active.
- `acf_flat_q::Float64=0.20`: keep rows where both adjacent |Δs| are ≤ this quantile.
- `acf_trimfrac::Float64=0.0`: fraction of extreme residual pairs to trim when computing ACF.

# Returns
A named tuple with fields:
- `params::QuadVarParams`
- `s::Vector{Vector{Float64}}`
- `o::Vector{Matrix{Float64}}`
- `g::Vector{Vector{Float64}}`
- `c::Float64`
- `stats::NamedTuple` with keys `z_rms`, `cov68`, `cov95`, `s_min`, `s_max`
- `acf::Float64`
- `acf_lag::Int`
- `n_pairs::Int`

# Notes
- Weights use `1 ./ max(varpred(s, p; varfloor=varpred_floor), varfloor_eff)`
  with `varfloor_eff = max(qc_varfloor, τ_floor * median(per-batch medians))`.
- At least one variance update is guaranteed even if the main loop exits early.
"""
function fitquadvarmodel(
    Ms::Vector{<:AbstractMatrix{<:Real}};
    # signal/offset/gain
    λ_signal::Float64=0.1,
    λ_offset::Float64=100.0,
    gain_ridge::Float64=0.05,
    damping::Float64=0.5,
    qc_varfloor::Float64=1e-12,
    maxiter::Int=80,
    var_update_every::Int=5,
    obj_reltol::Float64=1e-5,
    # variance update knobs
    quantlevel::Float64=0.90,
    bins::Int=20,
    var_stat::Symbol=:var,
    σ₀²_min::Float64=0.0,
    ϕ_min::Float64=1.0,
    κ_min::Float64=0.0,
    ϕ_ridge::Float64=2e-4,
    κ_ridge::Float64=2e-4,
    z_max::Float64=5.0,
    # guards / robust
    c_eps::Float64=10.0,
    v_clip_q::Float64=0.97,
    # adaptive variance floor
    τ_floor::Float64=1e-3,
    # global variance-scale calibration
    calibrate_scale::Bool=true,
    cal_trim_lo::Float64=0.20,
    cal_trim_hi::Float64=0.80,
    cal_clip_lo::Float64=0.7,
    cal_clip_hi::Float64=1.6,
    # residual autocorrelation (single integer lag)
    ac_lag::Int=1,
    feas_tol::Float64=1e-12,
    acf_restrict::Symbol=:flat,         # :none | :low | :flat | :lowflat
    acf_z_quant::Float64=quantlevel,
    acf_z_max_keep::Float64=0.10,            # keep rows with s/c ≤ this
    acf_flat_q::Float64=0.10,           # keep rows with |Δs| ≤ this quantile
    acf_trimfrac::Float64=0.0,
    # numerical stability
    chol_jitter::Float64=1e-8,
    varpred_floor::Float64=0.0,
    )

    G = length(Ms)
    Ts = map(M -> size(M,1), Ms)
    Rs = map(M -> size(M,2), Ms)
    all(Ts .≥ 3) || error("each batch needs ≥3 scans")
    all(Rs .≥ 2) || error("each batch needs ≥2 replicates")

    # per-batch state
    Y  = [Float64.(M) for M in Ms]
    s  = [vec(mean(Y[g], dims=2)) for g in 1:G]   # init shared signal per batch
    o  = [zeros(Float64, Ts[g], Rs[g]) for g in 1:G]
    gn = [ones(Float64, Rs[g]) for g in 1:G]
    DT = [makeDTD(Ts[g]) for g in 1:G]
    D  = [makediff2(Ts[g]) for g in 1:G]

    # initial variance params and scale
    p = QuadVarParams(1.0, 0.0, 0.0)
    c = max(quantile(vcat(s...), quantlevel), eps(Float64))  # for logging only

    # initialize adaptive variance floor
    v_meds = Float64[]
    for g in 1:G
        push!(v_meds, median(varpred(s[g], p; varfloor=varpred_floor)))
    end
    varfloor_eff = max(qc_varfloor, τ_floor * median(v_meds))

    Jprev = Inf
    did_var_update = false

    for it in 1:maxiter
        # weights from current variance (physical) with adaptive floor
        w = [1.0 ./ max.(varpred(s[g], p; varfloor=varpred_floor), varfloor_eff) for g in 1:G]

        # offsets
        F = [cholridge(w[g], DT[g], λ_offset; jitter=chol_jitter) for g in 1:G]
        for g in 1:G
            fitoffsets!(o[g], Y[g], gn[g], s[g], F[g], w[g])
        end

        # gains (with optional ridge to 1) + rescale s
        for g in 1:G
            fitgains!(gn[g], Y[g], o[g], s[g], w[g]; ridge_to_1=gain_ridge)
        end

        # signal update + soft smoothing via (I + λ D'D)⁻¹
        for g in 1:G
            updatesignal!(s[g], Y[g], o[g], gn[g], w[g]; damping=damping)
            A = I + λ_signal * (D[g]' * D[g])
            s[g] .= A \ s[g]
            @. s[g] = max(s[g], 0.0)
        end

        # pooled variance update
        if (it == 1) || (var_update_every ≤ 1) || (it % var_update_every == 0)
            Rm = [residuals(Y[g], o[g], s[g], gn[g]) for g in 1:G]
            p, c = fitquadvar(
                s, Rm;
                quantlevel=quantlevel,
                bins=bins,
                var_stat=var_stat,
                σ₀²_min=σ₀²_min,
                ϕ_min=ϕ_min,
                κ_min=κ_min,
                ϕ_ridge=ϕ_ridge,
                κ_ridge=κ_ridge,
                z_max=z_max,
                c_eps=c_eps,
                v_clip_q=v_clip_q,
                feas_tol=feas_tol
            )

            if calibrate_scale
                p = rescalevarparams(p, s, Rm; 
                                     varfloor=varfloor_eff,
                                     trim_lo=cal_trim_lo, trim_hi=cal_trim_hi,
                                     clip_lo=cal_clip_lo, clip_hi=cal_clip_hi)
            end

            # hard bounds after any modification of p
            p = clampvarparams(p; σ₀²_min=σ₀²_min, ϕ_min=ϕ_min, κ_min=κ_min)

            did_var_update = true

            # refresh adaptive floor from current predicted variances
            empty!(v_meds)
            for g in 1:G
                push!(v_meds, median(varpred(s[g], p; varfloor=varpred_floor)))
            end
            varfloor_eff = max(qc_varfloor, τ_floor * median(v_meds))
        end

        # monitor objective
        J = 0.0
        for g in 1:G
            v  = max.(varpred(s[g], p; varfloor=varpred_floor), varfloor_eff)
            wᵢ = 1.0 ./ v
            r  = residuals(Y[g], o[g], s[g], gn[g])
            r2 = vec(sum(abs2, r; dims=2))
            J += dot(wᵢ, r2) + λ_signal * sum(abs2, D[g]*s[g])
        end
        if isfinite(Jprev) && abs(Jprev - J) / max(Jprev, 1.0) < obj_reltol
            break
        end
        Jprev = J
    end

    # final safety: guarantee at least one variance update
    if !did_var_update
        Rm = [residuals(Y[g], o[g], s[g], gn[g]) for g in 1:G]
        p, c = fitquadvar(
            s, Rm;
            quantlevel=quantlevel,
            bins=bins,
            var_stat=var_stat,
            σ₀²_min=σ₀²_min,
            ϕ_min=ϕ_min,
            κ_min=κ_min,
            ϕ_ridge=ϕ_ridge,
            κ_ridge=κ_ridge,
            z_max=z_max,
            c_eps=c_eps,
            v_clip_q=v_clip_q,
            feas_tol=feas_tol
        )

        if calibrate_scale
            p = rescalevarparams(p, s, Rm;
                                 varfloor=varfloor_eff,
                                 trim_lo=cal_trim_lo, trim_hi=cal_trim_hi,
                                 clip_lo=cal_clip_lo, clip_hi=cal_clip_hi)
        end

        # hard bounds again
        p = clampvarparams(p; σ₀²_min=σ₀²_min, ϕ_min=ϕ_min, κ_min=κ_min)

        empty!(v_meds)
        for g in 1:G
            push!(v_meds, median(varpred(s[g], p; varfloor=varpred_floor)))
        end
        varfloor_eff = max(qc_varfloor, τ_floor * median(v_meds))
    end

    # Pooled fit statistics
    stats = pooledfitstats(Ms, o, gn, s, p; qc_varfloor=qc_varfloor)

    # Residual autocorrelation
    acfinfo = pooledresidacf(Y, o, gn, s, p;
        lag=ac_lag,
        qc_varfloor=qc_varfloor, 
        restrict=acf_restrict,
        z_quant=acf_z_quant,
        z_max_keep=acf_z_max_keep,
        flat_q=acf_flat_q,
        trimfrac=acf_trimfrac)

    return (params=p, s=s, o=o, g=gn, c=c,
            stats=stats,
            acf=acfinfo.acf, acf_lag=acfinfo.lag, n_pairs=acfinfo.n_pairs)
end

# """
#     boundedlsβ₂(A::AbstractMatrix{<:Real}, b::AbstractVector{<:Real},
#                 lb::NTuple{2,<:Real}; feas_tol::Real=1e-12)

# Solve a 2-parameter **bounded least squares**:

#     minimize ‖A*β − b‖²  subject to  β ≥ lb

# by enumerating the 4 active-set patterns (free/fixed) and picking the best
# feasible solution.

# # Arguments
# - `A`: design matrix with 2 columns.
# - `b`: target vector (length matches `size(A,1)`).
# - `lb`: elementwise lower bounds `(β₁_min, β₂_min)`.

# # Keyword arguments
# - `feas_tol=1e-12`: feasibility slack; accept `β` if `β ≥ lb − feas_tol`.

# # Returns
# A tuple `(β1, β2)::Tuple{Float64,Float64}`.

# # Notes
# - Falls back to the unconstrained `A \\ b` if no enumerated pattern is feasible
#   (then clamps to the bounds).
# """
function boundedlsβ₂(A::AbstractMatrix{<:Real}, b::AbstractVector{<:Real},
                     lb::NTuple{2,<:Real}; feas_tol::Real=1e-12)

    lbv = lb isa Tuple ? Float64[lb...] : Float64.(lb)

    @assert size(A,2) == 2
    @assert length(lbv) == 2

    bestβ = nothing; bestJ = Inf
    for mask in 0:3
        fixed = Int[]; free = Int[]
        for j in 1:2
            (((mask >> (j-1)) & 1) == 1 ? push!(fixed,j) : push!(free,j))
        end
        β = zeros(2)
        if isempty(free)
            β .= lbv
        else
            r = b - A[:,fixed] * lbv[fixed]
            βfree = A[:,free] \ r
            β[fixed] = lbv[fixed]
            β[free] = βfree
        end
        if all(β .≥ lbv .- feas_tol)
            J = sum(abs2, A * β .- b)
            if J < bestJ
                bestJ = J; bestβ = copy(β)
            end
        end
    end

    if isnothing(bestβ)
        β = collect(A \ b)
        for j in 1:2
            β[j] = max(β[j], lbv[j])
        end
        return (β[1], β[2])
    else
        return (bestβ[1], bestβ[2])
    end
end

# """
#     cholridge(w::AbstractVector{<:Real}, DT_D, λ::Real; jitter::Real=1e-8)

# Form a **ridge-regularized SPD system**  

#     A = Diagonal(w .+ j) + λ * DT_D

# and return its Cholesky factorization.

# # Arguments
# - `w`: nonnegative weights (length `n`).
# - `DT_D`: smoothing matrix `D' * D` (size `n × n`).
# - `λ`: ridge/smoothing strength.

# # Keyword arguments
# - `jitter=1e-8`: small relative diagonal shift multiplier, prevents
#   singularity if `w` has zeros.

# # Returns
# A `Cholesky` factorization (`LinearAlgebra.Cholesky` for dense,  
# `SparseArrays.CHOLMOD.Factor` for sparse input).
# """
@inline function cholridge(w::AbstractVector{<:Real}, DT_D, λ::Real; jitter::Real=1e-8)
    # Compute diagonal jitter: scale by mean(w) to keep magnitude consistent,
    # fall back to machine epsilon if w is all zeros
    j = jitter * max(mean(w), eps(Float64))

    # Construct SPD system:
    # - diagonal part encodes observation weights (+ jitter)
    # - λ * DT_D adds smoothness/regularization

    A = Diagonal(w .+ j) + λ * DT_D

    # Return Cholesky factorization (lower-triangular).
    # `check=false` skips SPD validation (assumed guaranteed by construction).
    cholesky(Symmetric(A, :L); check=false)
end

# """
#     clampvarparams(p::QuadVarParams; σ₀²_min::Float64=0.0, ϕ_min::Float64=0.0, κ_min::Float64=0.0)

# Clamp variance parameters to elementwise lower bounds.

# # Keyword arguments
# - `σ₀²_min`, `ϕ_min`, `κ_min`: lower bounds for each parameter.

# # Returns
# A `QuadVarParams` with each field replaced by `max(field, bound)`.
# """
@inline function clampvarparams(p::QuadVarParams; σ₀²_min::Number=0, ϕ_min::Number=0,
    κ_min::Number=0)
    σ_floor = coerce_like_reference(σ₀²_min, p.σ₀², "σ₀²_min")
    ϕ_floor = coerce_like_reference(ϕ_min, p.ϕ, "ϕ_min")
    κ_floor = coerce_like_reference(κ_min, p.κ, "κ_min")
    QuadVarParams(max(p.σ₀², σ_floor), max(p.ϕ, ϕ_floor), max(p.κ, κ_floor))
end

# """
#     fitgains!(g::Vector{Float64}, Y::Matrix{Float64}, o::Matrix{Float64},
#               s::Vector{Float64}, w::Vector{Float64}; ridge_to_1::Float64=0.0)

# Weighted least-squares estimate of per-replicate gains, with optional
# **ridge shrinkage toward 1** and subsequent rescaling to keep `mean(g) = 1`
# (by compensating `s`).

# # Arguments
# - `g`: gains (length R), updated in place.
# - `Y`: data matrix (N × R).
# - `o`: offsets (N × R).
# - `s`: shared signal (length N).
# - `w`: weights (length N).

# # Keyword arguments
# - `ridge_to_1=0.0`: shrinkage strength toward 1 (`0`=none, `1`=force all gains to 1).

# # Returns
# `nothing` (mutates `g` and rescales `s` to preserve overall level).

# # Notes
# - If the normalizer `∑ w .* s.^2` is nonpositive or non-finite, falls back to `g .= 1`.
# """
function fitgains!(g::Vector{Float64}, Y::Matrix{Float64}, o::Matrix{Float64},
                    s::Vector{Float64}, w::Vector{Float64}; ridge_to_1::Float64=0.0)
    den = sum(@. w * s * s)
    if !(isfinite(den) && den > 0)
        fill!(g, 1.0); return
    end
    @inbounds for r in eachindex(g)
        num = sum(@. w * s * (Y[:,r] - o[:,r]))
        g[r] = num / den
    end
    if ridge_to_1 > 0
        @. g = g - ridge_to_1 * (g - 1.0)
    end
    gm = mean(g)
    if isfinite(gm) && gm > 0 && gm != 1
        @. g = g / gm
        @. s = s * gm
    end
    nothing
end

# """
#     fitoffsets!(o::Matrix{Float64}, Y::Matrix{Float64}, g::Vector{Float64},
#                 s::Vector{Float64}, F, w::Vector{Float64})

# Estimate per-replicate offsets by solving

#     (diag(w) + λ D'D) * o_r = w .* (Y[:,r] − g[r] * s)

# for each column `r`, using the supplied factor `F`, then effect-code rows to
# zero mean.

# # Arguments
# - `o`: offsets (N × R), updated in place.
# - `Y`: data matrix (N × R).
# - `g`: gains vector (length R).
# - `s`: shared signal (length N).
# - `F`: linear solver/factor supporting `\` for the system above (e.g. from [`cholridge`]).
# - `w`: weights (length N).

# # Returns
# `nothing` (mutates `o`).
# """
function fitoffsets!(o::Matrix{Float64}, Y::Matrix{Float64}, g::Vector{Float64},
                     s::Vector{Float64}, F, w::Vector{Float64})
    R = size(Y, 2)
    @inbounds for r in 1:R
        rhs = @. w * (Y[:,r] - g[r] * s)
        o[:,r] .= F \ rhs
    end
    o .-= mean(o, dims=2)  # effect coding
    nothing
end

# """
#     fitquadvar(s_list::Vector{Vector{Float64}},
#                R_list::Vector{Matrix{Float64}}; kwargs...)

# Two-stage, **scale-free** pooled fit of quadratic variance parameters

#     σ²(s) = σ₀² + ϕ*s + κ*s²

# across batches:
# 1) estimate `σ₀²` robustly from low-leverage (`z = s/c`) rows,  
# 2) fit `(ϕ*c, κ*c²)` by nonnegative LS with optional ridge and Huber IRLS.

# # Arguments
# - `s_list`: list of signal vectors per batch.
# - `R_list`: list of residual matrices per batch (N_b × R_b).

# # Keyword arguments
# - `quantlevel=0.97`: quantile for pooled scale `c` (defines `z = s/c`).
# - `bins=0`: if ≥3, equal-frequency binning of `(z, v)` for robustness.
# - `var_stat=:var`: per-row variance statistic (`:var` or `:mad2`).
# - `σ₀²_min`, `ϕ_min`, `κ_min`: nonnegativity lower bounds.
# - `ϕ_ridge`, `κ_ridge`: ridge strengths on `(ϕ*c, κ*c²)`.
# - `z_max=6.0`: clamp `z` to `[0, z_max]`.
# - `c_eps=1e-12`: minimum pooled scale.
# - `v_clip_q=0.99`: winsorize `v` at this upper quantile if `0.5<v_clip_q<1`.
# - `z_low_q=0.10`: low-`z` fraction for robust `σ₀²` estimate.
# - `feas_tol=1e-12`: feasibility slack for the bounded LS.

# # Returns
# `(params::QuadVarParams, c::Float64)` where `c` is the pooled scale.
# """
function fitquadvar(
    s_list::Vector{Vector{Float64}},
    R_list::Vector{Matrix{Float64}};
    quantlevel::Float64=0.97,
    bins::Int=0,
    var_stat::Symbol=:var,  # :mad2 or :var
    σ₀²_min::Float64=0.0,
    ϕ_min::Float64=0.0,
    κ_min::Float64=0.0,
    ϕ_ridge::Float64=0.0,
    κ_ridge::Float64=0.0,
    z_max::Float64=6.0,
    c_eps::Float64=1e-12,
    v_clip_q::Float64=0.99,
    z_low_q::Float64=0.10,
    feas_tol::Float64=1e-12       
)
    # pooled scale c (for z = s/c)
    s_all = vcat(s_list...)
    c = max(quantile(s_all, quantlevel), c_eps)

    # collect z and per-row variances
    z_all = Float64[]; v_all = Float64[]
    @inbounds for (s, R) in zip(s_list, R_list)
        v = if var_stat === :mad2
            med = vec(median(R, dims=2))
            mad = vec(median(abs.(R .- med), dims=2))
            (1.4826 .* mad) .^ 2
        elseif var_stat === :var
            vec(var(R; dims=2))
        else
            error("var_stat must be :mad2 or :var")
        end
        append!(z_all, s ./ c)
        append!(v_all, v)
    end

    # optional equal-frequency binning (robust compression)
    if bins ≥ 3
        qs = collect(range(0.0, 1.0; length=bins+1))
        edges = quantile(z_all, qs)
        zb = Float64[]; vb = Float64[]
        for j in 1:bins
            lo, hi = edges[j], edges[j+1]
            keep = (j == bins) ? ((z_all .≥ lo) .& (z_all .≤ hi)) : ((z_all .≥ lo) .& (z_all .< hi))
            idx = findall(keep)
            if !isempty(idx)
                push!(zb, median(z_all[idx]))
                push!(vb, median(v_all[idx]))
            end
        end
        z_all, v_all = zb, vb
    end

    N = length(v_all)
    if N == 0
        return QuadVarParams(1.0, 0.0, 0.0), c
    end

    # clamp leverage and winsorize v (scale-free)
    zc = clamp.(z_all, 0.0, z_max)
    if 0.5 < v_clip_q < 1.0
        v_hi = quantile(v_all, v_clip_q)
        @. v_all = min(v_all, v_hi)
    end

    # ---- Stage 1: σ₀² from low-z only (robust, data-driven) ----
    z_cut = quantile(z_all, z_low_q)
    idx_lo = findall(z_all .≤ z_cut)
    σ₀²_hat = isempty(idx_lo) ? median(v_all) : median(view(v_all, idx_lo))
    σ₀²_hat = max(σ₀²_hat, σ₀²_min)  # 0 allowed if supported

    # ---- Stage 2: fit (ϕ * c, κ * c^2) on y = v − σ₀²_hat with bounds ----
    y = @. max(v_all - σ₀²_hat, 0.0)     # keep nonnegative target
    X = hcat(zc, zc.^2)                  # columns correspond to (ϕ * c, κ * c^2)

    # optional tiny ridge on (ϕ, κ) -> append rows
    if (ϕ_ridge > 0.0) || (κ_ridge > 0.0)
        X = vcat(X, [sqrt(ϕ_ridge) 0.0; 0.0 sqrt(κ_ridge)])
        y = vcat(y, 0.0, 0.0)
    end

    # robust IRLS (Huber) for the 2-parameter NNLS
    βϕκ = zeros(2)
    lb = (ϕ_min * c, κ_min * c * c)
    w = ones(size(X,1))
    for iter in 1:3
        sw = sqrt.(w)
        Xw = X .* sw
        yw = y .* sw
        β1, β2 = boundedlsβ₂(Xw, yw, lb, feas_tol=feas_tol)
        βϕκ .= (β1, β2)

        # residuals on real rows only (exclude ridge rows)
        n_real = N
        r = @view(y[1:n_real]) .- (@view(X[1:n_real, :]) * βϕκ)
        r̃ = r .- median(r)
        s_r = 1.4826 * median(abs.(r̃))
        s_r = isfinite(s_r) && s_r > 0 ? s_r : sqrt(mean(r.^2) + eps())
        k = 1.345 * s_r
        @inbounds for i in 1:n_real
            ai = abs(r[i])
            w[i] = (ai ≤ k) ? 1.0 : (k / ai)
        end
        iter == 3 && break
    end

    # back to physical parameters
    σ₀² = σ₀²_hat
    ϕ = max(βϕκ[1] / c, ϕ_min)
    κ = max(βϕκ[2] / (c * c), κ_min)

    QuadVarParams(σ₀², ϕ, κ), c
end

# """
#     makediff2(n::Integer, ::Type{T}=Float64) where {T}

# Construct the second-difference operator `D` with Dirichlet-like interior rows:

# Main diagonal = −2, first off-diagonals = +1. The returned matrix has size
# `(n-2) × n` (interior rows only).

# # Arguments
# - `n`: length of the 1D grid (must be ≥ 3).
# - `T=Float64`: element type of the sparse operator.

# # Returns
# A sparse matrix `D::SparseMatrixCSC{T}` of size `(n-2) × n`.

# # Notes
# - Useful as a smoothness operator; `D'*D` is a discrete 1D Laplacian.
# - Endpoints are not included (no boundary rows).
# """
function makediff2(n::Integer, ::Type{T}=Float64) where {T}
    # Require at least 3 grid points; otherwise no interior second differences exist
    n ≥ 3 || throw(ArgumentError("need n ≥ 3"))
    
    # Build diagonals of the tridiagonal stencil [1, -2, 1]
    main = fill(T(-2), n)
    off  = fill(T( 1), n - 1)
    
    # Assemble full n × n banded matrix with the stencil
    Dloc = spdiagm(-1 => off, 0 => main, 1 => off)
    
    # Keep only interior rows (2 to n - 1), giving shape (n - 2) × n
    Dloc[2:n-1, :]
end

# """
#     makeDTD(n::Integer, ::Type{T}=Float64) where {T}

# Build the discrete 1D smoothing matrix `D' * D` based on the second-difference
# operator from [`makediff2`].

# # Arguments
# - `n`: length of the 1D grid (must be ≥ 3).
# - `T=Float64`: element type of the operator.

# # Returns
# A symmetric positive semi-definite sparse matrix `D' * D :: SparseMatrixCSC{T}`
# of size `n × n`.
# """
makeDTD(n::Integer, ::Type{T}=Float64) where {T} = (D = makediff2(n, T); D' * D)

# """
#     pooledfitstats(Ms, o, g, s, p; qc_varfloor=1e-12)

# Compute pooled per-observation z-scores across batches and summarize fit quality.

# # Arguments
# - `Ms`: list of data matrices (N_b × R_b).
# - `o`: list of offsets (N_b × R_b).
# - `g`: list of gains vectors (length R_b).
# - `s`: list of signal vectors (length N_b).
# - `p`: variance parameters.

# # Keyword arguments
# - `qc_varfloor=1e-12`: variance floor for QC standardization.

# # Returns
# Named tuple `(z_rms, cov68, cov95, s_min, s_max, nz)` where:
# - `z_rms` is the pooled RMS of z-scores
# - `cov68`/`cov95` are coverage rates for |z| ≤ 1 and ≤ 1.96
# - `s_min`/`s_max` bound the shared signal across batches
# - `nz` is the total number of **finite** z-scores used
# """
function pooledfitstats(Ms::Vector{<:AbstractMatrix{<:Real}},
                        o::Vector{<:AbstractMatrix{<:Real}},
                        g::Vector{<:AbstractVector{<:Real}},
                        s::Vector{<:AbstractVector{<:Real}},
                        p::QuadVarParams; qc_varfloor::Float64=1e-12)

    sumsq = 0.0
    nfin  = 0
    n68   = 0
    n95   = 0
    s_min = +Inf
    s_max = -Inf

    for b in eachindex(Ms)
        v_b   = varpred(s[b], p; varfloor=qc_varfloor)
        sd_b  = sqrt.(v_b)
        denom = @. max(sd_b, eps(Float64))
        Rb    = size(Ms[b], 2)

        @inbounds for r in 1:Rb
            μr = @. o[b][:, r] + g[b][r] * s[b]
            z  = (Ms[b][:, r] .- μr) ./ denom

            # finite mask
            m = isfinite.(z)
            if any(m)
                zv = @view z[m]
                n  = length(zv)
                nfin += n
                sumsq += sum(abs2, zv)
                az = abs.(zv)
                n68 += count(<=(1.0), az)
                n95 += count(<=(1.96), az)
            end
        end

        s_min = min(s_min, minimum(s[b]))
        s_max = max(s_max, maximum(s[b]))
    end

    z_rms = (nfin == 0) ? NaN : sqrt(sumsq / nfin)
    cov68 = (nfin == 0) ? 0.0 : n68 / nfin
    cov95 = (nfin == 0) ? 0.0 : n95 / nfin
    (z_rms=z_rms, cov68=cov68, cov95=cov95, s_min=s_min, s_max=s_max, nz=nfin)
end

# """
#     pooledresidacf(Y, o, g, s, p;
#                    lag=1, qc_varfloor=1e-12,
#                    restrict=:none, z_quant=0.97, z_max_keep=0.30,
#                    flat_q=0.20, trimfrac=0.0)
#
# Pool per-series residual autocorrelations via Fisher–z weighting.
#
# # Arguments
# - `Y`: list of data matrices (N_b × R_b).
# - `o`: list of offsets (N_b × R_b).
# - `g`: list of gains vectors (length R_b).
# - `s`: list of signal vectors (length N_b).
# - `p`: variance parameters.
#
# # Keyword arguments
# - `lag=1`: autocorrelation lag (≥ 1).
# - `qc_varfloor=1e-12`: variance floor for standardization inside ACF.
# - `restrict=:none`: how to mask observations before ACF.
#     - `:none`     → use all rows.
#     - `:low`      → keep rows with low normalized signal z = s/c′.
#     - `:flat`     → keep rows with small local slope (|Δs| below a quantile).
#     - `:lowflat`  → apply both low and flat masks.
# - `z_quant=0.97`: quantile used to form the pooled scale `c′` from all `s`.
#     Used **only** when `restrict ∈ (:low, :lowflat)`.
# - `z_max_keep=0.30`: keep rows with `s/c′ ≤ z_max_keep` (when low is active).
# - `flat_q=0.20`: keep rows where both adjacent |Δs| are ≤ the `flat_q` quantile.
# - `trimfrac=0.0`: pair-level trimming fraction applied inside `residautocorr`.
#
# # Returns
# Named tuple `(acf, lag, n_pairs)`; `acf` is `NaN` if no series contributes.
# """
function pooledresidacf(
    Y::Vector{<:AbstractMatrix{<:Real}},
    o::Vector{<:AbstractMatrix{<:Real}},
    g::Vector{<:AbstractVector{<:Real}},
    s::Vector{<:AbstractVector{<:Real}},
    p::QuadVarParams;
    lag::Int=1,
    qc_varfloor::Float64=1e-12,
    restrict::Symbol=:none,
    z_quant::Float64=0.97,
    z_max_keep::Float64=0.30,
    flat_q::Float64=0.20,
    trimfrac::Float64=0.0
)
    Σz_w = 0.0
    Σw = 0.0
    total_pairs = 0

    # ----- Build pooled scale c′ only if needed (low-intensity modes)
    cprime = 1.0
    if restrict === :low || restrict === :lowflat
        sall = vcat(s...)
        cprime = max(quantile(sall, z_quant), eps(Float64))
    end

    for b in eachindex(Y)
        varmodel = x -> varpred(x, p)

        # ----- Optional per-batch mask from s[b]
        mask_b = nothing
        if restrict != :none
            sb = s[b]
            keep = trues(length(sb))

            if restrict == :low || restrict == :lowflat
                z = sb ./ cprime
                keep .&= (z .≤ z_max_keep)
            end
            if restrict == :flat || restrict == :lowflat
                ds  = abs.(diff(sb))
                thr = quantile(ds, flat_q)
                flat = trues(length(sb))
                flat[2:end] .&= (ds .≤ thr)
                flat[1:end-1] .&= (ds .≤ thr)
                keep .&= flat
            end
            mask_b = keep
        end

        # compute per-series ACF with optional mask + trimming
        _, ρ_series_b, n_pairs_b = residautocorr(
            Y[b], o[b], g[b], s[b], varmodel;
            lag=lag, mask=mask_b, varfloor=qc_varfloor, trimfrac=trimfrac
        )

        @inbounds for j in eachindex(ρ_series_b)
            n = n_pairs_b[j]; ρ = ρ_series_b[j]
            if isfinite(ρ) && n ≥ 4
                ρc = clamp(ρ, -prevfloat(1.0), prevfloat(1.0))
                z = atanh(ρc); w = n - 3
                Σz_w += w * z; Σw += w; total_pairs += n
            end
        end
    end

    acf = (Σw > 0) ? tanh(Σz_w / Σw) : NaN
    (acf=acf, lag=lag, n_pairs=total_pairs)
end

# """
#     rescalevarparams(p::QuadVarParams,
#                      s_list::Vector{Vector{Float64}},
#                      R_list::Vector{Matrix{Float64}}; kwargs...)

# Post-hoc **global variance scale** calibration. Computes pooled standardized
# RMS residuals, then rescales all parameters by a single factor `β` so that
# typical z-scores are near 1, preserving the `(σ₀² : ϕ : κ)` shape.

# # Keyword arguments
# - `varfloor=1e-12`: minimum variance when forming z-scores.
# - `trim_lo=0.20`, `trim_hi=0.80`: central fraction used for the robust RMS.
# - `clip_lo=0.7`, `clip_hi=1.4`: clamp on the scale factor `β`.

# # Returns
# A new `QuadVarParams(β*σ₀², β*ϕ, β*κ)`.
# """
function rescalevarparams(
    p::QuadVarParams,
    s_list::Vector{Vector{Float64}},
    R_list::Vector{Matrix{Float64}};
    varfloor::Float64=1e-12,
    trim_lo::Float64=0.20,
    trim_hi::Float64=0.80,
    clip_lo::Float64=0.7,
    clip_hi::Float64=1.4)

    zr = Float64[]
    @inbounds for (s, R) in zip(s_list, R_list)
        v  = varpred(s, p; varfloor=varfloor)          # length T
        r2 = vec(sum(abs2, R; dims=2))                 # length T (sum over replicates)
        z_row = sqrt.(r2 ./ size(R,2)) ./ sqrt.(v)     # per-row z_rms
        append!(zr, z_row)                             # append vector
    end
    isempty(zr) && return p

    sort!(zr)
    n  = length(zr)
    i1 = clamp(floor(Int, trim_lo * n), 1, n)
    i2 = clamp(ceil(Int, trim_hi * n), 1, n)
    # ensure i1 ≤ i2 (write plainly to avoid parser weirdness)
    i1, i2 = min(i1, i2), max(i1, i2)

    m = median(view(zr, i1:i2))
    (isfinite(m) && m > 0) || return p

    β = clamp(1.0 / (m * m), clip_lo, clip_hi)          # scale factor on variance
    QuadVarParams(β * p.σ₀², β * p.ϕ, β * p.κ)
end

# """
#     residuals(Y, o, s, g)

# Convenience: compute residuals

#     R = Y − o − s * g'

# # Arguments
# - `Y`: data matrix (N × R).
# - `o`: offsets (N × R).
# - `s`: shared signal (length N).
# - `g`: gains (length R).

# # Returns
# Residual matrix `R` of size (N × R).
# """
@inline residuals(Y, o, s, g) = (Y .- o .- s*g')

# """
#     updatesignal!(s::Vector{Float64}, Y::Matrix{Float64}, o::Matrix{Float64},
#                   g::Vector{Float64}, w::Vector{Float64}; damping::Float64=0.5)

# Perform one **damped weighted LS** update of the shared signal:

#     s ← (1 − damping) * s + damping * (num ./ den)

# with `num = w .* (Y*g − o*g)` and `den = (g⋅g) * w`, followed by a
# nonnegativity projection.

# # Arguments
# - `s`: shared signal (length N), updated in place.
# - `Y`: data matrix (N × R).
# - `o`: offsets (N × R).
# - `g`: gains (length R).
# - `w`: weights (length N).

# # Keyword arguments
# - `damping=0.5`: step size in `[0, 1]`.

# # Returns
# `nothing` (mutates `s`).
# """
function updatesignal!(s::Vector{Float64}, Y::Matrix{Float64}, o::Matrix{Float64},
                       g::Vector{Float64}, w::Vector{Float64}; damping::Float64=0.5)
    T = length(s)

    # Allocate workspace for numerator and denominator
    num = zeros(Float64, T); den = zeros(Float64, T)

    # Compute Y * g (projection of data onto gains)
    mul!(num, Y, g)

    # Compute o * g (projection of offsets onto gains)
    mul!(den, o, g)

    # Weighted numerator: w .* (Y * g − o * g)
    @. num = w * (num - den)

    # Weighted denominator: (g⋅g) * w
    g2 = dot(g, g)
    @. den = g2 * w

    # Damped update: convex combination of old s and LS update (num ./ den)
    s_new = similar(s)
    @. s_new = (1 - damping) * s + damping * (num / max(den, eps(Float64)))

    # Nonnegativity projection: ensure signal is nonnegative
    @. s = max(s_new, 0.0)
    nothing
end

# ── residautocorr ─────────────────────────────────────────────────────────────────────────

"""
    residautocorr(Y, o, g, s, varmodel; lag=1, mask=nothing,
                  varfloor=1e-12, trimfrac=0.0)

Compute the lag-`lag` residual autocorrelation for each replicate series.

# Arguments
- `Y`: data matrix (N × R).
- `o`: offset matrix (N × R).
- `g`: gains vector (length R).
- `s`: shared signal vector (length N).
- `varmodel`: function mapping `s` → variance.

# Keyword arguments
- `lag::Int=1`: autocorrelation lag (≥ 1).
- `mask::Union{Nothing, AbstractVector{Bool}}=nothing`: optional observation
  mask (length N); rows with `mask[i] == false` are excluded.
- `varfloor::Float64=1e-12`: variance floor to avoid division by zero.
- `trimfrac::Float64=0.0`: fraction (0 ≤ trimfrac < 0.5) of largest standardized
  residual pairs to trim before computing autocorrelation.

# Returns
Tuple `(ρ_all, ρ_series, n_pairs)`:
- `ρ_all`: overall pooled ACF across replicates (Fisher–z weighted).
- `ρ_series`: vector of per-replicate ACF values.
- `n_pairs`: vector of counts of contributing lagged pairs per replicate.
"""
function residautocorr(
    Y::AbstractMatrix{<:Real},
    offsets::AbstractMatrix{<:Real},
    gains::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    varmodel::Function;
    lag::T1=1,
    mask::T2=nothing,
    varfloor::T3=1e-12,
    trimfrac::T4=0.0
    ) where {
        T1<:Integer,
        T2<:Union{Nothing, <:AbstractVector{Bool}},
        T3<:Real,
        T4<:Real
    }

    n_obs, n_series = size(Y)
    
    # Validate inputs / shapes / ranges
    lag ≥ 1 || throw(ArgumentError("lag must be ≥ 1"))
    size(offsets) == (n_obs, n_series) || throw(
        ArgumentError("offsets must be N×S matrix, got $(size(offsets))"))
    length(gains) == n_series || throw(
        ArgumentError("gains must be S vector, got length $(length(gains))"))
    length(x) == n_obs || throw(
        ArgumentError("covariate vector x must have length N, got $(length(x))"))
    0.0 ≤ trimfrac < 0.5 || throw(
        ArgumentError("trimfrac must be in [0, 0.5)"))  

    # Variance per observation, then SD (heteroskedastic scale)
    σ² = max.(varmodel.(x), varfloor)   # floor prevents zero/neg variances
    σ  = sqrt.(σ²)

    # Optional observation mask: true => include observation
    valid = isnothing(mask) ? trues(n_obs) : BitVector(mask)
    @assert length(valid) == n_obs

    # Outputs per series and pooled stats across series
    ρ_series = fill(Float64(NaN), n_series)  # per-series lag-`lag` autocorr
    n_pairs_series = fill(0, n_series)       # number of usable pairs per series
    Σz_w = 0.0                               # sum of w_i * z_i (Fisher-z)
    Σw   = 0.0                               # sum of weights w_i = n_i - 3

    for s_idx in 1:n_series
        # Standardize residuals using model mean and variance model
        μ = offsets[:, s_idx] .+ gains[s_idx] .* x
        ε = Y[:, s_idx] .- μ
        z = ε ./ σ  # standardized residuals (unitless)

        # Base keep mask: honor `valid` and drop non-finite values early
        basekeep = valid .& isfinite.(z)
        z_used = z[basekeep]
        length(z_used) ≤ lag && continue  # cannot form any lag-`lag` pair

        # Optional trimming by absolute magnitude (pair-level semantics)
        keep = trues(length(z_used))
        if trimfrac > 0
            az = abs.(z_used)
            _, qhi = quantile(az, (trimfrac, 1 - trimfrac))
            keep = az .≤ qhi
        end

        # Pair-level mask: both ends of a lagged pair must survive trimming
        pairkeep = keep[1:end-lag] .& keep[1+lag:end]
        n = count(pairkeep)
        n < 4 && continue  # Fisher weight (n-3) would be ≤ 0

        # Construct lagged vectors with pair-level filtering
        idx = findall(pairkeep)
        @views z_prev = z_used[1:end-lag][idx]
        @views z_next = z_used[1+lag:end][idx]

        # Correlation at lag `lag`; then Fisher-z for pooling across series
        ρ = cor(z_prev, z_next)
        ρ_series[s_idx] = ρ
        n_pairs_series[s_idx] = n

        # Clip away from ±1 to avoid atanh blow-up
        ρ_clipped = clamp(ρ, -prevfloat(1.0), prevfloat(1.0))
        z_fisher = atanh(ρ_clipped)

        # Inverse-variance weight: Var(z) ≈ 1/(n-3) ⇒ w = n-3
        w = n - 3
        Σz_w += w * z_fisher
        Σw += w
    end

    # Pooled correlation across series: mean Fisher-z back-transformed
    ρ̄ = (Σw > 0) ? tanh(Σz_w / Σw) : NaN
    
    # Return results
    ρ̄, ρ_series, n_pairs_series
end

# ── vif ───────────────────────────────────────────────────────────────────────────────────

"""
    vif(ρ::Real, n::Integer; ρmax=0.8, nonnegative=true, nmin=1)

Return the variance inflation factor (VIF) associated with correlation ρ and sample
size n. The factor is defined as

    VIF = 1 + 2 * ρ * (n_eff - 1) / n_eff

where n_eff = max(n, nmin) ensures a minimum effective sample size. Non-finite ρ
values are treated as 0, and correlations are sanitized (clipped and capped) before
use. The result is always ≥ 1 (variance cannot be deflated). The keyword arguments are
`ρmax` (maximum absolute correlation after capping), `nonnegative` (clip negative
correlations to 0 when true; otherwise cap symmetrically), and `nmin` (minimum sample
size for n_eff). Returns a floating-point variance inflation factor of the same type
as ρ.

See also
[`binretentions`](@ref)
"""
vif(ρ::Real, n::Integer; kwargs...) = vif(float(ρ), n; kwargs...)

function vif(ρ::T1, n::Integer; ρmax::Real=0.8, nonnegative::Bool=true, nmin::Integer=1
    ) where {T1<:AbstractFloat}

    # Validate simple invariants
    0 ≤ ρmax ≤ 1 || throw(ArgumentError("ρmax must be in [0, 1]"))
    nmin ≥ 1 || throw(ArgumentError("nmin must be ≥ 1"))

    # Promote scalars to the working float type T1
    ρcap = T1(ρmax)
    neff = T1(max(n, nmin))

    # Sanitize correlation input
    ρ₀ = isfinite(ρ) ? ρ : zero(T1)  # drop NaN/Inf → 0
    ρ₁ = nonnegative ? max(ρ₀, zero(T1)) : ρ₀  # optionally clip negatives
    ρ₂ = nonnegative ? min(ρ₁, ρcap) : clamp(ρ₁, -ρcap, ρcap)  # cap extreme correlations

    # Variance inflation factor formula
    VIF = one(T1) + T1(2) * ρ₂ * (neff - one(T1)) / neff

    # Enforce VIF ≥ 1 (cannot reduce variance)
    max(VIF, one(T1))
end
