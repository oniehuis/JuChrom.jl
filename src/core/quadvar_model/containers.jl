# ── QuadVarParams ──────────────────────────────────────────────────────────────────────────

"""
    QuadVarParams(σ₀²::Number, ϕ::Number, κ::Number)

Parameter container for a quadratic variance model of the form

    σ²(y) = σ₀² + ϕ * y + κ * y^2

# Fields
- `σ₀²`: baseline variance (intercept term, shares the variance units).
- `ϕ`: linear coefficient in `y` (shares the intensity units).
- `κ`: quadratic coefficient in `y` (dimensionless; multiplies `y^2`).

# Usage
`QuadVarParams` is typically used with [`varpred`] and [`varpredbias`] to
evaluate the variance model given a true or estimated input `y`. All fields accept
plain numbers or `Unitful` quantities; when intensities have units, set `σ₀²` and `ϕ`
using the corresponding squared or linear intensity units so predicted variances
retain their physical dimensions. The constructor enforces this relationship and
also requires `κ` to remain unitless.
"""
struct QuadVarParams{T1<:Number, T2<:Number, T3<:Number}
    σ₀²::T1
    ϕ::T2
    κ::T3

    function QuadVarParams(σ₀²::T1, ϕ::T2, κ::T3) where {T1<:Number, T2<:Number, T3<:Number}
        σ_has_units = isunitful(σ₀²)
        ϕ_has_units = isunitful(ϕ)
        κ_has_units = isunitful(κ)

        if σ_has_units
            intensity_ref = sqrt(Unitful.oneunit(σ₀²))
            ϕ_has_units || throw(Unitful.DimensionError(ϕ, intensity_ref))
            Unitful.dimension(ϕ) == Unitful.dimension(intensity_ref) ||
                throw(Unitful.DimensionError(ϕ, intensity_ref))
        elseif ϕ_has_units
            throw(Unitful.DimensionError(ϕ, one(Float64)))
        end

        κ_has_units && throw(Unitful.DimensionError(κ, one(Float64)))

        new{T1, T2, T3}(σ₀², ϕ, κ)
    end
end

function QuadVarParams{T1,T2,T3}(σ₀²::Number, ϕ::Number, κ::Number
    ) where {T1<:Number, T2<:Number, T3<:Number}
    QuadVarParams(convert(T1, σ₀²), convert(T2, ϕ), convert(T3, κ))
end

Base.broadcastable(p::QuadVarParams) = Base.RefValue(p)

# ── QuadVarFit ──────────────────────────────────────────────────────────────

"""
    QuadVarFit()

Container for fitted **signal/offset/gain + quadratic variance** results
across a selected set of m/z values.

# Fields
- `mz_unit::Union{Unitful.Units, Nothing}`: m/z unit tracked from the input grid
  (or `nothing` if unitless).
- `mz_ref::Vector{<:Real}`: numeric m/z grid present in the input data (values
  expressed in `mz_unit` when that is set).
- `mz_idx::Vector{Int}`: indices into `mz_ref` that were fitted.
- `mz_values::Vector{<:Real}`: selected m/z values (`mz_ref[mz_idx]`, numeric).
- `batchcount::Int`: number of batches.
- `n_reps_per_batch::Vector{Int}`: replicates per batch.
- `n_scans_per_batch::Vector{Int}`: scans per batch.
- `params::Vector{QuadVarParams}`: variance parameters per selected m/z.
- `intensity_unit::Union{Unitful.Units, Nothing}`: intensity unit tracked from input data
  (or `nothing` if unitless).
- `signal::Vector{Vector{Vector{Float64}}}`: fitted signal per m/z, per batch.
- `offsets::Vector{Vector{Matrix{Float64}}}`: offsets per m/z, per batch.
- `gains::Vector{Vector{Vector{Float64}}}`: gains per m/z, per batch.
- `scale_c::Vector{Float64}`: pooled scale `c` per m/z.
- `acf::Vector{Float64}`: pooled residual autocorrelation per m/z.
- `acf_lag::Vector{Int}`: lag used for `acf` (echoed) per m/z.
- `n_pairs::Vector{Int}`: total pairs contributing to `acf` per m/z.
- `qc_z_rms::Vector{Float64}`: pooled z-RMS per m/z.
- `qc_cov68::Vector{Float64}`: |z|≤1 coverage per m/z.
- `qc_cov95::Vector{Float64}`: |z|≤1.96 coverage per m/z.
- `qc_s_min::Vector{Float64}` / `qc_s_max::Vector{Float64}`: signal range per m/z.
- `qc_nz::Vector{Int}`: count of finite standardized residuals per m/z.
- `observed::Vector{Vector{Matrix{Float64}}}`: raw data **restricted to selected m/z**,
  shaped as `observed[b][r] :: n_scans_b × n_mz_select`.

# Notes
- Use `mz_values[j]` for labels (multiply by `mz_unit` to recover Unitful values);
  use `mz_idx[j]` if you need to index into a full grid.
- `observed[b][r][:, j]` matches the **j-th selected** m/z (i.e., `mz_values[j]`).
"""
struct QuadVarFit{T1<:AbstractVector{<:Real}, T2<:QuadVarParams, T3<:Union{Nothing, Unitful.Units}, T4<:Union{Nothing, Unitful.Units}}
    mz_unit::T4
    mz_ref::T1
    mz_idx::Vector{Int}
    mz_values::T1
    batchcount::Int
    n_reps_per_batch::Vector{Int}
    n_scans_per_batch::Vector{Int}
    params::Vector{T2}
    intensity_unit::T3
    signal::Vector{Vector{Vector{Float64}}}
    offsets::Vector{Vector{Matrix{Float64}}}
    gains::Vector{Vector{Vector{Float64}}}
    scale_c::Vector{Float64}
    acf::Vector{Float64}
    acf_lag::Vector{Int}
    n_pairs::Vector{Int}
    qc_z_rms::Vector{Float64}
    qc_cov68::Vector{Float64}
    qc_cov95::Vector{Float64}
    qc_s_min::Vector{Float64}
    qc_s_max::Vector{Float64}
    qc_nz::Vector{Int}
    observed::Vector{Vector{Matrix{Float64}}}
end

# Treat `QuadVarFit` as a scalar in broadcasting contexts
Base.Broadcast.broadcastable(f::QuadVarFit) = Base.RefValue(f)

@inline _with_mz_unit(val, mz_unit::Union{Nothing, Unitful.Units}) =
    isnothing(mz_unit) ? val : val * mz_unit

function Base.summary(q::QuadVarFit)
    n_batches = let x = q.n_reps_per_batch; x === nothing ? 0 : length(x) end
    n_sel     = let x = q.mz_idx; x === nothing ? 0 : length(x) end
    "QuadVarFit (batches=$(n_batches), ions=$(n_sel))"
end

Base.show(io::IO, q::QuadVarFit) = print(io, summary(q))

function Base.show(io::IO, ::MIME"text/plain", q::QuadVarFit)
    println(io, summary(q))

    nreps  = q.n_reps_per_batch
    nscans = q.n_scans_per_batch
    mzref  = q.mz_ref
    mzidx  = q.mz_idx
    mzvals = q.mz_values
    mzunit = q.mz_unit
    unit_suffix = isnothing(mzunit) ? "" : " (" * string(mzunit) * ")"
    mzvals_disp = _with_mz_unit.(mzvals, Base.RefValue(mzunit))

    println(io, "  batches          : ", length(nreps))
    !isempty(nreps)  && println(io, "    reps per batch : ", nreps)
    !isempty(nscans) && println(io, "    scans per batch: ", nscans)

    println(io, "  ion grid         : ", length(mzref), " total; selected: ", length(mzidx))
    if !isempty(mzvals)
        if length(mzvals) ≤ 6
            println(io, "    m/z values", unit_suffix, "     : ", round.(mzvals_disp; sigdigits=6))
        else
            head = join(round.(mzvals_disp[1:3]; sigdigits=6), ", ")
            tail = join(round.(mzvals_disp[end-2:end]; sigdigits=6), ", ")
            println(io, "    m/z values", unit_suffix, "     : [", head, ", …, ", tail, "]")
        end
    end

    # Parameter ranges (from q.params)
    if !isempty(q.params)
        σ0 = getfield.(q.params, :σ₀²)
        ϕ  = getfield.(q.params, :ϕ)
        κ  = getfield.(q.params, :κ)
        println(io, "  params (ranges)  : ",
                "σ₀²∈[", round(minimum(σ0); sigdigits=6), ", ", round(maximum(σ0); sigdigits=6), "], ",
                "ϕ∈[",   round(minimum(ϕ);  sigdigits=6), ", ", round(maximum(ϕ);  sigdigits=6), "], ",
                "κ∈[",   round(minimum(κ);  sigdigits=6), ", ", round(maximum(κ);  sigdigits=6), "]")
    end

    # QC summaries (from q.qc_z_rms, q.qc_cov95)
    fzr  = filter(isfinite, q.qc_z_rms)
    fc95 = filter(isfinite, q.qc_cov95)
    !isempty(fzr)  && println(io, "  z_rms (median)   : ", round(median(fzr);  sigdigits=6))
    !isempty(fc95) && println(io, "  cov95 (median)   : ", round(median(fc95); sigdigits=6))

    # ACF summary (from q.acf)
    acfs = filter(isfinite, q.acf)
    !isempty(acfs) && println(io, "  acf (median)     : ", round(median(acfs); sigdigits=6))
end

# ── varpred ───────────────────────────────────────────────────────────────────────────────

"""
    varpred(y, p::QuadVarParams; varfloor::Number=0)

Evaluate the quadratic variance model **elementwise** using a parameter bundle:

    σ²(y) = p.σ₀² + p.ϕ * y + p.κ * y^2

Works for scalars, vectors, and matrices via broadcasting.

# Arguments
- `y`: true input value(s) at which to evaluate the variance model (number or array).
- `p`: `QuadVarParams` with fields `σ₀²`, `ϕ`, `κ`.

# Keyword arguments
- `varfloor=0`: elementwise lower bound applied to predicted variances; plain numbers
  are automatically converted to the variance units when `σ₀²` carries units. Also
  replaces non-finite results.

# Returns
`σ²(y)` with the same shape as `y`. The element type is the promotion of
`eltype(y)`, the parameter field types in `p`, and `typeof(varfloor)`.

# Notes
- Non-finite intermediate results are replaced by `varfloor`.
- Predicted values retain the units implied by `σ₀²`/`ϕ`.
- Use this when `y` is a true (non-estimated) input. For estimated inputs, see [`varpredbias`].
"""
@inline function varpred(y, p::QuadVarParams; varfloor::Number=0)
    @. apply_variance_floor(muladd(p.κ, y * y, muladd(p.ϕ, y, p.σ₀²)), varfloor)
end

@inline function coerce_like_reference(value, reference, label::AbstractString)
    if isunitful(reference)
        if isunitful(value)
            Unitful.dimension(value) == Unitful.dimension(reference) ||
                throw(Unitful.DimensionError(value, reference))
            value
        else
            value * Unitful.oneunit(reference)
        end
    else
        isunitful(value) && throw(ArgumentError(
            "$label has units but the reference quantity is unitless"))
        value
    end
end

@inline function floorfinite(x::T1, floor::T2) where {T1<:Number, T2<:Number}
    T3 = promote_type(T1, T2)
    x_p = T3(x)
    floor_p = T3(floor)
    (isfinite(x) && x_p ≥ floor_p) ? x_p : floor_p
end

@inline function apply_variance_floor(val, varfloor)
    floor_adj = coerce_like_reference(varfloor, val, "varfloor")
    floorfinite(val, floor_adj)
end

"""
    varpred(y, σ₀²::Number, ϕ::Number, κ::Number; varfloor::Number=0)

Evaluate the quadratic variance model **elementwise** with scalar parameters:

    σ²(y) = σ₀² + ϕ * y + κ * y^2

Works for scalars, vectors, and matrices via broadcasting.

# Arguments
- `y`: true input value(s) (number or array).
- `σ₀²`, `ϕ`, `κ`: scalar parameters of the quadratic variance model.

# Keyword arguments
- `varfloor=0`: elementwise lower bound applied to predicted variances; coerced to the
  correct units when `σ₀²` is unitful. Also replaces non-finite results.

# Returns
`σ²(y)` with the same shape as `y`. The element type is the promotion of
`eltype(y)`, the parameter types, and `typeof(varfloor)`.

# Notes
- Non-finite intermediate results are replaced by `varfloor`.
- Predicted values inherit the units of `σ₀²`.
- If you have a parameter bundle, use [`varpred(y, p::QuadVarParams; varfloor=...)`].
"""
@inline function varpred(y, σ₀²::Number, ϕ::Number, κ::Number; varfloor::Number=0)
    @. apply_variance_floor(muladd(κ, y * y, muladd(ϕ, y, σ₀²)), varfloor)
end

# ── varpredbias ───────────────────────────────────────────────────────────────────────────

"""
    varpredbias(y, p::QuadVarParams; varfloor::Number=0)

Evaluate the quadratic variance model with a **bias correction** (for estimated `y`):

    σ²(y) = (p.σ₀² + p.ϕ * y + p.κ * y^2) / (1 + p.κ)

Works for scalars, vectors, and matrices via broadcasting.

# Arguments
- `y`: estimated input value(s) (number or array).
- `p`: `QuadVarParams` with fields `σ₀²`, `ϕ`, `κ`.

# Keyword arguments
- `varfloor=0`: elementwise lower bound applied to predicted variances; coercion to the
  variance units happens automatically when parameters are unitful. Also replaces
  non-finite results.

# Returns
`σ²(y)` with the same shape as `y`. The element type is the promotion of
`eltype(y)`, the parameter field types in `p`, and `typeof(varfloor)`.

# Notes
- Non-finite intermediate results are replaced by `varfloor`.
- Predicted values retain the units implied by `σ₀²`/`ϕ`.
- This form assumes `1 + p.κ > 0` (typically ensured by enforcing `κ_min ≥ 0`).
- For true (non-estimated) inputs, see [`varpred`].
"""
@inline function varpredbias(y, p::QuadVarParams; varfloor::Number=0)
    den = one(p.κ) + p.κ
    @. apply_variance_floor(muladd(p.κ, y*y, muladd(p.ϕ, y, p.σ₀²)) / den, varfloor)
end
