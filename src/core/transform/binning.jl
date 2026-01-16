# ── binretentions ─────────────────────────────────────────────────────────────────────────

"""
    binretentions(
        retentions::AbstractVector{<:Number},
        intensities::AbstractVector{<:Number}, 
        bin_edges::AbstractVector{<:Number}, 
        variances::AbstractVector{<:Number};
        zero_threshold::Number=1e-8, 
        rho_lag1::Real=0.0, 
        rho_max::Real=0.8
    ) -> (retention_centers::Vector, bin_intensities::Vector, bin_vars::Vector)

Return bin centers, binned intensities, and propagated variances by aggregating
`retentions` and their `corresponding` intensity pairs into the bins defined by 
`bin_edges` using inverse-variance weighted means.

Each observation contributes to the bin covering its retention value. Variances are
clamped to `zero_threshold` before weighting, the bin mean is computed as
`Σ wᵢ yᵢ / Σ wᵢ` with `wᵢ = 1 / Varᵢ`, and reported variances are inflated by a lag-1
correlation term controlled by `rho_lag1` and `rho_max`. Retentions and `bin_edges`
must both be unitless or share the same physical dimension. Bin edges must be a
monotonic vector defining left-closed, right-open intervals `[eᵢ, eᵢ₊₁)`; the final edge
is exclusive. Unitful intensities require variances with squared intensity units, and 
negative intensities are preserved so baseline-corrected signals can contribute. The 
function returns `(centers, bin_means, bin_vars)` and throws `ArgumentError` for empty 
inputs, mismatched lengths, negative variances, or invalid bin edges, and 
`Unitful.DimensionError` for incompatible units.

See also
[`vif`](@ref JuChrom.vif),
[`varpredbias`](@ref JuChrom.varpredbias(::Number, ::QuadVarParams)),
[`QuadVarParams`](@ref JuChrom.QuadVarParams).

# Examples
```jldoctest
julia> rets = collect(0.0:0.5:4.5);
       ints = [10, 12, 9, 4, 6, 7, 3, 2, 1, 0];
       edges = 0.0:1.0:5.0;
       vars = fill(0.25, length(ints));

julia> centers, means, varŝ = binretentions(rets, ints, edges, vars);

julia> centers
5-element Vector{Float64}:
 0.5
 1.5
 2.5
 3.5
 4.5

julia> means[1] ≈ (10 + 12) / 2
true
```
"""
function binretentions(
    retentions::AbstractVector{<:Number},
    intensities::AbstractVector{<:Number},
    bin_edges::AbstractVector{<:Number},
    variances::AbstractVector{<:Number};
    zero_threshold::Number=1e-8,
    rho_lag1::Real=0.0,
    rho_max::Real=0.8)

    # Sanity checks: make sure all arrays exist and align before any computation.
    isempty(retentions) && throw(ArgumentError("retentions empty"))
    isempty(intensities) && throw(ArgumentError("intensities empty"))
    isempty(bin_edges) && throw(ArgumentError("bin_edges empty"))
    isempty(variances) && throw(ArgumentError("variances empty"))

    length(retentions) == length(intensities) || throw(
        ArgumentError("retentions and intensities must have same length"))
    length(retentions) == length(variances)   || throw(
        ArgumentError("retentions and variances must have same length"))

    any(v -> v < zero(v), variances) && throw(ArgumentError("All variances must be ≥ 0"))

    n_bins = length(bin_edges) - 1
    n_bins > 0 || throw(ArgumentError("Need at least two bin edges"))

    r₁ = first(retentions)
    e₁ = first(bin_edges)
    v₁ = first(variances)
    i₁ = first(intensities)

    # Track whether representative samples carry units so we can enforce consistency.
    r_unitful = isunitful(r₁)
    e_unitful = isunitful(e₁)
    v_unitful = isunitful(v₁)
    z_unitful = isunitful(zero_threshold)

    # Retention / edge consistency: coordinates must both be unitless or share units.
    if r_unitful ≠ e_unitful throw(ArgumentError(
            "retentions and bin_edges must both have units or both be unitless"))
    elseif r_unitful
        dimension(r₁) == dimension(e₁) || throw(Unitful.DimensionError(r₁, e₁))
    end

    # Variance / threshold consistency mirrors the retention logic.
    if v_unitful
        if z_unitful
            dimension(v₁) == dimension(zero_threshold) || throw(
                Unitful.DimensionError(zero_threshold, v₁))
        else
            zero_threshold = zero_threshold * oneunit(v₁)
        end
    else
        z_unitful && throw(
            ArgumentError("zero_threshold has units but variances are unitless"))
    end

    # Variance units must match squared intensity units when intensities carry units.
    if isunitful(i₁)
        if v_unitful
            dim_int = Unitful.dimension(i₁)
            dim_var = Unitful.dimension(v₁)
            dim_int_sq = dim_int * dim_int
            dim_var == dim_int_sq || throw(Unitful.DimensionError(v₁, i₁ * i₁))
        else
            throw(ArgumentError(
                "variances must have squared intensity units when intensities are unitful"))
        end
    else
        v_unitful && throw(
            ArgumentError("variances have units but intensities are unitless"))
    end

    # Promote representative samples to Float64 numeric part (keep units if present)
    i₁_promoted = isunitful(i₁) ? 1.0 * i₁ : float(i₁)
    v₁_promoted = isunitful(v₁) ? 1.0 * v₁ : float(v₁)

    weight_type = typeof(inv(v₁_promoted))
    term_type = typeof(inv(v₁_promoted) * i₁_promoted)

    # avg_type / var_type capture the eventual element type of the binned results.
    avg_type = isunitful(i₁) ? typeof(i₁_promoted) : Float64
    var_type = isunitful(v₁) ? typeof(v₁_promoted) : Float64

    # Pre-compute bin midpoints reported alongside averages.
    bin_centers = similar(bin_edges, n_bins)
    @inbounds for j in 1:n_bins
        bin_centers[j] = (bin_edges[j] + bin_edges[j+1]) / 2
    end

    # Weighted accumulators implement Σ wᵢ yᵢ / Σ wᵢ, plus counts for variance inflation.
    weighted_sum = Vector{term_type}(undef, n_bins)
    weight_sum = Vector{weight_type}(undef, n_bins)
    point_count = Vector{Int}(undef, n_bins)
    @inbounds for j in 1:n_bins
        weighted_sum[j] = zero(term_type)
        weight_sum[j] = zero(weight_type)
        point_count[j] = 0
    end

    # Pre-compute the target bin index for every retention via vectorized search.
    bin_indices = searchsortedlast.(Base.RefValue(bin_edges), retentions)

    # Sweep all observations once, accumulating weighted moments per bin.
    @inbounds for k in eachindex(retentions)
        intensity = intensities[k]

        variance = variances[k]
        b = bin_indices[k]
        (1 ≤ b ≤ n_bins) || continue
        retentions[k] > bin_edges[end] && continue

        # Clamp tiny variances before inverting to avoid runaway weights.
        v_clamped = variance < zero_threshold ? zero_threshold : variance
        v_clamped_f = isunitful(v_clamped) ? 1.0 * v_clamped : float(v_clamped)
        w = inv(v_clamped_f)

        intensity_f = isunitful(intensity) ? 1.0 * intensity : float(intensity)

        # Accumulate weighted contributions for this bin.
        weighted_sum[b] += w * intensity_f
        weight_sum[b] += w
        point_count[b] += 1
    end

    avg_intensities = Vector{avg_type}(undef, n_bins)
    bin_variances = Vector{var_type}(undef, n_bins)

    # Finalize each bin: compute weighted mean and inflated variance from the moments.
    @inbounds for j in 1:n_bins
        if point_count[j] > 0 && !iszero(weight_sum[j])
            avg_intensities[j] = (weighted_sum[j] / weight_sum[j]) :: avg_type
            # base (iid) variance of weighted mean
            vbar = inv(weight_sum[j]) :: var_type

            # Apply a variance inflation factor (VIF) to account for serial correlation.
            nbin = point_count[j]
            # ρp = rho_lag1 < 0 ? 0.0 : min(rho_lag1, rho_max)
            VIF = vif(float(rho_lag1), nbin; ρmax=rho_max, nonnegative=true, nmin=1)
            #F  = 1 + 2*ρp*(nbin - 1)/nbin
            #bin_variances[j] = vbar * F
            bin_variances[j] = vbar * VIF
        else
            avg_intensities[j] = zero(avg_type)
            if var_type <: Unitful.AbstractQuantity
                bin_variances[j] = Inf * oneunit(v₁_promoted)
            else
                bin_variances[j] = Inf
            end
        end
    end

    # Return midpoints, inverse-variance means, and inflated variance estimates.
    bin_centers, avg_intensities, bin_variances
end

"""
    binretentions(
        msm::MassScanMatrix, 
        bin_edges::AbstractVector{<:Number}, 
        quadvar_params::Union{AbstractVector{<:QuadVarParams}, QuadVarParams}, 
        rho_lag1::Union{AbstractVector{<:Real}, Real};          
        jacobian_scale=nothing::Union{Nothing, AbstractVector{<:Real}, Function}, 
        zero_threshold::Number=1e-8, 
        rho_max::Real=0.8
    ) -> (msm_binned::MassScanMatrix, var_matrix::Matrix)

Return a binned `MassScanMatrix` and per-bin variances by applying
[`binretentions`](@ref JuChrom.binretentions(::AbstractVector{<:Number}, ::AbstractVector{<:Number}, ::AbstractVector{<:Number}, ::AbstractVector{<:Number}))
to each m/z trace in `msm`.

This method differs from the scalar `binretentions` in that it predicts per-scan variances 
via `QuadVarParams`, optionally applies a Jacobian scale `(f'(t))⁻²`, and processes all 
m/z columns in a matrix. If the matrix stores unitful intensities, per-scan and binned 
variances retain squared intensity units, while `msm_binned` stores stripped numeric 
intensities and carries the original metadata. Scalar `quadvar_params` and `rho_lag1` 
inputs are broadcast across m/z columns, and `zero_threshold` is promoted to squared 
intensity units when needed.

See also
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix),
[`QuadVarParams`](@ref JuChrom.QuadVarParams),
[`binretentions`](@ref JuChrom.binretentions(::AbstractVector{<:Number}, ::AbstractVector{<:Number}, ::AbstractVector{<:Number}, ::AbstractVector{<:Number})),
[`varpredbias`](@ref JuChrom.varpredbias(::Number, ::QuadVarParams)),
[`vif`](@ref JuChrom.vif).

# Example
```jldoctest
julia> msm = MassScanMatrix(collect(0.0:0.5:2.0), nothing, [100.0, 101.0], nothing,
                            [10 12; 9 11; 8 10; 7 9; 6 8], nothing);
       quad_params = fill(QuadVarParams(0.1, 0.01, 0.0), 2);
       lag1 = fill(0.2, 2);
       edges = 0:1:3.0;

julia> msm_binned, vars = binretentions(msm, edges, quad_params, lag1);

julia> size(msm_binned.intensities)
(3, 2)
```
"""
function binretentions(
    msm::MassScanMatrix,
    bin_edges::AbstractVector{<:Number},
    quadvar_params::Union{AbstractVector{<:QuadVarParams}, QuadVarParams},
    rho_lag1::Union{AbstractVector{<:Real}, Real};
    jacobian_scale::Union{Nothing, AbstractVector{<:Real}, Function}=nothing,
    zero_threshold::Number=1e-8,
    rho_max::Real=0.8)
    
    # Extract shared views to avoid recomputing getters in loops.
    rets = retentions(msm)  # RT values (unitful or unitless)
    ints = rawintensities(msm)
    n_scans, n_mz = size(ints)

    intensity_unit = intensityunit(msm)
    has_intensity_unit = intensity_unit !== nothing
    intensity_quantity = has_intensity_unit ? (one(Float64) * intensity_unit) : nothing
    variance_quantity = has_intensity_unit ? intensity_quantity^2 : nothing
    variance_dimension = has_intensity_unit ? Unitful.dimension(variance_quantity) : nothing

    function attach_variance_units(val)
        if has_intensity_unit
            if isunitful(val)
                Unitful.dimension(val) == variance_dimension || throw(
                    Unitful.DimensionError(val, variance_quantity))
                val
            else
                val * variance_quantity
            end
        else
            isunitful(val) && throw(ArgumentError(
                "msm intensities are unitless but derived variances have units"))
            val
        end
    end

    attach_intensity_units(vec) = has_intensity_unit ? (vec .* intensity_quantity) : vec
    strip_intensity_units(vec) = has_intensity_unit ?
        Unitful.ustrip.(Base.RefValue(intensity_unit), vec) : vec

    param_requires_units(p) = isunitful(p.σ₀²) || isunitful(p.ϕ)

    function promote_for_varpred(val, p)
        if param_requires_units(p)
            has_intensity_unit || throw(ArgumentError(
                "QuadVarParams carry units but msm intensities are unitless"))
            val * intensity_quantity
        else
            val
        end
    end

    zero_threshold_eff = zero_threshold
    if has_intensity_unit
        if isunitful(zero_threshold)
            Unitful.dimension(zero_threshold) == variance_dimension || throw(
                Unitful.DimensionError(zero_threshold, variance_quantity))
        else
            zero_threshold_eff = zero_threshold * variance_quantity
        end
    else
        isunitful(zero_threshold) && throw(ArgumentError(
            "zero_threshold has units but msm intensities are unitless"))
    end

    # Normalize parameter inputs so downstream code can assume one entry per m/z.
    params_vec = quadvar_params isa AbstractVector ? quadvar_params :
                 fill(quadvar_params, n_mz)
    length(params_vec) == n_mz || throw(ArgumentError(
        "quadvar_params length $(length(params_vec)) ≠ n_mz $n_mz"))

    # Broadcast scalar lag-1 correlations across columns if needed.
    rhos = rho_lag1 isa AbstractVector ? rho_lag1 : fill(rho_lag1, n_mz)
    length(rhos) == n_mz || throw(ArgumentError(
        "rho_lag1 length $(length(rhos)) ≠ n_mz $n_mz"))

    # Build the Jacobian scaling factor J^{-2} once (vector or callable input).
    if jacobian_scale === nothing
        Jinv² = fill(1.0, n_scans)
    elseif jacobian_scale isa AbstractVector
        length(jacobian_scale) == n_scans || throw(ArgumentError(
            "jacobian_scale vector must match number of scans"))
        Tj = typeof(inv(jacobian_scale[1])^2)
        Jinv² = Vector{Tj}(undef, n_scans)
        @inbounds for k in 1:n_scans
            J = jacobian_scale[k]
            abs(J) > 0 || throw(ArgumentError("jacobian_scale contains zero at scan $k"))
            Jinv²[k] = inv(J)^2
        end
    else
        J₁ = jacobian_scale(rets[1])
        abs(J₁) > 0 || throw(ArgumentError("jacobian_scale(t) returned zero at scan 1"))
        Tj = typeof(inv(J₁)^2)
        Jinv² = Vector{Tj}(undef, n_scans)
        Jinv²[1] = inv(J₁)^2
        @inbounds for k in 2:n_scans
            J = jacobian_scale(rets[k])
            abs(J) > 0 || throw(ArgumentError("jacobian_scale(t) returned zero at scan $k"))
            Jinv²[k] = inv(J)^2
        end
    end

    # First column: compute per-scan variances from model, then bin with serial inflation.
    first_ints_view = @view ints[:, 1]
    params₁ = params_vec[1]
    ρ₁ = rhos[1]
    first_var_base = varpredbias(promote_for_varpred(first_ints_view[1], params₁), params₁) * Jinv²[1]
    first_var = attach_variance_units(first_var_base)

    # Allocate vars₁ with correct type
    Tvar₁ = typeof(first_var)
    vars₁ = Vector{Tvar₁}(undef, n_scans)
    vars₁[1] = first_var
    @inbounds for k in 2:n_scans
        Y = promote_for_varpred(first_ints_view[k], params₁)
        vars₁[k] = attach_variance_units(varpredbias(Y, params₁) * Jinv²[k])
    end

    first_ints_for_scalar = attach_intensity_units(first_ints_view)
    bin_centers, first_bin_ints, first_bin_vars = binretentions(retentions(msm), 
        first_ints_for_scalar, bin_edges, vars₁; zero_threshold=zero_threshold_eff,
        rho_lag1=ρ₁, rho_max=rho_max)

    n_bins = length(bin_centers)
    raw_first_bin_ints = strip_intensity_units(first_bin_ints)
    im = Matrix{eltype(raw_first_bin_ints)}(undef, n_bins, n_mz)
    vars = similar(first_bin_vars, n_bins, n_mz)
    @inbounds begin
        im[:, 1] .= raw_first_bin_ints  # seed matrices with first m/z results
        vars[:, 1] .= first_bin_vars
    end

    # Remaining m/z columns
    @inbounds for col in 2:n_mz
        mzints = @view ints[:, col]
        p = params_vec[col]
        ρ = rhos[col]

        # Allocate mzvars with correct type
        first_Y = promote_for_varpred(mzints[1], p)
        first_var_mz = attach_variance_units(varpredbias(first_Y, p) * Jinv²[1])
        Tvar = typeof(first_var_mz)
        mzvars = Vector{Tvar}(undef, n_scans)
        mzvars[1] = first_var_mz
        for k in 2:n_scans
            Y = promote_for_varpred(mzints[k], p)
            # Variance propagation mirrors the first-column logic.
            mzvars[k] = attach_variance_units(varpredbias(Y, p) * Jinv²[k])
        end

        # Reuse scalar helper to bin this column with its own variance model + ρ.
        mzints_for_scalar = attach_intensity_units(mzints)
        _, bin_ints, bin_vars = binretentions(retentions(msm), mzints_for_scalar, bin_edges, mzvars;
            zero_threshold=zero_threshold_eff, rho_lag1=ρ, rho_max=rho_max)

        raw_bin_ints = strip_intensity_units(bin_ints)
        im[:, col]   .= raw_bin_ints
        vars[:, col] .= bin_vars
    end

    # Convert bin centers back to raw numeric retentions + unit metadata.
    raw_retentions, retentionunit = (isunitful(first(bin_centers)) ? 
        (Unitful.ustrip.(bin_centers), Unitful.unit(first(bin_centers))) : (bin_centers, 
        nothing))

    # Build the output matrix, cloning metadata to avoid side effects.
    msm_out = MassScanMatrix(
        raw_retentions,
        retentionunit,
        deepcopy(rawmzvalues(msm)),
        mzunit(msm),
        im,
        intensityunit(msm),
        level=deepcopy(level(msm)),
        instrument=deepcopy(instrument(msm)),
        acquisition=deepcopy(acquisition(msm)),
        user=deepcopy(user(msm)),
        sample=deepcopy(sample(msm)),
        extras=deepcopy(extras(msm))
    )

    # Return both the binned matrix and the per-bin variance estimates.
    msm_out, vars
end

# ── binmzvalues ───────────────────────────────────────────────────────────────────────────

"""
    binmzvalues(series::MassScanSeries, mzbin_fn::Function=integer) -> MassScanSeries

Return a new `MassScanSeries` with per-scan m/z values grouped by `mzbin_fn` and
intensities summed within each bin. All other fields are copied from the input series.

The output series uses concrete element types inferred from `mzbin_fn`.

See also
[`AbstractMassScanSeries`](@ref JuChrom.AbstractMassScanSeries),
[`MassScanSeries`](@ref JuChrom.MassScanSeries),
[`integer`](@ref),
[`intensities`](@ref JuChrom.intensities(::AbstractMassScanSeries, ::Integer)),
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassScanSeries, ::Integer)),
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScanSeries, ::Integer)),
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanSeries)),
[`uniquemzvalues`](@ref JuChrom.uniquemzvalues(::AbstractMassScanSeries)),
[`scan`](@ref JuChrom.scan(::AbstractScanSeries, ::Integer)),
[`scans`](@ref  JuChrom.scan(::AbstractScanSeries)).

```jldoctest
julia> scan1 = MassScan(1.0u"s", [101.2, 101.6, 102.1], [10, 20, 30])
       scan2 = MassScan(1.0u"s", [100.9, 101.2, 102.6], [20, 5, 30])
       series = MassScanSeries([scan1, scan2]);

julia> binned = binmzvalues(series);

julia> mzvalues(binned, 1)
2-element Vector{Int64}:
 101
 102

julia> intensities(binned, 1)
2-element Vector{Int64}:
 30
 30
```
"""
function binmzvalues(series::MassScanSeries, mzbin_fn::F=integer) where {F<:Function}
    @assert scancount(series) > 0 "Cannot bin m/z values in an empty series."

    # Use the first scan as a reference for type inference
    reference_scan = first(scans(series))
    reference_mz = first(rawmzvalues(reference_scan))
    binned_mz_type = typeof(mzbin_fn(reference_mz))

    # Infer types from the reference scan
    intensity_type        = eltype(rawintensities(reference_scan))
    retention_type        = typeof(reference_scan.retention)
    mz_unit_type          = typeof(reference_scan.mzunit)
    retention_unit_type   = typeof(reference_scan.retentionunit)
    intensity_unit_type   = typeof(reference_scan.intensityunit)
    metadata_type         = typeof(reference_scan.attrs)

    # Define concrete MassScan type with updated m/z value type
    ScanT = MassScan{
        retention_type,
        retention_unit_type,
        Vector{binned_mz_type},
        mz_unit_type,
        Vector{intensity_type},
        intensity_unit_type,
        Int,
        metadata_type
    }

    # Preallocate new scan array
    new_scans = Vector{ScanT}(undef, scancount(series))

    @inbounds for (i, scan) in enumerate(series)
        raw_mzs  = rawmzvalues(scan)
        raw_ints = rawintensities(scan)

        binned_mzs = mzbin_fn.(raw_mzs)
        unique_mzs = sort(unique(binned_mzs))

        bin_index_map = Dict(bin => idx for (idx, bin) in enumerate(unique_mzs))
        binned_ints = zeros(eltype(raw_ints), length(unique_mzs))

        @inbounds @simd for j in eachindex(raw_mzs)
            bin = binned_mzs[j]
            binned_ints[bin_index_map[bin]] += raw_ints[j]
        end

        new_scans[i] = ScanT(
            scan.retention,
            scan.retentionunit,
            unique_mzs,
            scan.mzunit,
            binned_ints,
            scan.intensityunit,
            scan.level,
            deepcopy(scan.attrs),
        )
    end

    # Infer full concrete MassScanSeries type
    SeriesT = MassScanSeries{
        ScanT,
        retention_unit_type,
        mz_unit_type,
        intensity_unit_type,
        typeof(new_scans),
        typeof(series.instrument),
        typeof(series.acquisition),
        typeof(series.user),
        typeof(series.sample),
        typeof(series.extras)
    }

    SeriesT(
        new_scans,
        deepcopy(series.instrument),
        deepcopy(series.acquisition),
        deepcopy(series.user),
        deepcopy(series.sample),
        deepcopy(series.extras)
    )
end

# ── integer ───────────────────────────────────────────────────────────────────────────────

"""
    integer(value::Real; offset::Real=0.3) -> Int

Assigns a real number `value` to an integer bin such that the integer bin `n` includes
the interval [n - offset, n + (1 - offset)).

See also [`binmzvalues`](@ref).

```julia
julia> integer(3.2)
3

julia> integer(3.69)
3

julia> integer(3.7)
4

julia> integer(3.7, 0.999)
4

julia> integer(3.7, 1)
ERROR: ArgumentError: `offset` must be in the range [0, 1).
```
"""
function integer(value::Real, offset::Real=0.3)
    0 ≤ offset < 1 || throw(ArgumentError("Offset must be in the range [0, 1)."))
    floor(Int, value + offset)
end
