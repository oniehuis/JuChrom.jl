# ── binretentions ─────────────────────────────────────────────────────────────────────────

function _binretention_trace(
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

function _validate_retention_grid(msm::AbstractMassScanMatrix, rgrid::RetentionGrid)
    data_unit = retentionunit(msm)
    grid_unit = retentionunit(rgrid)

    if isnothing(data_unit)
        isnothing(grid_unit) || throw(ArgumentError(
            "RetentionGrid has a unit but input matrix retentions are unitless"))
    elseif isnothing(grid_unit)
        throw(ArgumentError(
            "RetentionGrid is unitless but input matrix retentions have a unit"))
    else
        data_quantity = one(Float64) * data_unit
        grid_quantity = one(Float64) * grid_unit
        Unitful.dimension(data_quantity) == Unitful.dimension(grid_quantity) ||
            throw(Unitful.DimensionError(data_quantity, grid_quantity))
    end

    nothing
end

function _binretentions_with_grid(
    vmsm::AbstractVarianceMassScanMatrix,
    rgrid::RetentionGrid,
    rho_lag1::Union{AbstractVector{<:Real}, Real};
    zero_threshold::Number=1e-8,
    rho_max::Real=0.8)

    msm = parent(vmsm)
    _validate_retention_grid(msm, rgrid)

    grid_unit = retentionunit(rgrid)
    rets = isnothing(grid_unit) ? rawretentions(msm) : rawretentions(msm; unit=grid_unit)
    bin_edges = rawbinedges(rgrid)
    ints = rawintensities(msm)
    vars_in = variances(vmsm)
    n_scans, n_mz = size(ints)

    size(vars_in) == (n_scans, n_mz) || throw(ArgumentError(
        "variances size $(size(vars_in)) ≠ ($(n_scans), $(n_mz))"))

    intensity_unit = intensityunit(msm)
    has_intensity_unit = intensity_unit !== nothing
    intensity_quantity = has_intensity_unit ? (one(Float64) * intensity_unit) : nothing
    variance_quantity = has_intensity_unit ? intensity_quantity^2 : nothing
    variance_dimension = has_intensity_unit ? Unitful.dimension(variance_quantity) : nothing

    attach_intensity_units(vec) = has_intensity_unit ? (vec .* intensity_quantity) : vec
    strip_intensity_units(vec) = has_intensity_unit ?
        Unitful.ustrip.(Base.RefValue(intensity_unit), vec) : vec

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
            "zero_threshold has units but input matrix intensities are unitless"))
    end

    rhos = rho_lag1 isa AbstractVector ? rho_lag1 : fill(rho_lag1, n_mz)
    length(rhos) == n_mz || throw(ArgumentError(
        "rho_lag1 length $(length(rhos)) ≠ n_mz $n_mz"))

    function prepare_variance_column(col_view)
        first_val = col_view[1]
        if has_intensity_unit
            if isunitful(first_val)
                Unitful.dimension(first_val) == variance_dimension || throw(
                    Unitful.DimensionError(first_val, variance_quantity))
                return col_view
            end

            return col_view .* variance_quantity
        end

        isunitful(first_val) && throw(ArgumentError(
            "input matrix intensities are unitless but variances have units"))
        col_view
    end

    first_ints_view = @view ints[:, 1]
    first_vars = prepare_variance_column(@view vars_in[:, 1])
    first_ints_for_trace = attach_intensity_units(first_ints_view)
    bin_centers, first_bin_ints, first_bin_vars = _binretention_trace(
        rets, first_ints_for_trace, bin_edges, first_vars;
        zero_threshold=zero_threshold_eff, rho_lag1=rhos[1], rho_max=rho_max)

    n_bins = length(bin_centers)
    raw_first_bin_ints = strip_intensity_units(first_bin_ints)
    im = Matrix{eltype(raw_first_bin_ints)}(undef, n_bins, n_mz)
    vars = similar(first_bin_vars, n_bins, n_mz)
    @inbounds begin
        im[:, 1] .= raw_first_bin_ints
        vars[:, 1] .= first_bin_vars
    end

    @inbounds for col in 2:n_mz
        mzints = @view ints[:, col]
        mzvars = prepare_variance_column(@view vars_in[:, col])
        mzints_for_trace = attach_intensity_units(mzints)
        _, bin_ints, bin_vars = _binretention_trace(
            rets, mzints_for_trace, bin_edges, mzvars;
            zero_threshold=zero_threshold_eff, rho_lag1=rhos[col], rho_max=rho_max)

        im[:, col] .= strip_intensity_units(bin_ints)
        vars[:, col] .= bin_vars
    end

    msm_out = MassScanMatrix(
        bin_centers,
        grid_unit,
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

    msm_out, vars
end

"""
    binretentions(
        vmsm::AbstractVarianceMassScanMatrix,
        rgrid::RetentionGrid,
        rho_lag1::Union{AbstractVector{<:Real}, Real};
        zero_threshold::Number=1e-8,
        rho_max::Real=0.8
    ) -> VarianceMassScanMatrix

    binretentions(vmsm::AbstractVarianceMassScanMatrix, rgrid::RetentionGrid; kwargs...)

Return a binned `VarianceMassScanMatrix` by applying `rgrid` to the parent mass-scan
matrix and propagating the stored per-cell variances.

Retentions are converted to the grid's retention unit before bin assignment. The returned
matrix stores bin centers in the grid's raw retention scale and uses the grid's
`retentionunit`. Intensities are inverse-variance weighted within each bin. The stored
variances are propagated as variances of the weighted means and inflated with
[`vif`](@ref JuChrom.vif) using `rho_lag1`.

Scalar `rho_lag1` inputs are broadcast across m/z columns; vector inputs must have one
lag-1 correlation per m/z channel. No Jacobian scaling is applied. If the data were
retention-mapped before binning, apply the mapping to the `VarianceMassScanMatrix` first
so the stored intensities and variances are already on the mapped scale.

See also
[`VarianceMassScanMatrix`](@ref JuChrom.VarianceMassScanMatrix),
[`RetentionGrid`](@ref JuChrom.RetentionGrid),
[`densestgrid`](@ref JuChrom.densestgrid),
[`vif`](@ref JuChrom.vif).
"""
function binretentions(
    vmsm::AbstractVarianceMassScanMatrix,
    rgrid::RetentionGrid,
    rho_lag1::Union{AbstractVector{<:Real}, Real};
    zero_threshold::Number=1e-8,
    rho_max::Real=0.8)

    binned_msm, binned_variances = _binretentions_with_grid(
        vmsm,
        rgrid,
        rho_lag1;
        zero_threshold=zero_threshold,
        rho_max=rho_max)

    VarianceMassScanMatrix(binned_msm, binned_variances)
end

function binretentions(
    vmsm::AbstractVarianceMassScanMatrix,
    rgrid::RetentionGrid;
    zero_threshold::Number=1e-8,
    rho_max::Real=0.8)

    binretentions(
        vmsm,
        rgrid,
        0.0;
        zero_threshold=zero_threshold,
        rho_max=rho_max)
end

# ── binmzvalues ───────────────────────────────────────────────────────────────────────────

"""
    binmzvalues(
        series::MassScanSeries,
        mzbin_fn::Function=integer;
        validmzvalues::Union{Nothing, AbstractVector}=nothing
    ) -> MassScanSeries

Return a new `MassScanSeries` with per-scan m/z values grouped by `mzbin_fn` and
intensities summed within each bin. All other fields are copied from the input series.
If `validmzvalues` is provided, only binned m/z values present in `validmzvalues` are 
retained; values that bin outside the accepted set are discarded.

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
function binmzvalues(
        series::MassScanSeries,
        mzbin_fn::F=integer;
        validmzvalues::Union{Nothing, AbstractVector}=nothing
    ) where {F<:Function}

    @assert scancount(series) > 0 "Cannot bin m/z values in an empty series."

    # Use the first scan as a reference for type inference
    reference_scan = first(scans(series))
    reference_mz = first(rawmzvalues(reference_scan))
    binned_mz_type = typeof(mzbin_fn(reference_mz))
    accepted_mzs = validmzvalues ≡ nothing ? nothing : Set(validmzvalues)

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
        if !isnothing(accepted_mzs)
            unique_mzs = [mz for mz in unique_mzs if mz ∈ accepted_mzs]
        end

        bin_index_map = Dict(bin => idx for (idx, bin) in enumerate(unique_mzs))
        binned_ints = zeros(eltype(raw_ints), length(unique_mzs))

        @inbounds for j in eachindex(raw_mzs)
            bin = binned_mzs[j]
            if !isnothing(accepted_mzs) && !(bin ∈ accepted_mzs)
                continue
            end
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
