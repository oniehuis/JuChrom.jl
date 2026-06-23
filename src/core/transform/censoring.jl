# ── replacecensored ───────────────────────────────────────────────────────────────────────

"""
    CensoredReplacementInfo

Diagnostic summary returned by [`replacecensored`](@ref) when `returninfo=true`.

The object records the replacement parameters, empirical floors, and replacement counts
without storing full-size masks or replacement matrices.
"""
struct CensoredReplacementInfo
    k::Float64
    q::Float64
    floor_rule::Symbol
    global_floor::Float64
    column_floors::Vector{Float64}
    minpositives::Int
    n_replaced::Int
    fraction_replaced::Float64
    n_nonpositive::Int
    n_below_limit_positive::Int
    replaced_by_column::Vector{Int}
end

function Base.show(io::IO, info::CensoredReplacementInfo)
    print(io, "CensoredReplacementInfo(")
    print(io, info.n_replaced, " replaced")
    print(io, ", fraction=", round(info.fraction_replaced; sigdigits=4))
    print(io, ", floor_rule=:", info.floor_rule)
    print(io, ", k=", info.k)
    print(io, ", q=", info.q)
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", info::CensoredReplacementInfo)
    println(io, "CensoredReplacementInfo")
    println(io, "  Replaced cells: ", info.n_replaced,
        " (", round(100 * info.fraction_replaced; sigdigits=4), "%)")
    println(io, "  Non-positive cells replaced: ", info.n_nonpositive)
    println(io, "  Positive below-limit cells replaced: ", info.n_below_limit_positive)
    println(io, "  Floor rule: :", info.floor_rule)
    println(io, "  Global floor: ", info.global_floor)
    println(io, "  Column floors: ", length(info.column_floors), " values")
    println(io, "  k: ", info.k)
    println(io, "  q: ", info.q)
    print(io, "  minpositives: ", info.minpositives)
end

function validatereplacecensoredkwargs(k, q, minpositives)
    k >= 0 || throw(ArgumentError("k must be nonnegative."))
    0 <= q <= 1 || throw(ArgumentError("q must be between 0 and 1."))
    minpositives >= 1 || throw(ArgumentError("minpositives must be at least 1."))
    nothing
end

function positivequantile(vals::AbstractVector{<:Real}, q::Real)
    positives = vals[vals .> zero(eltype(vals))]
    isempty(positives) && return nothing
    Float64(quantile(positives, q))
end

function replacecensoredfloors(
    x::AbstractMatrix{<:Real},
    floor,
    q::Real,
    minpositives::Integer,
)
    n_mz = size(x, 2)
    global_floor = positivequantile(vec(x), q)
    isnothing(global_floor) && throw(ArgumentError(
        "replacecensored requires at least one positive intensity value."))

    if floor === :global_quantile
        floors = fill(global_floor, n_mz)
        return floors, global_floor, :global_quantile
    elseif floor === :column_quantile || floor === :feature_quantile
        floors = Vector{Float64}(undef, n_mz)
        @inbounds for col in 1:n_mz
            column_floor = positivequantile(@view(x[:, col]), q)
            floors[col] = isnothing(column_floor) ||
                count(>(0), @view(x[:, col])) < minpositives ?
                global_floor : column_floor
        end
        return floors, global_floor, :column_quantile
    elseif floor isa Real
        floor > 0 || throw(ArgumentError("floor must be positive."))
        isfinite(floor) || throw(ArgumentError("floor must be finite."))
        floors = fill(Float64(floor), n_mz)
        return floors, global_floor, :scalar
    elseif floor isa AbstractVector{<:Real}
        length(floor) == n_mz || throw(ArgumentError(
            "floor vector length $(length(floor)) must match m/z count $n_mz."))
        all(isfinite, floor) || throw(ArgumentError("floor vector must be finite."))
        all(>(0), floor) || throw(ArgumentError("floor vector must be positive."))
        return Float64.(floor), global_floor, :vector
    end

    throw(ArgumentError(
        "floor must be :column_quantile, :feature_quantile, :global_quantile, " *
        "a positive real value, or a positive vector."))
end

function copyreplacecensoredparent(
    vmsm::AbstractVarianceMassScanMatrix,
    intensities_out::AbstractMatrix{<:Real},
    extras_out::Dict{String,Any},
)
    MassScanMatrix(
        copy(rawretentions(vmsm)),
        retentionunit(vmsm),
        copy(rawmzvalues(vmsm)),
        mzunit(vmsm),
        intensities_out,
        intensityunit(vmsm);
        level=level(vmsm),
        instrument=deepcopy(instrument(vmsm)),
        acquisition=deepcopy(acquisition(vmsm)),
        user=deepcopy(user(vmsm)),
        sample=deepcopy(sample(vmsm)),
        extras=extras_out,
    )
end

"""
    replacecensored(vmsm::AbstractVarianceMassScanMatrix; kwargs...) -> VarianceMassScanMatrix
    replacecensored(vmsm::AbstractVarianceMassScanMatrix; returninfo=true, kwargs...) ->
        (VarianceMassScanMatrix, CensoredReplacementInfo)

Replace non-positive and below-limit intensities in `vmsm` by positive censored-value
estimates and update their variances.

For each intensity cell `x[i, j]` with variance `v[i, j]`, a reliability limit is computed
as `L[i, j] = max(floor[j], k * sqrt(v[i, j]))`. Values with `x[i, j] <= L[i, j]` are
treated as censored and replaced by `L[i, j] / 2`. Their variances are replaced by
`max(v[i, j], L[i, j]^2 / 12)`, corresponding to the uncertainty of an unknown value in
`[0, L[i, j]]` under a uniform censoring model. Values above the limit are copied
unchanged.

By default, `floor=:column_quantile`, so `floor[j]` is the `q` quantile of positive
intensities in m/z column `j`. Columns with fewer than `minpositives` positive values use
the global positive-intensity quantile. Use `floor=:global_quantile` for one shared floor,
a scalar for a fixed floor, or a vector with one floor per m/z column.

Keyword arguments:
- `k::Real=3.0`: multiplier for the propagated standard deviation.
- `q::Real=0.01`: positive-intensity quantile used for empirical floors.
- `floor=:column_quantile`: floor rule, scalar floor, or per-column floor vector.
- `minpositives::Integer=10`: minimum positive values required for a column floor.
- `returninfo::Bool=false`: return a [`CensoredReplacementInfo`](@ref) when `true`.

The returned matrix preserves retention, m/z, units, metadata, extras, and variance units.
"""
function replacecensored(
    vmsm::AbstractVarianceMassScanMatrix;
    k::Real=3.0,
    q::Real=0.01,
    floor=:column_quantile,
    minpositives::Integer=10,
    returninfo::Bool=false,
)
    validatereplacecensoredkwargs(k, q, minpositives)

    x = rawintensities(vmsm)
    v = rawvariances(vmsm)
    size(x) == size(v) || throw(DimensionMismatch(
        "intensity and variance matrices must have identical sizes."))

    floors, global_floor, floor_rule = replacecensoredfloors(x, floor, q, minpositives)

    x_out = Matrix{Float64}(x)
    v_out = Matrix{Float64}(v)
    n_scans, n_mz = size(x)

    n_replaced = 0
    n_nonpositive = 0
    n_below_limit_positive = 0
    replaced_by_column = zeros(Int, n_mz)

    @inbounds for col in 1:n_mz
        empirical_floor = floors[col]
        for row in 1:n_scans
            variance = v_out[row, col]
            limit = max(empirical_floor, k * sqrt(variance))
            if x_out[row, col] <= limit
                limit > 0 || throw(ArgumentError(
                    "censored replacement limit must be positive."))
                n_replaced += 1
                replaced_by_column[col] += 1
                if x_out[row, col] <= 0
                    n_nonpositive += 1
                else
                    n_below_limit_positive += 1
                end
                x_out[row, col] = limit / 2
                v_out[row, col] = max(variance, abs2(limit) / 12)
            end
        end
    end

    info = CensoredReplacementInfo(
        k,
        q,
        floor_rule,
        global_floor,
        floors,
        minpositives,
        n_replaced,
        n_replaced / length(x),
        n_nonpositive,
        n_below_limit_positive,
        replaced_by_column,
    )

    msm_out = copyreplacecensoredparent(vmsm, x_out, deepcopy(extras(vmsm)))
    vmsm_out = VarianceMassScanMatrix(msm_out, v_out, varianceunit(vmsm))

    returninfo ? (vmsm_out, info) : vmsm_out
end
