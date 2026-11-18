# ── mscanmatrix ───────────────────────────────────────────────────────────────────────────

abstract type AbstractMatrixFormat end
struct DenseMatrix <: AbstractMatrixFormat end
struct SparseMatrix <: AbstractMatrixFormat end

const DENSE = DenseMatrix()
const SPARSE = SparseMatrix()

"""
    mscanmatrix(series::MassScanSeries, mzgrid;
                format=DENSE, target_level::Integer=1, threshold::Real=0.0,
                unmatched::Symbol=:error) -> MassScanMatrix

Convert scans at `target_level` to a retention × m/z matrix on the explicit `mzgrid`.

`mzgrid` must be a nonempty vector, range, or tuple of finite, positive, strictly
increasing m/z values. Spacing may be arbitrary. Every requested value gets a column,
including values absent from every scan. The grid is copied, never inferred from the
recorded ions, sorted, or deduplicated. Use `format=DENSE` (the default) or
`format=SPARSE` to choose storage.

Matching is exact, with no rounding, binning, or tolerance. For example, `57.1f0` and
`57.1` do not match. Apply [`binmzvalues`](@ref) beforehand when needed. With
`unmatched=:error`, any recorded m/z absent from the grid raises `ArgumentError`.
Use `unmatched=:drop` to deliberately select a subset. Only scans at `target_level`
are checked, but all their m/z values are checked **before** intensity filtering.
Only intensities strictly greater than `threshold` are retained.

Grid entries must be either all plain real numbers or all real-valued Unitful
quantities. Plain values are interpreted in the series' m/z unit, if any. Unitful
values are converted to that unit before exact matching; compatible units may differ
between entries. A unitful grid requires a unitful series. Mixed plain/unitful values
and incompatible units raise `ArgumentError`. The result preserves the series' m/z,
retention, and intensity units, numeric intensity type, and metadata.

An absent entry is represented by zero; this does not establish whether the ion was
acquired or whether its signal was suppressed during acquisition. An entirely zero
matrix is valid. No scans at the requested MS level raises `ArgumentError`.

# Examples
```jldoctest
julia> series = MassScanSeries([
           MassScan(1.0u"s", [57, 71], [10, 20]),
           MassScan(2.0u"s", [57, 85], [5, 30])]);

julia> matrix = mscanmatrix(series, 57:14:99);

julia> rawmzvalues(matrix)
4-element Vector{Int64}:
 57
 71
 85
 99

julia> rawintensities(matrix) == [10 20 0 0; 5 0 30 0]
true

julia> sparsematrix = mscanmatrix(series, (57, 71, 85, 99); format=SPARSE);

julia> rawintensities(sparsematrix) isa SparseMatrixCSC
true

julia> rawintensities(sparsematrix) == rawintensities(matrix)
true

julia> selected = mscanmatrix(series, (57, 71); unmatched=:drop, threshold=5);

julia> rawintensities(selected) == [10 20; 0 0]
true

julia> unitful = MassScanSeries([MassScan(1u"s", [57, 71]u"Th", [10, 20]u"pA")]);

julia> converted = mscanmatrix(unitful, (57000u"mTh", 71000u"mTh", 85000u"mTh"));

julia> mzvalues(converted) == [57, 71, 85]u"Th"
true

julia> intensities(converted) == [10 20 0]u"pA"
true
```
"""
function mscanmatrix(series::MassScanSeries, mzgrid::Union{AbstractVector, Tuple};
                     format::AbstractMatrixFormat=DENSE,
                     target_level::Integer=1,
                     threshold::Real=0.0,
                     unmatched::Symbol=:error)
    unmatched in (:error, :drop) || throw(ArgumentError(
        "unmatched must be :error or :drop."))
    filtered_series = levelscans(series, target_level)
    grid = mscanmatrix_grid(mzgrid, mzunit(filtered_series))
    mz_to_index = Dict(mz => i for (i, mz) in enumerate(grid))

    if unmatched === :error
        for (i, scan) in enumerate(scans(filtered_series)), mz in rawmzvalues(scan)
            haskey(mz_to_index, mz) || throw(ArgumentError(
                "Recorded m/z $mz at selected scan $i is absent from mzgrid; " *
                "use unmatched=:drop to intentionally exclude it."))
        end
    end

    values = create_matrix(format, filtered_series, mz_to_index,
        scancount(filtered_series), length(grid), threshold)
    MassScanMatrix(copy(rawretentions(filtered_series)), retentionunit(filtered_series),
        grid, mzunit(filtered_series), values, intensityunit(filtered_series);
        level=target_level, instrument=instrument(series), acquisition=acquisition(series),
        user=user(series), sample=sample(series), extras=extras(series))
end

function mscanmatrix_grid(mzgrid, mz_unit)
    isempty(mzgrid) && throw(ArgumentError("mzgrid must not be empty."))
    if all(x -> x isa Real, mzgrid)
        # Narrow abstract element types without rounding through numeric promotion.
        values = map(identity, collect(mzgrid))
    elseif all(x -> x isa AbstractQuantity{<:Real}, mzgrid)
        isnothing(mz_unit) && throw(ArgumentError(
            "A unitful mzgrid requires a series with an m/z unit."))
        all(x -> Unitful.dimension(x) == Unitful.dimension(mz_unit), mzgrid) ||
            throw(ArgumentError("mzgrid units must be compatible with the series' m/z unit."))
        values = map(x -> ustrip(mz_unit, x), collect(mzgrid))
    else
        throw(ArgumentError(
            "mzgrid must contain either only real numbers or only real-valued quantities."))
    end
    all(isfinite, values) || throw(ArgumentError("mzgrid values must be finite."))
    all(>(0), values) || throw(ArgumentError("mzgrid values must be positive."))
    all(i -> values[i - 1] < values[i], 2:length(values)) || throw(ArgumentError(
        "mzgrid values must be strictly increasing (no duplicates allowed)."))
    values
end

function create_matrix(::DenseMatrix, series, mz_to_index, n_scans, n_ions, threshold)
    values = zeros(eltype(rawintensities(first(scans(series)))), n_scans, n_ions)
    Threads.@threads for i in 1:n_scans
        scan = series[i]
        for (mz, intensity) in zip(rawmzvalues(scan), rawintensities(scan))
            j = get(mz_to_index, mz, 0)
            if j != 0 && intensity > threshold
                values[i, j] = intensity
            end
        end
    end
    values
end

function create_matrix(::SparseMatrix, series, mz_to_index, n_scans, n_ions, threshold)
    rows = Int[]
    cols = Int[]
    values = eltype(rawintensities(first(scans(series))))[]
    for (i, scan) in enumerate(scans(series))
        for (mz, intensity) in zip(rawmzvalues(scan), rawintensities(scan))
            j = get(mz_to_index, mz, 0)
            if j != 0 && intensity > threshold
                push!(rows, i)
                push!(cols, j)
                push!(values, intensity)
            end
        end
    end
    sparse(rows, cols, values, n_scans, n_ions)
end
