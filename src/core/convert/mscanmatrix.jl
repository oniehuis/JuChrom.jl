# ── mscanmatrix ───────────────────────────────────────────────────────────────────────────

# Define matrix format types for dispatch
abstract type AbstractMatrixFormat end
struct DenseMatrix <: AbstractMatrixFormat end
struct SparseMatrix <: AbstractMatrixFormat end

# Convenience constants for easy use
const DENSE = DenseMatrix()
const SPARSE = SparseMatrix()

"""
    mscanmatrix(series::MassScanSeries, format::AbstractMatrixFormat; 
                target_level::Integer=1, threshold::Real=0.0)

Convert a MassScanSeries to a MassScanMatrix with the specified matrix format.

# Arguments
- `series`: The input MassScanSeries
- `format`: Matrix format - `DENSE` or `SPARSE` (or `DenseMatrix()`, `SparseMatrix()`); 
  default is `DENSE`
- `target_level`: Scan level to filter for (default: 1)
- `threshold`: Minimum intensity threshold for inclusion (default: 0.0)

# Examples
```jldoctest
julia> scan1 = MassScan(1.0u"s", [101.2, 101.6, 102.1], [10, 20, 30])
       scan2 = MassScan(2.0u"s", [100.9, 101.2, 102.6], [20, 5, 30])
       series = MassScanSeries([scan1, scan2]);

julia> matrix_dense = mscanmatrix(series);  # Create dense matrix (default)

julia> intensities(matrix_dense) isa Matrix
true

julia> matrix_sparse = mscanmatrix(series, SPARSE);  # Create sparse matrix

julia> intensities(matrix_sparse) isa SparseArrays.SparseMatrixCSC
true

julia> intensities(matrix_dense) ≈ [0 10 20 30 0; 20 5 0 0 30]
true

julia> intensities(matrix_dense) ≈ intensities(matrix_sparse)
true

julia> retentions(matrix_dense) ≈ [1.0, 2.0]u"s"
true

julia> retentionunit(matrix_dense)
s

julia> mzvalues(matrix_dense) ≈ [100.9, 101.2, 101.6, 102.1, 102.6]
true

julia> mzunit(matrix_dense) === nothing
true

julia> matrix_filtered = mscanmatrix(series, threshold=7.5);  # Apply intensity threshold

julia> intensities(matrix_filtered) ≈ [0 10 20 30 0; 20 0 0 0 30]
true

julia> scan1 = MassScan(1.0u"s", [101.2, 101.6, 102.1], [10, 20, 30], level=1)
       scan2 = MassScan(2.0u"s", [100.9, 101.2, 102.6], [20, 5, 30], level=2)
       series = MassScanSeries([scan1, scan2]);

julia> matrix_level_2 = mscanmatrix(series, target_level=2);

julia> level(matrix_level_2)
2

julia> intensities(matrix_level_2) ≈ [20 5 30]
true

julia> series_annotated = MassScanSeries([scan1, scan2];
                          instrument=(detector="Orbitrap", manufacturer="Thermo"),
                          acquisition=(mode="FullScan", method="DDA", polarity="positive"),
                          user=(operator="Alice", project="QC-2025"),
                          sample=(ID="sample_001", matrix="plasma", prep="SPE"),
                          extras=Dict("injection_volume" => 5.0, "comment" => "QC run"));


julia> matrix_annotated = mscanmatrix(series_annotated);

julia> instrument(matrix_annotated)
(detector = "Orbitrap", manufacturer = "Thermo")

julia> acquisition(matrix_annotated)
(mode = "FullScan", method = "DDA", polarity = "positive")

julia> user(matrix_annotated)
(operator = "Alice", project = "QC-2025")

julia> sample(matrix_annotated)
(ID = "sample_001", matrix = "plasma", prep = "SPE")

julia> extras(matrix_annotated)
Dict{String, Any} with 2 entries:
  "injection_volume" => 5.0
  "comment"          => "QC run"
```
"""
function mscanmatrix(series::MassScanSeries, format::DenseMatrix;
                     target_level::Integer=1,
                     threshold::Real=0.0)

    filtered_series = levelscans(series, target_level)
    n_scans = scancount(filtered_series)
    n_scans > 0 || throw(ArgumentError("No scans found at target level $target_level"))

    uniquemzs = uniquemzvalues(filtered_series, target_level)
    n_ions = length(uniquemzs)

    mz_to_index = Dict{eltype(uniquemzs), Int}(ion => i for (i, ion) in enumerate(uniquemzs))

    # Build the dense intensity matrix per your existing logic
    intensity_matrix = create_matrix(format, filtered_series, mz_to_index, n_scans, n_ions, 
                                     threshold)

    # Normalize axes/matrix and propagate units (fallback when arrays are unitless)
    rt_vals, rt_unit = norm_axis(rawretentions(filtered_series), 
                                  retentionunit(filtered_series))
    mz_vals, mz_unit = norm_axis(uniquemzs, mzunit(filtered_series))
    I_vals, iu_unit = norm_matrix(intensity_matrix, intensityunit(filtered_series))

    MassScanMatrix(rt_vals, rt_unit, mz_vals, mz_unit, I_vals, iu_unit;
                   level=target_level,
                   instrument=instrument(series),
                   acquisition=acquisition(series),
                   user=user(series),
                   sample=sample(series),
                   extras=extras(series))
end


function mscanmatrix(series::MassScanSeries, format::SparseMatrix;
                     target_level::Integer=1,
                     threshold::Real=0.0)

    filtered_series = levelscans(series, target_level)
    n_scans = scancount(filtered_series)
    n_scans > 0 || throw(ArgumentError("No scans found at target level $target_level"))

    uniquemzs = uniquemzvalues(filtered_series, target_level)
    n_ions = length(uniquemzs)

    mz_to_index = Dict{eltype(uniquemzs), Int}(ion => i for (i, ion) in enumerate(uniquemzs))

    # Build the sparse intensity matrix per your existing logic
    intensity_matrix = create_matrix(format, filtered_series, mz_to_index, n_scans, n_ions, 
                                     threshold)

    # Normalize axes/matrix and propagate units (fallback when arrays are unitless)
    rt_vals, rt_unit = norm_axis(rawretentions(filtered_series), 
                                  retentionunit(filtered_series))
    mz_vals, mz_unit = norm_axis(uniquemzs, mzunit(filtered_series))
    I_vals, iu_unit = norm_matrix(intensity_matrix, intensityunit(filtered_series))

    MassScanMatrix(rt_vals, rt_unit, mz_vals, mz_unit, I_vals,  iu_unit;
                   level=target_level,
                   instrument=instrument(series),
                   acquisition=acquisition(series),
                   user=user(series),
                   sample=sample(series),
                   extras=extras(series))
end


# Default to dense matrix
function mscanmatrix(series::MassScanSeries; target_level::Integer=1, threshold::Real=0.0)
    mscanmatrix(series, DenseMatrix(); target_level=target_level, threshold=threshold)
end

# Format-specific matrix creation methods
function create_matrix(::DenseMatrix, filtered_series, mz_to_index, n_scans, n_ions, 
                       threshold)
    # Get the intensity type from first scan
    IntensityType = eltype(intensities(filtered_series[1]))
    
    # Pre-allocate matrix with correct type
    intensity_matrix = zeros(IntensityType, n_scans, n_ions)
    
    # Fill matrix
    Threads.@threads for i in 1:n_scans
        scan = @inbounds filtered_series[i]
        mz_vals = mzvalues(scan)
        intensity_vals = intensities(scan)
        
        @inbounds @simd for idx in eachindex(mz_vals, intensity_vals)
            mz = mz_vals[idx]
            intensity = intensity_vals[idx]
            if intensity > threshold
                j = mz_to_index[mz]
                intensity_matrix[i, j] = intensity
            end
        end
    end
    
    intensity_matrix
end

function create_matrix(::SparseMatrix, filtered_series, mz_to_index, n_scans, n_ions, 
                       threshold)
    # Get the intensity type from first scan
    IntensityType = eltype(intensities(filtered_series[1]))
    
    # Collect all non-zero entries
    rows = Int[]
    cols = Int[]
    vals = IntensityType[]
    
    # Estimate capacity to reduce allocations
    estimated_nnz = estimate_nonzeros(filtered_series, threshold)
    sizehint!(rows, estimated_nnz)
    sizehint!(cols, estimated_nnz)
    sizehint!(vals, estimated_nnz)
    
    # Sequential collection (avoid threading due to race conditions)
    @inbounds for i in 1:n_scans
        scan = filtered_series[i]
        mz_vals = mzvalues(scan)
        intensity_vals = intensities(scan)
        
        @simd for idx in eachindex(mz_vals, intensity_vals)
            mz = mz_vals[idx]
            intensity = intensity_vals[idx]
            if intensity > threshold
                j = mz_to_index[mz]
                push!(rows, i)
                push!(cols, j)
                push!(vals, intensity)
            end
        end
    end
    
    sparse(rows, cols, vals, n_scans, n_ions)
end

function estimate_nonzeros(filtered_series, threshold)
    # Simple heuristic: sample a few scans to estimate sparsity
    n_scans = scancount(filtered_series)
    sample_size = min(5, n_scans)
    
    total_nonzeros = 0
    @inbounds for i in 1:sample_size
        scan = filtered_series[i]
        total_nonzeros += count(x -> x > threshold, intensities(scan))
    end
    
    # Extrapolate to full dataset
    div(total_nonzeros * n_scans, sample_size)
end

# Normalize an axis vector to plain numbers + capture unit.
# Falls back to `fallback_unit` when the vector itself is unitless.
@inline function norm_axis(v::AbstractVector, fallback_unit::Union{Nothing, Unitful.Units})
    isempty(v) && return Float64[], fallback_unit
    x1 = first(v)
    if x1 isa Unitful.AbstractQuantity
        u = Unitful.unit(x1)
        return [float(ustrip(uconvert(u, x))) for x in v], u
    else
        return collect(float.(v)), fallback_unit
    end
end

# Normalize an intensity matrix + capture unit.
# If matrix elements are Quantities, strip to Float64 and capture their unit.
# If matrix is unitless, preserve eltype (e.g., Int, Float64) and use fallback_unit.
function norm_matrix(M::AbstractMatrix, fallback_unit::Union{Nothing, Unitful.Units})
    if isempty(M)
        return similar(M, eltype(M), size(M)), fallback_unit
    end
    x1 = M[begin]
    if x1 isa Unitful.AbstractQuantity
        u = Unitful.unit(x1)
        A = similar(M, Float64)
        @inbounds for j in axes(M,2), i in axes(M,1)
            A[i,j] = float(ustrip(uconvert(u, M[i,j])))
        end
        return A, u
    else
        # unitless: keep original eltype (don’t force Float64)
        if eltype(M) <: Number
            return copy(M), fallback_unit
        else
            # extremely rare: non-Number eltype; just convert to Float64
            A = similar(M, Float64)
            @inbounds for j in axes(M,2), i in axes(M,1)
                A[i,j] = float(M[i,j])
            end
            return A, fallback_unit
        end
    end
end
