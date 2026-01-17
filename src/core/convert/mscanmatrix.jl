# ── mscanmatrix ───────────────────────────────────────────────────────────────────────────

# Define matrix format types for dispatch
abstract type AbstractMatrixFormat end
struct DenseMatrix <: AbstractMatrixFormat end
struct SparseMatrix <: AbstractMatrixFormat end

# Convenience constants for easy use
const DENSE = DenseMatrix()
const SPARSE = SparseMatrix()

"""
    mscanmatrix(
        series::MassScanSeries,
        format::AbstractMatrixFormat; 
        target_level::Integer=1,
        threshold::Real=0.0
    ) -> MassScanMatrix

Convert a MassScanSeries to a MassScanMatrix with the specified matrix format.

`series` supplies the input scans, and `format` selects `DENSE` or `SPARSE` storage (or
`DenseMatrix()` / `SparseMatrix()`), with `DENSE` as the default. `DENSE` stores the full
retention × m/z grid as a standard dense matrix, while `SPARSE` stores only the nonzero
entries in a sparse matrix to reduce memory for sparse spectra. `target_level` filters
scans by MS level, and `threshold` excludes intensities below the specified minimum.

See also
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix),
[`MassScanSeries`](@ref JuChrom.MassScanSeries).

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

    # Use raw values so units remain in metadata (no Unitful values in the grid).
    uniquemzs = sort(unique(vcat(rawmzvalues.(scans(filtered_series))...)))
    n_ions = length(uniquemzs)

    mz_to_index = Dict{eltype(uniquemzs), Int}(ion => i for (i, ion) in enumerate(uniquemzs))

    # Build the dense intensity matrix per your existing logic
    intensity_matrix = create_matrix(format, filtered_series, mz_to_index, n_scans, n_ions, 
                                     threshold)

    # Preserve original eltypes; units are stored separately in metadata.
    rt_vals = copy(rawretentions(filtered_series))
    rt_unit = retentionunit(filtered_series)
    mz_vals = copy(uniquemzs)
    mz_unit = mzunit(filtered_series)
    I_vals = intensity_matrix
    iu_unit = intensityunit(filtered_series)

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

    # Use raw values so units remain in metadata (no Unitful values in the grid).
    uniquemzs = sort(unique(vcat(rawmzvalues.(scans(filtered_series))...)))
    n_ions = length(uniquemzs)

    mz_to_index = Dict{eltype(uniquemzs), Int}(ion => i for (i, ion) in enumerate(uniquemzs))

    # Build the sparse intensity matrix per your existing logic
    intensity_matrix = create_matrix(format, filtered_series, mz_to_index, n_scans, n_ions, 
                                     threshold)

    # Preserve original eltypes; units are stored separately in metadata.
    rt_vals = copy(rawretentions(filtered_series))
    rt_unit = retentionunit(filtered_series)
    mz_vals = copy(uniquemzs)
    mz_unit = mzunit(filtered_series)
    I_vals = intensity_matrix
    iu_unit = intensityunit(filtered_series)

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
        mz_vals = rawmzvalues(scan)
        intensity_vals = rawintensities(scan)
        
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
        mz_vals = rawmzvalues(scan)
        intensity_vals = rawintensities(scan)
        
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
