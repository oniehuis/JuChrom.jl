# ── mzchrom ───────────────────────────────────────────────────────────────────────────────

# helpers

# Works for mz grids that are Vector{<:Real} *or* Vector{<:Unitful.AbstractQuantity}.
function resolve_mz_index(mz_values::AbstractVector,
                          sel::Number,
                          tol;
                          mz_unit=nothing,
                          warning::Bool=true)::Int
    # empty grid -> soft fail
    if isempty(mz_values)
        warning && @warn "m/z grid is empty; returning 0."
        return 0
    end

    # normalize grid to Float64
    v = if first(mz_values) isa Unitful.AbstractQuantity
        mz_unit === nothing && throw(ArgumentError("m/z grid has units but `mz_unit` was not provided."))
        @. float(ustrip(uconvert(mz_unit, mz_values)))
    else
        @. float(mz_values)
    end

    # normalize selection to Float64 in same unit system
    sel_val = if sel isa Unitful.AbstractQuantity
        mz_unit === nothing && throw(ArgumentError("unitful selection requires a unitful m/z grid"))
        try
            float(ustrip(uconvert(mz_unit, sel)))
        catch e
            if e isa Unitful.DimensionError
                throw(ArgumentError("incompatible units for selection: $(sel) vs $(mz_unit)"))
            else
                rethrow()
            end
        end
    else
        float(sel)
    end

    # normalize tol to Float64 in same unit system
    tol_val = if tol isa Unitful.AbstractQuantity
        mz_unit === nothing && throw(
            ArgumentError("unitful tolerance requires a unitful m/z grid"))
        try
            float(ustrip(uconvert(mz_unit, tol)))
        catch e
            if e isa Unitful.DimensionError
                throw(ArgumentError("incompatible units for tol: $(tol) vs $(mz_unit)"))
            else
                rethrow()
            end
        end
    else
        float(tol)
    end

    Δ = @. abs(v - sel_val)
    i = findmin(Δ)[2]  # index of nearest
    if Δ[i] ≤ tol_val
        return i
    else
        warning && @warn ("No m/z match for $(sel) within tolerance $(tol) " 
                          * "$(mz_unit === nothing ? "" : string(mz_unit)); returning 0.")
        return 0
    end
end

# Resolve a 1-based index with optional soft failure
@inline function resolve_index(n::Int, k::Integer; warning::Bool)
    if 1 ≤ k ≤ n
        return k
    else
        warning && @warn "index $k out of range 1:$n; returning zero"
        return 0
    end
end

# Total Ion Current from a MassScan → scalar
@inline tic_value(scan::MassScan) = sum(intensities(scan))

# Build a ChromScan from (retention, scalar value), preserving attrs only.
# Uses the declared constructor:
#   ChromScan(::Unitful.AbstractQuantity, ::Real; attrs)
@inline function chrom_from_scalar(rt, y; attrs=NamedTuple())
    return ChromScan(rt, y; attrs=attrs)
end

# Pull one trace value out of a single MassScan by m/z selection (TIC-like scalar)
function extract_value_by_mz(scan::MassScan, sel::Number, tol::Number, warning::Bool)
    mzv   = mzvalues(scan)
    mz_u  = mzunit(scan)  # may be `nothing`

    # If m/z are unitless but tol is unitful, fail *here* (so tests see ArgumentError)
    if mz_u === nothing && tol isa Unitful.AbstractQuantity
        throw(ArgumentError("unitless series: `tol` must be unitless"))
    end

    j = resolve_mz_index(mzv, sel, tol; mz_unit=mz_u, warning=warning)
    if j == 0
        return zero(eltype(intensities(scan)))
    else
        return intensities(scan)[j]
    end
end

# ── mzchrom(series::MassScanSeries, ...) ──────────────────────────────────────────────────

"""
    mzchrom(series::MassScanSeries,
            selection::Union{Nothing, <:Number, AbstractVector{<:Number}}=nothing;
            tol::Number=3e-4,
            warning::Bool=true)

Extract a chromatogram scan series from a `MassScanSeries`.

- `selection`: either `nothing` (TIC) or a scalar/collection of m/z targets. When a
  scalar or vector is given, intensities at the selected m/z are summed per scan
  (XIC).  
- `tol`: the absolute tolerance used when matching `selection` m/z values to each scan’s
  discrete m/z grid. The nearest m/z value within `tol` is chosen for each target. `tol` 
  may be numeric or unitful. If the series’ m/z are unitless, `tol` must be unitless 
  (unitful `tol` errors). If the series’ m/z have a unit, numeric `tol` is interpreted in 
  that unit; unitful `tol` is converted to the m/z unit; incompatible units error.
- `warning`: whether to emit warnings if m/z values are not found or indices are out
  of bounds.

Returns a `ChromScanSeries` with retention time and intensity units preserved.
"""
function mzchrom(series::MassScanSeries, selection=nothing; tol=3e-4, warning::Bool=true)
    N = length(scans(series))
    N > 0 || throw(ArgumentError("Series has no scans"))

    compute_value = if selection === nothing
        scan -> tic_value(scan)
    else
        if selection isa AbstractVector
            scan -> sum(extract_value_by_mz(scan, s, tol, warning) for s in selection)
        else
            scan -> extract_value_by_mz(scan, selection, tol, warning)
        end
    end

    cs1 = chrom_from_scalar(retentions(series)[1], compute_value(scans(series)[1]); 
                            attrs=NamedTuple())
    chroms = Vector{typeof(cs1)}(undef, N)
    chroms[1] = cs1

    @inbounds for i in 2:N
        y = compute_value(scans(series)[i])
        chroms[i] = chrom_from_scalar(retentions(series)[i], y; attrs=NamedTuple())
    end

    ChromScanSeries(chroms;
        instrument=instrument(series),
        acquisition=acquisition(series),
        user=user(series),
        sample=sample(series),
        extras=extras(series))
end

# ── mzchrom(msm::MassScanMatrix, ...) ─────────────────────────────────────────────────────

"""
    mzchrom(msm::MassScanMatrix, 
            selection::Union{Nothing, <:Number, AbstractVector{<:Number}}=nothing;
            by::Symbol=:mz, tol::Number=3e-4, warning::Bool=true)

Extract a chromatogram scan series from a `MassScanMatrix`.

- `selection`: either `nothing` (TIC) or a scalar/collection of m/z or index targets.
  When given, intensities at the selected m/z or column indices are summed per scan.  
- `by`: specifies whether `selection` refers to m/z values (`:mz`, default) or to
  column indices (`:index`) in the matrix.  
- `tol`: the absolute tolerance used when matching `selection` m/z values to the discrete 
  m/z grid of the matrix. The nearest m/z value within `tol` is chosen for each target. 
  `tol` may be numeric or unitful. If the matrix m/z are unitless, `tol` must be unitless 
  (unitful `tol` errors). If the matrix m/z have a unit, numeric `tol` is interpreted in 
  that unit; unitful `tol` is converted to the m/z unit; incompatible units error.    
- `warning`: whether to emit warnings if m/z values are not found or indices are out
  of bounds.

Returns a `ChromScanSeries` with retention time and intensity units preserved.
"""
function mzchrom(msm::MassScanMatrix,
    selection::Union{Nothing, <:Number, AbstractVector{<:Number}}=nothing;
    by::Symbol=:mz, 
    tol::Number=3e-4, 
    warning::Bool=true)

    nrows = scancount(msm)
    nrows > 0 || throw(ArgumentError("Matrix has no rows (scans)"))

    if by === :index && selection isa Real && !(selection isa Integer)
        throw(ArgumentError(
            "`selection` must be an Integer or a collection of Integers when `by = :index`"))
    end
    if (by === :index && selection isa AbstractVector{<:Number} 
        && !(eltype(selection) <: Integer))

        throw(ArgumentError(
            "`selection` must be a collection of Integers or an Integer when `by = :index`"))
    end

    # Pre-resolve scalar selection to a single column index for speed
    scalar_j = nothing
    if selection !== nothing && !(selection isa AbstractVector)
        if by === :mz
            scalar_j = resolve_mz_index(mzvalues(msm), selection, tol; mz_unit=mzunit(msm), 
                                        warning=warning)
        else
            scalar_j = resolve_index(length(mzvalues(msm)), Int(selection); warning=warning)
        end
    end

    # Resolve vector selections to column indices (dedup + sort)
    resolved_sel = selection
    if selection isa AbstractVector
        if by === :mz
            resolve = sel -> resolve_mz_index(mzvalues(msm), sel, tol; mz_unit=mzunit(msm), 
                                              warning=warning)
            idx = [resolve(sel) for sel in selection]
            filter!(!=(0), idx)
            unique!(idx); sort!(idx)
            resolved_sel = idx
        elseif by === :index
            idx = [resolve_index(length(mzvalues(msm)), Int(k); warning=warning) 
                   for k in selection]
            filter!(!=(0), idx)
            unique!(idx); sort!(idx)
            resolved_sel = idx
        else
            throw(ArgumentError("Unsupported `by` = $by"))
        end
    end

    I   = intensities(msm)
    rts = retentions(msm)

    # Precompute TIC row sums when building a TIC
    row_sums = resolved_sel === nothing ? vec(sum(I, dims=2)) : nothing

    y1 = if resolved_sel === nothing
        row_sums[1]
    else
        if resolved_sel isa AbstractVector
            sum(@view I[1, resolved_sel])
        else
            j = scalar_j === nothing ? 0 : scalar_j
            j == 0 ? zero(eltype(I)) : I[1, j]
        end
    end

    cs1 = chrom_from_scalar(rts[1], y1; attrs=NamedTuple())
    chroms = Vector{typeof(cs1)}(undef, nrows)
    chroms[1] = cs1

    @inbounds for r in 2:nrows
        y = if resolved_sel === nothing
            row_sums[r]
        else
            if resolved_sel isa AbstractVector
                sum(@view I[r, resolved_sel])
            else
                j = scalar_j === nothing ? 0 : scalar_j
                j == 0 ? zero(eltype(I)) : I[r, j]
            end
        end
        chroms[r] = chrom_from_scalar(rts[r], y; attrs=NamedTuple())
    end

    ChromScanSeries(chroms;
        instrument=instrument(msm),
        acquisition=acquisition(msm),
        user=user(msm),
        sample=sample(msm),
        extras=extras(msm))
end
