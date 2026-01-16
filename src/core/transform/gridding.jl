# ── densestgrid ───────────────────────────────────────────────────────────────────────────

"""
    densestgrid(
        datavectors::AbstractVector{<:AbstractVector{<:Real}}; 
        minwidth::Union{Nothing,Real}=nothing, 
        maxwidth::Union{Nothing,Real}=nothing,
        tolerance::Real=1e-8,
        coarse_inflation::Real=1.001,
        primary_refine_iters::Int=35,
        secondary_refine_iters::Int=20,
        max_anchors_basic::Int=256,
        max_anchors_enriched::Int=1024
    )

Return `(edges, w)` where `edges` is a regular grid and `w` is the near‑minimal bin width 
that ensures every bin in `[global_min, global_max]` contains at least one observation from 
each input dataset in `datavectors`. The grid uses left‑closed, right‑open bins 
`[eᵢ, eᵢ₊₁)`, with the final bin including the right edge within `tolerance`.

The method begins from a rigorous lower bound on `w` and expands the candidate width by
`coarse_inflation` until a basic anchor phase yields full occupancy, then performs two 
rounds of binary refinement (`primary_refine_iters` with a basic anchor set and 
`secondary_refine_iters` with an enriched set) to tighten the admissible interval for `w`. 
The basic anchor set uses residue 0 together with residues `(global_min - x) % w` from 
early points of each dataset, capped at `max_anchors_basic` with de‑duplication at 
floating‑point resolution (unit in the last place); the enriched set adds residues from
terminal points and midpoints of large internal gaps up to `max_anchors_enriched`.

The function throws `ArgumentError` if any dataset is empty, the overlap is empty, or no 
admissible width is found within `maxwidth`. For a width closer to the theoretical minimum, 
decrease `coarse_inflation` and/or increase the refinement iteration counts and anchor caps, 
at the cost of additional runtime.
"""
function densestgrid(
    datavectors::AbstractVector{<:AbstractVector{<:Real}};
    minwidth::Union{Nothing,Real}=nothing,
    maxwidth::Union{Nothing,Real}=nothing,
    tolerance::Real=1e-8,
    coarse_inflation::Real=1.001,
    primary_refine_iters::Int=35,
    secondary_refine_iters::Int=20,
    max_anchors_basic::Int=256,
    max_anchors_enriched::Int=1024
)
    # Validation
    isempty(datavectors) && throw(ArgumentError("datavectors cannot be empty"))
    @inbounds for (i,v) in enumerate(datavectors)
        isempty(v) && throw(ArgumentError("datavectors[$i] cannot be empty"))
    end

    # Identify the overlap every dataset shares; bins outside this region would miss
    # at least one dataset regardless of width/phase.
    global_min = maximum(minimum(v) for v in datavectors)
    global_max = minimum(maximum(v) for v in datavectors)
    global_max ≤ global_min && throw(
        ArgumentError("No overlapping range across all datavectors"))

    # Strip values outside the common overlap so subsequent bin checks only consider
    # points that can actually contribute to the shared grid.
    datasets = [sort(filter(x -> global_min ≤ x ≤ global_max, v)) for v in datavectors]
    any(isempty, datasets) && throw(
        ArgumentError("A dataset lost all points inside common overlap"))

    # Positive internal gaps (keep only > 0)
    # These gaps deliver the tightest theoretical lower bound on the bin width.
    pos_gaps = [g for r in datasets for g in (length(r) ≥ 2 ? diff(r) : Float64[]) if g > 0.0]
    isempty(pos_gaps) && throw(
        ArgumentError("All datasets singleton or duplicates; cannot form bins"))
    computed_minwidth = minimum(pos_gaps)

    # Respect optional width constraints provided by the caller.
    effective_minwidth = isnothing(minwidth) ? computed_minwidth : minwidth
    effective_minwidth ≤ 0 && throw(ArgumentError("minwidth must be positive"))
    effective_maxwidth = isnothing(maxwidth) ? (global_max - global_min) : maxwidth
    effective_minwidth > effective_maxwidth &&
        throw(ArgumentError("minwidth ($effective_minwidth) > maxwidth ($effective_maxwidth)"))

    # Max gap including overlap edges (captures situations where a dataset starts later
    # than `global_min` or ends earlier than `global_max`).
    max_internal_gap = maximum((maximum(diff(r)) for r in datasets if length(r) ≥ 2); init=0.0)
    max_edge_gap = maximum(max(first(r) - global_min, global_max - last(r)) for r in datasets)
    max_gap = max(max_internal_gap, max_edge_gap)

    # Lower bound (strict > half any gap)
    # Any bin narrower than half the maximum gap would force a hole for some dataset.
    lb = max(effective_minwidth, (max_gap / 2) * (1 + 10 * eps(1.0)))
    lb > effective_maxwidth && throw(ArgumentError(
        "Lower bound width $lb exceeds maxwidth $effective_maxwidth"))

    # Quick rejector: widths that cannot close the largest gap no matter the phase.
    width_immediately_invalid(w) = (2w - max_gap) ≤ 0

    # Anchor generation (refactored earlier)
    # Build candidate phase offsets that we test for a given width. The enriched mode
    # incorporates residues from leading/trailing points and wide internal gaps.
    function anchors_for(w; enriched::Bool=false)
        capacity = enriched ? max_anchors_enriched : max_anchors_basic
        ulp_w = 2 * eps(w)
        A = Float64[0.0]

        # Attempt to insert a single residue unless it already exists within an ULP.
        function add_or_full!(x)::Bool
            @inbounds for z in A
                if abs(z - x) ≤ ulp_w
                    return length(A) ≥ capacity
                end
            end
            push!(A, x)
            length(A) ≥ capacity
        end

        # Batch insert helper; stops early if capacity cap is reached.
        function add_seq!(iter)::Bool
            @inbounds for x in iter
                add_or_full!(x) && return true
            end
            length(A) ≥ capacity
        end

        @inbounds for r in datasets
            limit = enriched ? length(r) : min(length(r), 8)
            add_seq!((mod(global_min - r[k], w) for k in 1:limit)) && break
        end
        if enriched
            @inbounds for r in datasets
                add_or_full!(mod(global_min - last(r), w)) && break
            end
            @inbounds for r in datasets
                if length(r) ≥ 2
                    diffs = diff(r)
                    kmax = min(5, length(diffs))
                    if kmax > 0
                        idxs = partialsortperm(diffs, rev=true, 1:kmax)
                        add_seq!((mod(global_min - (r[idx] + r[idx+1]) / 2, w) for idx in idxs)) && break
                    end
                end
                length(A) ≥ capacity && break
            end
        end

        sort!(A)
        A
    end

    # Helpers
    # Translate a width/offset pair into grid edges covering [global_min, global_max].
    function compute_bins(w, off)
        first_edge = global_min - off
        n_bins = ceil(Int, (global_max - first_edge)/w)
        needs_extra = first_edge + n_bins * w < global_max - 2 * eps(global_max)
        n_bins += needs_extra ? 1 : 0
        first_edge, n_bins
    end

    # Check whether every bin contains at least one observation from each dataset.
    function occupancy_ok(w, off)::Bool
        first_edge, n_bins = compute_bins(w, off)
        nb = n_bins
        @inbounds for r in datasets
            @inbounds for j in 0:nb-1
                a = first_edge + j*w
                b = a + w
                lastbin = (j == nb - 1)  # final bin gets tolerance inclusion on the right edge
                ia = searchsortedfirst(r, a)
                hit = false
                @inbounds for k in ia:length(r)
                    x = r[k]
                    if (!lastbin && x ≥ b) || (lastbin && x > b + tolerance)
                        break
                    end
                    if x ≥ a && ((lastbin && x ≤ b + tolerance) || x < b)
                        hit = true
                        break
                    end
                end
                hit || return false
            end
        end
        true
    end

    function build_edges(w, off)
        first_edge, n_bins = compute_bins(w, off)
        collect(first_edge + i * w for i in 0:n_bins)
    end

    # Try the available anchors for a width and return the first successful grid.
    function edges_for_width(w; enriched::Bool=false)
        acands = anchors_for(w; enriched=enriched)
        @inbounds for off in acands
            if occupancy_ok(w, off)
                return true, build_edges(w, off)
            end
        end
        false, nothing
    end

    # Binary refinement between failing/passing widths using either anchor set.
    function refine_stage(w_low, w_high, iters; enriched::Bool=false)
        edges_found = nothing
        for _ in 1:iters
            mid = (w_low + w_high) / 2
            if mid == w_low || mid == w_high
                break
            end
            ok, edges = edges_for_width(mid; enriched=enriched)
            if ok
                w_high = mid
                edges_found = edges
            else
                w_low = mid
            end
        end
        w_low, w_high, edges_found
    end

    # Exponentially inflate the width until a viable anchor is found.
    function coarse_stage(lb, wmax, inflation)
        edges_found = nothing
        w_succ = nothing
        widths = Iterators.takewhile(
            x -> x ≤ wmax,
            Iterators.flatten((
                (lb,),
                Iterators.accumulate(*, Iterators.repeated(inflation); init=lb),
            )),
        )
        w_fail = lb
        for w in widths
            if !width_immediately_invalid(w)
                ok, edges = edges_for_width(w; enriched=false)
                if ok
                    edges_found = edges
                    w_succ = w
                    break
                end
            end
            w_fail = w
        end
        w_succ === nothing && throw(
            ArgumentError("Could not find occupancy width up to maxwidth"))
        w_fail, w_succ, edges_found
    end

    # Stages
    # Stage 1: inflate width until some anchor yields full occupancy.
    w_fail, w_succ, edges_found = coarse_stage(lb, effective_maxwidth, coarse_inflation)
    w_low = w_fail
    w_high = w_succ
    # Stage 2: binary search with the basic anchor set.
    w_low, w_high, edgesB = refine_stage(w_low, w_high, primary_refine_iters; enriched=false)
    if edgesB !== nothing
        edges_found = edgesB
    end
    # Stage 3: binary search with the enriched anchor set for a final tightening.
    w_low, w_high, edgesC = refine_stage(w_low, w_high, secondary_refine_iters; enriched=true)
    if edgesC !== nothing
        edges_found = edgesC
    end

    if edges_found === nothing
        ok, edges_final = edges_for_width(w_high; enriched=true)
        ok || throw(ArgumentError("Failed to recover grid edges at final width"))
        edges_found = edges_final
    end

    edges_found, w_high
end


"""
    densestgrid(
        datavectors::AbstractVector{<:AbstractVector{<:Unitful.AbstractQuantity}};
        minwidth::Union{Nothing,Unitful.AbstractQuantity}=nothing,
        maxwidth::Union{Nothing,Unitful.AbstractQuantity}=nothing,
        tolerance::Union{Nothing,Unitful.AbstractQuantity}=nothing,
        coarse_inflation::Real=1.001,
        primary_refine_iters::Int=35,
        secondary_refine_iters::Int=20,
        max_anchors_basic::Int=256,
        max_anchors_enriched::Int=1024
    )

Unit‑aware wrapper of the numeric `densestgrid`.

Computes a near‑minimal bin width `w` and regular grid `edges` such that every
bin contains at least one observation from each dataset (three‑stage procedure:
coarse inflation, primary refinement, secondary refinement with enriched
anchors).

Requirements:
* All elements of all inner vectors must be `Unitful.AbstractQuantity` with the **same
  physical dimension** (units themselves may differ).
* `minwidth`, `maxwidth`, `tolerance` when provided must be quantities of that
  dimension (`nothing` allowed).
* If `tolerance === nothing`, it defaults to `1e-8 * ref_unit`.

Keyword arguments (forwarded to numeric core):
`coarse_inflation`, `primary_refine_iters`, `secondary_refine_iters`,
`max_anchors_basic`, `max_anchors_enriched`.

Returns `(edges, w)` with the original unit attached.
"""
function densestgrid(
    datavectors::AbstractVector{<:AbstractVector{<:Unitful.AbstractQuantity}};
    minwidth::Union{Nothing,Unitful.AbstractQuantity}=nothing,
    maxwidth::Union{Nothing,Unitful.AbstractQuantity}=nothing,
    tolerance::Union{Nothing,Unitful.AbstractQuantity}=nothing,
    coarse_inflation::Real=1.001,
    primary_refine_iters::Int=35,
    secondary_refine_iters::Int=20,
    max_anchors_basic::Int=256,
    max_anchors_enriched::Int=1024
)
    isempty(datavectors) && throw(ArgumentError("datavectors cannot be empty"))
    firstvec = datavectors[1]
    isempty(firstvec) && throw(ArgumentError("datavectors[1] cannot be empty"))

    # Use the first datapoint as a reference to enforce consistent physical dimension.
    ref_unit = unit(firstvec[1])
    ref_dim  = Unitful.dimension(ref_unit)

    @inbounds for (i, vec) in enumerate(datavectors)
        isempty(vec) && throw(ArgumentError("datavectors[$i] cannot be empty"))
        @inbounds for x in vec
            isa(x, Unitful.AbstractQuantity) ||
                throw(ArgumentError("datapoint without unit in datavectors[$i]"))
            Unitful.dimension(unit(x)) == ref_dim ||
                throw(ArgumentError("mixed physical dimensions in datavectors"))
        end
    end

    # Convert all values to the reference unit and strip them to plain Float64 so the
    # numeric core can operate without carrying AbstractQuantity types.
    strip = x -> Unitful.ustrip(uconvert(ref_unit, x))
    numeric_datasets = [strip.(v) for v in datavectors]

    # Helper that strips keyword quantities (if provided) to the reference unit.
    to_num(q, default) =
        q === nothing ? default : Unitful.ustrip(uconvert(ref_unit, q))

    min_num = to_num(minwidth, nothing)
    max_num = to_num(maxwidth, nothing)
    tol_num = tolerance === nothing ? 1e-8 :
              Unitful.ustrip(uconvert(ref_unit, tolerance))

    edges_num, width_num = densestgrid(
        numeric_datasets;
        minwidth=min_num,
        maxwidth=max_num,
        tolerance=tol_num,
        coarse_inflation=coarse_inflation,
        primary_refine_iters=primary_refine_iters,
        secondary_refine_iters=secondary_refine_iters,
        max_anchors_basic=max_anchors_basic,
        max_anchors_enriched=max_anchors_enriched
    )

    # Re-attach the reference unit to the edges/width before returning.
    edges_num .* ref_unit, width_num * ref_unit
end

"""
    densestgrid(msmatrices::AbstractDict{<:Any, <:AbstractMassScanMatrix};
                minwidth::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
                maxwidth::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
                tolerance::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
                coarse_inflation=1.001,
                primary_refine_iters=35,
                secondary_refine_iters=20,
                max_anchors_basic=256,
                max_anchors_enriched=1024)

Dictionary wrapper for [`densestgrid`](@ref JuChrom.densestgrid(::AbstractVector{<:AbstractMassScanMatrix})).
The keys are ignored;
the values are treated as the dataset list, and all keyword arguments and behavior
match the vector method.

Use this when your matrices are already keyed by sample or run identifiers. The
returned `(edges, w)` and error conditions are identical to the vector-based API.
"""
function densestgrid(
    msmatrices::AbstractDict{<:Any, <:AbstractMassScanMatrix};
    minwidth::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
    maxwidth::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
    tolerance::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
    coarse_inflation::Real=1.001,
    primary_refine_iters::Int=35,
    secondary_refine_iters::Int=20,
    max_anchors_basic::Int=256,
    max_anchors_enriched::Int=1024
)
    isempty(msmatrices) &&
        throw(ArgumentError("msmatrices dictionary cannot be empty"))

    # Pull out the retention vectors once so we can inspect units and forward them.
    retention_vectors = retentions.(values(msmatrices))

    # Track a sample point to infer whether the set is unitful or unitless.
    sample = nothing
    for v in retention_vectors
        isempty(v) && throw(ArgumentError("at least one MassScanMatrix has no data"))
        sample === nothing && (sample = first(v))
    end

    if isa(sample, Unitful.AbstractQuantity)
        # Unitful path: enforce consistent dimensions and adapt keyword quantities.
        ref_unit = unit(sample)
        ref_dim  = Unitful.dimension(ref_unit)
        for vec in retention_vectors, x in vec
            isa(x, Unitful.AbstractQuantity) ||
                throw(ArgumentError("mixture of unitless and unitful retention data"))
            Unitful.dimension(unit(x)) == ref_dim ||
                throw(ArgumentError("mixed physical dimensions in retention data"))
        end
        kw_ok(q) = q === nothing || isa(q, Unitful.AbstractQuantity) ||
                   throw(ArgumentError("keyword must be AbstractQuantity or nothing for unitful data"))
        kw_ok(minwidth); kw_ok(maxwidth)
        tol_qty = tolerance === nothing ? 1e-8 * ref_unit :
                  isa(tolerance, Unitful.AbstractQuantity) ? tolerance :
                  throw(ArgumentError("tolerance must be AbstractQuantity or nothing for unitful data"))

        # Delegate to the unit-aware vector overload with minimally adjusted keywords.
        return densestgrid(
            retention_vectors;
            minwidth=minwidth,
            maxwidth=maxwidth,
            tolerance=tol_qty,
            coarse_inflation=coarse_inflation,
            primary_refine_iters=primary_refine_iters,
            secondary_refine_iters=secondary_refine_iters,
            max_anchors_basic=max_anchors_basic,
            max_anchors_enriched=max_anchors_enriched
        )
    end

    # Unitless path: all entries must be plain reals and keywords remain numeric.
    for vec in retention_vectors, x in vec
        (isa(x, Real) && !isa(x, Unitful.AbstractQuantity)) ||
            throw(ArgumentError("mixture of unitful and unitless retention data"))
    end
    kw_real(q) = q === nothing || (isa(q, Real) && !isa(q, Unitful.AbstractQuantity)) ||
                 throw(ArgumentError("keyword must be Real or nothing for unitless data"))
    kw_real(minwidth)
    kw_real(maxwidth)
    tol_num = tolerance === nothing ? 1e-8 :
              (isa(tolerance, Real) && !isa(tolerance, Unitful.AbstractQuantity)) ?
                  tolerance :
                  throw(ArgumentError("tolerance must be Real or nothing for unitless data"))

    # Unitless vectors can go straight into the numeric solver.
    densestgrid(
        retention_vectors;
        minwidth=minwidth,
        maxwidth=maxwidth,
        tolerance=tol_num,
        coarse_inflation=coarse_inflation,
        primary_refine_iters=primary_refine_iters,
        secondary_refine_iters=secondary_refine_iters,
        max_anchors_basic=max_anchors_basic,
        max_anchors_enriched=max_anchors_enriched
    )
end

"""
    densestgrid(msmatrices::AbstractVector{<:AbstractMassScanMatrix};
                minwidth::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
                maxwidth::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
                tolerance::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
                coarse_inflation=1.001,
                primary_refine_iters=35,
                secondary_refine_iters=20,
                max_anchors_basic=256,
                max_anchors_enriched=1024)

Determine a near-minimal bin width `w` and matching retention grid `edges` so that
every bin in the shared overlap contains at least one observation from each matrix.

The method extracts retention vectors from each matrix, trims them to their common
overlap, and searches for the smallest bin width that admits a fully occupied grid.
The search uses a coarse inflation step followed by two rounds of binary refinement
with a capped set of anchor phases (`max_anchors_basic` and `max_anchors_enriched`).
The inputs are a non-empty vector of `MassScanMatrix` values and optional bounds
(`minwidth`, `maxwidth`) plus a right-edge inclusion `tolerance`. Retentions must be
either unitless or `Unitful.AbstractQuantity` values sharing the same physical
dimension. If retentions are unitful, `minwidth`, `maxwidth`, and `tolerance` must be
`nothing` or dimensionally compatible quantities and the default `tolerance` is
`1e-8 * unit`; if retentions are unitless, these keywords must be `nothing` or plain
reals and the default `tolerance` is `1e-8`. The function returns `(edges, w)` (with
units attached when appropriate) and throws `ArgumentError` if inputs are empty, units
are inconsistent, the overlap is empty, or no admissible width is found up to
`maxwidth`. See also [`densestgrid`](@ref JuChrom.densestgrid(::AbstractDict{<:Any,<:AbstractMassScanMatrix})).
"""
function densestgrid(
    msmatrices::AbstractVector{<:AbstractMassScanMatrix};
    minwidth::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
    maxwidth::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
    tolerance::Union{Nothing,Real,Unitful.AbstractQuantity}=nothing,
    coarse_inflation::Real=1.001,
    primary_refine_iters::Int=35,
    secondary_refine_iters::Int=20,
    max_anchors_basic::Int=256,
    max_anchors_enriched::Int=1024
)
    isempty(msmatrices) &&
        throw(ArgumentError("msmatrices vector cannot be empty"))

    # Work with the raw retention vectors directly so we can validate and pass them
    # into the numeric/unitful cores without repeatedly calling accessors.
    retention_vectors = retentions.(msmatrices)
    isempty(first(retention_vectors)) &&
        throw(ArgumentError("at least one MassScanMatrix has no data"))

    # Inspect the first available datapoint to determine whether we're in the unitful
    # or unitless branch.
    sample = first(first(retention_vectors))

    if isa(sample, Unitful.AbstractQuantity)
        # Unitful retentions: ensure all datasets share the same physical dimension and
        # adapt keyword arguments/tolerance to the reference unit.
        ref_unit = unit(sample)
        ref_dim  = Unitful.dimension(ref_unit)
        for v in retention_vectors, x in v
            isa(x, Unitful.AbstractQuantity) ||
                throw(ArgumentError("mixture of unitless and unitful retention data"))
            Unitful.dimension(unit(x)) == ref_dim ||
                throw(ArgumentError("mixed physical dimensions in retention data"))
        end
        kw = q -> q === nothing || isa(q, Unitful.AbstractQuantity) ||
                  throw(ArgumentError("keyword must be AbstractQuantity or nothing"))
        kw(minwidth)
        kw(maxwidth)
        tol_qty = tolerance === nothing ? 1e-8 * ref_unit :
                  isa(tolerance, Unitful.AbstractQuantity) ? tolerance :
                  throw(ArgumentError("tolerance must be AbstractQuantity or nothing"))
        # Send everything to the unit-aware vector overload, reusing the extracted datasets.
        return densestgrid(
            retention_vectors;
            minwidth=minwidth,
            maxwidth=maxwidth,
            tolerance=tol_qty,
            coarse_inflation=coarse_inflation,
            primary_refine_iters=primary_refine_iters,
            secondary_refine_iters=secondary_refine_iters,
            max_anchors_basic=max_anchors_basic,
            max_anchors_enriched=max_anchors_enriched
        )
    end

    # Unitless retentions: reject any quantities and leave keywords as numeric values.
    for v in retention_vectors, x in v
        (isa(x, Real) && !isa(x, Unitful.AbstractQuantity)) ||
            throw(ArgumentError("mixture of unitful and unitless retention data"))
    end
    kw = q -> q === nothing || (isa(q, Real) && !isa(q, Unitful.AbstractQuantity)) ||
                 throw(ArgumentError("keyword must be Real or nothing"))
    kw(minwidth)
    kw(maxwidth)
    tol_num = tolerance === nothing ? 1e-8 :
              (isa(tolerance, Real) && !isa(tolerance, Unitful.AbstractQuantity)) ?
                  tolerance :
                  throw(ArgumentError("tolerance must be Real or nothing"))
    # Unitless retention vectors can be forwarded directly to the numeric solver.
    densestgrid(
        retention_vectors;
        minwidth=minwidth,
        maxwidth=maxwidth,
        tolerance=tol_num,
        coarse_inflation=coarse_inflation,
        primary_refine_iters=primary_refine_iters,
        secondary_refine_iters=secondary_refine_iters,
        max_anchors_basic=max_anchors_basic,
        max_anchors_enriched=max_anchors_enriched
    )
end
