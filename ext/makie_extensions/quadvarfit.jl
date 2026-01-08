using Printf
using Unitful

const DEFAULT_FIGSIZE = (1600, 800)

"""
    plot(q::JuChrom.QuadVarFit; mzi=1, batch=1,
         mode=:align, scan_range=nothing, varpred_floor=0.0,
         figsize=DEFAULT_FIGSIZE, legend=true)

Plot a single batch view from a self-contained `QuadVarFit`.

- `mode = :align`  → aligned replicates ((Y - o)/g), plus **one** CI band around `s` and the `s` line
- `mode = :offsets` → offsets per replicate

`mzi` indexes into the selected m/z (`q.mz_idx`). `batch` selects which batch to plot.
Works whether `q.observed` stores the **full** m/z grid or only the **selected** m/z columns.
"""
function Makie.plot(q::JuChrom.QuadVarFit;
              mzi::Integer = 1,
              batch::Integer = 1,
              mode::Symbol = :align,   # :align (with CI) or :offsets
              scan_range::Union{Nothing,AbstractUnitRange{<:Integer}} = nothing,
              varpred_floor::Real = 0.0,   # only used for :align
              figsize::Tuple{<:Integer,<:Integer} = DEFAULT_FIGSIZE,
              legend::Bool = true)

    mode in (:align, :offsets) ||
        throw(ArgumentError("mode must be :align or :offsets"))
    1 ≤ mzi ≤ length(q.mz_idx) ||
        throw(ArgumentError("mzi out of range: 1:$(length(q.mz_idx))"))
    1 ≤ batch ≤ q.batchcount ||
        throw(ArgumentError("batch out of range: 1:$(q.batchcount)"))

    hasfield(typeof(q), :observed) ||
        throw(ArgumentError("QuadVarFit has no field `observed`. Store observed matrices in the fit to plot alignments."))

    # ---- bookkeeping ----
    # Detect whether observed matrices are full grid or reduced-to-selection
    obs_batch = q.observed[batch]
    n_reps_b  = q.n_reps_per_batch[batch]
    n_scans_b = q.n_scans_per_batch[batch]
    length(obs_batch) == n_reps_b ||
        throw(ArgumentError("batch $batch: observed has $(length(obs_batch)) reps; expected $n_reps_b"))

    obs_ncols = size(obs_batch[1], 2)
    full_ncols = length(q.mz_ref)
    sel_ncols  = length(q.mz_values)

    observed_is_full = obs_ncols == full_ncols
    observed_is_sel  = obs_ncols == sel_ncols

    (observed_is_full || observed_is_sel) ||
        throw(ArgumentError("observed[$batch][1] has $obs_ncols m/z; expected $full_ncols (full grid) or $sel_ncols (selection)."))

    # Column index into observed matrices
    col_idx  = observed_is_full ? q.mz_idx[mzi] : mzi
    mzval    = q.mz_values[mzi]
    mzunit   = q.mz_unit
    p        = q.params[mzi]
    acf_val  = q.acf[mzi]

    # Shapes per replicate
    for r in 1:n_reps_b
        size(obs_batch[r], 1) == n_scans_b ||
            throw(ArgumentError("observed[$batch][$r] has $(size(obs_batch[r],1)) scans; expected $n_scans_b"))
        size(obs_batch[r], 2) == obs_ncols ||
            throw(ArgumentError("observed[$batch][$r] has inconsistent column count"))
    end

    s_b = q.signal[mzi][batch]   # n_scans_b
    o_b = q.offsets[mzi][batch]  # n_scans_b × n_reps_b
    g_b = q.gains[mzi][batch]    # n_reps_b

    fig = Makie.Figure(; size=figsize)
    ylabel = mode === :offsets ? "offset" : "aligned intensity ((Y−o)/g)"
    mzunit_str = isnothing(mzunit) ? "" : " " * string(mzunit)
    title  = @sprintf("m/z %.3f%s — batch %d; σ₀²=%g  ϕ=%g  κ=%.4f  acf=%+.2f",
                      mzval, mzunit_str, batch,
                      round(p.σ₀², sigdigits=3),
                      round(p.ϕ,  sigdigits=4),
                      p.κ,
                      acf_val)

    ax = Makie.Axis(fig[1, 1];
                    xlabel = "scan",
                    ylabel = ylabel,
                    title  = title)

    palette = Makie.wong_colors()

    if mode === :align
        # One CI band around the shared signal + the shared signal + aligned replicates
        x     = 1:n_scans_b
        sd_b  = sqrt.(JuChrom.varpred(s_b, p; varfloor = varpred_floor))
        lower = clamp.(s_b .- 1.96 .* sd_b, 0.0, Inf)
        upper = s_b .+ 1.96 .* sd_b

        for r in 1:n_reps_b
            y_obs     = @view obs_batch[r][:, col_idx]
            den       = (g_b[r] > 0) ? g_b[r] : eps(Float64)
            y_aligned = (y_obs .- @view o_b[:, r]) ./ den
            col       = palette[(r - 1) % length(palette) + 1]
            Makie.lines!(ax, x, y_aligned; color=col, linewidth=1, alpha=1.0, label = legend ? "rep $r" : nothing)
        end

        Makie.lines!(ax, x, s_b; color=:black, linewidth=0.5, label="shared signal")
        Makie.band!(ax, x, lower, upper; transparency=true, alpha=0.75)

    else  # :offsets
        for r in 1:n_reps_b
            col = palette[(r - 1) % length(palette) + 1]
            Makie.lines!(ax, 1:n_scans_b, @view o_b[:, r]; color=col, linewidth=1, label = legend ? "rep $r" : nothing)
        end
    end

    legend && Makie.axislegend(ax; position=:rt, framevisible=false)

    if scan_range !== nothing
        lo = clamp(first(scan_range), 1, n_scans_b)
        hi = clamp(last(scan_range),  1, n_scans_b)
        Makie.xlims!(ax, lo, hi)
    end

    fig
end
