import JuChrom: tictrace, tictrace!
using Unitful

@recipe(
    TicTrace,
    retentions,
    intensities,
    baseline_intensities,
    step_retentions,
    step_numbers
) do scene
    Theme(
        linecolor=Makie.wong_colors()[1],
        linewidth=1.0,
        baselinecolor=:green,
        baselinealpha=0.8,
        baselinelinewidth=3.0,
        stepcolor=:orange,
        stepalpha=0.8,
        steplinewidth=2.0,
        labelcolor=:red,
        labels=true,
        labelfont=:bold,
        labelfontsize=11,
        labelpadding=1,
        labeltopgappixels=3,
        labelspacingfactor=1.15,
        yheadroom=1.05,
        steplabelgappixels=3,
        ticstepgappixels=6
    )
end

argument_names(::Type{<:TicTrace}) = (
    :retentions,
    :intensities,
    :baseline_intensities,
    :step_retentions,
    :step_numbers
)

function Makie.plot!(trace::TicTrace)
    ax = current_axis()

    Makie.add_input!(trace.attributes, :viewport_obs, ax.scene.viewport)
    Makie.add_input!(trace.attributes, :limits_obs, ax.finallimits)

    tictrace_register_limit_updates!(ax, trace)

    map!(
        trace.attributes,
        [
            :step_retentions,
            :step_numbers,
            :labels,
            :labelfont,
            :labelfontsize,
            :labelpadding,
            :labeltopgappixels,
            :labelspacingfactor,
            :steplabelgappixels,
            :viewport_obs,
            :limits_obs
        ],
        [:label_positions, :label_texts, :step_ymax]
    ) do step_retentions, step_numbers, labels, labelfont, labelfontsize,
            labelpadding, labeltopgappixels, labelspacingfactor, steplabelgappixels,
            viewport, limits
        label_positions, label_texts = tictrace_visible_label_data(
            ax,
            step_retentions,
            step_numbers,
            labels,
            labelfont,
            labelfontsize,
            labelpadding,
            labeltopgappixels,
            labelspacingfactor,
            steplabelgappixels,
            viewport,
            limits
        )
        layout = tictrace_ladder_y_layout(
            step_numbers,
            labels,
            labelfont,
            labelfontsize,
            labelpadding,
            labeltopgappixels,
            steplabelgappixels,
            0,
            viewport
        )
        step_ymax = layout.step_ymax
        label_positions, label_texts, step_ymax
    end

    lines!(
        trace,
        trace.retentions,
        trace.baseline_intensities;
        color=lift((color, alpha) -> (color, alpha), trace.baselinecolor, trace.baselinealpha),
        linewidth=trace.baselinelinewidth
    )
    vlines!(
        trace,
        trace.step_retentions;
        ymax=trace.step_ymax,
        color=lift((color, alpha) -> (color, alpha), trace.stepcolor, trace.stepalpha),
        linewidth=trace.steplinewidth
    )
    lines!(
        trace,
        trace.retentions,
        trace.intensities;
        color=trace.linecolor,
        linewidth=trace.linewidth
    )
    text!(
        trace,
        trace.label_positions;
        text=trace.label_texts,
        rotation=pi / 2,
        align=(:right, :center),
        color=trace.labelcolor,
        font=trace.labelfont,
        fontsize=trace.labelfontsize
    )

    trace
end

function tictrace_register_limit_updates!(ax::Axis, trace::TicTrace)
    Makie.Observables.onany(
        trace,
        trace.retentions,
        trace.intensities,
        trace.baseline_intensities,
        trace.yheadroom,
        trace.step_numbers,
        trace.labels,
        trace.labelfont,
        trace.labelfontsize,
        trace.labelpadding,
        trace.labeltopgappixels,
        trace.steplabelgappixels,
        trace.ticstepgappixels,
        ax.scene.viewport;
        update=true
    ) do retentions, intensities, baseline_intensities, yheadroom, step_numbers, labels,
            labelfont, labelfontsize, labelpadding, labeltopgappixels, steplabelgappixels,
            ticstepgappixels, viewport
        tictrace_set_limits!(
            ax,
            retentions,
            yheadroom,
            intensities,
            baseline_intensities;
            step_numbers=step_numbers,
            labels=labels,
            labelfont=labelfont,
            labelfontsize=labelfontsize,
            labelpadding=labelpadding,
            labeltopgappixels=labeltopgappixels,
            steplabelgappixels=steplabelgappixels,
            ticstepgappixels=ticstepgappixels,
            viewport=viewport
        )
        nothing
    end

    trace
end

function tictrace(
    msm::JuChrom.AbstractMassScanMatrix,
    result::JuChrom.AlkaneSeriesResult;
    size=(900, 450),
    figure::NamedTuple=NamedTuple(),
    axis::NamedTuple=NamedTuple(),
    title::Union{Nothing, AbstractString}=nothing,
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    color=Makie.wong_colors()[1],
    linewidth::Real=1.0,
    baseline::Bool=false,
    baselinecolor=:green,
    baselinealpha::Real=0.8,
    baselinelinewidth::Real=3.0 * linewidth,
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true,
    stepcolor=:orange,
    stepalpha::Real=0.8,
    steplinewidth::Real=2.0,
    labelcolor=:red,
    labels::Bool=true,
    labelfont=:bold,
    labelfontsize::Real=11,
    labelpadding::Real=1,
    labeltopgappixels::Real=3,
    labelspacingfactor::Real=1.15,
    yheadroom::Real=1.05,
    steplabelgappixels::Real=3,
    ticstepgappixels::Real=6
)
    fig = Figure(; merge((; size=size), figure)...)
    ax = Axis(
        fig[1, 1];
        tictrace_ladder_axis_attributes(
            msm,
            retentionunit,
            intensityunit,
            axis;
            title=title
        )...
    )
    tictrace!(
        ax,
        msm,
        result;
        axis=axis,
        title=title,
        retentionunit=retentionunit,
        intensityunit=intensityunit,
        color=color,
        linewidth=linewidth,
        baseline=baseline,
        baselinecolor=baselinecolor,
        baselinealpha=baselinealpha,
        baselinelinewidth=baselinelinewidth,
        molecularion=molecularion,
        gapfilled=gapfilled,
        edgeextended=edgeextended,
        stepcolor=stepcolor,
        stepalpha=stepalpha,
        steplinewidth=steplinewidth,
        labelcolor=labelcolor,
        labels=labels,
        labelfont=labelfont,
        labelfontsize=labelfontsize,
        labelpadding=labelpadding,
        labeltopgappixels=labeltopgappixels,
        labelspacingfactor=labelspacingfactor,
        yheadroom=yheadroom,
        steplabelgappixels=steplabelgappixels,
        ticstepgappixels=ticstepgappixels
    )

    fig
end

function tictrace!(
    ax::Axis,
    msm::JuChrom.AbstractMassScanMatrix,
    result::JuChrom.AlkaneSeriesResult;
    axis::NamedTuple=NamedTuple(),
    title::Union{Nothing, AbstractString}=nothing,
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    color=Makie.wong_colors()[1],
    linewidth::Real=1.0,
    baseline::Bool=false,
    baselinecolor=:green,
    baselinealpha::Real=0.8,
    baselinelinewidth::Real=3.0 * linewidth,
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true,
    stepcolor=:orange,
    stepalpha::Real=0.8,
    steplinewidth::Real=2.0,
    labelcolor=:red,
    labels::Bool=true,
    labelfont=:bold,
    labelfontsize::Real=11,
    labelpadding::Real=1,
    labeltopgappixels::Real=3,
    labelspacingfactor::Real=1.15,
    yheadroom::Real=1.05,
    steplabelgappixels::Real=3,
    ticstepgappixels::Real=6
)
    tictrace_validate_checksum(msm, result)
    tictrace_set_axis_attributes!(
        ax,
        tictrace_ladder_axis_attributes(
            msm,
            retentionunit,
            intensityunit,
            axis;
            title=title
        )
    )
    x, y = tictrace_xy(msm, retentionunit, intensityunit)
    baseline_y = tictrace_baseline_values(
        result,
        baseline,
        length(x);
        intensityunit=intensityunit
    )
    step_retentions, step_numbers = tictrace_ladder_step_data(
        result,
        retentionunit,
        molecularion,
        gapfilled,
        edgeextended
    )
    tictrace_set_limits!(
        ax,
        x,
        yheadroom,
        y,
        baseline_y;
        step_numbers=step_numbers,
        labels=labels,
        labelfont=labelfont,
        labelfontsize=labelfontsize,
        labelpadding=labelpadding,
        labeltopgappixels=labeltopgappixels,
        steplabelgappixels=steplabelgappixels,
        ticstepgappixels=ticstepgappixels,
        viewport=ax.scene.viewport[]
    )

    tictrace!(
        ax,
        x,
        y,
        baseline_y,
        step_retentions,
        step_numbers;
        linecolor=color,
        linewidth=linewidth,
        baselinecolor=baselinecolor,
        baselinealpha=baselinealpha,
        baselinelinewidth=baselinelinewidth,
        stepcolor=stepcolor,
        stepalpha=stepalpha,
        steplinewidth=steplinewidth,
        labelcolor=labelcolor,
        labels=labels,
        labelfont=labelfont,
        labelfontsize=labelfontsize,
        labelpadding=labelpadding,
        labeltopgappixels=labeltopgappixels,
        labelspacingfactor=labelspacingfactor,
        yheadroom=yheadroom,
        steplabelgappixels=steplabelgappixels,
        ticstepgappixels=ticstepgappixels
    )

    ax
end

const TicTraceLineData = Union{
    JuChrom.AbstractMassScanMatrix,
    JuChrom.AbstractChromScanSeries,
    JuChrom.AbstractMassScanSeries
}

"""
    tictrace(data; retentionunit=nothing, intensityunit=nothing,
             figure=(;), axis=(;), title=nothing, kwargs...) -> Figure
    tictrace(position, data; retentionunit=nothing, intensityunit=nothing,
             axis=(;), title=nothing, kwargs...) -> AxisPlot
    tictrace!(ax, data; retentionunit=nothing, intensityunit=nothing,
              axis=(;), title=nothing, kwargs...) -> Lines

Plot the total ion chromatogram for a JuChrom scan container with Makie.

`data` may be an `AbstractMassScanMatrix`, `AbstractVarianceMassScanMatrix`,
`AbstractMassScanSeries`, or `AbstractChromScanSeries`. For mass-scan data, the y values
are the per-scan sums of intensities. For chromatographic scan series, the y values are
the scan intensities.

`retentionunit` and `intensityunit` are converted before numeric values are stripped. The
default `nothing` strips the stored unit when present and leaves unitless data unchanged.
Line keywords are forwarded to Makie's `lines!`. Axis keywords are passed via `axis`, and
the title is empty unless `title` or `axis=(; title=...)` is supplied. Tick intervals can
be controlled with Makie's `xticks`, for example `axis=(; xticks=0:1:40)`.
"""
function tictrace(
    data::TicTraceLineData;
    figure::NamedTuple=NamedTuple(),
    axis::NamedTuple=NamedTuple(),
    title::Union{Nothing, AbstractString}=nothing,
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    kwargs...
)
    fig = Figure(; figure...)
    ax = Axis(
        fig[1, 1];
        tictrace_axis_attributes(
            data,
            retentionunit,
            intensityunit,
            axis;
            title=title,
            defaulttitle=true
        )...
    )
    tictrace!(
        ax,
        data;
        axis=axis,
        title=title,
        retentionunit=retentionunit,
        intensityunit=intensityunit,
        kwargs...
    )

    fig
end

function tictrace(
    position::Union{Makie.GridPosition, Makie.GridSubposition},
    data::TicTraceLineData;
    figure::NamedTuple=NamedTuple(),
    axis::NamedTuple=NamedTuple(),
    title::Union{Nothing, AbstractString}=nothing,
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    kwargs...
)
    isempty(pairs(figure)) || throw(ArgumentError(
        "figure keyword is not supported when plotting into a grid position"))
    isempty(Makie.contents(position; exact=true)) || error(
        "non-mutating tictrace requires an empty grid position; use tictrace! for existing axes")
    ax = Axis(
        position;
        tictrace_axis_attributes(
            data,
            retentionunit,
            intensityunit,
            axis;
            title=title,
            defaulttitle=true
        )...
    )
    plt = tictrace!(
        ax,
        data;
        axis=axis,
        title=title,
        retentionunit=retentionunit,
        intensityunit=intensityunit,
        kwargs...
    )

    Makie.AxisPlot(ax, plt)
end

function tictrace!(
    ax::Axis,
    data::TicTraceLineData;
    axis::NamedTuple=NamedTuple(),
    title::Union{Nothing, AbstractString}=nothing,
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    kwargs...
)
    tictrace_set_axis_attributes!(
        ax,
        tictrace_axis_attributes(
            data,
            retentionunit,
            intensityunit,
            axis;
            title=title,
            defaulttitle=false
        )
    )
    x, y = tictrace_xy(data, retentionunit, intensityunit)
    plt = Makie.lines!(ax, x, y; kwargs...)
    tictrace_set_xlimits!(ax, x)

    plt
end

"""
    lines(data; retentionunit=nothing, intensityunit=nothing, figure=(;), axis=(;), kwargs...)
    lines(position, data; retentionunit=nothing, intensityunit=nothing, axis=(;), kwargs...)
    lines!(ax, data; retentionunit=nothing, intensityunit=nothing, kwargs...)

Lower-level Makie line methods for plotting the total ion chromatogram of a JuChrom scan
container. These methods use the same TIC extraction and unit-conversion rules as
[`tictrace`](@ref), but do not set axis labels.
"""

function Makie.lines(
    data::TicTraceLineData;
    figure::NamedTuple=NamedTuple(),
    axis::NamedTuple=NamedTuple(),
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    kwargs...
)
    tictrace_lines(data, figure, axis, retentionunit, intensityunit; kwargs...)
end

function Makie.lines(
    position::Union{Makie.GridPosition, Makie.GridSubposition},
    data::TicTraceLineData;
    figure::NamedTuple=NamedTuple(),
    axis::NamedTuple=NamedTuple(),
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    kwargs...
)
    isempty(pairs(figure)) || throw(ArgumentError(
        "figure keyword is not supported when plotting into a grid position"))
    tictrace_lines(position, axis, data, retentionunit, intensityunit; kwargs...)
end

function Makie.lines!(
    ax::Axis,
    data::TicTraceLineData;
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    kwargs...
)
    tictrace_lines!(ax, data, retentionunit, intensityunit; kwargs...)
end

function tictrace_lines(
    data,
    figure::NamedTuple,
    axis::NamedTuple,
    retentionunit::Union{Nothing, Unitful.Units},
    intensityunit::Union{Nothing, Unitful.Units};
    kwargs...
)
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; axis...)
    plt = Makie.lines!(
        ax,
        data;
        retentionunit=retentionunit,
        intensityunit=intensityunit,
        kwargs...
    )

    Makie.FigureAxisPlot(fig, ax, plt)
end

function tictrace_lines(
    position::Union{Makie.GridPosition, Makie.GridSubposition},
    axis::NamedTuple,
    data,
    retentionunit::Union{Nothing, Unitful.Units},
    intensityunit::Union{Nothing, Unitful.Units};
    kwargs...
)
    isempty(Makie.contents(position; exact=true)) || error(
        "non-mutating lines requires an empty grid position; use lines! for existing axes")
    ax = Axis(position; axis...)
    plt = Makie.lines!(
        ax,
        data;
        retentionunit=retentionunit,
        intensityunit=intensityunit,
        kwargs...
    )

    Makie.AxisPlot(ax, plt)
end

function tictrace_lines!(
    ax::Axis,
    data,
    retentionunit::Union{Nothing, Unitful.Units},
    intensityunit::Union{Nothing, Unitful.Units};
    kwargs...
)
    x, y = tictrace_xy(data, retentionunit, intensityunit)

    Makie.lines!(ax, x, y; kwargs...)
end

function tictrace_xy(
    data,
    retentionunit::Union{Nothing, Unitful.Units},
    intensityunit::Union{Nothing, Unitful.Units}
)
    (
        JuChrom.rawretentions(data; unit=retentionunit),
        tictrace_values(data; intensityunit=intensityunit)
    )
end

function tictrace_validate_checksum(
    msm::JuChrom.AbstractMassScanMatrix,
    result::JuChrom.AlkaneSeriesResult
)
    JuChrom.alkane_validate_raw_msm_checksum(msm, result)
end

function tictrace_values(
    msm::JuChrom.AbstractMassScanMatrix;
    intensityunit::Union{Nothing, Unitful.Units}=nothing
)
    vec(sum(JuChrom.rawintensities(msm; unit=intensityunit); dims=2))
end

function tictrace_values(
    series::JuChrom.AbstractChromScanSeries;
    intensityunit::Union{Nothing, Unitful.Units}=nothing
)
    JuChrom.rawintensities(series; unit=intensityunit)
end

function tictrace_values(
    series::JuChrom.AbstractMassScanSeries;
    intensityunit::Union{Nothing, Unitful.Units}=nothing
)
    [
        sum(JuChrom.rawintensities(series, scanindex; unit=intensityunit))
        for scanindex in 1:JuChrom.scancount(series)
    ]
end

function tictrace_baseline_values(
    result::JuChrom.AlkaneSeriesResult,
    baseline::Bool,
    nscans::Integer;
    intensityunit::Union{Nothing, Unitful.Units}=nothing
)
    baseline && !isnothing(result.baselineinfo) ||
        return fill(NaN, nscans)

    tictrace_values(result.baselineinfo.baselines; intensityunit=intensityunit)
end

function tictrace_unit_suffix(unit; unitless::Bool=false)
    isnothing(unit) ? (unitless ? " [unitless]" : "") : " [$unit]"
end

function tictrace_axis_unit_label(::Nothing)
    return "unitless"
end

tictrace_axis_unit_label(unit::Unitful.Units) = string(unit)

function tictrace_axis_label(label::AbstractString, unit)
    "$label [$(tictrace_axis_unit_label(unit))]"
end

function tictrace_axis_attributes(
    data,
    retentionunit::Union{Nothing, Unitful.Units},
    intensityunit::Union{Nothing, Unitful.Units},
    axis::NamedTuple;
    title::Union{Nothing, AbstractString}=nothing,
    defaulttitle::Bool=false
)
    defaults = (
        xlabel=tictrace_axis_label(
            "Retention",
            tictrace_display_retentionunit(data, retentionunit)
        ),
        ylabel=tictrace_axis_label(
            "Intensity",
            tictrace_display_intensityunit(data, intensityunit)
        )
    )
    defaulttitle && (defaults = merge(defaults, (; title="")))
    attributes = merge(defaults, axis)
    isnothing(title) ? attributes : merge(attributes, (; title=title))
end

function tictrace_ladder_axis_attributes(
    data,
    retentionunit::Union{Nothing, Unitful.Units},
    intensityunit::Union{Nothing, Unitful.Units},
    axis::NamedTuple;
    title::Union{Nothing, AbstractString}=nothing
)
    defaults = (
        xlabel="Retention" * tictrace_unit_suffix(
            tictrace_display_retentionunit(data, retentionunit)
        ),
        ylabel="TIC" * tictrace_unit_suffix(
            tictrace_display_intensityunit(data, intensityunit);
            unitless=true
        ),
        xgridvisible=false,
        ygridvisible=false
    )
    attributes = merge(defaults, axis)
    isnothing(title) ? attributes : merge(attributes, (; title=title))
end

function tictrace_set_axis_attributes!(ax::Axis, attributes::NamedTuple)
    for (key, value) in pairs(attributes)
        setproperty!(ax, key, value)
    end

    ax
end

function tictrace_display_retentionunit(data, retentionunit::Nothing)
    JuChrom.retentionunit(data)
end

function tictrace_display_retentionunit(data, retentionunit::Unitful.Units)
    retentionunit
end

function tictrace_display_intensityunit(data, intensityunit::Nothing)
    JuChrom.intensityunit(data)
end

function tictrace_display_intensityunit(data, intensityunit::Unitful.Units)
    intensityunit
end

function tictrace_ymax(yseries)
    ymax = 0.0
    for y in yseries
        isempty(y) && continue
        localmax = maximum(y)
        isfinite(localmax) && (ymax = max(ymax, localmax))
    end

    ymax > 0 ? ymax : 1.0
end

function tictrace_set_limits!(
    ax::Axis,
    x::AbstractVector,
    yheadroom::Real,
    yseries::AbstractVector...;
    step_numbers::AbstractVector{<:Integer}=Int[],
    labels::Bool=false,
    labelfont=:bold,
    labelfontsize::Real=11,
    labelpadding::Real=1,
    labeltopgappixels::Real=3,
    steplabelgappixels::Real=3,
    ticstepgappixels::Real=6,
    viewport=nothing
)
    ymax = tictrace_ymax(yseries)
    ylimit = yheadroom * ymax
    if labels && !isempty(step_numbers) && !isnothing(viewport)
        layout = tictrace_ladder_y_layout(
            step_numbers,
            labels,
            labelfont,
            labelfontsize,
            labelpadding,
            labeltopgappixels,
            steplabelgappixels,
            ticstepgappixels,
            viewport
        )
        layout.labelbandpixels > 0 && layout.data_ymax > 0 &&
            (ylimit = ymax / layout.data_ymax)
    end

    ylims!(ax, 0, ylimit)
    tictrace_set_xlimits!(ax, x)

    ax
end

function tictrace_set_xlimits!(ax::Axis, x::AbstractVector)
    isempty(x) && return ax

    xmin, xmax = extrema(x)
    xmin == xmax ? xlims!(ax, xmin - 0.5, xmax + 0.5) : xlims!(ax, xmin, xmax)

    ax
end

function tictrace_ladder_step_data(
    result::JuChrom.AlkaneSeriesResult,
    unit::Union{Nothing, Unitful.Units},
    molecularion::Bool,
    gapfilled::Bool,
    edgeextended::Bool
)
    steps = JuChrom.alkaneladdersteps(
        result;
        molecularion=molecularion,
        gapfilled=gapfilled,
        edgeextended=edgeextended
    )

    step_retentions = [
        tictrace_step_retention(step, result.retentionunit, unit) for step in steps
    ]
    step_numbers = [step.ladderstep for step in steps]

    step_retentions, step_numbers
end

function tictrace_visible_label_data(
    ax::Axis,
    step_retentions::AbstractVector{<:Real},
    step_numbers::AbstractVector{<:Integer},
    labels::Bool,
    labelfont,
    labelfontsize::Real,
    labelpadding::Real,
    labeltopgappixels::Real,
    labelspacingfactor::Real,
    steplabelgappixels::Real,
    viewport,
    limits
)
    labels || return Point2f[], String[]
    isempty(step_retentions) && return Point2f[], String[]
    tictrace_limits_are_usable(limits) || return Point2f[], String[]

    layout = tictrace_ladder_y_layout(
        step_numbers,
        labels,
        labelfont,
        labelfontsize,
        labelpadding,
        labeltopgappixels,
        steplabelgappixels,
        0,
        viewport
    )
    y = limits.origin[2] + layout.label_yfraction * limits.widths[2]
    positions = Point2f[]
    labeltexts = String[]
    occupied = Tuple{Float64, Float64}[]

    for index in sortperm(step_retentions)
        x = step_retentions[index]
        limits.origin[1] ≤ x ≤ limits.origin[1] + limits.widths[1] || continue
        text = string(step_numbers[index])
        x0, x1 = tictrace_label_pixel_xspan(
            ax,
            viewport,
            limits,
            text,
            x,
            y,
            labelfont,
            labelfontsize,
            labelpadding,
            labelspacingfactor
        )
        tictrace_xspan_inside_viewport(x0, x1, viewport) || continue
        any(region -> tictrace_xspans_overlap(x0, x1, region...), occupied) &&
            continue

        push!(positions, Point2f(x, y))
        push!(labeltexts, text)
        push!(occupied, (x0, x1))
    end

    positions, labeltexts
end

function tictrace_limits_are_usable(limits)
    all(isfinite, limits.origin) && all(isfinite, limits.widths) &&
        all(>(0), limits.widths)
end

function tictrace_ladder_y_layout(
    step_numbers::AbstractVector{<:Integer},
    labels::Bool,
    labelfont,
    labelfontsize::Real,
    labelpadding::Real,
    labeltopgappixels::Real,
    steplabelgappixels::Real,
    ticstepgappixels::Real,
    viewport
)
    fallback = (
        label_yfraction=1.0,
        step_ymax=1.0,
        data_ymax=1.0,
        labelbandpixels=0.0,
        labelreservepixels=0.0
    )
    labels || return fallback
    isempty(step_numbers) && return fallback
    isnothing(viewport) && return fallback
    viewportheight = Float64(viewport.widths[2])
    isfinite(viewportheight) && viewportheight > 0 || return fallback

    labelpadding = max(0.0, Float64(labelpadding))
    labeltopgappixels = max(0.0, Float64(labeltopgappixels))
    steplabelgappixels = max(0.0, Float64(steplabelgappixels))
    ticstepgappixels = max(0.0, Float64(ticstepgappixels))
    labelbandpixels = tictrace_label_band_pixels(
        step_numbers,
        labelfont,
        labelfontsize,
        labelpadding
    )
    labelreservepixels = max(0.0, labeltopgappixels + labelbandpixels - labelpadding)
    label_yfraction = clamp(
        (viewportheight - labeltopgappixels) / viewportheight,
        0.0,
        1.0
    )
    step_ymax = clamp(
        (viewportheight - labelreservepixels - steplabelgappixels) / viewportheight,
        0.0,
        1.0
    )
    data_ymax = clamp(
        (viewportheight - labelreservepixels - steplabelgappixels - ticstepgappixels) /
            viewportheight,
        0.0,
        1.0
    )

    (
        label_yfraction=label_yfraction,
        step_ymax=step_ymax,
        data_ymax=data_ymax,
        labelbandpixels=labelbandpixels,
        labelreservepixels=labelreservepixels
    )
end

function tictrace_label_band_pixels(
    step_numbers::AbstractVector{<:Integer},
    labelfont,
    labelfontsize::Real,
    labelpadding::Real
)
    font = tictrace_measurement_font(labelfont)
    labelheight = 0.0
    for step_number in step_numbers
        box = Makie.text_bb(string(step_number), font, labelfontsize)
        labelheight = max(labelheight, Float64(box.widths[1]))
    end

    labelheight + 2 * max(0.0, Float64(labelpadding))
end

function tictrace_label_pixel_xspan(
    ax::Axis,
    viewport,
    limits,
    text::AbstractString,
    x::Real,
    y::Real,
    labelfont,
    labelfontsize::Real,
    labelpadding::Real,
    labelspacingfactor::Real
)
    xpixel, _ = data2pixel(ax, viewport, limits, x, y)
    box = Makie.text_bb(text, tictrace_measurement_font(labelfont), labelfontsize)
    labelwidth = labelspacingfactor * (box.widths[2] + 2 * labelpadding)
    xpixel - labelwidth / 2, xpixel + labelwidth / 2
end

function tictrace_measurement_font(font::Symbol)
    Makie.to_font(Makie.MAKIE_DEFAULT_THEME.fonts, font)
end

function tictrace_measurement_font(font)
    Makie.to_font(font)
end

function tictrace_xspans_overlap(x0::Real, x1::Real, other0::Real, other1::Real)
    max(x0, other0) ≤ min(x1, other1)
end

function tictrace_xspan_inside_viewport(x0::Real, x1::Real, viewport)
    viewportwidth = Float64(viewport.widths[1])
    isfinite(viewportwidth) && viewportwidth > 0 || return false
    0 ≤ x0 && x1 ≤ viewportwidth
end

function tictrace_step_retention(
    step::JuChrom.AlkaneLadderStep,
    resultunit::Nothing,
    unit::Nothing
)
    step.apexretention
end

function tictrace_step_retention(
    step::JuChrom.AlkaneLadderStep,
    resultunit::Nothing,
    unit::Unitful.Units
)
    throw(ArgumentError("cannot convert unitless ladder step retentions to a unit"))
end

function tictrace_step_retention(
    step::JuChrom.AlkaneLadderStep,
    resultunit::Unitful.Units,
    unit::Nothing
)
    step.apexretention
end

function tictrace_step_retention(
    step::JuChrom.AlkaneLadderStep,
    resultunit::Unitful.Units,
    unit::Unitful.Units
)
    ustrip(unit, step.apexretention * resultunit)
end
