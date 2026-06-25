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
        linecolor=:black,
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
        labelspacingfactor=1.15,
        yheadroom=1.05,
        labelyfraction=0.99
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

    map!(
        trace.attributes,
        [:retentions, :intensities, :baseline_intensities, :yheadroom],
        :limits_update
    ) do retentions, intensities, baseline_intensities, yheadroom
        tictrace_set_limits!(ax, retentions, yheadroom, intensities, baseline_intensities)
        nothing
    end

    map!(
        trace.attributes,
        [
            :step_retentions,
            :step_numbers,
            :labels,
            :labelfont,
            :labelfontsize,
            :labelpadding,
            :labelspacingfactor,
            :labelyfraction,
            :viewport_obs,
            :limits_obs
        ],
        [:label_positions, :label_texts]
    ) do step_retentions, step_numbers, labels, labelfont, labelfontsize,
            labelpadding, labelspacingfactor, labelyfraction, viewport, limits
        tictrace_visible_label_data(
            ax,
            step_retentions,
            step_numbers,
            labels,
            labelfont,
            labelfontsize,
            labelpadding,
            labelspacingfactor,
            labelyfraction,
            viewport,
            limits
        )
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
    textlabel!(
        trace,
        trace.label_positions;
        text=trace.label_texts,
        text_rotation=pi / 2,
        text_align=(:right, :center),
        text_color=trace.labelcolor,
        background_color=:white,
        strokewidth=0,
        padding=trace.labelpadding,
        cornerradius=0,
        font=trace.labelfont,
        fontsize=trace.labelfontsize
    )

    trace
end

function tictrace(
    msm::JuChrom.AbstractMassScanMatrix,
    result::JuChrom.AlkaneSeriesResult;
    size=(900, 450),
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    color=:black,
    linewidth::Real=1.0,
    baseline::Bool=true,
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
    labelspacingfactor::Real=1.15,
    yheadroom::Real=1.05,
    labelyfraction::Real=0.99
)
    displayretentionunit = tictrace_display_retentionunit(msm, retentionunit)
    displayintensityunit = tictrace_display_intensityunit(msm, intensityunit)

    fig = Figure(; size=size)
    ax = Axis(
        fig[1, 1],
        xlabel="Retention" * tictrace_unit_suffix(displayretentionunit),
        ylabel="TIC" * tictrace_unit_suffix(displayintensityunit; unitless=true),
        xgridvisible=false,
        ygridvisible=false
    )
    tictrace!(
        ax,
        msm,
        result;
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
        labelspacingfactor=labelspacingfactor,
        yheadroom=yheadroom,
        labelyfraction=labelyfraction
    )

    fig
end

function tictrace!(
    ax::Axis,
    msm::JuChrom.AbstractMassScanMatrix,
    result::JuChrom.AlkaneSeriesResult;
    retentionunit::Union{Nothing, Unitful.Units}=nothing,
    intensityunit::Union{Nothing, Unitful.Units}=nothing,
    color=:black,
    linewidth::Real=1.0,
    baseline::Bool=true,
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
    labelspacingfactor::Real=1.15,
    yheadroom::Real=1.05,
    labelyfraction::Real=0.99
)
    tictrace_validate_checksum(msm, result)
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
    tictrace_set_limits!(ax, x, yheadroom, y, baseline_y)

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
        labelspacingfactor=labelspacingfactor,
        yheadroom=yheadroom,
        labelyfraction=labelyfraction
    )

    ax
end

"""
    lines(data; retentionunit=nothing, intensityunit=nothing, figure=(;), axis=(;), kwargs...)
    lines(position, data; retentionunit=nothing, intensityunit=nothing, axis=(;), kwargs...)
    lines!(ax, data; retentionunit=nothing, intensityunit=nothing, kwargs...)

Plot the total ion chromatogram for a JuChrom scan container with Makie.

`data` may be an `AbstractMassScanMatrix`, `AbstractVarianceMassScanMatrix`,
`AbstractMassScanSeries`, or `AbstractChromScanSeries`. For mass-scan data, the y values
are the per-scan sums of intensities. For chromatographic scan series, the y values are
the scan intensities.

`retentionunit` and `intensityunit` are converted before numeric values are stripped. The
default `nothing` strips the stored unit when present and leaves unitless data unchanged.
All other keywords are forwarded to Makie's `lines` plot.
"""
const TicTraceLineData = Union{
    JuChrom.AbstractMassScanMatrix,
    JuChrom.AbstractChromScanSeries,
    JuChrom.AbstractMassScanSeries
}

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
    yseries::AbstractVector...
)
    ymax = tictrace_ymax(yseries)
    ylims!(ax, 0, yheadroom * ymax)
    isempty(x) || xlims!(ax, minimum(x), maximum(x))

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
    labelspacingfactor::Real,
    labelyfraction::Real,
    viewport,
    limits
)
    labels || return Point2f[], String[]
    isempty(step_retentions) && return Point2f[], String[]
    tictrace_limits_are_usable(limits) || return Point2f[], String[]

    y = limits.origin[2] + labelyfraction * limits.widths[2]
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
