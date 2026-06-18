import JuChrom: tictrace, tictrace!
using Unitful

const DEFAULT_TICTRACE_SIZE = (900, 450)
const DEFAULT_TICTRACE_COLOR = :black
const DEFAULT_TICTRACE_BASELINE_COLOR = :green
const DEFAULT_TICTRACE_BASELINE_ALPHA = 0.8
const DEFAULT_TICTRACE_STEP_COLOR = :orange
const DEFAULT_TICTRACE_STEP_ALPHA = 0.8
const DEFAULT_TICTRACE_LABEL_COLOR = :red
const DEFAULT_TICTRACE_LINEWIDTH = 1.0
const DEFAULT_TICTRACE_BASELINE_LINEWIDTH_FACTOR = 3.0
const DEFAULT_TICTRACE_STEP_LINEWIDTH = 2.0
const DEFAULT_TICTRACE_LABEL_FONTSIZE = 11
const DEFAULT_TICTRACE_LABEL_SPACING_FACTOR = 1.15
const DEFAULT_TICTRACE_Y_HEADROOM = 1.05
const DEFAULT_TICTRACE_LABEL_Y_FRACTION = 0.99

function tictrace(
    msm::JuChrom.AbstractMassScanMatrix;
    size=DEFAULT_TICTRACE_SIZE,
    unit::Union{Nothing, Unitful.Units}=JuChrom.retentionunit(msm),
    color=DEFAULT_TICTRACE_COLOR,
    linewidth::Real=DEFAULT_TICTRACE_LINEWIDTH
)
    fig = Figure(; size=size)
    ax = Axis(
        fig[1, 1],
        xlabel="Retention" * tictrace_unit_suffix(unit),
        ylabel="TIC" * tictrace_unit_suffix(JuChrom.intensityunit(msm); unitless=true),
        xgridvisible=false,
        ygridvisible=false
    )
    tictrace!(ax, msm; unit=unit, color=color, linewidth=linewidth)

    fig
end

function tictrace(
    msm::JuChrom.AbstractMassScanMatrix,
    result::JuChrom.AlkaneSeriesResult;
    size=DEFAULT_TICTRACE_SIZE,
    unit::Union{Nothing, Unitful.Units}=JuChrom.retentionunit(msm),
    color=DEFAULT_TICTRACE_COLOR,
    linewidth::Real=DEFAULT_TICTRACE_LINEWIDTH,
    baseline::Bool=true,
    baselinecolor=DEFAULT_TICTRACE_BASELINE_COLOR,
    baselinealpha::Real=DEFAULT_TICTRACE_BASELINE_ALPHA,
    baselinelinewidth::Real=DEFAULT_TICTRACE_BASELINE_LINEWIDTH_FACTOR * linewidth,
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true,
    stepcolor=DEFAULT_TICTRACE_STEP_COLOR,
    stepalpha::Real=DEFAULT_TICTRACE_STEP_ALPHA,
    steplinewidth::Real=DEFAULT_TICTRACE_STEP_LINEWIDTH,
    labelcolor=DEFAULT_TICTRACE_LABEL_COLOR,
    labels::Bool=false,
    labelfontsize::Real=DEFAULT_TICTRACE_LABEL_FONTSIZE,
    labelspacingfactor::Real=DEFAULT_TICTRACE_LABEL_SPACING_FACTOR
)
    fig = Figure(; size=size)
    ax = Axis(
        fig[1, 1],
        xlabel="Retention" * tictrace_unit_suffix(unit),
        ylabel="TIC" * tictrace_unit_suffix(JuChrom.intensityunit(msm); unitless=true),
        xgridvisible=false,
        ygridvisible=false
    )
    tictrace!(
        ax,
        msm,
        result;
        unit=unit,
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
        labelfontsize=labelfontsize,
        labelspacingfactor=labelspacingfactor
    )

    fig
end

function tictrace!(
    ax::Axis,
    msm::JuChrom.AbstractMassScanMatrix;
    unit::Union{Nothing, Unitful.Units}=JuChrom.retentionunit(msm),
    color=DEFAULT_TICTRACE_COLOR,
    linewidth::Real=DEFAULT_TICTRACE_LINEWIDTH
)
    x = JuChrom.rawretentions(msm; unit=unit)
    y = tictrace_values(msm)

    lines!(ax, x, y; color=color, linewidth=linewidth)
    tictrace_set_limits!(ax, x, y)

    ax
end

function tictrace!(
    ax::Axis,
    msm::JuChrom.AbstractMassScanMatrix,
    result::JuChrom.AlkaneSeriesResult;
    unit::Union{Nothing, Unitful.Units}=JuChrom.retentionunit(msm),
    color=DEFAULT_TICTRACE_COLOR,
    linewidth::Real=DEFAULT_TICTRACE_LINEWIDTH,
    baseline::Bool=true,
    baselinecolor=DEFAULT_TICTRACE_BASELINE_COLOR,
    baselinealpha::Real=DEFAULT_TICTRACE_BASELINE_ALPHA,
    baselinelinewidth::Real=DEFAULT_TICTRACE_BASELINE_LINEWIDTH_FACTOR * linewidth,
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true,
    stepcolor=DEFAULT_TICTRACE_STEP_COLOR,
    stepalpha::Real=DEFAULT_TICTRACE_STEP_ALPHA,
    steplinewidth::Real=DEFAULT_TICTRACE_STEP_LINEWIDTH,
    labelcolor=DEFAULT_TICTRACE_LABEL_COLOR,
    labels::Bool=false,
    labelfontsize::Real=DEFAULT_TICTRACE_LABEL_FONTSIZE,
    labelspacingfactor::Real=DEFAULT_TICTRACE_LABEL_SPACING_FACTOR
)
    x = JuChrom.rawretentions(msm; unit=unit)
    y = tictrace_values(msm)

    yseries = if baseline && !isnothing(result.baselineinfo)
        baseline_y = tictrace_values(result.baselineinfo.baselines)
        (y, baseline_y)
    else
        (y,)
    end

    ymax = tictrace_ymax(yseries)
    ylimit = DEFAULT_TICTRACE_Y_HEADROOM * ymax
    tictrace_add_ladder_steps!(
        ax,
        result,
        ylimit;
        unit=unit,
        molecularion=molecularion,
        gapfilled=gapfilled,
        edgeextended=edgeextended,
        stepcolor=stepcolor,
        stepalpha=stepalpha,
        steplinewidth=steplinewidth,
        labelcolor=labelcolor,
        labels=labels,
        labelfontsize=labelfontsize,
        labelspacingfactor=labelspacingfactor,
        xrange=(minimum(x), maximum(x))
    )

    if length(yseries) == 2
        lines!(
            ax,
            x,
            yseries[2];
            color=(baselinecolor, baselinealpha),
            linewidth=baselinelinewidth
        )
    end

    lines!(ax, x, y; color=color, linewidth=linewidth)
    tictrace_set_limits!(ax, x, yseries...)

    ax
end

function tictrace_values(msm::JuChrom.AbstractMassScanMatrix)
    vec(sum(JuChrom.rawintensities(msm); dims=2))
end

function tictrace_unit_suffix(unit; unitless::Bool=false)
    isnothing(unit) ? (unitless ? " [unitless]" : "") : " [$unit]"
end

function tictrace_ymax(yseries)
    ymax = 0.0
    for y in yseries
        isempty(y) && continue
        localmax = maximum(y)
        isfinite(localmax) && (ymax = max(ymax, Float64(localmax)))
    end

    ymax > 0 ? ymax : 1.0
end

function tictrace_set_limits!(ax::Axis, x::AbstractVector, yseries::AbstractVector...)
    ymax = tictrace_ymax(yseries)
    ylims!(ax, 0, DEFAULT_TICTRACE_Y_HEADROOM * ymax)
    isempty(x) || xlims!(ax, minimum(x), maximum(x))

    ax
end

function tictrace_add_ladder_steps!(
    ax::Axis,
    result::JuChrom.AlkaneSeriesResult,
    ylimit::Real;
    unit::Union{Nothing, Unitful.Units},
    molecularion::Bool,
    gapfilled::Bool,
    edgeextended::Bool,
    stepcolor,
    stepalpha::Real,
    steplinewidth::Real,
    labelcolor,
    labels::Bool,
    labelfontsize::Real,
    labelspacingfactor::Real,
    xrange
)
    steps = JuChrom.alkaneladdersteps(
        result;
        molecularion=molecularion,
        gapfilled=gapfilled,
        edgeextended=edgeextended
    )
    isempty(steps) && return ax

    xpositions = [
        tictrace_step_retention(step, result.retentionunit, unit) for step in steps
    ]
    vlines!(ax, xpositions; color=(stepcolor, stepalpha), linewidth=steplinewidth)

    if labels
        labelindices = tictrace_nonoverlapping_label_indices(
            ax,
            xpositions,
            xrange,
            labelfontsize,
            labelspacingfactor
        )
        isempty(labelindices) && return ax
        labelpositions = xpositions[labelindices]
        labelsteps = steps[labelindices]

        textlabel!(
            ax,
            labelpositions,
            fill(
                DEFAULT_TICTRACE_LABEL_Y_FRACTION * Float64(ylimit),
                length(labelpositions)
            );
            text=[string(step.ladderstep) for step in labelsteps],
            text_rotation=pi / 2,
            text_align=(:right, :center),
            text_color=labelcolor,
            background_color=:white,
            strokewidth=0,
            padding=1,
            cornerradius=0,
            font=:bold,
            fontsize=labelfontsize
        )
    end

    ax
end

function tictrace_nonoverlapping_label_indices(
    ax::Axis,
    xpositions::AbstractVector{<:Real},
    xrange,
    labelfontsize::Real,
    labelspacingfactor::Real
)
    xmin, xmax = Float64(xrange[1]), Float64(xrange[2])
    isfinite(xmin) && isfinite(xmax) && xmax > xmin || return collect(eachindex(xpositions))

    axiswidth = tictrace_axis_pixel_width(ax)
    minpixels = max(1.0, Float64(labelfontsize) * Float64(labelspacingfactor))
    mindata = minpixels / axiswidth * (xmax - xmin)

    indices = Int[]
    lastx = -Inf
    for index in sortperm(xpositions)
        x = Float64(xpositions[index])
        isfinite(x) || continue
        if isempty(indices) || x - lastx ≥ mindata
            push!(indices, index)
            lastx = x
        end
    end

    sort!(indices)
    indices
end

function tictrace_axis_pixel_width(ax::Axis)
    try
        width = Float64(ax.scene.viewport[].widths[1])
        isfinite(width) && width > 0 && return width
    catch
    end

    Float64(DEFAULT_TICTRACE_SIZE[1])
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
    Float64(ustrip(unit, step.apexretention * resultunit))
end
