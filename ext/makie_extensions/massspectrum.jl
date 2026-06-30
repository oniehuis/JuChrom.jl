import JuChrom: massspectrum, massspectrum!
using Unitful

@recipe(MassSpectrum, mzvalues, intensities) do scene
    Theme(
        mz_label_font = Makie.defaultfont(),
        mz_label_fontsize = 10,
        mz_label_offset = (0, 5),
        mz_label_padding = 1,
        mz_label_top_padding_factor = 1.5,
        mz_label_color = :black,
        linecolor = :blue,
        linewidth = 1.0,
        intensity_threshold = 0.0
    )
end

argument_names(::Type{<:MassSpectrum}) = (:mzvalues, :intensities)

function massspectrum(
    mzvalues::AbstractVector{<:Unitful.AbstractQuantity},
    intensities::AbstractVector{<:Unitful.AbstractQuantity};
    figure::NamedTuple=NamedTuple(),
    axis::NamedTuple=NamedTuple(),
    kwargs...
)
    massspectrum_unitful_vectors(mzvalues, intensities, figure, axis; kwargs...)
end

function massspectrum(
    mzvalues::AbstractVector{<:Unitful.AbstractQuantity},
    intensities::AbstractVector;
    figure::NamedTuple=NamedTuple(),
    axis::NamedTuple=NamedTuple(),
    kwargs...
)
    massspectrum_unitful_vectors(mzvalues, intensities, figure, axis; kwargs...)
end

function massspectrum(
    mzvalues::AbstractVector,
    intensities::AbstractVector{<:Unitful.AbstractQuantity};
    figure::NamedTuple=NamedTuple(),
    axis::NamedTuple=NamedTuple(),
    kwargs...
)
    massspectrum_unitful_vectors(mzvalues, intensities, figure, axis; kwargs...)
end

function massspectrum!(
    ax::Axis,
    mzvalues::AbstractVector{<:Unitful.AbstractQuantity},
    intensities::AbstractVector{<:Unitful.AbstractQuantity};
    kwargs...
)
    massspectrum_unitful_vectors!(ax, mzvalues, intensities; kwargs...)
end

function massspectrum!(
    ax::Axis,
    mzvalues::AbstractVector{<:Unitful.AbstractQuantity},
    intensities::AbstractVector;
    kwargs...
)
    massspectrum_unitful_vectors!(ax, mzvalues, intensities; kwargs...)
end

function massspectrum!(
    ax::Axis,
    mzvalues::AbstractVector,
    intensities::AbstractVector{<:Unitful.AbstractQuantity};
    kwargs...
)
    massspectrum_unitful_vectors!(ax, mzvalues, intensities; kwargs...)
end

function Makie.convert_arguments(
    ::Type{<:MassSpectrum},
    mzvalues::AbstractVector,
    intensities::AbstractVector
)
    (mzvalues, intensities)
end

function Makie.convert_arguments(
    ::Type{<:MassSpectrum},
    spectrum::JuChrom.AbstractMassSpectrum
)
    (JuChrom.mzvalues(spectrum), JuChrom.intensities(spectrum))
end

function Makie.plot!(ms::MassSpectrum)
    ax = current_axis()

    Makie.add_input!(ms.attributes, :viewport_obs, ax.scene.viewport)
    Makie.add_input!(ms.attributes, :limits_obs, ax.finallimits)

    massspectrum_register_limit_updates!(ax, ms)

    input_nodes = [
        :mzvalues,
        :intensities,
        :mz_label_font,
        :mz_label_fontsize,
        :mz_label_offset,
        :mz_label_padding,
        :mz_label_top_padding_factor,
        :intensity_threshold,
        :viewport_obs,
        :limits_obs
    ]
    output_nodes = [:linesegments, :mz_label_positions, :mz_labels]

    map!(ms.attributes, input_nodes, output_nodes) do mzvalues, intensities, _label_font,
            _label_fontsize, _label_offset, _label_padding, _label_top_padding_factor,
            _intensity_threshold, viewport, limits
        mz_values, mz_unit = massspectrum_numeric_values(mzvalues)
        intensity_values, intensity_unit = massspectrum_numeric_values(intensities)
        massspectrum_set_axis_labels!(ax, mz_unit, intensity_unit)
        plotlimits = ax.finallimits[]

        mz_label_positions, mz_labels = calculate_visible_labels(
            ax,
            ms,
            mz_values,
            intensity_values,
            intensity_unit,
            viewport,
            plotlimits
        )

        linesegments = Point2f[]
        for (mz, intensity) in zip(mz_values, intensity_values)
            push!(linesegments, Point2f(mz, 0.0))
            push!(linesegments, Point2f(mz, intensity))
        end
        return (linesegments, mz_label_positions, mz_labels)
    end

    linesegments!(ms, ms.linesegments; color=ms.linecolor, linewidth=ms.linewidth)
    text!(ms, ms.mz_label_positions;
            text = ms.mz_labels,
            font = ms.mz_label_font,
            fontsize = ms.mz_label_fontsize,
            color = ms.mz_label_color,
            align = (:center, :bottom),
            offset = ms.mz_label_offset
        )

    return ms
end

function massspectrum_register_limit_updates!(ax::Axis, ms::MassSpectrum)
    Makie.Observables.onany(
        ms,
        ms.mzvalues,
        ms.intensities,
        ms.mz_label_font,
        ms.mz_label_fontsize,
        ms.mz_label_offset,
        ms.mz_label_padding,
        ms.mz_label_top_padding_factor,
        ms.intensity_threshold,
        ax.scene.viewport;
        update=true
    ) do mzvalues, intensities, _label_font, _label_fontsize, _label_offset,
            _label_padding, _label_top_padding_factor, _intensity_threshold, viewport
        mz_values, _mz_unit = massspectrum_numeric_values(mzvalues)
        intensity_values, intensity_unit = massspectrum_numeric_values(intensities)
        massspectrum_set_limits!(
            ax,
            mz_values,
            intensity_values,
            ms,
            intensity_unit,
            viewport
        )
        nothing
    end

    ms
end

function calculate_visible_labels(
    ax,
    ms,
    mzvalues,
    intensities,
    intensity_unit,
    viewport,
    limits
)
    threshold = massspectrum_threshold_value(ms.intensity_threshold[], intensity_unit)
    # Calculate label candidates
    label_candidates = [(mz=mz, intensity=intensity, label=string(mz)) for (mz, intensity) 
        in zip(mzvalues, intensities) if intensity ≥ threshold]
    
    sort!(label_candidates, by=x->x.intensity, rev=true)

    # Find non-overlapping labels
    visible_labels = String[]
    visible_positions = Point2f[]
    occupied_regions = []
    
    # Add vertical lines as occupied regions
    for (mz, int) in zip(mzvalues, intensities)
        x0, y0 = data2pixel(ax, viewport, limits, mz, 0)
        x1, y1 = data2pixel(ax, viewport, limits, mz, int)
        push!(occupied_regions, (x0, x1, y0, y1, :line))
    end

    minint = onepixelhight(viewport, limits)
    for label in label_candidates
        if label.intensity > minint
            x0p, x1p, y0p, y1p = covered_pixel_area(ax, ms, viewport, limits, label.label, 
                label.mz, label.intensity)

            overlaps = any(occupied_regions) do region
                x0, x1, y0, y1, type = region
                x_overlap = (x0p ≤ x0 ≤ x1p || x0p ≤ x1 ≤ x1p)
                
                if type ≡ :line
                    x_overlap && (y0 ≤ y0p ≤ y1)  # For lines: label must be within line's Y range
                else
                    x_overlap && (y0p ≤ y0 ≤ y1p || y0p ≤ y1 ≤ y1p)  # For labels: normal overlap
                end
            end
            
            if !overlaps
                massspectrum_yspan_inside_viewport(y0p, y1p, viewport) || continue

                push!(visible_labels, label.label)
                push!(visible_positions, Point2f(label.mz, label.intensity))  # Use intensity, not text_y
                push!(occupied_regions, (x0p, x1p, y0p, y1p, :label))
            end
        end
    end
    
    visible_positions, visible_labels
end

function data2pixel(ax::Axis, viewport, limits, x::Real, y::Real)
    # Calculate pixels per data unit
    ppd_x = viewport.widths[1] / limits.widths[1]
    ppd_y = viewport.widths[2] / limits.widths[2]
    
    # Transform: (data - data_origin) * scale + pixel_origin
    x_pixel = (x - limits.origin[1]) * ppd_x
    y_pixel = (y - limits.origin[2]) * ppd_y
    
    # Return the pixel coordinates
    (x_pixel, y_pixel)
end

function covered_pixel_area(ax, ms, viewport, limits, label, x_data, y_data)
    border = ms.mz_label_padding[]
    offset = ms.mz_label_offset[]
    box = Makie.text_bb(label, ms.mz_label_font[], ms.mz_label_fontsize[])
    width, height = box.widths[1], box.widths[2]
    x_pixel, y_pixel = data2pixel(ax, viewport, limits, x_data, y_data)
    x0 = x_pixel + offset[1] - border - width/2
    x1 = x_pixel + offset[1] + border + width/2
    y0 = y_pixel + offset[2] - border
    y1 = y_pixel + offset[2] + border + height
    (x0, x1, y0, y1)
end

onepixelhight(viewport, limits) = limits.widths[2] / viewport.widths[2]

function massspectrum_set_limits!(
    ax::Axis,
    mzvalues,
    intensities,
    ms,
    intensity_unit,
    viewport
)
    ymax = massspectrum_ymax(intensities)
    if ymax ≤ 0
        ylims!(ax, 0, 1)
        return ax
    end

    headroompixels = massspectrum_label_headroom_pixels(
        ms,
        mzvalues,
        intensities,
        intensity_unit
    )
    ydatafraction = massspectrum_data_yfraction(viewport, headroompixels)
    ylimit = ydatafraction > 0 ? ymax / ydatafraction : 1.1 * ymax
    ylims!(ax, 0, ylimit)

    ax
end

function massspectrum_ymax(intensities)
    isempty(intensities) && return 0.0
    ymax = maximum(intensities)

    isfinite(ymax) && ymax > 0 ? ymax : 0.0
end

function massspectrum_data_yfraction(viewport, headroompixels::Real)
    isnothing(viewport) && return 0.0
    viewportheight = Float64(viewport.widths[2])
    isfinite(viewportheight) && viewportheight > 0 || return 0.0

    clamp(
        (viewportheight - max(0.0, Float64(headroompixels))) / viewportheight,
        0.0,
        1.0
    )
end

function massspectrum_label_headroom_pixels(ms, mzvalues, intensities, intensity_unit)
    isempty(mzvalues) && return 0.0
    isempty(intensities) && return 0.0

    threshold = massspectrum_threshold_value(ms.intensity_threshold[], intensity_unit)
    labelindex = massspectrum_highest_label_index(intensities, threshold)
    isnothing(labelindex) && return 0.0

    offset = ms.mz_label_offset[]
    offset_y = length(offset) ≥ 2 ? Float64(offset[2]) : 0.0
    padding = max(0.0, Float64(ms.mz_label_padding[]))
    top_gap = massspectrum_label_top_gap_pixels(ms, offset_y, padding)
    box = Makie.text_bb(
        string(mzvalues[labelindex]),
        massspectrum_measurement_font(ms.mz_label_font[]),
        ms.mz_label_fontsize[]
    )
    labelheight = Float64(box.widths[2])

    max(0.0, offset_y + labelheight + top_gap)
end

function massspectrum_label_top_gap_pixels(ms, offset_y::Real, padding::Real)
    top_gap_factor = max(0.0, Float64(ms.mz_label_top_padding_factor[]))
    visible_lower_gap = max(0.0, Float64(offset_y))

    max(Float64(padding), top_gap_factor * visible_lower_gap)
end

function massspectrum_highest_label_index(intensities, threshold)
    bestindex = nothing
    bestintensity = -Inf
    for index in eachindex(intensities)
        intensity = intensities[index]
        intensity ≥ threshold || continue
        intensity > bestintensity || continue
        bestindex = index
        bestintensity = intensity
    end

    bestindex
end

function massspectrum_yspan_inside_viewport(y0::Real, y1::Real, viewport)
    viewportheight = Float64(viewport.widths[2])
    isfinite(viewportheight) && viewportheight > 0 || return false
    0 ≤ y0 && y1 ≤ viewportheight
end

function massspectrum_unitful_vectors(
    mzvalues::AbstractVector,
    intensities::AbstractVector,
    figure::NamedTuple,
    axis::NamedTuple;
    kwargs...
)
    mz_values, mz_unit = massspectrum_numeric_values(mzvalues)
    intensity_values, intensity_unit = massspectrum_numeric_values(intensities)
    plot_kwargs = massspectrum_numeric_plot_kwargs((; kwargs...), intensity_unit)
    axis_attributes = merge(massspectrum_axis_attributes(mz_unit, intensity_unit), axis)

    massspectrum(
        mz_values,
        intensity_values;
        figure=figure,
        axis=axis_attributes,
        plot_kwargs...
    )
end

function massspectrum_unitful_vectors!(
    ax::Axis,
    mzvalues::AbstractVector,
    intensities::AbstractVector;
    kwargs...
)
    mz_values, mz_unit = massspectrum_numeric_values(mzvalues)
    intensity_values, intensity_unit = massspectrum_numeric_values(intensities)
    plot_kwargs = massspectrum_numeric_plot_kwargs((; kwargs...), intensity_unit)
    massspectrum_set_axis_labels!(ax, mz_unit, intensity_unit)

    massspectrum!(ax, mz_values, intensity_values; plot_kwargs...)
end

function massspectrum_numeric_values(values)
    isempty(values) && return values, nothing
    firstvalue = first(values)
    firstvalue isa Unitful.AbstractQuantity || return values, nothing

    value_unit = Unitful.unit(firstvalue)
    Unitful.ustrip.(Ref(value_unit), values), value_unit
end

function massspectrum_threshold_value(threshold, intensity_unit)
    threshold isa Unitful.AbstractQuantity || return threshold
    isnothing(intensity_unit) && throw(ArgumentError(
        "unitful intensity_threshold requires unitful intensities"))

    Unitful.ustrip(intensity_unit, threshold)
end

function massspectrum_numeric_plot_kwargs(kwargs::NamedTuple, intensity_unit)
    :intensity_threshold in keys(kwargs) || return kwargs
    threshold = massspectrum_threshold_value(kwargs.intensity_threshold, intensity_unit)

    merge(kwargs, (; intensity_threshold=threshold))
end

function massspectrum_unit_suffix(unit; unitless::Bool=false)
    isnothing(unit) ? (unitless ? " [unitless]" : "") : " [$unit]"
end

function massspectrum_axis_attributes(mz_unit, intensity_unit)
    (
        xlabel="m/z" * massspectrum_unit_suffix(mz_unit),
        ylabel="Intensity" * massspectrum_unit_suffix(intensity_unit; unitless=true)
    )
end

function massspectrum_set_axis_labels!(ax::Axis, mz_unit, intensity_unit)
    attributes = massspectrum_axis_attributes(mz_unit, intensity_unit)
    massspectrum_set_axis_label!(ax.xlabel, attributes.xlabel)
    massspectrum_set_axis_label!(ax.ylabel, attributes.ylabel)

    ax
end

function massspectrum_set_axis_label!(label, default::AbstractString)
    isempty(label[]) && (label[] = default)

    label
end

function massspectrum_measurement_font(font::Symbol)
    Makie.to_font(Makie.MAKIE_DEFAULT_THEME.fonts, font)
end

function massspectrum_measurement_font(font)
    Makie.to_font(font)
end
