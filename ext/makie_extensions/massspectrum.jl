import JuChrom: massspectrum, massspectrum!

@recipe(MassSpectrum, mzvalues, intensities) do scene
    Theme(
        mz_label_font = Makie.defaultfont(),
        mz_label_fontsize = 10,
        mz_label_offset = (0, 5),
        mz_label_padding = 1,
        mz_label_color = :black,
        linecolor = :blue,
        linewidth = 1.0,
        intensity_threshold = 0.0
    )
end

argument_names(::Type{<:MassSpectrum}) = (:mzvalues, :intensities)

function Makie.plot!(ms::MassSpectrum)
    ax = current_axis()

    Makie.add_input!(ms.attributes, :viewport_obs, ax.scene.viewport)
    Makie.add_input!(ms.attributes, :limits_obs, ax.finallimits)

    input_nodes = [:mzvalues, :intensities, :viewport_obs, :limits_obs]
    output_nodes = [:linesegments, :mz_label_positions, :mz_labels]

    map!(ms.attributes, input_nodes, output_nodes) do mzvalues, intensities, viewport, limits

        # Auto-adjust y-axis limits when data changes
        if !isempty(intensities)
            max_intensity = maximum(intensities)
            ylims!(ax, 0, max_intensity * 1.1)  # Add 10% padding
        end

        mz_label_positions, mz_labels = calculate_visible_labels(ax, ms, mzvalues, intensities, viewport, limits)

        linesegments = Point2f[]
        for (mz, intensity) in zip(mzvalues, intensities)
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

function calculate_visible_labels(ax, ms, mzvalues, intensities, viewport, limits)
    threshold = ms.intensity_threshold[]
    # Calculate label candidates
    label_candidates = [(mz=mz, intensity=intensity, label=string(mz)) for (mz, intensity) 
        in zip(mzvalues, intensities) if intensity ≥ threshold]
    
    sort!(label_candidates, by=x->x.intensity, rev=true)

    # Find non-overlapping labels
    visible_labels = []
    visible_positions = []
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
                
                if type == :line
                    x_overlap && (y0 ≤ y0p ≤ y1)  # For lines: label must be within line's Y range
                else
                    x_overlap && (y0p ≤ y0 ≤ y1p || y0p ≤ y1 ≤ y1p)  # For labels: normal overlap
                end
            end
            
            if !overlaps
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
