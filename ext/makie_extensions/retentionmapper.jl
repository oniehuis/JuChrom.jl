# Constants
using Unitful
const DEFAULT_PLOT_SIZE = (1600, 800)
const DEFAULT_LINE_COLOR = :blue
const DEFAULT_LINE_WIDTH = 1.0
const DEFAULT_MARKER_SIZE = 10.0
const DEFAULT_MARKER_COLOR = :orange
const RANGE_EXTENSION_FACTOR = 0.05
const FINE_GRID_POINTS = 1000
const DEFAULT_DIGITS = 1
const DEFAULT_FONT_SIZE = 16

# Helper functions
function format_unit_string(unit)
    unit === nothing ? "" : " [$unit]"
end

function create_fine_range(min_val, max_val, extension_factor=RANGE_EXTENSION_FACTOR, points=FINE_GRID_POINTS)
    range_ext = (max_val - min_val) * extension_factor
    LinRange(min_val - range_ext, max_val + range_ext, points)
end

function format_derivative_unit_string(rA_unit, rB_unit)
    if rA_unit === nothing && rB_unit === nothing
        return ""
    elseif rA_unit === nothing
        return " [$rB_unit]"
    elseif rB_unit === nothing
        return " [$(JuChrom.inverse(rA_unit))]"
    else
        return " [$rB_unit × $(JuChrom.inverse(rA_unit))]"
    end
end

function format_inverse_derivative_unit_string(rA_unit, rB_unit)
    if rA_unit === nothing && rB_unit === nothing
        return ""
    elseif rB_unit === nothing
        return " [$rA_unit]"
    elseif rA_unit === nothing
        return " [$(JuChrom.inverse(rB_unit))]"
    else
        return " [$rA_unit × $(JuChrom.inverse(rB_unit))]"
    end
end

function create_retention_axis(fig, position, title, xlabel, ylabel)
    Axis(fig[position...], title=title, xlabel=xlabel, ylabel=ylabel)
end

function add_retention_plot!(ax, x_data, y_data, x_scatter, y_scatter, linecolor, linewidth, markercolor, markersize)
    scatter!(ax, x_scatter, y_scatter, color=markercolor, markersize=markersize, strokecolor=:black, strokewidth=0.5)
    lines!(ax, x_data, y_data, color=linecolor, linewidth=linewidth)
end

function add_derivative_plot!(ax, x_data, y_data, linecolor, linewidth)
    lines!(ax, x_data, y_data, color=linecolor, linewidth=linewidth)
end

function calculate_residuals(rm, retentions_B_pred, rB_unit, digits)
    actual_retentions_B = retentions_B(rm, unit=rB_unit)
    diff = retentions_B_pred .- actual_retentions_B
    # Unitful residuals: strip to numeric before rounding; unitless passes through.
    first_diff = first(diff)
    numeric_diff = first_diff isa Unitful.AbstractQuantity ? Unitful.ustrip.(diff) : diff
    round.(numeric_diff; digits=digits)
end

function determine_text_position(predicted_value, y_threshold, markersize)
    if predicted_value < y_threshold
        return ((:left, :center), (0, markersize))
    else
        return ((:right, :center), (0, -markersize))
    end
end

function add_residual_annotations!(ax, rm, retentions_B_pred, rA_unit, rB_unit, markersize, digits)
    residuals = calculate_residuals(rm, retentions_B_pred, rB_unit, digits)
    rA_raw = rawretentions_A(rm, unit=rA_unit)
    
    # Calculate y-threshold for text positioning
    y_values = [retentions_B_pred..., retentions_B(rm, unit=rB_unit)...]
    ymin, ymax = extrema(y_values)
    y_threshold = ymin + (ymax - ymin) / 2
    
    for (i, residual) in enumerate(residuals)
        residual == 0 && continue
        add_single_residual_annotation!(ax, rA_raw[i], retentions_B_pred[i], residual, y_threshold, markersize)
    end
end

function add_single_residual_annotation!(ax, x_pos, y_pos, residual, y_threshold, markersize)
    color = residual > 0 ? :black : :red
    align, offset = determine_text_position(y_pos, y_threshold, markersize)
    
    text!(ax, x_pos, y_pos,
          text="$residual",
          color=color,
          fontsize=DEFAULT_FONT_SIZE,
          rotation=π/2,
          align=align,
          offset=offset)
end

function validate_retention_mapper(rm)
    if isempty(retentions_A(rm))
        throw(ArgumentError("RetentionMapper has no data"))
    end
    
    if length(retentions_A(rm)) == 1
        @warn "Only one data point available, plot may not be meaningful"
    end
end

function calculate_plot_ranges(rm, rA_unit, rB_unit)
    r1_fine = create_fine_range(mapmin(rm, unit=rA_unit), mapmax(rm, unit=rA_unit))
    r2_fine = create_fine_range(invmapmin(rm, unit=rB_unit), invmapmax(rm, unit=rB_unit))
    return r1_fine, r2_fine
end

function add_forward_mapping_plot!(fig, rm, r1_fine, rA_unit, rB_unit, linecolor, linewidth, markercolor, markersize, digits)
    # Create axis
    A_unit_string = format_unit_string(rA_unit)
    B_unit_string = format_unit_string(rB_unit)
    
    ax = create_retention_axis(fig, (1, 1), 
                              "Retention A → Retention B",
                              "Retention A" * A_unit_string, 
                              "Retention B" * B_unit_string)
    
    # Add main plot
    y_data = rawapplymap.(rm, r1_fine, unit=rB_unit, warn=false)
    x_scatter = rawretentions_A(rm, unit=rA_unit)
    y_scatter = rawretentions_B(rm, unit=rB_unit)
    
    add_retention_plot!(ax, ustrip.(r1_fine), y_data, x_scatter, y_scatter, 
                       linecolor, linewidth, markercolor, markersize)
    
    # Add residual annotations
    retentions_B_pred = applymap.(rm, retentions_A(rm, unit=rA_unit), warn=false)
    add_residual_annotations!(ax, rm, retentions_B_pred, rA_unit, rB_unit, markersize, digits)
    
    return ax
end

function add_forward_derivative_plot!(fig, rm, r1_fine, rA_unit, rB_unit, linecolor, linewidth)
    # Create axis
    A_unit_string = format_unit_string(rA_unit)
    dBdA_unit_string = format_derivative_unit_string(rA_unit, rB_unit)
    
    ax = create_retention_axis(fig, (2, 1),
                              "Retention A → dB/dA",
                              "Retention A" * A_unit_string, 
                              "dB/dA" * dBdA_unit_string)
    
    # Add derivative plot
    y_data = rawderivmap.(rm, r1_fine, rA_unit=rA_unit, rB_unit=rB_unit, warn=false)
    add_derivative_plot!(ax, ustrip.(r1_fine), y_data, linecolor, linewidth)
    
    return ax
end

function add_reverse_mapping_plot!(fig, rm, r2_fine, rA_unit, rB_unit, linecolor, linewidth)
    # Create axis
    A_unit_string = format_unit_string(rA_unit)
    B_unit_string = format_unit_string(rB_unit)
    
    ax = create_retention_axis(fig, (1, 2), 
                              "Retention B → Retention A",
                              "Retention B" * B_unit_string, 
                              "Retention A" * A_unit_string)
    
    # Add reverse mapping plot
    y_data = rawinvmap.(rm, r2_fine, unit=rA_unit, warn=false)
    add_derivative_plot!(ax, ustrip.(r2_fine), y_data, linecolor, linewidth)
    
    return ax
end

function add_reverse_derivative_plot!(fig, rm, r2_fine, rA_unit, rB_unit, linecolor, linewidth)
    # Create axis
    B_unit_string = format_unit_string(rB_unit)
    dAdB_unit_string = format_inverse_derivative_unit_string(rA_unit, rB_unit)
    
    ax = create_retention_axis(fig, (2, 2),
                              "Retention B → dA/dB",
                              "Retention B" * B_unit_string, 
                              "dA/dB" * dAdB_unit_string)
    
    # Add inverse derivative plot
    y_data = rawderivinvmap.(rm, r2_fine)
    add_derivative_plot!(ax, ustrip.(r2_fine), y_data, linecolor, linewidth)
    
    return ax
end

function Makie.plot(rm::RetentionMapper;
    size=DEFAULT_PLOT_SIZE,
    linecolor=DEFAULT_LINE_COLOR,
    linewidth=DEFAULT_LINE_WIDTH, 
    markersize=DEFAULT_MARKER_SIZE,
    markercolor=DEFAULT_MARKER_COLOR,
    rA_unit::T1=retentionunit_A(rm),
    rB_unit::T2=retentionunit_B(rm),
    reverse::Bool=false,
    digits::Int=DEFAULT_DIGITS
    ) where {
    T1<:Union{Nothing, Unitful.Units},
    T2<:Union{Nothing, Unitful.Units}
    }
    
    # Validate inputs
    validate_retention_mapper(rm)
    
    # Create figure
    fig = Figure(; size=size)
    
    # Calculate plot ranges
    r1_fine, r2_fine = calculate_plot_ranges(rm, rA_unit, rB_unit)
    
    # Add forward mapping and derivative plots
    add_forward_mapping_plot!(fig, rm, r1_fine, rA_unit, rB_unit, 
                             linecolor, linewidth, markercolor, markersize, digits)
    add_forward_derivative_plot!(fig, rm, r1_fine, rA_unit, rB_unit, linecolor, linewidth)
    
    # Add reverse plots if requested
    if reverse
        add_reverse_mapping_plot!(fig, rm, r2_fine, rA_unit, rB_unit, linecolor, linewidth)
        add_reverse_derivative_plot!(fig, rm, r2_fine, rA_unit, rB_unit, linecolor, linewidth)
    end
    
    return fig
end
