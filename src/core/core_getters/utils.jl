# Core unit conversion helper functions

function _handle_unitless(vals, unit, field_display_name)
    if !isnothing(unit)
        throw(ArgumentError("Cannot convert unitless $(field_display_name) to a unit"))
    end
    vals
end

function _handle_unitful_convert(vals, stored_unit, unit)
    u = isnothing(unit) ? stored_unit : unit
    uconvert.(u, vals .* stored_unit)
end

function _handle_unitful_strip(vals, stored_unit, unit)
    u = isnothing(unit) ? stored_unit : unit
    ustrip.(u, vals .* stored_unit)
end

function _handle_unitful_convert_scalar(val, stored_unit, unit)
    u = isnothing(unit) ? stored_unit : unit
    uconvert(u, val * stored_unit)
end

function _handle_unitful_strip_scalar(val, stored_unit, unit)
    u = isnothing(unit) ? stored_unit : unit
    ustrip(u, val * stored_unit)
end