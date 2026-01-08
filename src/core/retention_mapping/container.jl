abstract type AbstractRetentionMapper{A, B} end

struct RetentionMapper{
    T1<:AbstractVector{<:Real},
    T2<:Union{Nothing, Unitful.Units},
    T3<:Real,
    T4<:Real,
    T5<:Real,
    T6<:Real,
    T7<:AbstractVector{<:Real},
    T8<:Union{Nothing, Unitful.Units},
    T9<:Real,
    T10<:Real,
    T11<:Real,
    T12<:Real,
    T13<:Real,
    T14<:Real,
    T15<:AbstractVector{<:Real},
    T16<:AbstractVector{<:Real},
    T17<:Spline,
    } <: AbstractRetentionMapper{T2, T8}

    rA::T1
    rA_unit::T2
    rA_min::T3
    rA_max::T4
    rA_norm_min::T5
    rA_norm_max::T6
    rB::T7
    rB_unit::T8
    rB_min::T9
    rB_max::T10
    rB_pred_min::T11
    rB_pred_max::T12
    rB_pred_norm_min::T13
    rB_pred_norm_max::T14
    knots::T15
    coefs::T16
    spline::T17
    extras::Dict{String, Any}
end

Base.broadcastable(rm::RetentionMapper) = Base.RefValue(rm)

function Base.show(io::IO, rm::RetentionMapper)
    # Get basic information
    n_points = length(rm.rA)
    
    # Determine decimals for formatting
    decimals_rA = maxdecimals(rm.rA)
    decimals_rB = maxdecimals([rm.rB_min, rm.rB_max]) + 1
    
    # Format domain bounds
    fmt_rA_min = @sprintf("%.*f", decimals_rA, rm.rA_min)
    fmt_rA_max = @sprintf("%.*f", decimals_rA, rm.rA_max)
    fmt_rB_min = @sprintf("%.*f", decimals_rB, rm.rB_pred_min)
    fmt_rB_max = @sprintf("%.*f", decimals_rB, rm.rB_pred_max)
    
    # Extract units for display
    rA_unit_str = rm.rA_unit === nothing ? "unitless" : string(rm.rA_unit)
    rB_unit_str = rm.rB_unit === nothing ? "unitless" : string(rm.rB_unit)
    
    # Calculate residuals with proper unit handling
    if rm.rA_unit !== nothing
        rA_with_units = rm.rA * rm.rA_unit
    else
        rA_with_units = rm.rA
    end
    vals_pred = applymap.(rm, rA_with_units)
    if rm.rB_unit !== nothing
        vals_true = rm.rB * rm.rB_unit
    else
        vals_true = rm.rB
    end
    
    residuals = vals_pred .- vals_true
    min_residual, max_residual = extrema(abs.(residuals))
    avg_residual = mean(abs.(residuals))
    
    # Format residuals with units
    residual_fmt(val) = @sprintf("%.*f %s", 3, Unitful.ustrip(unit(val), val), string(unit(val)))
    
    # Check if metadata exists
    has_metadata = !isempty(rm.extras)
    
    # Print header
    point_word = n_points == 1 ? "point" : "points"
    println(io, "RetentionMapper with $n_points calibration $point_word")
    
    # Print domain information (input)
    println(io, "├─ Domain A (input):")
    println(io, "│  ├─ Range: $fmt_rA_min to $fmt_rA_max ($rA_unit_str)")
    println(io, "│  └─ Type: $(eltype(rm.rA))")
    
    # Print codomain information (output)
    println(io, "├─ Domain B (output):")
    println(io, "│  ├─ Range: $fmt_rB_min to $fmt_rB_max ($rB_unit_str)")
    println(io, "│  └─ Type: $(eltype(rm.rB))")
    
    # Print spline information
    println(io, "├─ Spline:")
    println(io, "│  ├─ Order: 4 (cubic)")
    println(io, "│  ├─ Knots: $(length(rm.knots))")
    println(io, "│  └─ Coefficients: $(length(rm.coefs))")
    
    # Print fit quality (adjust connector based on whether metadata follows)
    quality_prefix = has_metadata ? "├─" : "└─"
    println(io, "$quality_prefix Fit quality:")
    quality_sub_prefix = has_metadata ? "│" : " "
    println(io, "$quality_sub_prefix  ├─ Min residual: ", residual_fmt(min_residual))
    println(io, "$quality_sub_prefix  ├─ Avg residual: ", residual_fmt(avg_residual))
    println(io, "$quality_sub_prefix  └─ Max residual: ", residual_fmt(max_residual))
    
    # Print metadata using similar approach to ScanSeries
    if has_metadata
        println(io, "└─ Metadata:")
        metadata_pairs = collect(rm.extras)
        for (i, (key, value)) in enumerate(metadata_pairs)
            is_last = i == length(metadata_pairs)
            prefix = is_last ? "   └─" : "   ├─"
            continuation = is_last ? "    " : "   │"
            
            if ScanSeriesDisplay.is_complex_value(value)
                println(io, "$prefix $key:")
                ScanSeriesDisplay.print_complex_value(io, value, "$(continuation)  ", "$(continuation)  ", true)
            else
                value_str = ScanSeriesDisplay.format_value_with_wrap(string(value), 50)
                if '\n' in value_str
                    ScanSeriesDisplay.print_wrapped_value(io, "$prefix $key =", value_str, "$(continuation)     ")
                else
                    if is_last
                        print(io, "$prefix $key = $value_str")
                    else
                        println(io, "$prefix $key = $value_str")
                    end
                end
            end
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", mappers::Vector{<:RetentionMapper})
    n = length(mappers)
    println(io, "$n RetentionMapper$(n == 1 ? "" : "s"):")
    
    for (i, mapper) in enumerate(mappers)
        # Extract key info
        n_points = length(mapper.rA)
        rA_unit_str = mapper.rA_unit === nothing ? "unitless" : string(mapper.rA_unit)
        rB_unit_str = mapper.rB_unit === nothing ? "unitless" : string(mapper.rB_unit)
        
        # Calculate average residual quickly
        if mapper.rA_unit !== nothing
            rA_with_units = mapper.rA * mapper.rA_unit
        else
            rA_with_units = mapper.rA
        end
        vals_pred = applymap.(mapper, rA_with_units)
        vals_true = mapper.rB_unit !== nothing ? mapper.rB * mapper.rB_unit : mapper.rB
        avg_residual = mean(abs.(vals_pred .- vals_true))
        
        # Format output
        prefix = i == n ? "└─" : "├─"
        residual_str = @sprintf("%.3f", Unitful.ustrip(unit(avg_residual), avg_residual))
        
        println(io, "$prefix [$i] $n_points pts, $rA_unit_str→$rB_unit_str, " *
                   "avg residual: $residual_str")
    end
end

# Also add the MIME show method for consistency
function Base.show(io::IO, ::MIME"text/plain", rm::RetentionMapper)
    show(io, rm)
end

# ── extras ────────────────────────────────────────────────────────────────────────────────

"""
    extras(rm::AbstractRetentionMapper)

Retrieve the metadata dictionary associated with the retention mapper.

The metadata can store additional information about the mapper, such as the calibration 
file name, creation date, experimental conditions, or other user-defined properties 
that were stored during mapper creation or added subsequently.

## Arguments

- `rm::AbstractRetentionMapper`: The retention mapper from which to retrieve metadata.

## Returns

A `Dict{String, Any}` containing the metadata key-value pairs. Returns an empty 
dictionary if no metadata was stored.

## Notes

- Metadata is stored as a mutable dictionary and can be modified after mapper creation.
- Common metadata keys include "cal_file", "creation_date", "instrument", etc.
- The metadata does not affect the mathematical mapping operations.

See also [`fitmap`](@ref).

## Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]
       mapper = fitmap(retention_times, retention_indices, 
                       extras=Dict("cal_file" => "test.cal", "instrument" => "GC-MS"));

julia> extras(mapper)
Dict{String, Any} with 2 entries:
  "cal_file"   => "test.cal"
  "instrument" => "GC-MS"

julia> extras(mapper)["cal_file"]
"test.cal"

julia> extras(mapper)["analysis_date"] = "2025-06-28";

julia> extras(mapper)
Dict{String, Any} with 3 entries:
  "cal_file"      => "test.cal"
  "instrument"    => "GC-MS"
  "analysis_date" => "2025-06-28"
```
"""
extras(rm::AbstractRetentionMapper) = rm.extras

# ── rawretentions_A ───────────────────────────────────────────────────────────────────────

"""
    rawretentions_A(rm::AbstractRetentionMapper;
                    unit::Union{Nothing, Unitful.Units}=nothing)

Return the raw numeric calibration data values from the input domain (domain A), stripped 
of units.

This function always returns the numeric values without units, regardless of whether 
the mapper was created with unitful or unitless data.

## Arguments

- `rm::AbstractRetentionMapper`: The fitted retention mapper.
- `unit`: For unitful mappers, specifies the unit to convert to before stripping. 
  For unitless mappers, this parameter is ignored.

## Returns

A vector of raw numeric calibration input values (no units), optionally converted to the 
specified unit first.

## Throws

- `ArgumentError`: When attempting to specify a `unit` for a unitless mapper. The error
  message will be "Cannot convert unitless retention A to a unit".

See also [`retentions_A`](@ref), [`retentions_B`](@ref), [`rawretentions_B`](@ref), 
[`mapmin`](@ref), [`mapmax`](@ref), [`fitmap`](@ref).

## Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> rawretentions_A(mapper) ≈ [1.2, 2.5, 4.1, 6.8, 9.3]
true

julia> rawretentions_A(mapper, unit=u"s") ≈ [72.0, 150.0, 246.0, 408.0, 558.0]
true

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> rawretentions_A(mapper_unitless) ≈ [1.2, 2.5, 4.1, 6.8, 9.3]
true

julia> rawretentions_A(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A to a unit
```
"""
@inline function rawretentions_A(
    rm::AbstractRetentionMapper{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rA, unit, "retention A")
end

@inline function rawretentions_A(
    rm::AbstractRetentionMapper{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip(rm.rA, rm.rA_unit, unit)
end

# ── rawretentions_B ───────────────────────────────────────────────────────────────────────

"""
    rawretentions_B(rm::AbstractRetentionMapper;
                    unit::Union{Nothing, Unitful.Units}=nothing)

Return the raw numeric calibration data values from the output domain (domain B), 
stripped of units.

This function always returns the numeric values without units, regardless of whether 
the mapper was created with unitful or unitless data.

## Arguments

- `rm::AbstractRetentionMapper`: The fitted retention mapper.
- `unit`: For unitful mappers, specifies the unit to convert to before stripping. 
  For unitless mappers, this parameter is ignored.

## Returns

A vector of raw numeric calibration output values (no units), optionally converted to the 
specified unit first.

## Throws

- `ArgumentError`: When attempting to specify a `unit` for a unitless mapper. The error
  message will be "Cannot convert unitless retention B to a unit".

See also [`retentions_B`](@ref), [`retentions_A`](@ref), [`rawretentions_A`](@ref), 
[`mapmin`](@ref), [`mapmax`](@ref), [`fitmap`](@ref).

## Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> rawretentions_B(mapper) ≈ [100.0, 200.0, 300.0, 400.0, 500.0]
true

julia> rawretentions_B(mapper, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention B to a unit
```
"""
@inline function rawretentions_B(
    rm::AbstractRetentionMapper{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rB, unit, "retention B")
end

@inline function rawretentions_B(
    rm::AbstractRetentionMapper{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip(rm.rB, rm.rB_unit, unit)
end

# ── retentions_A ──────────────────────────────────────────────────────────────────────────

"""
    retentions_A(rm::AbstractRetentionMapper; unit::Union{Nothing, Unitful.Units}=nothing)

Return the calibration data values from the input domain (domain A) used to fit the 
retention mapper.

This function provides access to the original input values that were used during 
calibration. For unitful mappers, returns the values with units; for unitless mappers, 
returns the raw numeric values.

## Arguments

- `rm::AbstractRetentionMapper`: The fitted retention mapper.
- `unit`: Desired output unit (optional). If specified for unitful mappers, converts to 
  the requested unit.

## Returns

A vector of calibration input values. For unitful mappers, returns a vector of `Quantity` 
objects with appropriate units (converted if `unit` is specified). For unitless mappers, 
returns a vector of raw numeric values.

## Throws

- `ArgumentError`: When attempting to specify a `unit` for a unitless mapper. The error
  message will be "Cannot convert unitless retention A to a unit".

See also [`rawretentions_A`](@ref), [`retentions_B`](@ref), [`rawretentions_B`](@ref), 
[`mapmin`](@ref), [`mapmax`](@ref), [`fitmap`](@ref).

## Examples

```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> retentions_A(mapper) ≈ [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
true

julia> retentions_A(mapper, unit=u"s") ≈ [72.0, 150.0, 246.0, 408.0, 558.0]u"s"
true

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> retentions_A(mapper_unitless) ≈ [1.2, 2.5, 4.1, 6.8, 9.3]
true

julia> retentions_A(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A to a unit
```
"""
@inline function retentions_A(
    rm::AbstractRetentionMapper{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rA, unit, "retention A")
end

@inline function retentions_A(
    rm::AbstractRetentionMapper{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(rm.rA, rm.rA_unit, unit)
end

# ── retentions_B ──────────────────────────────────────────────────────────────────────────

"""
    retentions_B(rm::AbstractRetentionMapper;
    unit::Union{Nothing, Unitful.Units}=nothing)

Return the calibration data values from the output domain (domain B) used to fit the 
retention mapper.

This function provides access to the original output values that were used during 
calibration. For unitful mappers, returns the values with units; for unitless mappers, 
returns the raw numeric values.

## Arguments

- `rm::AbstractRetentionMapper`: The fitted retention mapper.
- `unit`: Desired output unit (optional). If specified for unitful mappers, converts to 
  the requested unit.

## Returns

A vector of calibration output values. For unitful mappers, returns a vector of `Quantity` 
objects with appropriate units (converted if `unit` is specified). For unitless mappers, 
returns a vector of raw numeric values.

## Throws

- `ArgumentError`: When attempting to specify a `unit` for a unitless mapper. The error
  message will be "Cannot convert unitless retention B to a unit".

See also [`rawretentions_B`](@ref), [`retentions_A`](@ref), [`rawretentions_A`](@ref), 
[`mapmin`](@ref), [`mapmax`](@ref), [`fitmap`](@ref).

## Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> retentions_B(mapper) ≈ [100.0, 200.0, 300.0, 400.0, 500.0]
true

julia> retentions_B(mapper, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention B to a unit
```
"""
@inline function retentions_B(
    rm::AbstractRetentionMapper{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rB, unit, "retention B")
end

@inline function retentions_B(
    rm::AbstractRetentionMapper{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(rm.rB, rm.rB_unit, unit)
end

# ── retentionunit_A ───────────────────────────────────────────────────────────────────────

"""
    retentionunit_A(mapper::AbstractRetentionMapper)

Return the unit associated with the input domain (domain A) of the retention mapper.

For unitful mappers, this returns the `Unitful.Units` object representing the unit 
of the input values. For unitless mappers, this returns `Nothing`.

## Arguments

- `mapper::AbstractRetentionMapper`: The fitted retention mapper.

## Returns

The unit of the input domain. Returns a `Unitful.Units` object for unitful mappers, 
or `Nothing` for unitless mappers.

See also [`retentions_A`](@ref), [`rawretentions_A`](@ref), [`fitmap`](@ref).

## Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> retentionunit_A(mapper)
minute

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> retentionunit_A(mapper_unitless) === nothing
true
```
"""
retentionunit_A(mapper::AbstractRetentionMapper) = mapper.rA_unit

# ── retentionunit_B ───────────────────────────────────────────────────────────────────────

"""
    retentionunit_B(mapper::AbstractRetentionMapper)

Return the unit associated with the output domain (domain B) of the retention mapper.

For unitful mappers, this returns the `Unitful.Units` object representing the unit 
of the output values. For unitless mappers, this returns `Nothing`.

## Arguments

- `mapper::AbstractRetentionMapper`: The fitted retention mapper.

## Returns

The unit of the output domain. Returns a `Unitful.Units` object for unitful mappers, 
or `Nothing` for unitless mappers.

See also [`retentions_B`](@ref), [`rawretentions_B`](@ref), [`fitmap`](@ref).

## Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> retentionunit_B(mapper) === nothing
true
```
"""
retentionunit_B(mapper::AbstractRetentionMapper) = mapper.rB_unit
