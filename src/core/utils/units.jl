# ── consistentunits ───────────────────────────────────────────────────────────────────────

"""
    consistentunits(vals::AbstractArray{<:Union{Real, AbstractQuantity{<:Real}}}) -> Bool

Check whether all elements in `vals` have consistent Unitful units.

This function determines if an array contains only plain numbers, or only quantities with 
the same unit. Mixed arrays of numbers and quantities, or arrays of quantities with 
differing units, are considered inconsistent.

# Arguments
- `vals`: An array of values, which may be `Unitful.AbstractQuantity` or plain numbers.

# Return
- `Bool`: `true` if
    - `vals` is empty (trivially consistent),
    - all elements are plain numbers (no units), or
    - all elements are quantities with the **same** unit.
  Returns `false` if
    - there is a mixture of quantities and plain numbers, or
    - quantities have different units.

See also: [`strip_units_checked`](@ref)

# Examples
```jldoctest
julia> JuChrom.consistentunits(Int64[])
true

julia> JuChrom.consistentunits([1.0, 2.0, 3.0])
true

julia> JuChrom.consistentunits([1, 2, 3]u"s")
true

julia> ret = [1u"s", 2u"minute", 3u"s"];

julia> ret ≈ [1, 120, 3]u"s"
true

julia> JuChrom.consistentunits(ret)
true

julia> JuChrom.consistentunits([1u"m", 2.0u"minute", 3.0u"kg"])
false

julia> JuChrom.consistentunits([1u"s", 2.0, 3.0])
false
```
"""
function consistentunits(vals::AbstractArray{<:Union{Real, AbstractQuantity{<:Real}}})
    if isempty(vals)
        return true  # empty collection is trivially consistent
    elseif all(x -> x isa AbstractQuantity, vals)
        first_unit = unit(first(vals))
        return all(u -> u == first_unit, unit.(vals))
    elseif any(x -> x isa AbstractQuantity, vals)
        return false  # mixed units and plain numbers
    else
        return true   # all plain numbers
    end
end

# ── inverse ───────────────────────────────────────────────────────────────────────────────

"""
    inverse(u::Unitful.Units) -> Unitful.Units

Return the reciprocal (multiplicative inverse) of the given unit.

This function computes the inverse of a unit, returning a new unit representing `1/u`. 
This is useful for expressing rates or reciprocals of physical quantities (e.g. from 
meters to 1/meter).

# Arguments
- `u`: A unit from `Unitful`.

# Return
- `Unitful.Units`: The inverse unit, such as `1/u`.

# Example
```jldoctest
julia> JuChrom.inverse(u"m")
m⁻¹
julia> JuChrom.inverse(u"kg")
kg⁻¹
julia> JuChrom.inverse(u"s")
s⁻¹
```
"""
inverse(u::Unitful.Units) = unit(1.0 / u)

# ── isunitful ─────────────────────────────────────────────────────────────────────────────

"""
    isunitful(x) -> Bool

Return `true` if `x` is a `Unitful.AbstractQuantity`, regardless of whether it carries
nontrivial units. This is a light-weight check used to branch between unitful and
unitless code paths; values with `Unitful.NoUnits` still return `true` because they are
wrapped quantities.

# Arguments
- `x`: Any value that may or may not be a Unitful quantity.

# Return
- `Bool`: `true` for Unitful quantities, `false` otherwise.

# Example
```jldoctest
julia> JuChrom.isunitful(1.0u"s")
true

julia> JuChrom.isunitful(2.0)
false

julia> JuChrom.isunitful(ustrip(1.0u"s"))
false
```
"""
isunitful(x) = x isa Unitful.AbstractQuantity

# ── strip_units_checked ───────────────────────────────────────────────────────────────────

"""
    strip_units_checked(arr::AbstractArray, name::String) -> Tuple{AbstractArray, Union{Unitful.Units, Nothing}}

Strip units from an array, returning both the unitless values and the associated unit.

This function checks that all elements of `arr` have consistent units (if any), strips the 
units if present, and returns a tuple `(values, unit)`. If the array is empty, it returns 
the array and `nothing`. If the array contains quantities with inconsistent units, an 
`ArgumentError` is thrown. If the array contains only numbers (no units), the original 
array and `nothing` are returned.

# Arguments
- `arr`: An array of numbers or `Unitful.AbstractQuantity` values.
- `name`: A string used for error messages to identify the data source.

# Return
- A tuple `(values, unit)`, where `values` is an array of numbers (units stripped if 
  present), and `unit` is the associated `Unitful.Units` object or `nothing` if not 
  applicable.

# Throws
- `ArgumentError` if the array contains quantities with inconsistent units.

See also: [`consistentunits`](@ref)

# Examples
julia> JuChrom.strip_units_checked([1, 2]u"s", "times")
([1, 2], s)

julia> JuChrom.strip_units_checked([1.0, 2.0], "values")
([1.0, 2.0], nothing)

julia> JuChrom.strip_units_checked([], "empty")
(Any[], nothing)

julia> JuChrom.strip_units_checked([1u"s", 2u"m"], "mixed")
ERROR: ArgumentError: Inconsistent units in mixed.
[...]
"""
function strip_units_checked(mat::AbstractArray, name::String)
    isempty(mat) && return mat, nothing

    if eltype(mat) <: AbstractQuantity
        consistentunits(mat) || throw(ArgumentError("Inconsistent units in $name."))
        return ustrip.(mat), unit(first(mat))
    else
        return mat, nothing
    end
end

# ── Thomson ───────────────────────────────────────────────────────────────────────────────

# # Julia ≥ 1.12 enforces stricter world-age semantics for globals created via Core.eval.
# # Predefine Unitful's hidden base-factor dictionary so Unitful.register can see this
# # binding without emitting warnings about accessing it before creation.
# if !isdefined(@__MODULE__, Symbol("#Unitful_basefactors"))
#     @eval const $(Symbol("#Unitful_basefactors")) =
#         Dict{Symbol, Tuple{Float64, Rational{Int}}}()
# end

@unit Th "Th" Thomson (Unitful.u / Unitful.q) true  # u = unified atomic mass unit, q = e = elementary charge, Th = u / e
