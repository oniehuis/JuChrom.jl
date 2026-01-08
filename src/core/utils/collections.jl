# ── copy_with ─────────────────────────────────────────────────────────────────────────────

"""
    copy_with(obj, overrides::NamedTuple) -> typeof(obj)

Create a deep copy of `obj`, replacing specified fields with new values.

This function returns a new object of the same type as `obj`, with selected fields 
replaced by values provided in the `overrides` named tuple. All other fields are 
deep-copied from the original object. The function assumes that all overridden values 
are type-compatible with the corresponding fields.

This function does **not** infer new field types and expects all overrides to match the 
original field types.

# Arguments
- `obj`: The original object to copy.
- `overrides`: A `NamedTuple` of field-value pairs to override in `obj`.

# Return
- A new object of the same type as `obj`, with specified fields updated.

# Example
```jldoctest
julia> struct MyStruct
         x::Int
         y::String
       end

julia> obj = MyStruct(1, "a")
MyStruct(1, "a")

julia> JuChrom.copy_with(obj, (y = "b",))
MyStruct(1, "b")

julia> JuChrom.copy_with(obj, (x = 42, y = "hello",))
MyStruct(42, "hello")
```
"""
function copy_with(obj, overrides::NamedTuple)
    # Get the concrete type of the input object
    T = typeof(obj)

    # Retrieve all field names of the object's type
    names = fieldnames(T)

    # Build a list of values for the new object:
    # For each field, if an override is provided, use it;
    # otherwise, use the original object's field value.
    values = [haskey(overrides, name) ? overrides[name] : deepcopy(getfield(obj, name))
              for name in names]

    # Construct a new instance of type T with the collected field values.
    # Note: This assumes all overridden values are type-compatible with the original fields.
    T(values...)
end

# ── typify ────────────────────────────────────────────────────────────────────────────────

"""
    typify(collection) -> new_collection

Promotes and converts the elements of a dictionary or vector to a common type.

For a dictionary, `typify(dict::Dict)` inspects the types of the values,
computes their common promoted type `T` using `promote_type`, and returns
a new `Dict{K, T}` with each value converted via `convert(T, v)`.

For a vector, `typify(vec::AbstractVector)` promotes the element types to
a common type `T` and returns a `Vector{T}` with all elements converted.

If the promoted type `T` is not concrete (e.g., a union or abstract type),
a warning is emitted, as this may lead to type instability.

# Arguments
- `dict::Dict`: A dictionary with arbitrary value types.
- `vec::AbstractVector`: A vector with potentially mixed element types.

# Returns
- A new container of the same kind (`Dict` or `Vector`) with unified and converted element types.

# Throws
- A `MethodError` if any element cannot be converted to the promoted type `T`.

# Examples
```julia
d = Dict("a" => 1, "b" => 2.5)
typify(d)
# => Dict{String, Float64}

v = [1, 2.5, 3]
typify(v)
# => Vector{Float64}

v = [1, "x"]
typify(v)
# => Vector{Union{Int, String}}, with warning
"""
function typify(dict::Dict)
    T = promote_type(typeof.(values(dict))...)
    if !isconcretetype(T)
        @warn "Type $T promoted by typify is not concrete."
    end
    Dict(k => convert(T, v) for (k, v) in dict)
end

function typify(vec::AbstractVector)
    T = promote_type(typeof.(vec)...)
    if !isconcretetype(T)
        @warn "Type $T promoted by typify is not concrete." maxlog=1
    end
    convert(Vector{T}, vec)
end
