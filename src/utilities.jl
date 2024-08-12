"""
    copy_with_eltype(array::AbstractArray, elementtype::Type)

Create a mutable copy of the `array` with the type of its elements converted to 
`elementtype`.

# Example
```jldoctest
julia> JuChrom.copy_with_eltype(Int[1, 2, 3, 4, 5, 6], Float64)
6-element Vector{Float64}:
 1.0
 2.0
 3.0
 4.0
 5.0
 6.0

julia> JuChrom.copy_with_eltype(Float64[1, 2, 3, 4, 5, 6], Int32)
6-element Vector{Int32}:
 1
 2
 3
 4
 5
 6

julia> JuChrom.copy_with_eltype(Float64[1.1, 2, 3, 4, 5, 6], Int32)
ERROR: InexactError: Int32(1.1)
[...]

```
"""
copy_with_eltype(array::AbstractArray, elementtype::Type) = copyto!(similar(array, 
    elementtype), array)