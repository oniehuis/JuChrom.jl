"""
    JuChrom.copy_with_eltype(array::AbstractArray, elementtype::Type)

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


"""
    JuChrom.findclosest(A::AbstractVector{<:Number}, x::Number) -> Int

Return the index of the number closest to number x in a list of numbers sorted in ascending 
order. If there is a tie, the index of the larger number is returned.

# Example
```jldoctest
julia> JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 0)
3

julia> JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 1.5)
5

julia> JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], -1.5)
2
```
"""
function findclosest(A::AbstractVector{<:Number}, x::Number)
    length(A) ≤ 1 && return firstindex(A)
    i = searchsortedfirst(A, x)
    if i == firstindex(A)
        return i
    elseif i > lastindex(A)
        return lastindex(A)
    else
        (x - A[i-1]) < (A[i] - x) ? (return i - 1) : (return i)
    end
end


"""
    integer(value:::Real; start::Real=0.7) -> Int

Return the integer for the given `value` that satisfies the following condition: 
``integer - 1 + start ≤ value < integer + start``, where ``0 ≤ start < 1``.

See also [`AbstractGCMS`](@ref), [`GCMS`](@ref), [`binions`](@ref), [`ions`](@ref).

# Example
```jldoctest
julia> integer(29.7)
30

julia> integer(30.0)
30

julia> integer(30.69)
30

julia> integer(29.7, start=0.8)
29
```
"""
function integer(value::Real; start::Real=0.7)
    0 ≤ start < 1 || throw(ArgumentError(string("fractional digits of binning interval ",
        "start is outside the interval [0,1) of allowed values: $start")))
    start == 0 && return trunc(Int, value)
    trunc(Int, value + (1 - start))
end
