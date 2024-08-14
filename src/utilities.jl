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
    cosine(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Return the angle between two non-zero vectors, which can be considered a measure of the
`similarity` (i.e., `cosine` similarity) between the two vectors.

# Example
```jldoctest
julia> cosine([100, 500, 250], [200, 1000, 0])
0.8978872704229618

julia> cosine([100, 0, 50], [0, 20, 0])
0.0

julia> cosine([100, 500, 250], [10, 50, 25])
1.0
```
"""
function cosine(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || throw(
        DimensionMismatch("vectors x and y have different lengths"))
    length(x) == 0 && throw(ArgumentError("vectors x and y are empty"))
    iszero(x) && throw(ArgumentError("vector x contains only zeros"))
    iszero(y) && throw(ArgumentError("vector y contains only zeros"))
    s = sum(x .* y) / (sqrt(sum(x.^2)) * sqrt(sum(y.^2)))
    0 ≤ s ≤ 1 && return s
    (s < 0 || isnan(s)) ? zero(typeof(s)) : one(typeof(s))
end
## See also [`similarity`](@ref).


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

Return the `integer` for the given `value` that satisfies the following condition: 
`integer` - 1 + `start` ≤ `value` < `integer` + `start`, where 0 ≤ `start` < 1.

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
