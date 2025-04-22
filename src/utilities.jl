"""
    JuChrom.copy_with_eltype(array::AbstractArray, elementtype::Type)

Create a mutable copy of the `array` with the type of its elements converted to 
`elementtype`.

# Examples
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

Return the index of the number closest to `x` in a list `A` of numbers sorted in ascending 
order. If case of a tie, the index of the larger number is returned.

# Examples
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
    JuChrom.invert(dictionary::Dict)

Return a dictionary where the values from the input dictionary become the keys. Each key 
in the returned dictionary maps to a list of all original keys from dictionary that were 
associated with that value.

# Example
```jldoctest
julia> d = Dict(:a => 1.0, :b => 2.0, :c => 2.0, :d => 1.0)
Dict{Symbol, Float64} with 4 entries:
  :a => 1.0
  :b => 2.0
  :d => 1.0
  :c => 2.0

julia> JuChrom.invert(d)
Dict{Float64, Vector{Symbol}} with 2 entries:
  2.0 => [:b, :c]
  1.0 => [:a, :d]
```
"""
function invert(dictionary::Dict{T1, T2}) where {T1, T2}
    inverted_dictionary = Dict{T2, Vector{T1}}()
    for (key, value) in dictionary
        push!(get!(() -> T1[], inverted_dictionary, value), key)
    end
    inverted_dictionary
end


"""
    JuChrom.name(::Type)

Return the name of the type.

# Example
```jldoctest
julia> chrom = ChromMS(Int32[1, 2, 3]u"s", Int64[85, 100], Int32[0 12; 34 956; 23 1]);

julia> JuChrom.name(typeof(chrom))
ChromMS
```
"""
name(::Type{T}) where {T} = (isempty(T.parameters) ? T : T.name.wrapper)
