# """
#     findclosest(A::AbstractArray{<:Number}, x::Number)

# Return the index of the number closest to number x in a list of numbers sorted in ascending 
# order. If there is a tie, the index of the larger number is returned.

# # Example
# ```jldoctest
# julia> findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 0)
# 3

# julia> findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 1.5)
# 5

# julia> findclosest([-2, -1, 0, 1, 2, 3, 4, 5], -1.5)
# 2
# ```
# """
function findclosest(A::AbstractArray{<:Number}, x::Number)
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

# """
#     JuChrom.copy_with_eltype(array::AbstractArray, elementtype::Type)

# Create a mutable copy of the `array` with the type of its elements converted to 
# `elementtype`.

# # Example
# ```jldoctest
# julia> JuChrom.copy_with_eltype(Int[1, 2, 3, 4, 5, 6], Float64)
# 6-element Vector{Float64}:
#  1.0
#  2.0
#  3.0
#  4.0
#  5.0
#  6.0

# julia> JuChrom.copy_with_eltype(Float64[1, 2, 3, 4, 5, 6], Int32)
# 6-element Vector{Int32}:
#  1
#  2
#  3
#  4
#  5
#  6

# julia> JuChrom.copy_with_eltype(Float64[1.1, 2, 3, 4, 5, 6], Int32)
# ERROR: InexactError: Int32(1.1)
# [...]
# ```
# """
copy_with_eltype(array::AbstractArray, elementtype::Type) = copyto!(similar(array, 
    elementtype), array)