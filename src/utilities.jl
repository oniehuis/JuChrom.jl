# function naturalbesplinepolynomials(xs, ys)
#     n = length(xs) - 1
#     h = [xs[k+1] - xs[k] for k in 1:n]

#     # Preliminary definitions
#     In = I(n)
#     E = In[1:n-1, :]
#     Z = zeros(n, n)
#     J = diagm(0 => ones(n), 1 => -ones(n-1))
#     H = diagm(0 => h)

#     # Continuity of first derivative at internal nodes
#     A1 = E * [Z J 2H 3H^2]
#     v1 = zeros(n-1)

#     # Continuity of second derivative at internal nodes
#     A2 = E * [Z Z J 3H]
#     v2 = zeros(n-1)

#     # Left endpoint interpolation
#     AL = [In Z Z Z]
#     vL = ys[1:n]

#     # Right endpoint interpolation
#     AR = [In H H^2 H^3];
#     vR = ys[2:n+1]

#     # Boundary conditions: natural spline
#     nsL = [zeros(1, 2n) [2 zeros(1, 2n - 1)] ]
#     nsR = [zeros(1, 2n) [zeros(1, n - 1)  1] [zeros(1, n - 1) 3H[end]]]

#     # Assemble and solve the full system
#     A = [AL; AR; A1; A2; nsL; nsR]
#     v = [vL; vR; v1; v2; 0; 0]
#     z = A\v

#     # Break the coefficients into separate vectors
#     rows = 1:n
#     a = z[rows]
#     b = z[n .+ rows];  c = z[2n .+ rows];  d = z[3n .+ rows]
#     [Polynomial([a[k], b[k], c[k], d[k]]) for k in 1:n]
# end


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


const metricprefixes = (
    ("q", -30),
    ("r", -27),
    ("y", -24),
    ("z", -21),
    ("a", -18),            
    ("f", -15), 
    ("p", -12),
    ("n",  -9),
    ("µ",  -6),
    ("m",  -3),
    ("" ,   0),
    ("K",  +3),
    ("M",  +6),
    ("G",  +9),
    ("T", +12),
    ("P", +15),
    ("E", +18),
    ("Z", +21),
    ("Y", +24),
    ("R", +27),
    ("Q", +30))


"""
    JuChrom.metricprefix(number::Real) -> Tuple{String, Int}

Return the metric prefix and the corresponding exponent for a given number.

# Examples
```jldoctest
julia> JuChrom.metricprefix(347 * 10^7)
("G", 9)

julia> JuChrom.metricprefix(-42558335)
("M", 6)

julia> JuChrom.metricprefix(23)
("", 0)

julia> JuChrom.metricprefix(-0.00003623)
("µ", -6)

julia> JuChrom.metricprefix(425 * 10^-17)
("f", -15)

julia> JuChrom.metricprefix(big"136455347543543453485634875672536725423547234")
("Q", 30)

julia> JuChrom.metricprefix(Inf)
("", 0)

julia> JuChrom.metricprefix(0)
("", 0)

julia> JuChrom.metricprefix(NaN)
("", 0)
```
"""
function metricprefix(number::Real)
    (number == -Inf || number == 0 || number == +Inf || number === NaN) && return "", 0
    exponent = log10(abs(number))
    i = floor(Int, exponent / 3) + 11
    10.0^(floor(Int, exponent / 3) * 3) > abs(number) && (i -= 1)
    if i < 1
        metricprefixes[begin]
    elseif i > 21
        metricprefixes[end]
    else
        metricprefixes[i]
    end
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
