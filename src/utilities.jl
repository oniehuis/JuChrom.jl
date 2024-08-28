import BSplineKit
import BSplineKit: extrapolate, Smooth
#import BSplineKit: Linear as pt


function rt2ri(ts, timeunit, S, E)
    rt::Unitful.Time -> begin
        t = ustrip(timeunit, rt)
        first(ts) ≤ t ≤ last(ts) ? S(t) : E(t)
    end
end


"""
    bsplineinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
    retentionindices::AbstractVector{<:Real};
    extrapolation::Union{Nothing, Symbol}=nothing, bsplineorder::Int=4)

Return a function that maps a retention time to a retention index using interpolation and, 
optionally, extrapolation, with a B-spline computed from a vector of retention times and a
corresponding vector of retention indices. For time values that fall outside the retention 
time range used to compute the B-spline, the function returns missing by default. The 
optional keyword argument `extrapolation` allows you to override this default behavior and 
specify an extrapolation method. One method is currently supported: :linear. The optional 
keyword argument `bsplineorder` allows you to change the default B-spline order value (4, 
which corresponds to cubic splines). For more details on the extrapolation method and the 
B-spline order, see the [BSplineKit.jl](https://jipolanco.github.io/BSplineKit.jl) 
documentation.

See also [`scantimes`](@ref), [`retentionindices`](@ref).

# Examples
```jldoctest
julia> rts, ris = [1, 2, 3, 4]u"s", [1000, 1900, 3100, 3900];

julia> rt2ri = bsplineinterpolation(rts, ris);

julia> rt2ri(1u"s") ≈ 1000.0
true

julia> rt2ri(1.5u"s") ≈ 1368.7499999999998
true

julia> rt2ri((1//30)u"minute") ≈ 1900
true

julia> rt2ri(5u"s")
missing

julia> rt2ri = bsplineinterpolation(rts, ris, extrapolation=:linear);

```
"""
function bsplineinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
    retentionindices::AbstractVector{<:Real};
    extrapolation::Union{Nothing, Symbol}=nothing, bsplineorder::Int=4)

    # Sanity checks
    bsplineorder > 1 || throw(ArgumentError(
        "B-spline order must be larger than 1"))
    length(retentiontimes) == length(retentionindices) || throw(ArgumentError(
        "number of retention times and number of retention indices are different"))
    length(retentiontimes) ≥ bsplineorder || throw(ArgumentError(string("the order of ", 
        "the B-spline cannot exceed the number of retention time-retention index pairs ", 
        "provided")))
    length(Set(retentiontimes)) == length(retentiontimes) || throw(ArgumentError(
        "retention times contain identical values"))
    issorted(retentiontimes) || throw(ArgumentError(
        "retention times not in ascending order"))
    length(Set(retentionindices)) == length(retentionindices) || throw(ArgumentError(
        "retention times contain identical values"))
    issorted(retentionindices) || throw(ArgumentError(
        "retention indices not in ascending order"))

    # Extract and store the retention time unit and strip it from the retention times
    timeunit = unit(eltype(retentiontimes))
    ts = ustrip.(retentiontimes)

    # Interpolation
    S = BSplineKit.interpolate(ts, retentionindices, BSplineKit.BSplineOrder(bsplineorder))

    # Extrapolation
    if isnothing(extrapolation)
       E = _ -> missing
    elseif extrapolation == :linear
       E = extrapolate(S, BSplineKit.SplineExtrapolations.Linear())
    else
       throw(ArgumentError("invalid extrapolation method: $extrapolation"))
    end

    rt2ri(ts, timeunit, S, E)
end


# julia> rt2ri = bsplineinterpolation(rts, ris, extrapolation=:linear);

# julia> rt2ri(5u"s") ≈ 4266.666666666666
# true


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
julia> gcms = GCMS(Int32[1, 2, 3]u"s", Int64[85, 100], Int32[0 12; 34 956; 23 1]);

julia> JuChrom.name(typeof(gcms))
GCMS
```
"""
name(::Type{T}) where {T} = (isempty(T.parameters) ? T : T.name.wrapper)


"""
    JuChrom.nextlocalmaximum(values::AbstractVector{<:Real}; 
    startindex::Integer=firstindex(values), stopindex::Integer=lastindex(values))

Return the index range of the next local maximum. If the maximum is a single peak, the 
range will have a length of 1; if a plateau maximum is found, the range will have a 
length greater than 1. If no maximum is found, nothing is returned. The optional keyword 
arguments `startindex` and `startindex` allow you to restrict the search range for the 
maximum. Note that `startindex` and `startindex` must be at least two indices apart (e.g., 
`startindex=1`, `stopindex=3`).

See also [`LocalMaxima`](@ref).

# Examples
```jldoctest
julia> values = [4, 3, 5, 3, 6, 6, 6, 4, 7];

julia> JuChrom.nextlocalmaximum(values)  # finds 5
3:3

julia> JuChrom.nextlocalmaximum(values, startindex=4)  # finds pleateau maximum 6, 6, 6
5:7

julia> JuChrom.nextlocalmaximum(values, startindex=4, stopindex=7)  # finds nothing

julia> length(values)  
9

julia> JuChrom.nextlocalmaximum(values, startindex=8)  # implicit stopindex=9
ERROR: ArgumentError: stopindex not greater than startindex + 1: 8 9
[...]

julia> JuChrom.nextlocalmaximum(values, startindex=1, stopindex=2)
ERROR: ArgumentError: stopindex not greater than startindex + 1: 1 2
[...]
```
"""
function nextlocalmaximum(values::AbstractVector{<:Real}; 
    startindex::Integer=firstindex(values), stopindex::Integer=lastindex(values))
    start = stop = 0
    firstindex(values) ≤ startindex ≤ lastindex(values) || throw(
        BoundsError(values, startindex))
    firstindex(values) ≤ stopindex ≤ lastindex(values) || throw(
        BoundsError(values, stopindex))
    startindex + 2 ≤ stopindex || throw(ArgumentError(string("stopindex not greater than ", 
        "startindex + 1: $startindex $stopindex")))
    startindex + 2 ≤ stopindex ≤ lastindex(values) || return nothing
    for i in startindex+1:stopindex
        if values[i-1] < values[i]
            start = stop = i
        elseif start != 0
            if values[i-1] > values[i]
                return start:stop
            elseif values[i-1] == values[i]
                stop = i
            end
        end
    end
end


struct LocalMaxima{T1<:AbstractVector{<:Real}, T2<:Integer, T3<:Integer, T4<:Integer}
    values::T1
    startindex::T2
    stopindex::T3
    length::T4
    function LocalMaxima{T1, T2, T3}(values::T1, startindex::T2, stopindex::T3
        ) where {T1<:AbstractVector{<:Real}, T2<:Integer, T3<:Integer}
        firstindex(values) ≤ startindex ≤ lastindex(values) || throw(
            BoundsError(values, startindex))
        firstindex(values) ≤ stopindex ≤ lastindex(values) || throw(
            BoundsError(values, stopindex))
        startindex + 2 ≤ stopindex || throw(ArgumentError(string("stopindex not greater ",
            "than startindex + 1: $startindex $stopindex")))
        new{T1, T2, T3, typeof(length(values))}(values, startindex, stopindex, 
            length(values))
    end
end


"""
    JuChrom.LocalMaxima(values::AbstractVector{<:Real}; 
    startindex::Integer=firstindex(values), stopindex::Integer=lastindex(values))

Return an iterator that yields the index range of each local maximum until no more local 
maxima are found. If the maximum is a single peak, the range will have a length of 1; if a 
plateau maximum is found, the range will have a length greater than 1. The optional keyword 
arguments `startindex` and `stopindex` allow you to restrict the search range for the 
maximum. Note that `startindex` and `stopindex` must be at least two indices apart (e.g., 
`startindex=1`, `stopindex=3`).

See also [`nextlocalmaximum`](@ref).

# Examples
```jldoctest
julia> values = [4, 3, 5, 3, 6, 6, 6, 4, 7];

julia> for lm in JuChrom.LocalMaxima(values); println(lm); end
3:3
5:7

julia> for lm in JuChrom.LocalMaxima(values, startindex=3); println(lm); end
5:7

julia> for lm in JuChrom.LocalMaxima(values, stopindex=7); println(lm); end
3:3

julia> for lm in JuChrom.LocalMaxima(values, startindex=4, stopindex=8); println(lm); end
5:7

julia> for lm in JuChrom.LocalMaxima(values, startindex=4, stopindex=5); println(lm); end
ERROR: ArgumentError: stopindex not greater than startindex + 1: 4 5
[...]

julia> length(values)
9

julia> for lm in JuChrom.LocalMaxima(values, startindex=8); println(lm); end
ERROR: ArgumentError: stopindex not greater than startindex + 1: 8 9
[...]
```
"""
function LocalMaxima(values::T1; startindex::T2=firstindex(values), 
    stopindex::T3=lastindex(values)) where {
    T1<:AbstractVector{<:Real},
    T2<:Integer, 
    T3<:Integer}
    LocalMaxima{T1, T2, T3}(values, startindex, stopindex)
end


function Base.iterate(LM::LocalMaxima, currentindex=LM.startindex)
    length(LM.values) == LM.length || throw(
        AssertionError("number of items changed since LocalMaxima was built"))
    if currentindex ≥ LM.stopindex - 1
        return nothing
    else
        lm = nextlocalmaximum(LM.values, startindex=currentindex, stopindex=LM.stopindex)
        isnothing(lm) ? (return nothing) : (return (lm, last(lm)+1))
    end
end
