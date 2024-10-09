module Deconvolution

include("../NNLS/NNLS.jl")
import .NNLS: nnls
using Unitful

export LocalMaxima
export lsfit
export nextlocalmaximum
export nnls


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


"""
    JuChrom.lsfit(xs::AbstractVector{<:Union{Real, Unitful.Quantity{<:Real}}}, 
    ys::AbstractVector{<:Union{Real, Unitful.Quantity{<:Real}}})

Return the intercept and slope from a simple linear regression, fitting the line of best 
fit (in the form `y = mx + b`) for the given data points `xs` and `ys` using the least 
squares method. The function expects two arguments: a vector of x-coordinates (independent 
variable), `xs`, and a vector of y-coordinates (dependent variable), `ys`. Both vectors 
must have the same length and contain at least two values. The function throws an 
`ArgumentError` if the lengths of `xs` and `ys` differ or if either vector contains fewer 
than two values.

The function calculates the slope (`m`) and intercept (`b`) for the line of best fit using 
the formulas

- `m = Σ((xᵢ - x̄)(yᵢ - ȳ)) / Σ((xᵢ - x̄)²)`
- `b = ȳ - m * x̄`

where `x̄` and `ȳ` are the means of the `xs` and `ys` vectors, respectively.

# Examples
```jldoctest
julia> xs, ys = [1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0];

julia> JuChrom.lsfit(xs, ys) .≈ (0.0, 1.0)
true
"""
function lsfit(xs::AbstractVector{<:Union{Real, Unitful.Quantity{<:Real}}}, 
    ys::AbstractVector{<:Union{Real, Unitful.Quantity{<:Real}}})
    n = length(xs)
    n ≤ 1 && throw(ArgumentError("each vector must consist of at least two values"))
    n ≠ length(ys) && throw(ArgumentError("vectors differ in size"))
    x̅, y̅ = (sum(xs), sum(ys)) ./ n
    numerator = zero(eltype(xs)) * zero(eltype(ys))
    denominator = zero(eltype(xs))^2
    for i in eachindex(xs)
        numerator += (xs[i] - x̅) * (ys[i] - y̅)
        denominator += (xs[i] - x̅)^2
    end
    slope = numerator / denominator
    y̅ - slope * x̅, slope
end


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

end  # module