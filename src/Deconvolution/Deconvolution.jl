module Deconvolution

import NNLS: nnls
using JuChrom
using Unitful
import Statistics

export LocalMaxima
export lsfit
export nextlocalmaximum
export stddev


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
(true, true)
```
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


"""
    JuChrom.stddev(chrom::AbstractChromMS; windowsize::Integer=13, threshold::Real=0,
    nthreads::Integer=Threads.nthreads())::Tuple{Union{Float64, Nothing}, Int}

Return an estimate of the standard deviation (σ) in the intensity measurements of the 
instrument used to infer the chromatographic data by analyzing intensity fluctuations 
within user-defined scan windows (default: 13 scans). The function also returns the number 
of windows considered in the computation as a second output. This approach is based on the 
method described by Stein (1999) for estimating noise levels in chromatographic data, with 
some modifications. In brief, only windows with intensity values above the specified 
threshold (default: 0) are included. Additionally, the total number of transitions must 
exceed half the window's scan count, and no more than two consecutive intensity values can 
fall on the same side of the mean intensity within any given window. For windows that meet 
these criteria, the absolute differences between intensity values and the median intensity 
are used to calculate the median absolute deviation (MAD). Since intensity fluctuations are 
proportional to the square root of the measured intensity, the MAD is normalized by 
dividing it by the square root of the median intensity. The MAD is calculated for all 
possible non-overlapping windows across the intensity values for all ions. The median of 
the MAD values is then multiplied by 1.4826 to estimate σ, assuming that the 
intensity-normalized fluctuations follow a normal distribution. If no suitable windows are 
found for calculating σ, the function returns (`nothing`, 0).

# Reference
Stein SE (1999): An integrated method for spectrum extraction and compound identification 
from gas chromatography/mass spectrometry data. J. Am. Soc. Mass. Spectrom. 10: 770–781.

# Examples
```jldoctest
julia> dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D");

julia> chrom = binions(importdata(dfolder, ChemStationMS()));

julia> JuChrom.stddev(chrom) .≈ (1.989116630064713, 9874)
(true, true)

julia> dfolder = joinpath(JuChrom.agilent, "C7-C40_MassHunterMS.D");

julia> chrom = binions(importdata(dfolder, MassHunterMS()));

julia> JuChrom.stddev(chrom) .≈ (1.951017182665152, 10979)
(true, true)
```
"""
function stddev(chrom::AbstractChromMS; windowsize::Integer=13, 
    threshold::Real=0, nthreads::Integer=Threads.nthreads()
    )::Tuple{Union{Float64, Nothing}, Int}
    nthreads < 1 && throw(ArgumentError("the number of threads must be one or more"))
    nthreads > Threads.nthreads() && throw(
        ArgumentError("the number of threads exceeds the maximum available"))
    n = ioncount(chrom) ≥ nthreads ? nthreads : ioncount(chrom)
    chunks = Iterators.partition(1:ioncount(chrom), ioncount(chrom) ÷ n)
    tasks = map(chunks) do ionindices
        Threads.@spawn stddev(chrom, ionindices, windowsize=windowsize, threshold=threshold)
    end
    chunk_lms = fetch.(tasks)
    mads = collect(Iterators.flatten(chunk_lms))
    length(mads) > 1 ? (1.4826 * Statistics.median(mads), length(mads)) : (nothing, 0)
end


function stddev(chrom::AbstractChromMS, ionindices; windowsize::Integer=13, 
    threshold::Real=0)
    windowsize > 4 || throw(ArgumentError("window size must be greater than four scans"))
    windowsize ≤ scancount(chrom) || throw(
        ArgumentError("window size must not exceed the scan count"))
    mads = Vector{Float64}()
    for iᵢ in ionindices
        istart = 1
        while istart ≤ scancount(chrom) - windowsize + 1
            istop = istart + windowsize - 1
            x̄::Float64 = Statistics.mean(@view intensities(chrom)[istart:istop, iᵢ])
            x̃::Float64 = Statistics.median(@view intensities(chrom)[istart:istop, iᵢ])
            level = :na
            next = false
            mintransitions = floor(windowsize / 2) + 1
            abovecount = belowcount = transitions = 0
            absdiff = Vector{Float64}()
            for (scancount, iₛ) in enumerate(istart:istop)
                x = intensities(chrom)[iₛ, iᵢ]

                if x ≤ threshold
                    istart = iₛ + 1
                    next = true
                    break
                end
                
                # To increase the likelihood that only random fluctuations are considered, 
                # transitions within the window are measured relative to the mean instead 
                # of the median. Additionally, no more than two consecutive signal 
                # intensities may occur on the same side of the mean.

                # transition :above -> :below
                if level == :above && x < x̄
                    abovecount = 0
                    transitions += 1
                    level = :below
                    belowcount += 1
                # transition :below -> :above
                elseif level == :below && x > x̄
                    belowcount = 0
                    transitions += 1
                    level = :above
                    abovecount += 1
                # remains :above
                elseif level == :above && x ≥ x̄
                    abovecount += 1
                # remains :below
                elseif level == :below && x ≤ x̄
                    belowcount += 1
                # starts :above
                elseif level == :na && x > x̄
                    level = :above
                    abovecount += 1
                    abovecount = 1
                # starts :below
                elseif level == :na && x < x̄
                    level = :below
                    belowcount += 1
                    belowcount = 1
                end

                # Ensure that no more than two consecutive signal intensities occur on the 
                # same side of the mean!
                if abovecount == 3 || belowcount == 3
                    istart = iₛ + 1
                    next = true
                    break
                end

                # Ensure that the transition count is at least half the window size
                if windowsize - scancount < mintransitions - transitions
                    istart = iₛ + 1
                    next = true
                    break
                end

                # Although transitions are measured relative to the mean, absolute 
                # differences are recorded based on the median to enable the calculation 
                # of the median absolute deviation (mad).
                push!(absdiff, abs(x - x̃))
            end
            next && continue

            istart = istop + 1

            push!(mads, Statistics.median(absdiff) / sqrt(x̃))
        end
    end
    mads
end


end  # module