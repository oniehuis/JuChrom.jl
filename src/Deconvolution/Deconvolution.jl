module Deconvolution

using JuChrom
using Unitful
using CairoMakie

import BasicInterpolators: CubicSplineInterpolator
import Distributions: Normal, cquantile
import .JuChrom: NNLS
import Optim: minimizer, optimize
import QuadGK: quadgk
import Statistics: mean, median

export candidatepeaks
export componentbins
export componentpeakmodels
export deconvolutedtic
export deconvolutionwindows
export LocalMaxima
export lsfit
export massspectra
export nextlocalmaximum
export paraamdis
export saveplot
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
    length(mads) > 1 ? (1.4826 * median(mads), length(mads)) : (nothing, 0)
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
            x̄::Float64 = mean(@view intensities(chrom)[istart:istop, iᵢ])
            x̃::Float64 = median(@view intensities(chrom)[istart:istop, iᵢ])
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

            push!(mads, median(absdiff) / sqrt(x̃))
        end
    end
    mads
end


# """
#     JuChrom.paraamdis(chrom::AbstractChromMS; σ_signal::Real=stddev(chrom), q::Real=0.05,
#     rejectionthreshold::Real=4, baselinefraction::Real=0.5, dropfactor::Real=0.0, 
#     windowsize::Integer=scancount(chrom), ionscanorder::IonScanOrder=LinearDescending(), 
#     proportionalityfactor=10000, tic::Bool=false)

# Deconvolute.

# # Reference
# Stein SE (1999): An integrated method for spectrum extraction and compound identification 
# from gas chromatography/mass spectrometry data. J. Am. Soc. Mass. Spectrom. 10: 770–781.

# # Examples
# ```jldoctest
# julia> dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D");

# julia> chrom = binions(importdata(dfolder, ChemStationMS()));

# julia> mss = JuChrom.paraamdis(chrom)

# ```
# """
# function paraamdis(chrom::AbstractChromMS; σ_signal::Real=first(stddev(chrom)), 
#     q1::Real=0.01, q2::Real=0.01, baselinefraction::Real=0.5, dropfactor::Real=0.0, 
#     windowsize::Integer=scancount(chrom), ionscanorder::IonScanOrder=LinearDescending(), 
#     proportionalityfactor=10000, tic::Bool=false)

#     println("Identification of candidate peaks")
#     if tic
#        tic = totalionchromatogram(chrom)
#        @time cps = candidatepeaks(ChromMS(scantimes(tic), [0], reshape(intensities(tic), 
#            (:, 1))), σ_signal=σ_signal, q1=q1, q2=q2, baselinefraction=baselinefraction, 
#            dropfactor=dropfactor, windowsize=windowsize, ionscanorder=ionscanorder)
#     else
#        @time cps = candidatepeaks(chrom, σ_signal=σ_signal, q1=q1, q2=q2,
#             baselinefraction=baselinefraction, dropfactor=dropfactor, 
#             windowsize=windowsize, ionscanorder=ionscanorder)
#     end

#     println("Calculation of the deconvoluted total ion chromatogram")
#     @time bintimes, sharpnesssums, binpeaks = deconvolutedtic(chrom, cps)

#     println("Identification of compounds")
#     @time cbins = componentbins(sharpnesssums, proportionalityfactor=proportionalityfactor)

#     println("Selection of the best peak model")
#     @time peakmodels = componentpeakmodels(cbins, binpeaks)

#     println("Extracting mass spectra")
#     @time mss = massspectra(chrom, ionscanorder, peakmodels)    
# end


# """
#     JuChrom.candidatepeaks(chrom::AbstractChromMS; σ_signal::Real=stddev(chrom), 
#     q1::Real=0.05, q2::Real=0.001, baselinefraction::Real=0.5, dropfactor::Real=0.0, 
#     windowsize::Integer=scancount(chrom), ionscanorder::IonScanOrder=LinearDescending(), 
#     nthreads::Integer=Threads.nthreads())

# Generate a list of candidate peaks using the method outlined by Stein (1999) for component 
# perception, prior to identifying individual components, with some modifications.

# Keyword arguments:

# `σ_signal`: An estimate of the standard deviation (σ) in the instrument's intensity 
# measurements.

# `q1`: A probability value (default 0.05) used to specify the maximum allowable intensity 
# difference on one side of a local maximum, between the lowest intensity value and the 
# currently evaluated intensity. This intensity difference is calculated using the 
# expression `cquantile(Normal(0, σ_signal), q1) * sqrt(evaluated_intensity)`, where 
# `σ_signal` represents the standard deviation of random signal fluctuations, assuming these 
# fluctuations are normally distributed and proportional to the square root of the measured 
# signal. Setting `q1` to 0.05 means that the deconvolution window extension stops when 
# encountering an intensity difference that suggests a new peak because the previously 
# encountered lowest intensity could only have been measured with a 0.01 probability if the 
# current intensity represents the true baseline signal. The window boundary is then set to 
# the preceding scan. 

# Note: The algorithm evaluates the probability of the lowest value by considering the 
# current value as the true signal, rather than the other way around. This approach prevents 
# a premature stop of the window extension at scans with zero intensity, which would render 
# any intensity increase significant due to the multiplication of the cquantile term with 
# sqrt(0) (equal to zero).

# `q2`: A probability value (default 0.001) used to determine the minimum height a peak must 
# exceed above the estimated baseline to be considered a candidate. This height is calculated 
# using the formula `cquantile(Normal(0, σ_signal), q2) * sqrt(peak_max_intensity)`, where 
# `σ_signal` denotes the standard deviation of random signal fluctuations, which are assumed 
# to be normally distributed and proportional to the square root of the measured signal. 
# Setting `q2` to 0.001 indicates that the baseline intensity at the time of the peak's 
# maximum would be measured by chance with a probability of 0.001, assuming that the peak's 
# maximum intensity reflects the true signal.

# `baselinefraction`: The proportion of scans with intensity values that deviate minimally 
# from a tentative baseline, defined by the intensity values immediately to the left and 
# right of the local maximum. This baseline ensures that all scan intensity values are above 
# those predicted by it. The fraction of scans is used to compute a baseline against which 
# the local maximum is evaluated to determine if it represents a candidate peak.

# `dropfactor`: The fraction of the local maximum's intensity that a scan's intensity must 
# fall below for the deconvolution window extension in a given direction to stop.

# `windowsize`: The maximum number of scans included in a deconvolution window.

# `ionscanorder`: Object that defines the sequence in which ions were scanned.

# `nthreads`: The number of threads used for computation, with a default set to the maximum 
# available.

# See also [`JuChrom.AbstractChromMS`](@ref), [`JuChrom.stddev`](@ref), 
# [`JuChrom.IonScanOrder`](@ref).

# # Reference
# Stein SE (1999): An integrated method for spectrum extraction and compound identification 
# from gas chromatography/mass spectrometry data. J. Am. Soc. Mass. Spectrom. 10: 770–781.

# # Examples
# ```jldoctest
# julia> dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D");

# julia> chrom = binions(importdata(dfolder, ChemStationMS()));

# julia> JuChrom.paraamdis(chrom);

# ```
# """
# function candidatepeaks(chrom::AbstractChromMS; σ_signal::Real=first(stddev(chrom)), 
#     q1::Real=0.01, q2::Real=0.01, baselinefraction::Real=0.5, dropfactor::Real=0.0, 
#     windowsize::Integer=scancount(chrom), ionscanorder::IonScanOrder=LinearDescending(), 
#     nthreads::Integer=Threads.nthreads())

#     nthreads > Threads.nthreads() && throw(
#         ArgumentError("the number of threads exceeds the maximum available"))
#     δtᵢ = ionscantimeshift(chrom, ionscanorder)
#     n = ioncount(chrom) ≥ nthreads ? nthreads : ioncount(chrom)
#     chunks = Iterators.partition(1:ioncount(chrom), ioncount(chrom) ÷ n)
#     tasks = map(chunks) do chunk
#         Threads.@spawn candidatepeaks(chrom, chunk, σ_signal, q1, q2, baselinefraction, 
#             dropfactor, windowsize, δtᵢ)
#     end
#     chunk_lms = fetch.(tasks)
#     collect(Iterators.flatten(chunk_lms))
# end


# function candidatepeaks(chrom::AbstractChromMS, indices::AbstractVector{<:Integer}, 
#     σ_signal::Real, q1::Real, q2::Real, baselinefraction::Real, dropfactor::Real, 
#     windowsize::Real, δtᵢ::Function)
#     cps = Vector{CandidatePeak}()
#     for iᵢ in indices
#         ioncps = Vector{CandidatePeak}()
#         for lm in LocalMaxima(intensities(chrom)[:, iᵢ])
#             dw = deconvolutionwindow(chrom, δtᵢ, iᵢ, lm, windowsize, σ_signal, q1, 
#                 dropfactor)
#             isnothing(dw) && continue
#             cp = evaluatepeak(dw, σ_signal, q2, baselinefraction)
#             isnothing(cp) && continue
#             peaksharpness!(cp)
#             push!(ioncps, cp)
#         end
#         append!(cps, removeoverlaps(ioncps))
#     end
#     cps
# end


# function deconvolutionwindows(chrom::AbstractChromMS; σ_signal::Real=first(stddev(chrom)), 
#     q1::Real=0.01, dropfactor::Real=0.0, windowsize::Integer=scancount(chrom), 
#     ionscanorder::IonScanOrder=LinearDescending(), nthreads::Integer=Threads.nthreads())

#     nthreads > Threads.nthreads() && throw(
#         ArgumentError("the number of threads exceeds the maximum available"))
#     δtᵢ = ionscantimeshift(chrom, ionscanorder)
#     n = ioncount(chrom) ≥ nthreads ? nthreads : ioncount(chrom)
#     chunks = Iterators.partition(1:ioncount(chrom), ioncount(chrom) ÷ n)
#     tasks = map(chunks) do chunk
#         Threads.@spawn deconvolutionwindows(chrom, chunk, σ_signal, q1, dropfactor,
#             windowsize, δtᵢ)
#     end
#     chunk_lms = fetch.(tasks)
#     collect(Iterators.flatten(chunk_lms))
# end


# function deconvolutionwindows(chrom::AbstractChromMS, indices::AbstractVector{<:Integer}, 
#     σ_signal::Real, q1::Real, dropfactor::Real, windowsize::Real, δtᵢ::Function)
#     dws = Vector{DeconvolutionWindow}()
#     for iᵢ in indices, lm in LocalMaxima(intensities(chrom)[:, iᵢ])
#         dw = deconvolutionwindow(chrom, δtᵢ, iᵢ, lm, windowsize, σ_signal, q1, 
#             dropfactor)
#         isnothing(dw) && continue
#         push!(dws, dw)
#     end
#     dws
# end


# ############################################################################################
# # Code for delimiting DeconvolutionWindow (dw)
# ############################################################################################
# mutable struct DeconvolutionWindow{T<:AbstractChromMS}
#     const chromatogram::T                # chromatogram that the dw refers to
#     const δtᵢ::Function                  # ion scan time shift function
#     const ionindex::Int                  # ion index that the dw refers to
#     const maxindexrange::UnitRange{Int}  # scan index range of maximum
#     leftindex::Int                       # leftmost scan index of dw
#     rightindex::Int                      # rightmost scan index of dw
#     leftlowindex::Int                    # left index for drawing low value tangent
#     rightlowindex::Int                   # right index for drawing low value tangent
# end


# abstract type Direction end
# struct Forward <: Direction end
# struct Reverse <: Direction end


# onestep(::Forward) = +1
# onestep(::Reverse) = -1


# chromatogram(dw::DeconvolutionWindow) = dw.chromatogram
# δtᵢ(dw::DeconvolutionWindow) = dw.δtᵢ
# ionindex(dw::DeconvolutionWindow) = dw.ionindex
# maxindexrange(dw::DeconvolutionWindow) = dw.maxindexrange
# maxstartindex(dw::DeconvolutionWindow) = first(dw.maxindexrange)
# maxstopindex(dw::DeconvolutionWindow) = last(dw.maxindexrange)
# maxintensity(dw::DeconvolutionWindow) = intensity(chromatogram(dw), maxstartindex(dw), 
#     ionindex(dw))
# maxsize(dw::DeconvolutionWindow) = maxstopindex(dw) - maxstartindex(dw) + 1
# leftindex(dw::DeconvolutionWindow) = dw.leftindex
# rightindex(dw::DeconvolutionWindow) = dw.rightindex
# windowsize(dw::DeconvolutionWindow) = rightindex(dw) - leftindex(dw) + 1
# windowindices(dw::DeconvolutionWindow) = leftindex(dw):rightindex(dw)
# currentindex(dw::DeconvolutionWindow, ::Reverse) = dw.leftindex
# currentindex(dw::DeconvolutionWindow, ::Forward) = dw.rightindex
# currentindex!(dw::DeconvolutionWindow, ::Reverse, i::Integer) = dw.leftindex = i
# currentindex!(dw::DeconvolutionWindow, ::Forward, i::Integer) = dw.rightindex = i
# leftlowindex(dw::DeconvolutionWindow) = dw.leftlowindex
# rightlowindex(dw::DeconvolutionWindow) = dw.rightlowindex
# leftlowindex!(dw::DeconvolutionWindow, i::Integer) = dw.leftlowindex = i
# rightlowindex!(dw::DeconvolutionWindow, i::Integer) = dw.rightlowindex = i
# lowindex(dw::DeconvolutionWindow, ::Reverse) = dw.leftlowindex
# lowindex(dw::DeconvolutionWindow, ::Forward) = dw.rightlowindex
# lowindex!(dw::DeconvolutionWindow, ::Reverse, i::Integer) = dw.leftlowindex = i
# lowindex!(dw::DeconvolutionWindow, ::Forward, i::Integer) = dw.rightlowindex = i
# lowintensity(dw::DeconvolutionWindow, direction::Direction) = intensity(chromatogram(dw), 
#     lowindex(dw, direction), ionindex(dw))


# function deconvolutionwindow(chrom::AbstractChromMS, δtᵢ::Function, ionindex::Int, 
#         maxindexrange::UnitRange, windowsize::Integer, σ_signal::Real, q1::Real, 
#         dropfactor::Real)

#     # Initialize deconvolution window
#     dw = DeconvolutionWindow(
#         chrom, δtᵢ,            # general data
#         ionindex,              # ion index of local maximum
#         maxindexrange,         # scan index range of local maximum
#         first(maxindexrange),  # start index initially set to the local maximum start
#         last(maxindexrange),   # stop index initially set to the local maximum stop
#         first(maxindexrange),  # left low index initially set to the local maximum start
#         last(maxindexrange))   # right low index initially set to the local maximum stop

#     # Component perception sensu Stein (1999) p. 773, step a
#     # Try expanding the deconvolution window to both sides until the window size is reached
#     expand_f = expand_r = :true
#     s = maxsize(dw)
#     while s < windowsize
#         if expand_f == :true
#             expand_f = expand!(dw, σ_signal, q1, dropfactor, Forward())
#             expand_f == :abord && return nothing 
#             expand_f == :true && (s += 1)
#             s < windowsize || break
#         end
#         if expand_r == :true
#             expand_r = expand!(dw, σ_signal, q1, dropfactor, Reverse())
#             expand_r == :abord && return nothing 
#             expand_r == :true && (s += 1)
#         end
#         expand_f == :true || expand_r == :true || break
#     end

#     # Component perception sensu Stein (1999) p. 773, step b
#     # Adjust the low indices on either side of the maximum so that no scan intensity in the 
#     # in the deconvolution falls below a straight line drawn through these two low points
#     adjustlowvalindices!(dw)

#     dw
# end


# function expand!(dw::DeconvolutionWindow, σ_signal::Real, q1::Real, dropfactor::Real, 
#     direction::Direction)
#     iₛ = currentindex(dw, direction)

#     # Stop expansion when scan count limits are reached
#     1 < iₛ < scancount(chromatogram(dw)) || return :false

#     # Move to next scan index
#     iₛ += onestep(direction)

#     # Evaluate the intensity of the current scan relative the low and max intensities
#     int = intensity(chromatogram(dw), iₛ, ionindex(dw))
#     lowint = lowintensity(dw, direction)
#     maxint = maxintensity(dw)
#     if int < lowint
#         lowindex!(dw, direction, iₛ)
#         currentindex!(dw, direction, iₛ)
#         # Stop expansion when intensity falls below threshold
#         return int < dropfactor * maxint ? :false : :true
#     elseif int - lowint > cquantile(Normal(0, σ_signal), q1) * sqrt(int)
#         iₛ -= onestep(direction)
#         currentindex!(dw, direction, iₛ)
#         # Stop expansion because intensity has risen above threshold
#         return :false
#     # elseif int - lowint > (maxint - lowint) * 0.1  # NEW Peak must show peak characteristics
#     #     return :abord
#     else
#         currentindex!(dw, direction, iₛ)
#         # Continue expansion
#         return :true
#     end
# end


# function adjustlowvalindices!(dw::DeconvolutionWindow, iteration::Integer=1)
#     f = lowvalline(dw)
#     v₁ = v₂ = 0.0
#     i₁ = i₂ = 0
#     for iₛ in windowindices(dw)
#         iₛ == leftlowindex(dw) && continue
#         iₛ == rightlowindex(dw) && continue
#         t = ionscantime(δtᵢ(dw), chromatogram(dw), iₛ, ionindex(dw))
#         Δint = intensity(chromatogram(dw), iₛ, ionindex(dw)) - f(t)
#         if Δint < 0 
#             if iₛ < maxstartindex(dw) && Δint < v₁
#                 v₁ = Δint
#                 i₁ = iₛ
#             elseif iₛ > maxstopindex(dw) && Δint < v₂
#                 v₂ = Δint
#                 i₂ = iₛ
#             end
#         end
#     end
#     repeat = false
#     if v₁ < 0
#         leftlowindex!(dw, i₁)
#         repeat = true
#     end
#     if v₂ < 0
#         rightlowindex!(dw, i₂)
#         repeat = true
#     end
#     if repeat && iteration < windowsize(dw)
#         adjustlowvalindices!(dw, iteration + 1)
#     end
# end


# function lowvalline(dw::DeconvolutionWindow)
#     t₁ = ionscantime(δtᵢ(dw), chromatogram(dw), leftlowindex(dw), ionindex(dw))
#     t₂ = ionscantime(δtᵢ(dw), chromatogram(dw), rightlowindex(dw), ionindex(dw))
#     int₁ = intensity(chromatogram(dw), leftlowindex(dw), ionindex(dw))
#     int₂ = intensity(chromatogram(dw), rightlowindex(dw), ionindex(dw))
#     t -> (int₂ - int₁) / (t₂ - t₁) * (t - t₁) + int₁
# end


# ############################################################################################
# # Code for obtaining CandidatePeaks (cp)
# ############################################################################################
# mutable struct CandidatePeak{T1<:Unitful.Time, T2<:AbstractVector{T1}}
#     const starttime::T1             # RT of peak emergence
#     const stoptime::T1              # RT of peak disappearance
#     const startindex::Int           # Scan at which the peak emerges
#     const stopindex::Int            # Scan at which the peak disappears
#     const nzstartindex::Int         # Scan at which the first peak intensity was measured
#     const nzstopindex::Int          # Scan at which the last peak intensity was measured
#     const maxtime::T1               # RT of peak mode
#     const maxheight::Float64        # Height of the peak mode above the baseline
#     const maxintensity::Float64     # Absolute intensity of the peak mode
#     const model::Function           # Function that returns the peak intensity for any RT
#     const nodes::T2                 # RTs used with peak heights to compute cublic B-spline
#     const heights::Vector{Float64}  # Peak heights above baseline at the nodes 
#     const baseline::Function        # Baseline function used to infer the peak hights
#     const deconvolutionwindow::DeconvolutionWindow  # Data from which the peak is derived
#     sharpness                # Peak sharpness sensu Stein (1999), slightly modified
#     leftsharpnesstime        # RT when the left sharpness tangent converges with the peak
#     rightsharpnesstime       # RT when the right sharpness tangent converges with the peak
#     area                     # Area of the peak in a time unit
# end

# chromatogram(cp::CandidatePeak) = chromatogram(cp.deconvolutionwindow)
# δtᵢ(cp::CandidatePeak) = δtᵢ(cp.deconvolutionwindow)
# ionindex(cp::CandidatePeak) = ionindex(cp.deconvolutionwindow)
# model(cp::CandidatePeak) = cp.model
# maxtime(cp::CandidatePeak) = cp.maxtime
# nodes(cp::CandidatePeak) = cp.nodes
# maxintensity(cp::CandidatePeak) = cp.maxintensity
# maxheight(cp::CandidatePeak) = cp.maxheight
# sharpness(cp::CandidatePeak) = cp.sharpness
# sharpness!(cp::CandidatePeak, val) = cp.sharpness = val
# leftsharpnesstime(cp::CandidatePeak) = cp.leftsharpnesstime
# leftsharpnesstime!(cp::CandidatePeak, val) = cp.leftsharpnesstime = val
# rightsharpnesstime(cp::CandidatePeak) = cp.rightsharpnesstime
# rightsharpnesstime!(cp::CandidatePeak, val) = cp.rightsharpnesstime = val
# heights(cp::CandidatePeak) = cp.heights

# maxstartindex(cp::CandidatePeak) = maxstartindex(cp.deconvolutionwindow)
# maxstopindex(cp::CandidatePeak) = maxstopindex(cp.deconvolutionwindow)
# leftindex(cp::CandidatePeak) = leftindex(cp.deconvolutionwindow)
# rightindex(cp::CandidatePeak) = rightindex(cp.deconvolutionwindow)
# leftlowindex(cp::CandidatePeak) = leftlowindex(cp.deconvolutionwindow)
# rightlowindex(cp::CandidatePeak) = rightlowindex(cp.deconvolutionwindow)
# lowvalline(cp::CandidatePeak) = lowvalline(cp.deconvolutionwindow) 
# windowindices(cp::CandidatePeak) = windowindices(cp.deconvolutionwindow)
# maxindexrange(cp::CandidatePeak) = maxindexrange(cp.deconvolutionwindow)
# startindex(cp::CandidatePeak) = cp.startindex
# stopindex(cp::CandidatePeak) = cp.stopindex
# nzstartindex(cp::CandidatePeak) = cp.nzstartindex
# nzstopindex(cp::CandidatePeak) = cp.nzstopindex
# starttime(cp::CandidatePeak) = cp.starttime
# stoptime(cp::CandidatePeak) = cp.stoptime
# area(cp::CandidatePeak) = cp.area
# area!(cp::CandidatePeak, val) = cp.area = val


# function evaluatepeak(dw::DeconvolutionWindow, σ_signal::Real, q2::Real, 
#     baselineindexfraction::Real)

#     # Component perception sensu Stein (1999) p. 773, step c
#     # Infer a baseline
#     bl = baseline(dw, baselineindexfraction)

#     # Component perception sensu Stein (1999) p. 773, step d (modified)
#     # Substract the least-squares baseline from all intensities in the deconvolution window
#     # Make sure intensities are not negative
#     heights = windowheights(dw, bl)

#     # Determine peak limits: first and last index of peak that is not zero, relative to 
#     # all scans and in the decomnvolution window only
#     nzstartindex, nzstopindex, peakheights = peaklimits(dw, heights)

#     # Infer candidate peak model
#     nodes, model, startindex, stopindex, starttime, stoptime = peakmodel(chromatogram(dw), 
#         ionindex(dw), δtᵢ(dw), nzstartindex, nzstopindex, peakheights)

#     # Assess the peak height at the peak maximum
#     maxtime, maxheight = peakmaximum(nodes, model)

#     # Peak cannot have multiple maxima of the same height
#     isnothing(maxtime) && return nothing
#     maxintensity = bl(maxtime) + maxheight

#     # Peak hight, in noise units, above baseline must be greater than the rejection 
#     # threshold
#     maxheight ≤ cquantile(Normal(0, σ_signal), q2) * sqrt(maxintensity) && return nothing

#     CandidatePeak(starttime, stoptime, startindex, stopindex, nzstartindex, nzstopindex, 
#         maxtime, maxheight, maxintensity, model, nodes, peakheights, bl, dw, nothing, 
#         nothing, nothing, nothing)
# end


# function baseline(dw::DeconvolutionWindow, baselineindexfraction::Real)
#     bindices = baselineindices(dw, baselineindexfraction)
#     xs = [ionscantime(δtᵢ(dw), chromatogram(dw), iₛ, ionindex(dw)) for iₛ in bindices]
#     ys = [intensity(chromatogram(dw), iₛ, ionindex(dw)) for iₛ in bindices]
#     intercept, slope = baselinefit(xs, ys, ionscantime(δtᵢ(dw), chromatogram(dw), 
#         leftindex(dw), ionindex(dw)), ionscantime(δtᵢ(dw), chromatogram(dw), 
#         rightindex(dw), ionindex(dw)))
#     t -> intercept + slope * t
# end


# function baselineindices(dw::DeconvolutionWindow, baselineindexfraction::Real)
#     tentativebaseline = lowvalline(dw)
#     idcs = Vector{Int}(undef, windowsize(dw))
#     ints = Vector{Float64}(undef, windowsize(dw))
#     for (i, iₛ) in enumerate(windowindices(dw))
#         idcs[i] = iₛ
#         ints[i] = (intensity(chromatogram(dw), iₛ, ionindex(dw))
#             - tentativebaseline(ionscantime(δtᵢ(dw), chromatogram(dw), iₛ, ionindex(dw))))
#     end
#     p = sortperm(ints)
#     n = ceil(Int, windowsize(dw) * baselineindexfraction)
#     @view idcs[@view p[1:n]]
# end


# function windowheights(dw::DeconvolutionWindow, baseline::Function)
#     heights = Vector{Float64}(undef, windowsize(dw))
#     for (i, iₛ) in enumerate(windowindices(dw))
#         t = ionscantime(δtᵢ(dw), chromatogram(dw), iₛ, ionindex(dw))
#         Δint = intensity(chromatogram(dw), iₛ, ionindex(dw)) - baseline(t)
#         heights[i] = Δint < 0 ? zero(eltype(Δint)) : Δint
#     end
#     heights
# end


# function nnls_fit(xs::AbstractVector{T}, ys::AbstractVector{<:Real}, x_start::T, x_stop::T,
#     reversed::Bool=false, firstpass::Bool=true) where {T<:Real}

#     windowsize = length(xs)

#     if !reversed
#         xs_shifted = xs .- x_start
#     else
#         xs_shifted = [abs(x - x_stop) for x in xs]
#     end
    
#     matrix = hcat(ones(T, windowsize), xs_shifted)
#     S = promote_type(T, eltype(ys))
#     intercept, slope = NNLS.nnls(convert(Matrix{S}, matrix), convert(Vector{S}, ys))  

#     # If the slope is zero, reverse the data as we (probably) have a negative slope. If 
#     # the slope is really zero, the result should not be different.
#     firstpass && slope == 0 && return nnls_fit(xs, ys, x_start, x_stop, !reversed, false)

#     if !reversed
#         intercept -= x_start * slope
#     else
#         slope = -slope
#         intercept -= x_stop * slope
#     end

#     intercept, slope
# end


# function baselinefit(xs::AbstractVector{T}, ys::AbstractVector{<:Real}, x_start::T, 
#     x_stop::T) where {T}

#     # See if regular linear regression predicts non-negative baseline values across the 
#     # window
#     intercept, slope = lsfit(xs, ys)
#     f(t::Unitful.Time) = intercept + slope * t
#     f(x_start) ≥ 0 && f(x_stop) ≥ 0 && return intercept, slope

#     # Since regular linear regression also predicted negative values, try nnls regression
#     timeunit = unit(eltype(xs))
#     intercept, slope = nnls_fit(ustrip.(timeunit, xs), ys, 
#         ustrip.(timeunit, (x_start, x_stop))..., slope * timeunit ≤ 0)
#     intercept, slope / timeunit
# end


# function peakmaximum(times::Vector{<:Unitful.Time}, itp::Function)

#     n = length(times) - 1
#     timeunit = unit(eltype(times))
#     xs = convert(Vector{Float64}, ustrip.(times))  # Need Float64 for optimizer

#     intensityₚ = 0.0
#     xₚ = Float64[]
#     for k in 1:n

#         x = minimizer(optimize(x -> -itp(x * timeunit), xs[k], xs[k+1]))
#         intensity = itp(x * timeunit)
        
#         if intensityₚ < intensity
#            intensityₚ = intensity
#            empty!(xₚ)
#            push!(xₚ, x)
#         elseif intensityₚ ≠ 0 && intensityₚ ≈ intensity
#             unique = true
#             for xᵢ in xₚ
#                 if x ≈ xᵢ
#                     unique = false
#                     break
#                 end
#             end
#             unique && push!(xₚ, x)
#         end
#     end

#     peak_count = length(xₚ)
#     if peak_count ≠ 1
#         # Check whether the optimization failed to converge at the same point
#         if peak_count == 2
#             avgnodedistance = ustrip(timeunit, (last(times) - first(times)) / n)
#             maxdistance = abs(last(xₚ) - first(xₚ))
#             if maxdistance < avgnodedistance
#                 x = minimizer(optimize(x -> -itp(x * timeunit), first(xₚ), last(xₚ)))
#                 intensity = itp(x * timeunit)
#                 if intensityₚ < intensity && !(x ≈ first(xₚ) || x ≈ last(xₚ))
#                     return x * timeunit, intensity
#                 end
#             end
#         end
#         # Comment: There are two reasons why multiple maxima may be detected:
#         # - Multiple local maxima may exist within the specified range of peak intensities.
#         # - In rare cases, the peak model may show two local maxima on either side of the 
#         #   true single maximum. This usually happens when intensity values are low, so 
#         #   these peaks are likely to be filtered out. Since the program relies on B-Spline 
#         #   peak models, these cases cannot be handled and are therefore discarded.
#         nothing, nothing
#     else
#         first(xₚ) * timeunit, intensityₚ
#     end
# end


# function peaklimits(dw::DeconvolutionWindow, heights::AbstractVector{<:Real})
#     start, stop = 1, length(heights)
#     for i in (maxstartindex(dw) - leftindex(dw) + 1):-1:1
#         if heights[i] == 0
#             start = i + 1
#             break
#         end
#     end
#     for i in (maxstopindex(dw) - leftindex(dw) + 1):stop
#         if heights[i] == 0
#             stop = i - 1
#             break
#         end
#     end
#     leftindex(dw) + start - 1, leftindex(dw) + stop - 1, heights[start:stop]
# end


# function peakmodel(chrom::AbstractChromMS, ionindex::Integer, δtᵢ::Function, 
#     nzstartindex::Integer, nzstopindex::Integer, peakheights::AbstractVector{<:Real})

#     timeunit = unit(eltype(ionscantimes(δtᵢ, chrom, ionindex)))
#     xs = ustrip.([ionscantime(δtᵢ, chrom, iₛ, ionindex) for iₛ in nzstartindex:nzstopindex])  

#     startindex, stopindex = nzstartindex, nzstopindex

#     if nzstartindex == 1
#         pushfirst!(xs, xs[begin] - (xs[begin+1] - xs[begin]))
#     else
#         startindex -= 1
#         pushfirst!(xs, ionscantime(δtᵢ, chrom, startindex, ionindex, ustripped=true))
#     end
#     if nzstopindex == scancount(chrom)
#         push!(xs, xs[end] + (xs[end] - xs[end-1]))
#     else
#         stopindex += 1
#         push!(xs, ionscantime(δtᵢ, chrom, stopindex, ionindex, ustripped=true))
#     end
    
#     ys = vcat(0.0, peakheights, 0.0)
#     itp = CubicSplineInterpolator(xs, sqrt.(ys))

#     starttime, stoptime = ionscantime.(δtᵢ, chrom, [startindex, stopindex], ionindex)

#     xs * timeunit, rt::Unitful.Time -> begin
#         t = ustrip(timeunit, rt);
#         first(xs) ≤ t ≤ last(xs) ? itp(t)^2 : 0.0
#     end, startindex, stopindex, starttime, stoptime 
# end


# function peaksharpness!(cp::CandidatePeak)
#     timeunit = unit(eltype(ionscantimes(δtᵢ(cp), chromatogram(cp), ionindex(cp))))
#     xs = convert(Vector{Float64}, ustrip.(timeunit, nodes(cp)))
#     xₚ = ustrip(timeunit, maxtime(cp))
#     denominator = sqrt(maxintensity(cp))
#     sₗ = sᵣ = tₗ = tᵣ = zero(Float64)
#     f = x -> (maxheight(cp) - cp.model(x * timeunit)) / abs(x - xₚ) / denominator
#     first_iteration = true
#     for i in eachindex(xs)
#         first_iteration && (first_iteration = false; continue) 
#         if xs[i] ≤ xₚ
#             t = minimizer(optimize(x -> -f(x), xs[i-1], xs[i]))
#             s = f(t)
#             if sₗ < s
#                 sₗ = s
#                 tₗ = t
#             end
#         elseif xs[i-1] ≥ xₚ
#             t = minimizer(optimize(x -> -f(x), xs[i-1], xs[i]))
#             s = f(t)
#             if sᵣ < s
#                 sᵣ = s
#                 tᵣ = t
#             end
#         else
#             t₁ = minimizer(optimize(x -> -f(x), xs[i-1], xₚ))
#             t₂ = minimizer(optimize(x -> -f(x), xₚ, xs[i]))
#             s₁ = f(t₁)
#             s₂ = f(t₂)
#             if sₗ < s₁
#                 sₗ = s₁
#                 tₗ = t₁
#             end
#             if sᵣ < s₂
#                 sᵣ = s₂
#                 tᵣ = t₂
#             end
#         end
#     end
    
#     if sₗ == 0
#         println("sₗ is zero!")
#     end
#     if sᵣ == 0
#         println("sᵣ is zero!")
#     end

#     sharpness!(cp, uconvert(u"s^-1", (sₗ + sᵣ / 2) / timeunit)) 
#     leftsharpnesstime!(cp, uconvert(u"minute", tₗ * timeunit))
#     rightsharpnesstime!(cp, uconvert(u"minute", tᵣ * timeunit))

#     nothing
# end


# function removeoverlaps(candidatepeaks::AbstractVector{CandidatePeak})
#     # Filter peaks so that they do not overlap
#     cp_f = Vector{CandidatePeak}()
#     for cp in candidatepeaks
#         # First peak cannot be problematic per se
#         if length(cp_f) == 0
#             push!(cp_f, cp)
#             continue
#         end

#         # CASE 1: complete redundancy (ony the DW local maxima differ)
#         cp.starttime == cp_f[end].starttime && cp.stoptime == cp_f[end].stoptime && continue

#         # CASE 2: partial overlap
#         if cp.starttime < cp_f[end].stoptime
#             if cp_f[end].sharpness < cp.sharpness
#                 cp_f[end] = cp
#             end
#             continue
#         end

#         push!(cp_f, cp)
#     end
#     cp_f
# end


# ############################################################################################
# # Code for computing deconvoluted total ion chromatogram (DTIC)
# ############################################################################################
# # Problem: What if someone create a vector with candidate peaks from different runs?
# # How could I prevent from analyzing in such a situation?
# function deconvolutedtic(chrom::AbstractChromMS, 
#     candidatepeaks::AbstractVector{CandidatePeak}; binsperscan::Int=10)
#     bin = binfunction(chrom, binsperscan)
#     bincount = binsperscan * scancount(chrom)
#     sharpnesssums = zeros(Float64, bincount)
#     binpeaks = Dict{Int, Vector{CandidatePeak}}()
#     sharpnessunit = unit(sharpness(first(candidatepeaks)))
#     for cp in candidatepeaks
#         b = bin(maxtime(cp))
#         sharpnesssums[b] += ustrip(sharpness(cp))
#         binpeaks[b] = push!(get(binpeaks, b, Vector{CandidatePeak}()), cp)
#     end

#     starttime = minscantime(chrom) - scanduration(chrom)
#     Δt_bin = runduration(chrom) / bincount
#     bintimes = [starttime + x * Δt_bin for x in 1:bincount]

#     bintimes, sharpnesssums * sharpnessunit, binpeaks
# end


# function binfunction(chrom::AbstractChromMS, binsperscan::Integer)
#     Δt = scanduration(chrom)
#     starttime = minscantime(chrom) - Δt
#     rt -> begin
#         dividend = (rt - starttime) * binsperscan 
#         q = dividend / Δt
#         qr = round(Int, q)
#         q ≈ qr && return qr
#         convert(Int, fld(dividend, Δt))
#     end
# end

# ############################################################################################
# # Code for extracting mass spectra from dtic data
# ############################################################################################
# function componentbins(sharpnesssums; proportionalityfactor=5000)
#     bincount = length(sharpnesssums)
#     componentbins = Vector{Int}()
#     for i in sortperm(sharpnesssums, rev=true)
#         sharpnesssums[i] == 0.0u"s^-1" && break
#         iₗ = i > 1 ? i - 1 : i
#         iᵣ = i < bincount ? i + 1 : i
#         sharpnesssum = sum(@view sharpnesssums[iₗ:iᵣ])
#         bins = convert(Int, ceil(proportionalityfactor / ustrip(sharpnesssum)))

#         iₗₗ = i - bins > 0 ? i - bins : 1
#         sharpnesssums[i] ≤ maximum(@view sharpnesssums[iₗₗ:iₗ]) && continue

#         iᵣᵣ = i + bins ≤ bincount ? i + bins : bincount
#         sharpnesssums[i] ≤ maximum(@view sharpnesssums[iᵣ:iᵣᵣ]) && continue

#         push!(componentbins, i)
#     end
#     componentbins
# end


# function componentpeakmodels(componentbins, binpeaks, fraction::Real=1)
#     componentpeakmodels = []
#     for bin in componentbins
#         maxsharpness = maximum(cp -> sharpness(cp), binpeaks[bin])
#         bestpeaks = filter(cp -> sharpness(cp) ≥ fraction * maxsharpness, binpeaks[bin])
#         map(cp -> area!(cp, peakarea(model(cp), nodes(cp))), bestpeaks)
#         push!(componentpeakmodels, first(bestpeaks))
#     end
#     componentpeakmodels
# end


# function peakarea(itp::Function, nodes::AbstractVector{<:Unitful.Time})
#     first(quadgk(t -> itp(t), first(nodes), last(nodes)))
# end


# function massspectra(chrom, ionscanorder::IonScanOrder, cpms)
#     δtᵢ = ionscantimeshift(chrom, ionscanorder)
#     if length(cpms) ≥ Threads.nthreads()
#         indexchunks = Iterators.partition(1:length(cpms), length(cpms) ÷ Threads.nthreads())
#         tasks = map(indexchunks) do indices
#             Threads.@spawn massspectra(chrom, δtᵢ, cpms, indices)
#         end
#         chunk_lms = fetch.(tasks)
#         collect(Iterators.flatten(chunk_lms))
#     else
#         massspectra(chrom, δtᵢ, cpms, eachindex(cpms))
#     end
# end


# function massspectra(chrom, δtᵢ, cpms, indices)
#     timeunit = unit(eltype(scantimes(chrom)))
#     mss = Vector{MassSpectrum}(undef, length(indices))
#     for (msindex, cpindex) in enumerate(indices)
#         cp = cpms[cpindex]
#         # Get indices of mass spectrum extraction window
#         extra::Int = ceil((stopindex(cp) - startindex(cp) + 1) / 2)
#         start = startindex(cp) - extra ≥ 1 ? startindex(cp) - extra : 1 
#         stop = (stopindex(cp) + extra ≤ scancount(chrom) ? 
#             stopindex(cp) + extra : scancount(chrom))
        
#         # Identify additional components in extraction window
#         cpmsₐ = []
#         for cpₐ in cpms
#             cp == cpₐ && continue
#             startₐ = startindex(cpₐ) - 1 ≥ 1 ? startindex(cpₐ) - 1 : startindex(cpₐ)
#             stopₐ = (stopindex(cpₐ) + 1 ≤ scancount(chrom) ? 
#                 stopindex(cpₐ) + 1 : scancount(chrom))
#             if start < startₐ < stop || start < stopₐ < stop
#                 push!(cpmsₐ, cpₐ)
#             end
#         end

#         # Extract mass spectrum
#         msintensities = zeros(Float64, ioncount(chrom))
#         peak_indices = startindex(cp):stopindex(cp)  # NEW
#         peak_tic = zeros(Float64, length(peak_indices))  # NEW
#         for iᵢ in eachindex(ions(chrom))

#             # compile ion scan times and associated intensities
#             xs = [ionscantime(δtᵢ, chrom, iₛ, iᵢ) for iₛ in start:stop]
#             ys = @view intensities(chrom)[start:stop, iᵢ]

#             # Compute matrix with model intensities at the ion scan times
#             models = Matrix{Float64}(undef, length(xs), length(cpmsₐ) + 1)
#             models[:, 1] = model(cp).(xs) / ustrip(timeunit, area(cp))  # peak of interest
#             for i in eachindex(cpmsₐ)  # additional components in window
#                 models[:, i + 1] = model(cpmsₐ[i]).(xs) / ustrip(timeunit, area(cpmsₐ[i]))
#             end

#             # Use nnls regression to fit the peak model data to the run intensity values
#             intercept, slope, scalefactor, scalefactorₐ... = peakintensity(ustrip.(xs), ys, 
#                 models)
#             msintensities[iᵢ] = scalefactor

#             # NEW
#             xs_peak = [ionscantime(δtᵢ, chrom, iₛ, iᵢ) for iₛ in peak_indices]
#             peak_tic .+= model(cp).(xs_peak) / ustrip(timeunit, area(cp)) * scalefactor
#             peak_tic .+= intercept .+ slope * ustrip.(xs_peak)
#         end
#         mss[msindex] = MassSpectrum(ions(chrom)[:], msintensities, 
#             retentiontime=maxtime(cp), metadata=Dict(:peak_tic => peak_tic, :peak_scanindices => peak_indices))
#     end
#     mss
# end


# function peakintensity(
#     times::AbstractVector{T1},
#     intensities::AbstractVector{<:Real},
#     modelintensities::Matrix{<:Real},
#     firstpass::Bool=true,
#     stoptime::T1=zero(T1)
#     ) where {T1<:Real}

#     windowsize = length(times)

#     if firstpass
#         starttime, stoptime = times[begin], times[end]
#         shiftedtimes = times .- starttime
#     else
#         shiftedtimes = times
#         modelintensities = @view modelintensities[windowsize:-1:1, :]
#     end
    
#     matrix = hcat(ones(T1, windowsize), shiftedtimes, modelintensities)
#     T = promote_type(eltype(matrix), eltype(intensities))
#     θ = NNLS.nnls(convert(Matrix{T}, matrix), convert(Vector{T}, intensities))  

#     if firstpass
#         θ[2] == 0 && return peakintensity(shiftedtimes, reverse(intensities), 
#             modelintensities, false, stoptime)
#         θ[1] -= starttime * θ[2]
#     else
#         θ[1] += stoptime * θ[2]
#         θ[2] = -θ[2]
#     end
#     θ
# end


# function saveplot(dw::DeconvolutionWindow, file::AbstractString)
#     CairoMakie.activate!()

#     f = Figure()

#     # AXIS 1: Details on deconvolution window
#     iᵢ = ionindex(dw)
#     max = string(maxstartindex(dw), ":", maxstopindex(dw))
#     dw_limits = string(leftindex(dw), ":", rightindex(dw))
#     ax1 = Axis(f[1,1], 
#         title="ion index: $iᵢ, maximum: $max, devonv. limits: $dw_limits", 
#         xlabel="Retention time [minute]", ylabel="Abundance")
#     start_idx_extended = leftindex(dw) > 1 ? leftindex(dw) - 1 : leftindex(dw)
#     stop_idx_extended = rightindex(dw) < scancount(chromatogram(dw)) ? rightindex(dw) + 1 : rightindex(dw)
    
#     Δt = scanduration(chromatogram(dw), timeunit=u"minute", ustripped=true) / 2

#     xmin = ionscantime(δtᵢ(dw), chromatogram(dw), start_idx_extended, ionindex(dw), timeunit=u"minute", ustripped=true) - Δt
#     xmax = ionscantime(δtᵢ(dw),chromatogram(dw), stop_idx_extended, ionindex(dw), timeunit=u"minute", ustripped=true) + Δt
#     xlims!(ax1, xmin, xmax)

#     # Intensities in deconvolutioin window
#     ts1 = ionscantime.(δtᵢ(dw), chromatogram(dw), start_idx_extended:stop_idx_extended, ionindex(dw), 
#         timeunit=u"minute", ustripped=true)
#     ints1 = intensities(chromatogram(dw))[start_idx_extended:stop_idx_extended, ionindex(dw)]
#     lines!(ax1, ts1, ints1, color=:blue)

#     # Indicate deconvolution window limits
#     ts2 = ionscantime.(δtᵢ(dw), chromatogram(dw), [leftindex(dw), rightindex(dw)], ionindex(dw), 
#         timeunit=u"minute", ustripped=true)
#     vlines!(ax1, ts2, color=:orange)

#     # Indicate baseline-adjusted minima
#     icdes = [leftlowindex(dw), rightlowindex(dw)]
#     ts3 = ionscantime.(δtᵢ(dw), chromatogram(dw), icdes, ionindex(dw), timeunit=u"minute", ustripped=true)
#     ints3 = intensity.(chromatogram(dw), icdes, ionindex(dw))
#     scatter!(ax1, ts3, ints3, color=:grey)

#     # Draw tentative baseline
#     tentativebaseline = lowvalline(dw)
#     ts4 = ionscantime.(δtᵢ(dw), chromatogram(dw), windowindices(dw), ionindex(dw), timeunit=u"minute",
#        ustripped=true)
#     ints4 = tentativebaseline.(ionscantime.(δtᵢ(dw), chromatogram(dw), windowindices(dw), ionindex(dw)))
#     lines!(ax1, ts4, ints4, color=:grey)

#     # # Indicate maximum
#     ts5 = ionscantime.(δtᵢ(dw), chromatogram(dw), maxindexrange(dw), ionindex(dw), timeunit=u"minute", 
#        ustripped=true)
#     ints5 = intensity.(chromatogram(dw), maxindexrange(dw), ionindex(dw))
#     scatter!(ax1, ts5, ints5, color=:red)

#     # Draw least-squares baseline
#     #ts6 = ionscantime.(δtᵢ(dw), chromatogram(dw), windowindices(dw), ionindex(dw), timeunit=u"minute",
#     #    ustripped=true)
#     #ints6 = dw.baseline.(ionscantime.(δtᵢ(dw), chromatogram(dw), windowindices(dw), ionindex(dw)))
#     #lines!(ax1, ts6, ints6, color=:red)

#     save(file, f)
# end


# function saveplot(cp::CandidatePeak, file::AbstractString)
#     CairoMakie.activate!()

#     f = Figure()

#     # AXIS 1: Details on deconvolution window
#     iᵢ = ionindex(cp)
#     max = string(maxstartindex(cp), ":", maxstopindex(cp))
#     dw_limits = string(leftindex(cp), ":", rightindex(cp))
#     ax1 = Axis(f[1,1], 
#         title="ion index: $iᵢ, maximum: $max, devonv. limits: $dw_limits", 
#         xlabel="Retention time [minute]", ylabel="Abundance")
#     start_idx_extended = leftindex(cp) > 1 ? leftindex(cp) - 1 : leftindex(cp)
#     stop_idx_extended = rightindex(cp) < scancount(chromatogram(cp)) ? rightindex(cp) + 1 : rightindex(cp)
    
#     Δt = scanduration(chromatogram(cp), timeunit=u"minute", ustripped=true) / 2

#     xmin = ionscantime(δtᵢ(cp), chromatogram(cp), start_idx_extended, ionindex(cp), timeunit=u"minute", ustripped=true) - Δt
#     xmax = ionscantime(δtᵢ(cp),chromatogram(cp), stop_idx_extended, ionindex(cp), timeunit=u"minute", ustripped=true) + Δt
#     xlims!(ax1, xmin, xmax)

#     # Intensities in deconvolutioin window
#     ts1 = ionscantime.(δtᵢ(cp), chromatogram(cp), start_idx_extended:stop_idx_extended, ionindex(cp), 
#         timeunit=u"minute", ustripped=true)
#     ints1 = intensities(chromatogram(cp))[start_idx_extended:stop_idx_extended, ionindex(cp)]
#     lines!(ax1, ts1, ints1, color=:blue)

#     # Indicate deconvolution window limits
#     ts2 = ionscantime.(δtᵢ(cp), chromatogram(cp), [leftindex(cp), rightindex(cp)], ionindex(cp), 
#         timeunit=u"minute", ustripped=true)
#     vlines!(ax1, ts2, color=:orange)

#     # Indicate baseline-adjusted minima
#     icdes = [leftlowindex(cp), rightlowindex(cp)]
#     ts3 = ionscantime.(δtᵢ(cp), chromatogram(cp), icdes, ionindex(cp), timeunit=u"minute", ustripped=true)
#     ints3 = intensity.(chromatogram(cp), icdes, ionindex(cp))
#     scatter!(ax1, ts3, ints3, color=:grey)

#     # Draw tentative baseline
#     tentativebaseline = lowvalline(cp)
#     ts4 = ionscantime.(δtᵢ(cp), chromatogram(cp), windowindices(cp), ionindex(cp), timeunit=u"minute",
#        ustripped=true)
#     ints4 = tentativebaseline.(ionscantime.(δtᵢ(cp), chromatogram(cp), windowindices(cp), ionindex(cp)))
#     lines!(ax1, ts4, ints4, color=:grey)

#     # # Indicate maximum
#     ts5 = ionscantime.(δtᵢ(cp), chromatogram(cp), maxindexrange(cp), ionindex(cp), timeunit=u"minute", 
#        ustripped=true)
#     ints5 = intensity.(chromatogram(cp), maxindexrange(cp), ionindex(cp))
#     scatter!(ax1, ts5, ints5, color=:red)

#     # Draw least-squares baseline
#     ts6 = ionscantime.(δtᵢ(cp), chromatogram(cp), windowindices(cp), ionindex(cp), timeunit=u"minute",
#         ustripped=true)
#     ints6 = cp.baseline.(ionscantime.(δtᵢ(cp), chromatogram(cp), windowindices(cp), ionindex(cp)))
#     lines!(ax1, ts6, ints6, color=:red)

#     ########################################################################################
#     # AXIS 2: Details on peak model
#     ########################################################################################
#     tₚ = round(u"minute", maxtime(cp), digits=4)
#     s = round(u"s^-1", sharpness(cp), digits=1)
#     ax2 = Axis(f[2,1], title="mode: $tₚ, sharpness: $s", 
#         xlabel="Retention time [minute]", ylabel="Abundance")
#     xlims!(ax2, ax1.xaxis.attributes.limits.val...)
    
#     ts7 = ionscantime.(δtᵢ(cp), chromatogram(cp), nzstartindex(cp):nzstopindex(cp), 
#         ionindex(cp), timeunit=u"minute", ustripped=true)
#     ints7 = heights(cp)
#     scatter!(ax2, ts7, ints7, color=:blue)

#     # Indicate scan time peak
#     scatter!(ax2, ustrip(u"minute", maxtime(cp)), maxheight(cp), color=:red)

#     # Draw peak b-spline
#     ts8 = LinRange(ionscantime(δtᵢ(cp), chromatogram(cp), start_idx_extended, ionindex(cp)), 
#         ionscantime(δtᵢ(cp), chromatogram(cp), stop_idx_extended, ionindex(cp)), 100)
#     ints8 = model(cp).(ts8)
#     lines!(ax2, ustrip.(u"minute", ts8), ints8, color=:green)

#     # Indicate peak limits
#     vlines!(ax2, ustrip.(u"minute", [starttime(cp), stoptime(cp)]), color=:orange)

#     # Indicate lines that were used to infer the peak sharpness
#     t₁ = leftsharpnesstime(cp)
#     t₂ = maxtime(cp)
#     t₃ = rightsharpnesstime(cp)
#     ts10 = ustrip.(u"minute", [t₁, t₂, t₃])

#     ints₁ = model(cp)(t₁)
#     ints₂ = maxheight(cp)
#     ints₃ = model(cp)(t₃)
#     ints10 = [ints₁, ints₂, ints₃]
#     lines!(ax2, ts10, ints10, color=:violet)

#     save(file, f)
# end

end  # module