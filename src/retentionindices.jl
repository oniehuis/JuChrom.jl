struct RiMapper{T1<:AbstractString, T2<:AbstractVector{<:Unitful.Time}, 
    T3<:AbstractVector{<:Real}, T4<:AbstractString, T5<:Union{AbstractString, Nothing}
    } <: AbstractRiMapper
    retentionindexname::T1
    retentiontimes::T2
    retentionindices::T3
    interpolationmethod::T4
    extrapolationmethod::T5
    rt2ri::Function
    metadata::Dict{Any, Any}
    function RiMapper{T1, T2, T3, T4, T5}(retentionindexname::T1, retentiontimes::T2, 
        retentionindices::T3, interpolationmethod::T4, extrapolationmethod::T5, 
        rt2ri::Function, metadata::Dict
        ) where {T1<:AbstractString, T2<:AbstractVector{<:Unitful.Time}, 
        T3<:AbstractVector{<:Real}, T4<:AbstractString, T5<:Union{AbstractString, Nothing}}

        length(retentionindexname) > 0 || throw(ArgumentError(
            "no retention index name provided"))
        length(interpolationmethod) > 0 || throw(ArgumentError(
            "interpolation method not specified"))
        !isnothing(extrapolationmethod) && length(extrapolationmethod) == 0 && throw(
            ArgumentError("extrapolation method not specified"))
        length(retentiontimes) == length(retentionindices) || throw(DimensionMismatch(
            "number of retention times and number of retention indices are different"))
        length(retentiontimes) > 1 || throw(
            ArgumentError("need at least two retention time-retention index pairs"))
        length(Set(retentiontimes)) == length(retentiontimes) || throw(ArgumentError(
            "retention times contain identical values"))
        issorted(retentiontimes) || throw(ArgumentError(
            "retention times not in ascending order"))
        length(Set(retentionindices)) == length(retentionindices) || throw(ArgumentError(
            "retention times contain identical values"))
        issorted(retentionindices) || throw(ArgumentError(
            "retention indices not in ascending order"))

        # testing inter und extrapolation?

        new(retentionindexname, retentiontimes, retentionindices, interpolationmethod, 
            extrapolationmethod, rt2ri, metadata)
    end
end


Base.broadcastable(mapper::RiMapper) = Ref(mapper)


function Base.show(io::IO, mapper::RiMapper)
    # Information about the element types
    type = typeof(mapper)
    typename = name(type)
    samplecount = length(retentiontimes(mapper))

    # Informan about the index and the number of ladder steps
    println(io, "$typename {index name: ", retentionindexname(mapper), 
        ", calibration points: ", samplecount, "}") 
    
    # Information about retention times
    if samplecount > 10
        println(io, "retention time range: ", minretentiontime(mapper), " - ", 
            maxretentiontime(mapper))
    else
        print(io, "retention times: ")
        for i in eachindex(retentiontimes(mapper))
            print(io, retentiontimes(mapper)[i])
            i < samplecount ? print(io, ", ") : println(io)
        end
    end

    # Information about retention indices
    if samplecount > 10
        println(io, "retention index range: ", minretentionindex(mapper), " - ", 
            maxretentionindex(mapper))
    else
        print(io, "retention indices: ")
        for i in eachindex(retentionindices(mapper))
            print(io, retentionindices(mapper)[i])
            i < samplecount ? print(io, ", ") : println(io)
        end
    end

    # Information about inter- and extrapolation
    println(io, "interpolation method: ", interpolationmethod(mapper))
    println(io, "extrapolation method: ", extrapolationmethod(mapper))

    # Information about the metadata
    n = length(metadata(mapper))
    print(io, "metadata: ", n, (n == 0 || n > 1) ? " entries" : " entry" )
end


"""
    RiMapper(retentionindexname::AbstractString, 
    retentiontimes::AbstractVector{<:Unitful.Time},
    retentionindices::AbstractVector{<:Real}; extrapolation::Bool=true, 
    metadata::Dict=Dict())

Construct an RiMapper object for mapping retention times to retention indices and vice 
versa using interpolation and, if enabled (default), extrapolation. The function utilizes 
a B-spline for interpolation, calculated from a vector of retention times and a 
corresponding vector of retention indices. For retention time values outside the range 
used to compute the B-spline, the function employs linear extrapolation to estimate a 
retention index. An optional keyword argument, `extrapolation`, can be used to disable 
extrapolation, in which case the function returns nothing for values outside the retention 
time range. The function will raise an error if the resulting mapping function does not 
produce continuously increasing values. Note that both retention times and 
retention indices must be in ascending order.

See also [`AbstractRiMapper`](@ref), [`retentionindexname`](@ref), 
[`retentionindices`](@ref), [`retentiontimes`](@ref), [`interpolationmethod`](@ref),
[`extrapolationmethod`](@ref), [`rt2ri`](@ref), [`metadata`](@ref), 
[`retentionindex`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000)
RiMapper {index name: Kovats, calibration points: 5}
retention times: 1 minute, 2 minute, 3 minute, 4 minute, 5 minute
retention indices: 1000, 2000, 3000, 4000, 5000
interpolation method: natural cubic b-spline
extrapolation method: linear
metadata: 0 entries

julia> retentionindex(ld, 1u"minute") ≈ 1000.0
true

julia> retentionindex(ld, 1.5u"minute") ≈ 1500.0
true

julia> retentionindex(ld, 11u"minute") ≈ 11000.0
[ Info: extrapolated value
true

julia> retentionindices(ld)
1000:1000:5000

julia> retentiontimes(ld)
(1:5) minute

julia> interpolationmethod(ld)
"natural cubic b-spline"

julia> extrapolationmethod(ld)
"linear"

julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000, extrapolation=false)
RiMapper {index name: Kovats, calibration points: 5}
retention times: 1 minute, 2 minute, 3 minute, 4 minute, 5 minute
retention indices: 1000, 2000, 3000, 4000, 5000
interpolation method: natural cubic b-spline
extrapolation method: nothing
metadata: 0 entries

julia> retentionindex(ld, 11u"minute") === nothing
true

julia> extrapolationmethod(ld) === nothing
true
```
"""
function RiMapper(retentionindexname::T1, retentiontimes::T2, retentionindices::T3;
    extrapolation::Bool=true, metadata::Dict=Dict{Any, Any}()) where {T1<:AbstractString, 
    T2<:AbstractVector{<:Unitful.Time}, T3<:AbstractVector{<:Real}}
    rt2ri = bsplineinterpolation(retentiontimes, retentionindices, 
        extrapolation=extrapolation)
    
    interpolationmethod = "natural cubic b-spline"
    extrapolationmethod = extrapolation ? "linear" : nothing
    T4 = typeof(interpolationmethod)
    T5 = typeof(extrapolationmethod)
    RiMapper{T1, T2, T3, T4, T5}(retentionindexname, retentiontimes, retentionindices, 
        interpolationmethod, extrapolationmethod, rt2ri, metadata)
end


"""
    interpolationmethod(mapper::RiMapper)

Return the name of the applied interpolation method.

See also [`RiMapper`](@ref), [`extrapolationmethod`](@ref), [`retentionindex`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);

julia> interpolationmethod(ld)
"natural cubic b-spline"
```
"""
interpolationmethod(mapper::RiMapper) = mapper.interpolationmethod


"""
    extrapolationmethod(mapper::RiMapper)

Return the name of the extrapolation method. If no extrapolation is applied, the function 
returns nothing.

See also [`RiMapper`](@ref), [`interpolationmethod`](@ref), [`retentionindex`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);

julia> extrapolationmethod(ld)
"linear"

julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"s", [100, 200, 300], extrapolation=false);

julia> extrapolationmethod(ld) === nothing
true
```
"""
extrapolationmethod(mapper::RiMapper) = mapper.extrapolationmethod


"""
    maxretentionindex(mapper::RiMapper)

Return the highest retention index used to construct the retention index mapper.

See also [`RiMapper`](@ref), [`minretentionindex`](@ref), [`retentionindices`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);

julia> maxretentionindex(ld)
3000
```
"""
maxretentionindex(mapper::RiMapper) = last(retentionindices(mapper))


"""
    maxretentiontime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the highest retention time used to construct the retention index mapper. The 
optional keyword argument `timeunit` lets you change the unit of the returned retention 
time. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether to include the 
unit in the returned value.

See also [`RiMapper`](@ref), [`maxretentiontime`](@ref), [`retentiontimes`](@ref).

# Examples
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);

julia> maxretentiontime(ld) ≈ 3.8u"minute"
true

julia> maxretentiontime(ld, timeunit=u"s") ≈ 228.0u"s"
true

julia> maxretentiontime(ld, timeunit=u"s", ustripped=true) ≈ 228.0
true
```
"""
function maxretentiontime(mapper::RiMapper;
    timeunit::Unitful.TimeUnits=unit(eltype(retentiontimes(mapper))), 
    ustripped::Bool=false)

    t = last(retentiontimes(mapper))
    ustripped ? ustrip(timeunit, t) : uconvert(timeunit, t)
end


"""
    metadata(mapper::AbstractRiMapper)

Return the metadata.

See also [`AbstractRiMapper`](@ref), [`RiMapper`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", (1:5)u"minute", 10:10:50, metadata=Dict(:id => 7))
RiMapper {index name: Kovats, calibration points: 5}
retention times: 1 minute, 2 minute, 3 minute, 4 minute, 5 minute
retention indices: 10, 20, 30, 40, 50
interpolation method: natural cubic b-spline
extrapolation method: linear
metadata: 1 entry

julia> metadata(ld)
Dict{Any, Any} with 1 entry:
  :id => 7
```
"""
metadata(mapper::AbstractRiMapper) = mapper.metadata


"""
    minretentionindex(mapper::RiMapper)

Return the lowest retention index used to construct the retention index mapper.

See also [`RiMapper`](@ref), [`maxretentionindex`](@ref), [`retentionindices`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);

julia> maxretentionindex(ld)
3000
```
"""
minretentionindex(mapper::RiMapper) = first(retentionindices(mapper))


"""
    minretentiontime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the lowest retention time used to construct the retention index mapper. The optional 
keyword argument `timeunit` lets you change the unit of the returned retention time. All 
time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether to include the 
unit in the returned value.

See also [`RiMapper`](@ref), [`maxretentiontime`](@ref), [`retentiontimes`](@ref).

# Examples
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);

julia> minretentiontime(ld) ≈ 1.2u"minute"
true

julia> minretentiontime(ld, timeunit=u"s") ≈ 72.0u"s"
true

julia> minretentiontime(ld, timeunit=u"s", ustripped=true) ≈ 72.0
true
```
"""
function minretentiontime(mapper::RiMapper;
    timeunit::Unitful.TimeUnits=unit(eltype(retentiontimes(mapper))), 
    ustripped::Bool=false)

    t = first(retentiontimes(mapper))
    ustripped ? ustrip(timeunit, t) : uconvert(timeunit, t)
end


"""
    retentionindex(mapper::RiMapper, retentiontime::Unitful.Time; info::Bool=false)

Return the retention index associated with a given retention time. If the return value is 
calculated through extrapolation, a message informing the user of this fact is displayed. 
The optional keyword argument `info` can be used to disable this message.

See also [`RiMapper`](@ref), [`retentionindices`](@ref), [`maxretentionindex`](@ref), 
[`minretentionindex`](@ref).

# Examples
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000])
RiMapper {index name: Kovats, calibration points: 3}
retention times: 1.2 minute, 2.4 minute, 3.8 minute
retention indices: 1000, 2000, 3000
interpolation method: natural cubic b-spline
extrapolation method: linear
metadata: 0 entries

julia> retentionindex(ld, 1.8u"minute") ≈ 1512.3626373626375
true

julia> retentionindex(ld, 1.1u"minute") ≈ 913.9194139194141
[ Info: extrapolated value
true

julia> retentionindex(ld, 1.1u"minute", info=false) ≈ 913.9194139194141
true

julia> ris = retentionindex.(ld, [1, 2]u"minute");  # broadcasting across multiple RT values
[ Info: extrapolated value

julia> ris ≈ [827.8388278388279, 1678.876678876679]
true
```
"""
function retentionindex(mapper::RiMapper, retentiontime::Unitful.Time; info::Bool=true)
    minretentiontime(mapper) ≤ retentiontime ≤ maxretentiontime(mapper) || (
        info && !isnothing(extrapolationmethod(mapper)) && @info "extrapolated value")
    rt2ri(mapper)(retentiontime)
end


"""
    retentionindex(chrom::AbstractChromatogram, retentiontime::Unitful.Time; 
    info::Bool=false)

Return the retention index corresponding to a given retention time. If the value is 
calculated through extrapolation, a message will notify the user of this. To suppress this 
message, use the optional keyword argument `info`. The function will raise an error if the 
chromatogram does not contain a retention index mapper with the required functionality 
implemented.

See also [`RiMapper`](@ref), [`retentionindices`](@ref), [`maxretentionindex`](@ref), 
[`minretentionindex`](@ref).

# Examples
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);

julia> chrom = Chrom(Int64[1, 2, 3]u"s", Int32[12, 956, 1], rimapper=ld);

julia> retentionindex(chrom, 1.8u"minute") ≈ 1512.3626373626375
true

julia> retentionindex(chrom, 1.1u"minute") ≈ 913.9194139194141
[ Info: extrapolated value
true

julia> retentionindex(chrom, 1.1u"minute", info=false) ≈ 913.9194139194141
true

julia> ris = retentionindex.(chrom, [1, 2]u"minute");  # broadcasting across multiple RT values
[ Info: extrapolated value

julia> ris ≈ [827.8388278388279, 1678.876678876679]
true

julia> chrom = Chrom(Int64[1, 2, 3]u"s", Int32[12, 956, 1]);

julia> retentionindex(chrom, 120u"s")
ERROR: ArgumentError: no retention index mapper implemented
[...]
```
"""
function retentionindex(chrom::AbstractChromatogram, 
    retentiontime::Unitful.Time; info::Bool=true)
    isnothing(rimapper(chrom)) && throw(ArgumentError(
        "no retention index mapper implemented"))
    retentionindex(rimapper(chrom), retentiontime, info=info)    
end


"""
    retentionindexname(mapper::RiMapper)

Return the retention index name.

See also [`AbstractRiMapper`](@ref), [`RiMapper`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", (1:10)u"minute", 1000:1000:10000);

julia> retentionindexname(ld)
"Kovats"
```
"""
retentionindexname(mapper::RiMapper) = mapper.retentionindexname


"""
    retentionindices(mapper::RiMapper)

Return the retention indices used to construct the retention index mapper. Note: the 
function returns a reference to the data structure.

See also [`RiMapper`](@ref), [`maxretentionindex`](@ref), [`minretentionindex`](@ref), 
[`retentiontimes`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);

julia> retentionindices(ld)  # reference to the data structure
3-element Vector{Int64}:
 1000
 2000
 3000

julia> retentionindices(ld)[:]  # a copy of these values
3-element Vector{Int64}:
 1000
 2000
 3000
```
"""
retentionindices(mapper::RiMapper) = mapper.retentionindices


# """
#     retentiontime(mapper::RiMapper, retentionindex::Real; timeunit::Unitful.TimeUnits, 
#     ustripped::Bool=false, info::Bool=false)

# Return the retention time associated with a given retention index. If the return value is 
# calculated through extrapolation, a message informing the user of this fact is displayed. 
# The optional keyword argument `info` can be used to disable this message. The optional 
# keyword parameter `timeunit` lets you specify the unit for the returned retention time. 
# All time units defined in the package 
# [Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
# supported. The optional keyword argument `ustripped` lets you choose whether the unit is 
# included in the returned value.

# See also [`RiMapper`](@ref), [`retentiontimes`](@ref), [`maxretentiontime`](@ref), 
# [`minretentiontime`](@ref).

# # Examples
# ```jldoctest
# julia> ld = RiMapper("Kovats", [1.2, 2.4, 4.3]u"minute", [1000, 2000, 3000]);

# julia> retentionindex(ld, 1.9u"minute") ≈ 1610.7750896057348
# true

# julia> retentiontime(ld, 1610.7750896057348) ≈ 1.8999999999999995u"minute"
# true

# julia> retentiontime(ld, 1610.7750896057348, ustripped=true) ≈ 1.8999999999999995

# julia> retentiontime(ld, 900) ≈ 1.087987321711569u"minute"
# [ Info: extrapolated value
# true

# julia> rts = retentiontime.(ld, [900, 1100]);  # broadcasting across multiple RI values
# [ Info: extrapolated value

# julia> rts ≈ [1.087987321711569, 1.3120777535153498]u"minute"
# true
# ```
# """
# function retentiontime(mapper::RiMapper, retentionindex::Real; 
#     timeunit::Unitful.TimeUnits=unit(eltype(retentiontimes(mapper))), 
#     ustripped::Bool=false, info::Bool=true)

#     rt = ri2rt(mapper)(retentionindex)
#     isnothing(rt) && (return rt)
#     minretentionindex(mapper) ≤ retentionindex ≤ maxretentionindex(mapper) || (info && 
#         @info "extrapolated value")
#     ustripped ? ustrip(timeunit, rt) : uconvert(timeunit, rt)
# end


"""
    retentiontimes(mapper::RiMapper; timeunit::Unitful.TimeUnits, ustripped::Bool=false)

Return the retention times used to construct the retention index mapper. The optional 
keyword argument `timeunit` lets you change the unit of the returned retention times. All 
time units defined in the package [Unitful.jl](https://painterqubits.github.io/Unitful.jl) 
(e.g., `u"s"`, `u"minute"`) are supported. The optional keyword argument `ustripped` lets 
you choose whether to include the unit in the returned values. Note: If no time unit 
conversion is applied and the unit is not stripped, the function returns a reference to the 
data structure.

See also [`RiMapper`](@ref), [`maxretentiontime`](@ref), [`minretentiontime`](@ref), 
[`retentionindices`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 3.8]u"minute", [1000, 2000, 3000]);

julia> retentiontimes(ld)  # reference to the data structure
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:
 1.2 minute
 2.4 minute
 3.8 minute

julia> retentiontimes(ld)[:]  # a copy of these values
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:
 1.2 minute
 2.4 minute
 3.8 minute

julia> retentiontimes(ld, timeunit=u"s")
3-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
  72.0 s
 144.0 s
 228.0 s

julia> retentiontimes(ld, timeunit=u"s", ustripped=true)
3-element Vector{Float64}:
  72.0
 144.0
 228.0
```
"""
function retentiontimes(mapper::RiMapper;
    timeunit::Unitful.TimeUnits=unit(eltype(mapper.retentiontimes)), ustripped::Bool=false)
    convert = unit(eltype(mapper.retentiontimes)) ≠ timeunit
    if !convert && !ustripped
        mapper.retentiontimes
    elseif !convert && ustripped
        ustrip.(mapper.retentiontimes)
    elseif convert && !ustripped
        uconvert.(timeunit, mapper.retentiontimes)
    else
        ustrip.(timeunit, mapper.retentiontimes)
    end
end


"""
    JuChrom.rt2ri(mapper::RiMapper)

Return the function that maps a retention time to its corresponding retention index. 
Note that direct use of this function is discouraged; users are encouraged to use the 
more versatile function `retentionindex` instead.

See also [`RiMapper`](@ref), [`retentionindex`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", [1.2, 2.4, 4.3]u"minute", [1000, 2000, 3000]);

julia> JuChrom.rt2ri(ld)(1.5u"minute") ≈ 1266.712648556876
true

julia> JuChrom.rt2ri(ld)(1u"minute") ≈ 821.4487832484438
true

julia> retentionindex(ld, 1u"minute") ≈ 821.44878324844389
[ Info: extrapolated value
true
```
"""
rt2ri(mapper::RiMapper) = mapper.rt2ri


# """
#     ri2rt(mapper::RiMapper)

# Return the function that maps a retention index to the corresponding retention time. Note 
# that direct use of this function is discouraged; users are encouraged to use the more 
# versatile function `retentiontime` instead.

# See also [`RiMapper`](@ref), [`rt2ri`](@ref), [`retentiontime`](@ref).

# # Example
# ```jldoctest
# julia> ld = RiMapper("Kovats", [1.2, 2.4, 4.3]u"minute", [1000, 2000, 3000]);

# julia> rt2ri(ld)(1.9u"minute") ≈ 1610.7750896057348
# true

# julia> ri2rt(ld)(1610.7750896057348) ≈ 1.8999999999999995u"minute"
# true
# ```
# """
# ri2rt(mapper::RiMapper) = mapper.ri2rt


function criticalpoints(polynomials::AbstractVector{<:Polynomial}, ts)
    cpoints = Vector{Tuple{<:Any, String}}()
    for (i, polynomial) in enumerate(polynomials)
        c, b, a = polynomial[1:3]
        Δ₀ = b^2 - 3 * a * c
        Δ₀ ≤ 0 && continue  # no critical point or an inflection point
        α, β = ts[i], ts[i+1]
        x₀, x₁ = (-b + sqrt(Δ₀)) / 3a, (-b - sqrt(Δ₀)) / 3a
        if polynomial(x₀) < polynomial(x₁)
            α < α + x₀ < β && push!(cpoints, (α + x₀, "local minimum"))
            α < α + x₁ < β && push!(cpoints, (α + x₁, "local maximum"))
        else
            # currently not tested in unit testing!
            α < α + x₀ < β && push!(cpoints, (α + x₀, "local maximum"))
            α < α + x₁ < β && push!(cpoints, (α + x₁, "local minimum"))
        end
    end
    cpoints
end


function interextrapolate(t, timeunit, y, n, polynomials, f₁, f₂)
    time::Unitful.Time -> begin
        x = ustrip(timeunit, time)
        if x < t[1]
            f₁(x)
        elseif x > t[n+1]
            f₂(x)
        elseif x == t[1]
            Float64(y[1])
        else
            k = findlast(x .> t)
            polynomials[k](x - t[k])
        end
    end
end


function interpolate(t, timeunit, y, n, polynomials)
    time::Unitful.Time -> begin
        x = ustrip(timeunit, time)
        if x < t[1] || x > t[n+1]
            nothing
        elseif x==t[1]
            Float64(y[1])
        else
            k = findlast(x .> t)
            polynomials[k](x - t[k])
        end
    end
end


# function interextrapolate_ri(rts, timeunit, ris, n, polynomials, f₁_inverse, f₂_inverse)
#     ri::Real -> begin
#         if ri < ris[1]
#             f₁_inverse(ri) * timeunit
#         elseif ri > ris[n+1]
#             f₂_inverse(ri) * timeunit
#         elseif ri == ris[1]
#             Float64(rts[1]) * timeunit
#         else
#             k = findlast(ri .> ris)
#             d, cba = coeffs(polynomials[k])
#             println(polynomials[k])
#             println("ri: ", ri)
#             println("ri interval: ", k, " ", ris[k])

#             # Determine the roots of the polynomial
#             # rs = PolynomialRoots.roots([d - ri, cba...])
#             rs = roots(Polynomial([d - ri, cba...]))
#             println(rs)

#             # Consider only the roots with a zero imaginary part
#             rs_im_zero = filter(x -> imag(x) == 0, rs)

#             # Convert complex numbers to a real numbers
#             rs_real = map(x -> real(x), rs_im_zero)

#             # Only consider roots which fall in the kth rt interval
#             Δt = rts[k+1] - rts[k]
#             println(rts[k], " ", Δt, " ", rts[k+1])
#             println(rts[k]+rs_real[1])
#             t = filter(x -> 0 < x ≤ Δt || x ≈ Δt, rs_real)
#             println(t)
#             @assert length(t) == 1 "more than one solution found"

#             (rts[k] + t[1]) * timeunit
#         end
#     end
# end


# function interpolate_ri(rts, timeunit, ris, n, polynomials)
#     ri::Real -> begin
#         if ri < ris[1]
#             nothing
#         elseif ri > ris[n+1]
#             nothing
#         elseif ri == ris[1]
#             Float64(rts[1]) * timeunit
#         else
#             k = findlast(ri .> ris)
#             d, cba = coeffs(polynomials[k])

#             # Determine the roots of the polynomial
#             rs = PolynomialRoots.roots([d - ri, cba...])

#             # Consider only the roots with a zero imaginary part
#             rs_im_zero = filter(x -> imag(x) == 0, rs)

#             # Convert complex numbers to a real numbers
#             rs_real = map(x -> real(x), rs_im_zero)

#             # Only consider roots which fall in the kth rt interval
#             Δt = rts[k+1] - rts[k]
#             t = filter(x -> 0 < x ≤ Δt || x ≈ Δt, rs_real)
#             @assert length(t) == 1 "more than one solution found"

#             (rts[k] + t[1]) * timeunit
#         end
#     end
# end



"""
    bsplineinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
    retentionindices::AbstractVector{<:Real}; extrapolation::Bool=true)

Return a function that maps retention time to a retention index. The function uses a 
B-spline for interpolation calculated from a vector of retention times and a corresponding 
vector of retention indices. For retention time values outside the range used to 
compute the B-spline, the function employs linear extrapolation to estimate a retention 
index. However, an optional keyword argument, `extrapolation`, can be used to disable 
extrapolation, in which case the function returns nothing for values outside the retention 
time range. The function will raise an error if the resulting mapping function does not 
yield continuously increasing values.

See also [`scantimes`](@ref), [`retentionindices`](@ref).

# Examples
```jldoctest
julia> rts = [1, 2, 3, 4, 5, 6, 7, 8]*u"s";

julia> ris = [1000, 1800, 3050, 3800, 5500, 6600, 6900, 7400];

julia> rt2ri = bsplineinterpolation(rts, ris);

julia> rt2ri(1u"s") ≈ 1000.0
true

julia> rt2ri(1.5u"s") ≈ 1333.7469941600825
true

julia> rt2ri((1//30)u"minute") ≈ 1800.0
true

julia> rt2ri(9.1u"s") ≈ 7950.0
true

julia> rt2ri = bsplineinterpolation(rts, ris; extrapolation=false);

julia> rt2ri(9.1u"s") === nothing
true
```
"""
function bsplineinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
    retentionindices::AbstractVector{<:Real}; extrapolation::Bool=true)

    # Sanity checks
    length(retentiontimes) == length(retentionindices) || throw(DimensionMismatch(
        "number of retention times and number of retention indices are different"))
    length(retentiontimes) > 1 || throw(
        ArgumentError("need at least two retention time-retention index pairs"))
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
    t = ustrip.(retentiontimes)
    y = retentionindices

    polynomials = naturalbesplinepolynomials(t, y)

    cpoint = criticalpoints(polynomials, t)
    length(cpoint) == 0 || @warn string("the derived function does not ", 
        "yield continuously increasing values: $cpoint")

    n = length(t) - 1

    # Return desired interpolation function
    if extrapolation

        slope₁ = first(polynomials)[1]
        f₁(rt) = first(polynomials)[0] - slope₁ * (t[1] - rt)

        Δrt = t[end]-t[end-1]
        slope₂ = (last(polynomials)[1] + last(polynomials)[2] * (Δrt) 
            + last(polynomials)[3] * (Δrt)^2)
        f₂(rt) = last(polynomials)(Δrt) + slope₂ * (rt - t[end])

        # The slopes of both extrapolated lines must be greater than zero
        slope₁ > 0 && slope₂ > 0 || @warn string("the derived function ", 
            "does not yield continuously increasing values")

        interextrapolate(t, timeunit, y, n, polynomials, f₁, f₂)
    else
        interpolate(t, timeunit, y, n, polynomials)
    end
end


# julia> rts = [1, 2, 3, 4, 5, 6, 7, 8]*u"s";

# julia> ris = [1000, 2000, 3050, 3360, 5500, 6600, 6900, 7400];

# julia> rt2ri = bsplineinterpolation(rts, ris);
# Warning: the derived function does not yield continuously increasing values: Tuple{Any, String}[(3.3848989851524127, "local minimum"), (3.3460161285522876, "local maximum")]
# [...]
