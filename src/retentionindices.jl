"""
    PolationMethod

Supertype for all data interpolation and extrapolation methods implemented for mapping 
retention times to retention indices, and potentially vice versa.

See also [`Linear`](@ref), [`NaturalCubicBSpline`](@ref), [`RiMapper`](@ref).
"""
abstract type PolationMethod end


"""
    Linear() <: PolationMethod

PolationMethod type that specifies data is linearly extrapolated using the slope 
calculated by the applied interpolator at the nearest retention time–retention index 
calibration point.

See also [`PolationMethod`](@ref), [`NaturalCubicBSpline`](@ref), [`RiMapper`](@ref).
"""
struct Linear <: PolationMethod end 


struct NaturalCubicBSpline <: PolationMethod
    force::Bool
    NaturalCubicBSpline(force::Bool) = new(force)
end


"""
    NaturalCubicBSpline(; force::Bool=false) <: PolationMethod

PolationMethod type that specifies data is interpolated using a natural cubic B-spline. 
Application of this type will raise an error if the resulting mapping function does not 
produce continuously increasing values. However, by setting the `force` keyword argument to 
true, a mapping function will be returned even if it does not yield continuously increasing 
values. This can be useful for identifying problematic or erroneous calibration points.

See also [`PolationMethod`](@ref), [`NaturalCubicBSpline`](@ref), [`RiMapper`](@ref).

# Examples
```jldoctest
julia> NaturalCubicBSpline()
NaturalCubicBSpline(false)

julia> NaturalCubicBSpline(; force=true)
NaturalCubicBSpline(true)
```
"""
NaturalCubicBSpline(; force::Bool=false) = NaturalCubicBSpline(force)


force(obj::NaturalCubicBSpline) = obj.force


"""
    PiecewiseLinear() <: PolationMethod

PolationMethod type that specifies data is interpolated using a piecewise linear approach.

See also [`PolationMethod`](@ref), [`NaturalCubicBSpline`](@ref), [`RiMapper`](@ref).
"""
struct PiecewiseLinear <: PolationMethod end


struct RiMapper{T1<:AbstractString, T2<:AbstractVector{<:Unitful.Time}, 
    T3<:AbstractVector{<:Real}, T4<:PolationMethod, T5<:Union{<:PolationMethod, Nothing}
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
        T3<:AbstractVector{<:Real}, T4<:PolationMethod, 
        T5<:Union{<:PolationMethod, Nothing}}

        length(retentionindexname) > 0 || throw(ArgumentError(
            "no retention index name provided"))
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

    # Information about the index and the number of ladder steps
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
    retentionindices::AbstractVector{<:Real};
    interpolationmethod::PolationMethod=NaturalCubicBSpline(), 
    extrapolationmethod::Union{Nothing, <:PolationMethod}=nothing,
    metadata::Dict=Dict())

Create an `RiMapper` object to map retention times to retention indices using 
interpolation, and extrapolation by default. The optional keyword arguments 
`interpolationmethod` and `extrapolationmethod` allow you to explicitly specify the 
methods used. Currently, the available interpolators are `NaturalCubicBSpline()` (default) 
and `PiecewiseLinear()`. The only available extrapolator is `Linear()`. If 
`extrapolationmethod` is set to nothing (default), the function will raise an error for 
retention time values outside the calibration range. Note that both retention times and 
retention indices must be provided in ascending order. Additionally, the optional 
`metadata` keyword argument allows you to associate metadata with the mapper.

See also [`AbstractRiMapper`](@ref), [`retentionindexname`](@ref), [`Linear`](@ref), 
[`NaturalCubicBSpline`](@ref), [`PiecewiseLinear`](@ref), [`retentionindices`](@ref), 
[`retentiontimes`](@ref), [`interpolationmethod`](@ref), [`extrapolationmethod`](@ref), 
[`rt2ri`](@ref), [`metadata`](@ref), [`retentionindex`](@ref).

# Examples
```jldoctest
julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000)
RiMapper {index name: Kovats, calibration points: 5}
retention times: 1 minute, 2 minute, 3 minute, 4 minute, 5 minute
retention indices: 1000, 2000, 3000, 4000, 5000
interpolation method: NaturalCubicBSpline(false)
extrapolation method: nothing
metadata: 0 entries

julia> retentionindex(ld, 1u"minute") ≈ 1000.0
true

julia> retentionindex(ld, 1.5u"minute") ≈ 1500.0
true

julia> retentionindex(ld, 11u"minute") ≈ 11000.0
ERROR: ArgumentError: retention time outside range for calculating retention index
[...]

julia> retentionindices(ld)
1000:1000:5000

julia> retentiontimes(ld)
(1:5) minute

julia> interpolationmethod(ld)
NaturalCubicBSpline(false)

julia> extrapolationmethod(ld) === nothing
true

julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000, extrapolationmethod=Linear())
RiMapper {index name: Kovats, calibration points: 5}
retention times: 1 minute, 2 minute, 3 minute, 4 minute, 5 minute
retention indices: 1000, 2000, 3000, 4000, 5000
interpolation method: NaturalCubicBSpline(false)
extrapolation method: Linear()
metadata: 0 entries

julia> retentionindex(ld, 11u"minute") === 11000.000000000011
true

julia> extrapolationmethod(ld)
Linear()

julia> ld = RiMapper("Kovats", (1:5)u"s", 10:10:50, interpolationmethod=PiecewiseLinear())
RiMapper {index name: Kovats, calibration points: 5}
retention times: 1 s, 2 s, 3 s, 4 s, 5 s
retention indices: 10, 20, 30, 40, 50
interpolation method: PiecewiseLinear()
extrapolation method: nothing
metadata: 0 entries
```
"""
function RiMapper(
    retentionindexname::AbstractString,
    retentiontimes::AbstractVector{<:Unitful.Time},
    retentionindices::AbstractVector{<:Real}; 
    interpolationmethod::PolationMethod=NaturalCubicBSpline(),
    extrapolationmethod::Union{Nothing, PolationMethod}=nothing,
    metadata::Dict=Dict{Any, Any}())

    RiMapper(retentionindexname, retentiontimes, retentionindices, interpolationmethod,
        extrapolationmethod, metadata)
end

function RiMapper(
    retentionindexname::T1,
    retentiontimes::T2,
    retentionindices::T3,
    interpolationmethod::T4,
    extrapolationmethod::T5,
    metadata::Dict=Dict{Any, Any}()
    ) where {
        T1<:AbstractString,
        T2<:AbstractVector{<:Unitful.Time}, 
        T3<:AbstractVector{<:Real},
        T4<:NaturalCubicBSpline,
        T5<:Linear}

    rt2ri = bsplineinterpolation(retentiontimes, retentionindices, extrapolation=true, 
        force=force(interpolationmethod))
    RiMapper{T1, T2, T3, T4, T5}(retentionindexname, retentiontimes, retentionindices, 
        interpolationmethod, extrapolationmethod, rt2ri, metadata)
end


function RiMapper(
    retentionindexname::T1,
    retentiontimes::T2,
    retentionindices::T3,
    interpolationmethod::T4,
    extrapolationmethod::T5,
    metadata::Dict=Dict{Any, Any}()
    ) where {
        T1<:AbstractString,
        T2<:AbstractVector{<:Unitful.Time}, 
        T3<:AbstractVector{<:Real},
        T4<:NaturalCubicBSpline,
        T5<:Nothing}

    rt2ri = bsplineinterpolation(retentiontimes, retentionindices, extrapolation=false, 
        force=force(interpolationmethod))
    RiMapper{T1, T2, T3, T4, T5}(retentionindexname, retentiontimes, retentionindices, 
        interpolationmethod, extrapolationmethod, rt2ri, metadata)
end


function RiMapper(
    retentionindexname::T1,
    retentiontimes::T2,
    retentionindices::T3,
    interpolationmethod::T4,
    extrapolationmethod::T5,
    metadata::Dict=Dict{Any, Any}()
    ) where {
        T1<:AbstractString,
        T2<:AbstractVector{<:Unitful.Time}, 
        T3<:AbstractVector{<:Real},
        T4<:PiecewiseLinear,
        T5<:Linear}

    rt2ri = piecewiselinearinterpolation(retentiontimes, retentionindices, 
        extrapolation=true)
    RiMapper{T1, T2, T3, T4, T5}(retentionindexname, retentiontimes, retentionindices, 
        interpolationmethod, extrapolationmethod, rt2ri, metadata)
end


function RiMapper(
    retentionindexname::T1,
    retentiontimes::T2,
    retentionindices::T3,
    interpolationmethod::T4,
    extrapolationmethod::T5,
    metadata::Dict=Dict{Any, Any}()
    ) where {
        T1<:AbstractString,
        T2<:AbstractVector{<:Unitful.Time}, 
        T3<:AbstractVector{<:Real},
        T4<:PiecewiseLinear,
        T5<:Nothing}

    rt2ri = piecewiselinearinterpolation(retentiontimes, retentionindices, 
        extrapolation=false)
    RiMapper{T1, T2, T3, T4, T5}(retentionindexname, retentiontimes, retentionindices, 
        interpolationmethod, extrapolationmethod, rt2ri, metadata)
end


"""
    interpolationmethod(mapper::RiMapper)

Return the name of the applied interpolation method.

See also [`RiMapper`](@ref), [`extrapolationmethod`](@ref), [`retentionindex`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000);

julia> interpolationmethod(ld)
NaturalCubicBSpline(false)
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
julia> rts, ris = (1:5)u"minute", 1000:1000:5000;

julia> ld = RiMapper("Kovats", rts, ris);

julia> extrapolationmethod(ld) === nothing
true

julia> ld = RiMapper("Kovats", rts, ris, extrapolationmethod=Linear());

julia> extrapolationmethod(ld)
Linear()
```
"""
extrapolationmethod(mapper::RiMapper) = mapper.extrapolationmethod


"""
    maxretentionindex(mapper::RiMapper)

Return the highest retention index used to construct the retention index mapper.

See also [`RiMapper`](@ref), [`minretentionindex`](@ref), [`retentionindices`](@ref).

# Example
```jldoctest
julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000);

julia> maxretentionindex(ld)
5000
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
julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000);

julia> maxretentiontime(ld) ≈ 5u"minute"
true

julia> maxretentiontime(ld, timeunit=u"s") ≈ 300u"s"
true

julia> maxretentiontime(ld, timeunit=u"s", ustripped=true) ≈ 300
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
julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000, metadata=Dict(:id => 7))
RiMapper {index name: Kovats, calibration points: 5}
retention times: 1 minute, 2 minute, 3 minute, 4 minute, 5 minute
retention indices: 1000, 2000, 3000, 4000, 5000
interpolation method: NaturalCubicBSpline(false)
extrapolation method: nothing
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
julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000);

julia> maxretentionindex(ld)
5000
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
julia> ld = RiMapper("Kovats", (1:5)u"minute", 1000:1000:5000);

julia> minretentiontime(ld) ≈ 1u"minute"
true

julia> minretentiontime(ld, timeunit=u"s") ≈ 60u"s"
true

julia> minretentiontime(ld, timeunit=u"s", ustripped=true) ≈ 60
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
    retentionindex(mapper::RiMapper, retentiontime::Unitful.Time)

Return the retention index associated with a given retention time.

See also [`RiMapper`](@ref), [`retentionindices`](@ref), [`maxretentionindex`](@ref), 
[`minretentionindex`](@ref).

# Examples
```jldoctest
julia> rts, ris = [1.2, 2.4, 3.8, 5.0]u"minute", [1000, 2000, 3000, 4000];

julia> ld = RiMapper("Kovats", rts, ris)
RiMapper {index name: Kovats, calibration points: 4}
retention times: 1.2 minute, 2.4 minute, 3.8 minute, 5.0 minute
retention indices: 1000, 2000, 3000, 4000
interpolation method: NaturalCubicBSpline(false)
extrapolation method: nothing
metadata: 0 entries

julia> retentionindex(ld, 1.8u"minute") ≈ 1516.9172932330828
true

julia> retentionindex(ld, 1.1u"minute") ≈ 913.9194139194141
ERROR: ArgumentError: retention time outside range for calculating retention index
[...]

julia> ris = retentionindex.(ld, [2, 3]u"minute");  # broadcasting across multiple RT values

julia> ris ≈ [1683.3751044277362, 2430.719656283566]
true
```
"""
retentionindex(mapper::RiMapper, retentiontime::Unitful.Time) = rt2ri(mapper)(retentiontime)


"""
    retentionindex(chrom::AbstractChromatogram, retentiontime::Unitful.Time; 
    info::Bool=false)

Return the retention index corresponding to a given retention time.

See also [`RiMapper`](@ref), [`retentionindices`](@ref), [`maxretentionindex`](@ref), 
[`minretentionindex`](@ref).

# Examples
```jldoctest
julia> rts, ris = [1.2, 2.4, 3.8, 5.0]u"minute", [1000, 2000, 3000, 4000];

julia> ld = RiMapper("Kovats", rts, ris);

julia> chrom = Chrom(Int64[1, 2, 3]u"s", Int32[12, 956, 1], rimapper=ld);

julia> retentionindex(chrom, 1.8u"minute") ≈ 1516.9172932330828
true

julia> retentionindex(chrom, 1.1u"minute") ≈ 913.9194139194141
ERROR: ArgumentError: retention time outside range for calculating retention index
[...]

julia> ris = retentionindex.(chrom, [2, 3]u"minute");  # broadcasting across multiple RT values

julia> ris ≈ [1683.3751044277362, 2430.719656283566]
true

julia> chrom = Chrom(Int64[1, 2, 3]u"s", Int32[12, 956, 1]);

julia> retentionindex(chrom, 120u"s")
ERROR: ArgumentError: no retention index mapper implemented
[...]
```
"""
function retentionindex(chrom::AbstractChromatogram, retentiontime::Unitful.Time)
    isnothing(rimapper(chrom)) && throw(ArgumentError(
        "no retention index mapper implemented"))
    retentionindex(rimapper(chrom), retentiontime)
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
julia> rts, ris = [1.2, 2.4, 3.8, 5.0]u"minute", [1000, 2000, 3000, 4000];

julia> ld = RiMapper("Kovats", rts, ris);

julia> retentionindices(ld)  # reference to the data structure
4-element Vector{Int64}:
 1000
 2000
 3000
 4000

julia> retentionindices(ld)[:]  # a copy of these values
4-element Vector{Int64}:
 1000
 2000
 3000
 4000
```
"""
retentionindices(mapper::RiMapper) = mapper.retentionindices


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
julia> rts, ris = [1.2, 2.4, 3.8, 5.0]u"minute", [1000, 2000, 3000, 4000];

julia> ld = RiMapper("Kovats", rts, ris);

julia> retentiontimes(ld)  # reference to the data structure
4-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:
 1.2 minute
 2.4 minute
 3.8 minute
 5.0 minute

julia> retentiontimes(ld)[:]  # a copy of these values
4-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:
 1.2 minute
 2.4 minute
 3.8 minute
 5.0 minute

julia> retentiontimes(ld, timeunit=u"s")
4-element Vector{Quantity{Float64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:
  72.0 s
 144.0 s
 228.0 s
 300.0 s

julia> retentiontimes(ld, timeunit=u"s", ustripped=true)
4-element Vector{Float64}:
  72.0
 144.0
 228.0
 300.0
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
julia> rts, ris = [1.2, 2.4, 3.8, 5.0]u"minute", [1000, 2000, 3000, 4000];

julia> ld = RiMapper("Kovats", rts, ris);

julia> JuChrom.rt2ri(ld)(1.5u"minute") ≈ 1260.5733082706768
true

julia> JuChrom.rt2ri(ld)(1u"minute") ≈ 821.4487832484438
ERROR: ArgumentError: retention time outside range for calculating retention index
[...]
```
"""
rt2ri(mapper::RiMapper) = mapper.rt2ri


function criticalpoints(spline, nodes, ϵ)
    S′ = BSplineKit.diff(spline, BSplineKit.Derivative(1))
    S″ = BSplineKit.diff(spline, BSplineKit.Derivative(2))
    critical_points = Vector{Float64}()
    for i in 1:(length(nodes)-1)
        result = Roots.find_zeros(S′, nodes[i], nodes[i+1])
        for x in result
            abs(S″(x)) > ϵ && push!(critical_points, x)
        end
    end
    critical_points
end


"""
    bsplineinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
    retentionindices::AbstractVector{<:Real}; extrapolation::Bool=false, force::Bool=false)

Return a function that maps retention time to a retention index. The function uses a 
B-spline for interpolation calculated from a vector of retention times and a corresponding 
vector of retention indices. For retention time values outside the range used to 
compute the B-spline, the function employs linear extrapolation to estimate a retention 
index. However, an optional keyword argument, `extrapolation`, can be used to disable 
extrapolation, in which case the function returns nothing for values outside the retention 
time range. The function will raise an error if the resulting mapping function does not 
produce continuously increasing values. However, by setting the force keyword argument to 
true, the function will return a mapping function even if it does not yield continuously 
increasing values. This can be useful for identifying problematic or erroneous calibration 
points.

See also [`scantimes`](@ref), [`retentionindices`](@ref).

# Examples
```jldoctest
julia> rts, ris = (1:8)*u"s", [1000, 1800, 3050, 3800, 5500, 6600, 6900, 7400];

julia> rt2ri = JuChrom.bsplineinterpolation(rts, ris);

julia> rt2ri(1u"s") ≈ 1000.0
true

julia> rt2ri(1.5u"s") ≈ 1333.7469941600825
true

julia> rt2ri((1//30)u"minute") ≈ 1800.0
true

julia> rt2ri(9.1u"s") ≈ 7950.0
ERROR: ArgumentError: retention time outside range for calculating retention index
[...]

julia> rt2ri = JuChrom.bsplineinterpolation(rts, ris; extrapolation=true);

julia> rt2ri(9.1u"s") ≈ 8053.1226382686355
true
```
"""
function bsplineinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
    retentionindices::AbstractVector{<:Real}; extrapolation::Bool=false, 
    force::Bool=false)

    # Sanity checks
    length(retentiontimes) == length(retentionindices) || throw(DimensionMismatch(
        "number of retention times and number of retention indices are different"))
    length(retentiontimes) ≥ 4 || throw(
        ArgumentError("need at least four retention time-retention index pairs"))
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
    rts_ustripped = convert(Vector{Float64}, ustrip.(retentiontimes))
    ris = retentionindices

    S = BSplineKit.interpolate(rts_ustripped, ris, BSplineKit.BSplineOrder(4),
        BSplineKit.Natural()) 

    # Check for critical points
    ϵ = (maximum(retentionindices) - minimum(retentionindices)) * 1e-10
    cpoints = criticalpoints(S, rts_ustripped, ϵ)
    if length(cpoints) > 0 && !force
        throw(ArgumentError(
            "computed Bspline contains critical point(s): $cpoints * $timeunit"))
    end

    if extrapolation
        return rt::Unitful.Time -> begin
            rt_ustripped = ustrip(timeunit, rt)
            BSplineKit.extrapolate(S, BSplineKit.Linear())(rt_ustripped)
        end
    else
        return rt::Unitful.Time -> begin
            rt_ustripped = ustrip(timeunit, rt)
            if first(retentiontimes) ≤ rt ≤ last(retentiontimes)
                S(rt_ustripped)
            else
                throw(ArgumentError("retention time outside range for calculating "
                    * "retention index"))
            end
        end
    end
end


"""
    piecewiselinearinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
    retentionindices::AbstractVector{<:Real}; extrapolation::Bool=false) -> Float64

Return a function that maps retention time to retention index using piecewise linear 
interpolation based on a vector of retention times and a corresponding vector of retention 
indices. For retention time values outside the interpolation range, the function applies 
linear extrapolation to estimate the retention index. However, an optional extrapolation 
keyword can disable extrapolation, in which case the function will raise an error for 
values outside the retention time range.

See also [`scantimes`](@ref), [`retentionindices`](@ref).

# Examples
```jldoctest
julia> rts, ris = (1:8)*u"s", [1000, 1800, 3050, 3800, 5500, 6600, 6900, 7400];

julia> rt2ri = JuChrom.piecewiselinearinterpolation(rts, ris);

julia> rt2ri(1u"s") ≈ 1000.0
true

julia> rt2ri(1.5u"s") ≈ 1400.0
true

julia> rt2ri((1//30)u"minute") ≈ 1800.0
true

julia> rt2ri(9.1u"s")
ERROR: ArgumentError: retention time outside range for calculating retention index
[...]

julia> rt2ri = JuChrom.piecewiselinearinterpolation(rts, ris; extrapolation=true);

julia> rt2ri(9.1u"s") ≈ 7950.0
true
```
"""
function piecewiselinearinterpolation(retentiontimes::AbstractVector{<:Unitful.Time}, 
    retentionindices::AbstractVector{<:Real}; extrapolation::Bool=false)

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
    rts_ustripped = ustrip.(retentiontimes)
    ris = retentionindices

    # Return desired interpolation function
    if extrapolation
        L = BasicInterpolators.LinearInterpolator(rts_ustripped, ris, 
            BasicInterpolators.NoBoundaries())
    else
        L = BasicInterpolators.LinearInterpolator(rts_ustripped, ris, 
            BasicInterpolators.StrictBoundaries())
    end

    if extrapolation
        return rt::Unitful.Time -> begin
            rt_ustripped = ustrip(timeunit, rt)
            L(rt_ustripped)
        end
    else
        return rt::Unitful.Time -> begin
            rt_ustripped = ustrip(timeunit, rt)
            if first(retentiontimes) ≤ rt ≤ last(retentiontimes)
                L(rt_ustripped)
            else
                throw(ArgumentError("retention time outside range for calculating "
                    * "retention index"))
            end
        end
    end
end
