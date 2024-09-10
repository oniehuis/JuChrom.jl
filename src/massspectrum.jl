"""
    AbstractMassSpectrum

Supertype for all mass spectrum implementations. All subtypes (e.g., `MassSpectrum`) 
include fields to store ions, their associated intensities, and optionally, retention time, 
retention index, and metadata.

See also [`MassSpectrum`](@ref), [`ChromMS`](@ref).
"""
abstract type AbstractMassSpectrum end


struct MassSpectrum{
    T1<:AbstractVector{<:Real},
    T2<:AbstractVector{<:Real},
    T3<:Union{Unitful.Time, Nothing},
    T4<:Union{AbstractString, Nothing},
    T5<:Union{<:Real, Nothing}} <: AbstractMassSpectrum
    ions::T1
    intensities::T2
    retentiontime::T3
    retentionindexname::T4
    retentionindex::T5
    metadata::Dict{Any, Any}

    function MassSpectrum{T1, T2, T3, T4, T5}(
        ions::T1,
        intensities::T2,
        retentiontime::T3,
        retentionindexname::T4,
        retentionindex::T5,
        metadata::Dict
        ) where {
            T1<:AbstractVector{<:Real},
            T2<:AbstractVector{<:Real},
            T3<:Union{Unitful.Time, Nothing},
            T4<:Union{AbstractString, Nothing},
            T5<:Union{<:Real, Nothing}}

        length(ions) > 0 || throw(ArgumentError("no ion(s) provided"))
        length(intensities) > 0 || throw(ArgumentError("no intensity value(s) provided"))
        length(ions) == length(intensities) || throw(DimensionMismatch(
            "ion count does not match intensity count"))
        length(Set(ions)) == length(ions) || throw(
            ArgumentError("ions contain identical values"))
        issorted(ions) || throw(ArgumentError("ions not in ascending order"))
        count(i -> i < 0, intensities) == 0 || throw(ArgumentError(
            "intensity values contain at least one value less than zero"))
        !isnothing(retentionindex) && isnothing(retentionindexname) && throw(ArgumentError(
            "retention index value is not paired with a retention index name"))
        !isnothing(retentionindexname) && isnothing(retentionindex) && throw(ArgumentError(
                "retention index name is not paired with a retention index value"))

        new(ions, intensities, retentiontime, retentionindexname, retentionindex, metadata)
    end
end


Base.broadcastable(ms::MassSpectrum) = Ref(ms)


function Base.show(io::IO, ms::AbstractMassSpectrum)
    # Information about the element types
    typename = name(typeof(ms))
    print(io, "$typename {ions: ", eltype(ions(ms))) 
    println(io, ", intensities: ", eltype(intensities(ms)), "}")

    # Information about the ions
    if ioncount(ms) > 10
        println(io, ioncount(ms), " ions; range: m/z ", minion(ms), " - ", maxion(ms))
    elseif 1 < ioncount(ms) ≤ 10
        print(io, ioncount(ms), " ions: m/z ")
        for (i, ion) in enumerate(ions(ms))
            print(io, ion)
            i < ioncount(ms) ? print(io, ", ") : println(io)
        end
    else
        println(io, ioncount(ms), " ion: m/z ", ion(ms, 1))
    end

    # Information about the intensities
    if length(intensities(ms)) > 1
        println(io, "intensity range: ", minintensity(ms), " - ", maxintensity(ms))
    else
        println(io, "intensity: ", first(intensities(ms)))
    end

    # Information about the retention time, if it is available
    if !isnothing(retentiontime(ms))
        println(io, "retention time: ", retentiontime(ms))
    end

    # Information about the retention index, if it is available
    if !isnothing(retentionindexname(ms))
        println(io, "retention index: ", retentionindex(ms), " (", retentionindexname(ms), 
            ")")
    end

    # Information about the metadata
    n = length(metadata(ms))
    print(io, "metadata: ", n, (n == 0 || n > 1) ? " entries" : " entry" )
end


"""
    MassSpectrum(MassSpectrum(ions::AbstractVector{<:Real}, 
    intensities::AbstractVector{<:Real};
    retentiontime::Union{Unitful.Time, Nothing}=nothing,
    retentionindexname::Union{AbstractString, Nothing}=nothing,
    retentionindex::Union{<:Real, Nothing}=nothing, metadata::Dict=Dict())

Construct a `ChromMS` object that includes `scantimes`, `ions`, `intensities`, and 
`metadata`. Note that `scantimes` and `ions` must be in ascending order, and `intensities` 
must not contain any values less than zero. The `retentiontime` must include a time unit. 
All time units supported by the 
[Unitful.jl package](https://painterqubits.github.io/Unitful.jl)(e.g., `u"s"`, `u"minute"`) 
are accepted.

See also [`AbstractMassSpectrum`](@ref).

# Examples
```jldoctest
julia> MassSpectrum([85.1, 112.2, 124.1], Int64[13, 0, 67])
MassSpectrum {ions: Float64, intensities: Int64}
3 ions: m/z 85.1, 112.2, 124.1
intensity range: 0 - 67
metadata: 0 entries

julia> MassSpectrum([85.1, 112.2, 124.1], Int64[13, 0, 67], retentiontime=3.2u"minute")
MassSpectrum {ions: Float64, intensities: Int64}
3 ions: m/z 85.1, 112.2, 124.1
intensity range: 0 - 67
retention time: 3.2 minute
metadata: 0 entries

julia> MassSpectrum([1.2, 2.1], [1.1, 8.1], retentionindex=131, retentionindexname="Kovats")
MassSpectrum {ions: Float64, intensities: Float64}
2 ions: m/z 1.2, 2.1
intensity range: 1.1 - 8.1
retention index: 131 (Kovats)
metadata: 0 entries

julia> MassSpectrum([85.1, 112.2, 124.1], [13.0, 0.0, 67.0], metadata=Dict(:name=>"Hexane"))
MassSpectrum {ions: Float64, intensities: Float64}
3 ions: m/z 85.1, 112.2, 124.1
intensity range: 0.0 - 67.0
metadata: 1 entry
```

"""
function MassSpectrum(
    ions::T1,
    intensities::T2;
    retentiontime::T3=nothing,
    retentionindexname::T4=nothing,
    retentionindex::T5=nothing,
    metadata::Dict=Dict()) where {
        T1<:AbstractVector{<:Real},
        T2<:AbstractVector{<:Real},
        T3<:Union{Unitful.Time, Nothing},
        T4<:Union{AbstractString, Nothing},
        T5<:Union{<:Real, Nothing}}
    MassSpectrum{T1, T2, T3, T4, T5}(ions, intensities, retentiontime, retentionindexname, 
        retentionindex, metadata)
end


"""
    intensities(ms::AbstractMassSpectrum)

Return the intensities.

See also [`AbstractMassSpectrum`](@ref).


# Example
```jldoctest
julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> intensities(ms)
3-element Vector{Int64}:
 13
  0
 67
```
"""
intensities(ms::AbstractMassSpectrum) = ms.intensities


"""
    intensity(ms::AbstractMassSpectrum, ionindex::Integer)

Return the intensity for an ion by specifying its index.

See also [`AbstractMassSpectrum`](@ref).

# Examples
```jldoctest
julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> intensity(ms, 1)
13

julia> intensity(ms, 2)
0

julia> intensity(ms, 4)
ERROR: BoundsError: attempt to access 3-element Vector{Int64} at index [4]
[...]
```
"""
intensity(ms::AbstractMassSpectrum, ionindex::Integer) = intensities(ms)[ionindex]



"""
    ion(gcms::AbstractMassSpectrum, ionindex::Integer)

Return the ion at the specified `ionindex`.

# Examples
```jldoctest
julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> ion(ms, 1)
85

julia> ion(ms, 2)
112

julia> ion(ms, 3)
124

julia> ion(ms, 4)
ERROR: BoundsError: attempt to access 3-element Vector{Int64} at index [4]
[...]

julia> ion.(ms, 1:2)  # broadcasting across multiple ion indices
2-element Vector{Int64}:
  85
 112
```

See also [`AbstractMassSpectrum`](@ref), [`ions`](@ref), [`ionindex`](@ref), [`minion`](@ref), [`maxion`](@ref), [`ioncount`](@ref).
"""
function ion(ms::AbstractMassSpectrum, ionindex::Integer)
    firstindex(ions(ms)) ≤ ionindex ≤ lastindex(ions(ms)) || throw(
        BoundsError(ions(ms), ionindex))
    ions(ms)[ionindex]
end


"""
    ionindex(ms::AbstractMassSpectrum, ion::Real) -> Int

Return the index of the `ion`. If the `ion` is not present in the mass spectrum, an error 
is raised.

# Example
```jldoctest
julia> ms = MassSpectrum([85.1, 112.0, 124.2], Int64[13, 0, 67])
MassSpectrum {ions: Float64, intensities: Int64}
3 ions: m/z 85.1, 112.0, 124.2
intensity range: 0 - 67
metadata: 0 entries

julia> ionindex(ms, 124.2)
3

julia> ionindex(ms, 124)
ERROR: ArgumentError: ion 124 does not exist
[...]
```

See also [`AbstractMassSpectrum`](@ref), [`ions`](@ref), [`ion`](@ref), [`minion`](@ref), [`maxion`](@ref), [`ioncount`](@ref).
"""
function ionindex(chrom::AbstractMassSpectrum, ion::Real)
    for (index, element) in enumerate(ions(chrom))
        element ≈ ion && return index
    end
    throw(ArgumentError("ion $ion does not exist"))
end


"""
    ions(ms::AbstractMassSpectrum)

Return the ions.

See also [`AbstractMassSpectrum`](@ref).

# Example
```jldoctest
julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> ions(ms)
3-element Vector{Int64}:
  85
 112
 124
```
"""
ions(ms::AbstractMassSpectrum) = ms.ions


"""
    ioncount(ms::AbstractMassSpectrum) -> Int

Return the number of ions.

See also [`AbstractMassSpectrum`](@ref).

# Example
```jldoctest
julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> ioncount(ms)
3
```
"""
ioncount(ms::AbstractMassSpectrum) = length(ions(ms))


"""
    massspectrum(gcms::AbstractChromMS, scanindex::Integer)

Return a MassSpectrum with the intensity values for all ions at the specified scan index.

See also [`AbstractChromMS`](@ref), [`MassSpectrum`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> massspectrum(chrom, 2)
MassSpectrum {ions: Int64, intensities: Int64}
2 ions: m/z 85, 100
intensity range: 34 - 956
retention time: 2 s
metadata: 0 entries

julia> massspectrum(chrom, 2, retentionindexname="Kovats", retentionindex=123.21)
MassSpectrum {ions: Int64, intensities: Int64}
2 ions: m/z 85, 100
intensity range: 34 - 956
retention time: 2 s
retention index: 123.21 (Kovats)
metadata: 0 entries

julia> ms = massspectrum(chrom, 2, metadata=Dict(:compound => "hexane"))
MassSpectrum {ions: Int64, intensities: Int64}
2 ions: m/z 85, 100
intensity range: 34 - 956
retention time: 2 s
metadata: 1 entry

julia> metadata(ms)
Dict{Any, Any} with 1 entry:
  :compound => "hexane"
```
"""
function massspectrum(chrom::AbstractChromMS, scanindex::Integer;  
    retentionindexname::Union{AbstractString, Nothing}=nothing,
    retentionindex::Union{Real, Nothing}=nothing, metadata::Dict=Dict(),)
    1 ≤ scanindex ≤ scancount(chrom) || throw(BoundsError(scantimes, scanindex))
    MassSpectrum(ions(chrom), intensities(chrom)[scanindex, :], 
        retentionindexname=retentionindexname, retentionindex=retentionindex,
        retentiontime=scantime(chrom, scanindex), metadata=metadata)
end


"""
    massspectrum(gcms::AbstractChromMS, time::Unitful.Time; precisetime::Bool=false)

Return a MassSpectrum containing the intensity values for all ions in the scan with a 
scan time closest to the specified `time``. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., u"s", u"minute") are 
supported. In case of a tie, the larger scan time is selected. If the optional parameter 
`precisetime` is set to true, the specified time must match exactly, otherwise an error 
will be thrown.

See also [`AbstractChromMS`](@ref), [`MassSpectrum`](@ref).

# Examples
```jldoctest
julia> chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
ChromMS {scan times: Int64, ions: Int64, intensities: Int64}
3 scans; scan times: 1 s, 2 s, 3 s
2 ions: m/z 85, 100
intensity range: 0 - 956
metadata: 0 entries

julia> massspectrum(chrom, 2u"s")
MassSpectrum {ions: Int64, intensities: Int64}
2 ions: m/z 85, 100
intensity range: 34 - 956
retention time: 2 s
metadata: 0 entries

julia> massspectrum(chrom, 2.5u"s")
MassSpectrum {ions: Int64, intensities: Int64}
2 ions: m/z 85, 100
intensity range: 1 - 23
retention time: 3 s
metadata: 0 entries

julia> massspectrum(chrom, 3.5u"s")
MassSpectrum {ions: Int64, intensities: Int64}
2 ions: m/z 85, 100
intensity range: 1 - 23
retention time: 3 s
metadata: 0 entries

julia> massspectrum(chrom, 2.1u"s", precisetime=true)
ERROR: ArgumentError: scantime 2.1 s does not exist
[...]
```
"""
function massspectrum(chrom::AbstractChromMS, time::Unitful.Time; precisetime::Bool=false)
    massspectrum(chrom, scantimeindex(chrom, time; precisetime=precisetime))
end


"""
    maxintensity(ms::AbstractMassSpectrum)

Return the maximum intensity.

See also [`AbstractMassSpectrum`](@ref).

# Example
```jldoctest
julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> maxintensity(ms)
67
```
"""
maxintensity(ms::AbstractMassSpectrum) = maximum(intensities(ms))


"""
    maxion(ms::AbstractMassSpectrum)

Return the largest ion.

See also [`AbstractMassSpectrum`](@ref).

# Example
```jldoctest
julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> maxion(ms)
124
```
"""
maxion(ms::AbstractMassSpectrum) = last(ions(ms))


"""
    metadata(ms::AbstractMassSpectrum) -> Dict{Any, Any}

Return the metadata.

See also [`AbstractMassSpectrum`](@ref).

# Examples
```jldoctest
julia> ms = MassSpectrum([85.1], [13.0], metadata=Dict(:name=>"unknown"))
MassSpectrum {ions: Float64, intensities: Float64}
1 ion: m/z 85.1
intensity: 13.0
metadata: 1 entry

julia> metadata(ms)
Dict{Any, Any} with 1 entry:
  :name => "unknown"

julia> ms = MassSpectrum([85.1], [13.0])
MassSpectrum {ions: Float64, intensities: Float64}
1 ion: m/z 85.1
intensity: 13.0
metadata: 0 entries

julia> metadata(ms)
Dict{Any, Any}()

julia> metadata(ms)[:species] = "Polistes dominula"
"Polistes dominula"

julia> metadata(ms)[:id] = 123
123

julia> metadata(ms)
Dict{Any, Any} with 2 entries:
  :species => "Polistes dominula"
  :id      => 123

julia> ms
MassSpectrum {ions: Float64, intensities: Float64}
1 ion: m/z 85.1
intensity: 13.0
metadata: 2 entries

julia> delete!(metadata(ms), :species)
Dict{Any, Any} with 1 entry:
  :id => 123

julia> ms
MassSpectrum {ions: Float64, intensities: Float64}
1 ion: m/z 85.1
intensity: 13.0
metadata: 1 entry
```
"""
metadata(ms::AbstractMassSpectrum) = ms.metadata


"""
    minintensity(ms::AbstractMassSpectrum)

Return the minimum intensity.

See also [`AbstractMassSpectrum`](@ref).

# Example
```jldoctest
julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> minintensity(ms)
0
```
"""
minintensity(ms::AbstractMassSpectrum) = minimum(intensities(ms))


"""
    minion(ms::AbstractMassSpectrum)

Return the smallest ion.

See also [`AbstractMassSpectrum`](@ref).

# Example
```jldoctest
julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> minion(ms)
85
```
"""
minion(ms::AbstractMassSpectrum) = first(ions(ms))


"""
    retentionindex(ms::AbstractMassSpectrum)

Return the retention index. If no retention index is associated with the mass spectrum, the 
function returns the value nothing.

See also [`AbstractMassSpectrum`](@ref).

# Example
```jldoctest
julia> ms = MassSpectrum([1.2], [1.1], retentionindex=131.1, retentionindexname="Kovats")
MassSpectrum {ions: Float64, intensities: Float64}
1 ion: m/z 1.2
intensity: 1.1
retention index: 131.1 (Kovats)
metadata: 0 entries

julia> retentionindex(ms) ≈ 131.1
true

julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> retentionindex(ms) === nothing
true
```
"""
retentionindex(ms::AbstractMassSpectrum) = ms.retentionindex


"""
    retentionindexname(ms::AbstractMassSpectrum)

Return the retention index name. If the mass spectrum has no associated retention index, 
the function returns the value nothing.

See also [`AbstractMassSpectrum`](@ref).

# Example
```jldoctest
julia> ms = MassSpectrum([1.2], [1.1], retentionindex=131.1, retentionindexname="Kovats")
MassSpectrum {ions: Float64, intensities: Float64}
1 ion: m/z 1.2
intensity: 1.1
retention index: 131.1 (Kovats)
metadata: 0 entries

julia> retentionindexname(ms)
"Kovats"

julia> ms = MassSpectrum(Int64[85, 112, 124], Int64[13, 0, 67])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 85, 112, 124
intensity range: 0 - 67
metadata: 0 entries

julia> retentionindexname(ms) === nothing
true
```
"""
retentionindexname(ms::AbstractMassSpectrum) = ms.retentionindexname


"""
    retentiontime(ms::AbstractMassSpectrum; timeunit::Unitful.TimeUnits, 
    ustripped::Bool=false)

Return the retention time. If no retention time is associated with the mass spectrum, the 
function returns the value nothing. The optional keyword argument `timeunit` lets you 
change the unit of the returned scan time. All time units defined in the package 
[Unitful.jl](https://painterqubits.github.io/Unitful.jl) (e.g., `u"s"`, `u"minute"`) are 
supported. The optional keyword argument `ustripped` lets you choose whether to include the 
unit in the returned value.

See also [`AbstractMassSpectrum`](@ref).

# Examples
```jldoctest
julia> ms = MassSpectrum([85.1, 112.2, 124.1], Int64[13, 0, 67], retentiontime=3.2u"minute")
MassSpectrum {ions: Float64, intensities: Int64}
3 ions: m/z 85.1, 112.2, 124.1
intensity range: 0 - 67
retention time: 3.2 minute
metadata: 0 entries

julia> retentiontime(ms) ≈ 3.2u"minute"
true

julia> retentiontime(ms, timeunit=u"s") ≈ 192.0u"s"
true

julia> retentiontime(ms, timeunit=u"s", ustripped=true) ≈ 192.0
true

julia> retentiontime(ms, ustripped=true) ≈ 3.2
true

julia> ms = MassSpectrum([85.1, 112.2, 124.1], Int64[13, 0, 67])
MassSpectrum {ions: Float64, intensities: Int64}
3 ions: m/z 85.1, 112.2, 124.1
intensity range: 0 - 67
metadata: 0 entries

julia> retentiontime(ms) === nothing
true

julia> retentiontime(ms, timeunit=u"ms", ustripped=true) === nothing
true
```
"""
function retentiontime(ms::AbstractMassSpectrum; timeunit::Union{Unitful.TimeUnits, 
    Nothing}=nothing, ustripped::Bool=false)
    isnothing(ms.retentiontime) && return nothing
    if isnothing(timeunit)
        ustripped ? ustrip(ms.retentiontime) : ms.retentiontime
    else
        ustripped ? ustrip(timeunit, ms.retentiontime) : uconvert(timeunit, 
            ms.retentiontime)
    end
end


"""
    sharedions(ms₁::AbstractMassSpectrum, ms₂::AbstractMassSpectrum)

Return the set of ions present in both mass spectra.

See also [`MassSpectrum`](@ref), [`ions`](@ref).

# Example
```jldoctest
julia> ms₁ = MassSpectrum([80, 85, 90], [100, 500, 250])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 80, 85, 90
intensity range: 100 - 500
metadata: 0 entries

julia> ms₂ = MassSpectrum([80, 85, 90], [100, 500, 250])
MassSpectrum {ions: Int64, intensities: Int64}
3 ions: m/z 80.0, 85.0, 90.1
intensity range: 100 - 500
metadata: 0 entries

julia> sharedions(ms₁, ms₂)
2-element Vector{Int64}:
 80
 85
```
"""
sharedions(ms₁::AbstractMassSpectrum, ms₂::AbstractMassSpectrum) = intersect(
    ions(ms₁), ions(ms₂))


"""
    similarity(ms₁::AbstractMassSpectrum, ms₂::AbstractMassSpectrum, f::Function)

Compute the similarity between the two mass spectra by applying the similarity function `f` 
(e.g., `cosine`) to the intensity values of the ions shared between `ms₁` and `ms₂`.

See also [`AbstractMassSpectrum`](@ref), [`cosine`](@ref).

# Examples
```jldoctest
julia> ms₁ = MassSpectrum(Int64[80, 85, 90], Int64[100, 500, 250]);

julia> ms₂ = MassSpectrum(Int32[80, 85], Int32[200, 1000]);

julia> similarity(ms₁, ms₂, cosine) ≈ 1.0

julia> ms₁ = MassSpectrum(Float64[80, 85, 90], Float32[100, 500, 250]);

julia> ms₂ = MassSpectrum(Int64[80, 85, 90], Int32[200, 1000, 0]);

julia> similarity(ms₁, ms₂, cosine) ≈ 0.8978873122816586
true
```
"""
function similarity(ms₁::AbstractMassSpectrum, ms₂::AbstractMassSpectrum, f::Function)
    ions = sharedions(ms₁, ms₂)
    f(intensity.(ms₁, ionindex.(ms₁, ions)), intensity.(ms₂, ionindex.(ms₂, ions)))
end
