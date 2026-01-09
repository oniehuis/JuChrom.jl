"""
    abstract type AbstractScanSeries{S, R, M, I} <: Any

Abstract supertype for objects representing a series of scan measurements, such as 
chromatograms or mass spectra.

`S` is the concrete scan type (a subtype of `AbstractScan`), `R` is the separation unit type 
(`Unitful.Units` subtype or `Nothing`), `M` is the m/z unit type (`Unitful.Units` subtype or 
`Nothing`), and `I` is the signal intensity unit type (`Unitful.Units` subtype or `Nothing`).

Concrete subtypes are expected to define `scans::AbstractVector{S}` along with
`instrument`, `acquisition`, `sample`, and `user` metadata as `NamedTuple`s, plus
`extras::Dict{String, Any}` for unstructured metadata. Subtypes may define additional
fields as needed.

See also [`AbstractChromScan`](@ref), [`AbstractMassScan`](@ref), [`scans`](@ref), 
[`instrument`](@ref), [`acquisition`](@ref), [`user`](@ref), [`sample`](@ref), 
[`extras`](@ref).
"""
abstract type AbstractScanSeries{S, R, M, I} end

Base.length(series::AbstractScanSeries) = length(scans(series))
Base.iterate(series::AbstractScanSeries) = iterate(scans(series))
Base.iterate(series::AbstractScanSeries, state) = iterate(scans(series), state)
Base.getindex(series::AbstractScanSeries, i::Integer) = scans(series)[i]
Base.eltype(series::AbstractScanSeries) = eltype(scans(series))
Base.first(series::AbstractScanSeries) = first(scans(series))
Base.last(series::AbstractScanSeries) = last(scans(series))
Base.getproperty(series::AbstractScanSeries, name::Symbol) = getfield(series, name)

"""
    abstract type AbstractChromScanSeries{S<:AbstractChromScan, R, I} 
        <: AbstractScanSeries{S, R, Nothing, I}

Abstract supertype for a series of chromatographic scans. This type represents a sequence 
of `AbstractChromScan` objects with associated metadata. It specializes `AbstractScanSeries` 
by setting the m/z unit type parameter `M` to `Nothing`, since chromatographic scans do not 
include m/z values. 

`S` is the concrete scan type (a subtype of `AbstractChromScan`), `R` is the separation unit 
type (`Unitful.Units` subtype or `Nothing`), and `I` is the signal intensity unit type 
(`Unitful.Units` subtype or `Nothing`). 

Concrete subtypes are expected to define `scans::AbstractVector{S}` along with `instrument`,
`acquisition`, `sample`, and `user` metadata as `NamedTuple`s, plus 
`extras::Dict{String, Any}`.

See also [`AbstractScanSeries`](@ref), [`AbstractChromScan`](@ref), [`scans`](@ref), 
[`instrument`](@ref), [`acquisition`](@ref), [`user`](@ref), [`sample`](@ref), 
[`extras`](@ref).
"""
abstract type AbstractChromScanSeries{S<:AbstractChromScan, R, I} <: AbstractScanSeries{S, R, Nothing, I} end

"""
    abstract type AbstractMassScanSeries{S<:AbstractMassScan, R, M, I} 
        <: AbstractScanSeries{S, R, M, I}

Abstract supertype for a series of mass spectrometry scans. This type represents a sequence 
of `AbstractMassScan` objects with associated metadata. It specializes `AbstractScanSeries` 
by including the m/z unit type parameter `M`, which is relevant for mass spectrometry data. 

`S` is the concrete scan type (a subtype of `AbstractMassScan`), `R` is the separation
unit type (`Unitful.Units` subtype or `Nothing`), `M` is the m/z unit type (`Unitful.Units` 
subtype or `Nothing`), and `I` is the signal intensity unit type (`Unitful.Units` subtype 
or `Nothing`). 

Concrete subtypes are expected to define `scans::AbstractVector{S}` along with `instrument`, 
`acquisition`, `sample`, and `user` metadata as `NamedTuple`s, plus 
`extras::Dict{String, Any}`.

See also [`AbstractScanSeries`](@ref), [`AbstractMassScan`](@ref), [`scans`](@ref), 
[`instrument`](@ref), [`acquisition`](@ref), [`user`](@ref), [`sample`](@ref), 
[`extras`](@ref).
"""
abstract type AbstractMassScanSeries{S<:AbstractMassScan, R, M, I} <: AbstractScanSeries{S, R, M, I} end

struct ChromScanSeries{
    S<:AbstractChromScan,
    R,
    I,
    T1<:AbstractVector{S},
    T2<:NamedTuple,
    T3<:NamedTuple,
    T4<:NamedTuple,
    T5<:NamedTuple,
    T6<:Dict{String, Any}
    } <: AbstractChromScanSeries{S, R, I}

    scans::T1
    instrument::T2
    acquisition::T3
    user::T4
    sample::T5
    extras::T6

    function ChromScanSeries{S, R, I, T1, T2, T3, T4, T5, T6}(
        scans::T1,
        instrument::T2,
        acquisition::T3,
        user::T4,
        sample::T5,
        extras::T6
    ) where {
        S<:AbstractChromScan,
        R,
        I,
        T1<:AbstractVector{S},
        T2<:NamedTuple,
        T3<:NamedTuple,
        T4<:NamedTuple,
        T5<:NamedTuple,
        T6<:Dict{String, Any}
    }

        length(scans) > 0 || throw(ArgumentError("no scan(s) provided"))

        E = eltype(scans)
        isconcretetype(E) || throw(ArgumentError(
            "scans must have a concrete element type, got $E"))

        new{S, R, I, T1, T2, T3, T4, T5, T6}(
            scans, instrument, acquisition, user, sample, extras)
    end
end

"""
    ChromScanSeries(
        scans::Vector{<:AbstractChromScan};
        instrument::NamedTuple=NamedTuple(),
        acquisition::NamedTuple=NamedTuple(),
        user::NamedTuple=NamedTuple(),
        sample::NamedTuple=NamedTuple(),
        extras::Dict{<:AbstractString, <:Any}=Dict()
    ) -> ChromScanSeries

Constructs a `ChromScanSeries`, representing an ordered collection of scans (subtypes of 
`AbstractChromScan`) with associated metadata.

`scans` must be a non-empty vector of chromatographic scan objects, all of the same
concrete subtype of `AbstractChromScan`. `instrument`, `acquisition`, `user`, and `sample`
accept optional `NamedTuple` metadata, and `extras` holds unstructured metadata. Returns a
`ChromScanSeries` containing the scans and metadata. Throws `ArgumentError` if `scans` is
empty or if the element type of `scans` is not concrete.

See also [`AbstractScanSeries`](@ref), [`AbstractChromScanSeries`](@ref), [`scans`](@ref), 
[`instrument`](@ref), [`acquisition`](@ref), [`user`](@ref), [`sample`](@ref), 
[`extras`](@ref).

# Examples
```jldoctest
julia> cs1 = ChromScan(1.0u"minute", 5.0u"pA");

julia> cs2 = ChromScan(2.0u"minute", 8.0u"pA");

julia> css = ChromScanSeries([cs1, cs2]; instrument=(detector="DAD",));

julia> scans(css) == [ChromScan(1.0u"minute", 5.0u"pA"), ChromScan(2.0u"minute", 8.0u"pA")]
true

julia> css.instrument
(detector = "DAD",)

julia> length(css)
2
```
"""
function ChromScanSeries(
    scans::T1;
    instrument::T2=NamedTuple(),
    acquisition::T3=NamedTuple(),
    user::T4=NamedTuple(),
    sample::T5=NamedTuple(),
    extras::T6=Dict{String, Any}()
    ) where {
    T1<:AbstractVector{<:AbstractChromScan},
    T2<:NamedTuple,
    T3<:NamedTuple,
    T4<:NamedTuple,
    T5<:NamedTuple,
    T6<:AbstractDict{<:AbstractString, <:Any}}

    length(scans) > 0 || throw(ArgumentError("no scan(s) provided"))
    S = eltype(scans)
    isconcretetype(S) || throw(ArgumentError(
        "scans must have a concrete element type, got $S"))

    R = typeof(first(scans).retention_unit)
    I = typeof(first(scans).intensity_unit)

    converted_extras = Dict{String, Any}(String(k) => v for (k, v) in extras)

    ChromScanSeries{S, R, I, T1, T2, T3, T4, T5, Dict{String, Any}}(
        scans, instrument, acquisition, user, sample, converted_extras)
end

Base.:(==)(a::ChromScanSeries, b::ChromScanSeries) =
    scancount(a) == scancount(b) &&
    all(x == y for (x, y) in zip(scans(a), scans(b))) &&
    instrument(a) == instrument(b) &&
    acquisition(a) == acquisition(b) &&
    user(a) == user(b) &&
    sample(a) == sample(b) &&
    extras(a) == extras(b)

struct MassScanSeries{
    S<:AbstractMassScan,
    R,
    M,
    I,
    T1<:AbstractVector{S},
    T2<:NamedTuple,
    T3<:NamedTuple,
    T4<:NamedTuple,
    T5<:NamedTuple,
    T6<:Dict{String, Any}
    } <: AbstractMassScanSeries{S, R, M, I}

    scans::T1
    instrument::T2
    acquisition::T3
    user::T4
    sample::T5
    extras::T6

    function MassScanSeries{S, R, M, I, T1, T2, T3, T4, T5, T6}(
        scans::T1,
        instrument::T2,
        acquisition::T3,
        user::T4,
        sample::T5,
        extras::T6
    ) where {
        S<:AbstractMassScan,
        R,
        M,
        I,
        T1<:AbstractVector{S},
        T2<:NamedTuple,
        T3<:NamedTuple,
        T4<:NamedTuple,
        T5<:NamedTuple,
        T6<:Dict{String, Any}
    }

        length(scans) > 0 || throw(ArgumentError("no scan(s) provided"))

        E = eltype(scans)
        isconcretetype(E) || throw(ArgumentError(
            "scans must have a concrete element type, got $E"))

        new{S, R, M, I, T1, T2, T3, T4, T5, T6}(
            scans, instrument, acquisition, user, sample, extras)
    end
end

"""
    MassScanSeries(
        scans::Vector{<:AbstractMassScan};
        instrument::NamedTuple=NamedTuple(),
        acquisition::NamedTuple=NamedTuple(),
        user::NamedTuple=NamedTuple(),
        sample::NamedTuple=NamedTuple(),
        extras::Dict{<:AbstractString, <:Any}=Dict()
    ) -> MassScanSeries

Constructs a `MassScanSeries`, representing an ordered collection of scans (subtypes of 
`AbstractMassScan`) with associated metadata.

`scans` must be a non-empty vector of mass scan objects, all of the same concrete subtype
of `AbstractMassScan`. `instrument`, `acquisition`, `user`, and `sample` accept optional
`NamedTuple` metadata, and `extras` holds unstructured metadata. Returns a `MassScanSeries`
containing the scans and metadata. Throws `ArgumentError` if `scans` is empty or if the
element type of `scans` is not concrete.

See also [`AbstractScanSeries`](@ref), [`AbstractMassScanSeries`](@ref), [`scans`](@ref), 
[`instrument`](@ref), [`acquisition`](@ref), [`user`](@ref), [`sample`](@ref), 
[`extras`](@ref).

# Examples
```jldoctest
julia> ms1 = MassScan(1.0u"s", [100.0, 150.0], [10.0, 20.0]);

julia> ms2 = MassScan(2.0u"s", [100.0, 150.0], [12.0, 22.0]);

julia> mss = MassScanSeries([ms1, ms2]; acquisition=(mode="FullScan",));

julia> mss.acquisition.mode
"FullScan"

julia> length(mss)
2
```
"""
function MassScanSeries(
    scans::T1;
    instrument::T2=NamedTuple(),
    acquisition::T3=NamedTuple(),
    user::T4=NamedTuple(),
    sample::T5=NamedTuple(),
    extras::T6=Dict{String, Any}()
    ) where {
    T1<:AbstractVector{<:AbstractMassScan},
    T2<:NamedTuple,
    T3<:NamedTuple,
    T4<:NamedTuple,
    T5<:NamedTuple,
    T6<:AbstractDict{<:AbstractString, <:Any}}

    length(scans) > 0 || throw(ArgumentError("no scan(s) provided"))
    S = eltype(scans)
    isconcretetype(S) || throw(ArgumentError(
        "scans must have a concrete element type, got $S"))

    R = typeof(first(scans).retention_unit)
    M = typeof(first(scans).mz_unit)
    I = typeof(first(scans).intensity_unit)

    converted_extras = Dict{String, Any}(String(k) => v for (k, v) in extras)

    MassScanSeries{S, R, M, I, T1, T2, T3, T4, T5, Dict{String, Any}}(
        scans, instrument, acquisition, user, sample, converted_extras)
end

Base.:(==)(a::MassScanSeries, b::MassScanSeries) =
    scancount(a) == scancount(b) &&
    all(x == y for (x, y) in zip(scans(a), scans(b))) &&
    instrument(a) == instrument(b) &&
    acquisition(a) == acquisition(b) &&
    user(a) == user(b) &&
    sample(a) == sample(b) &&
    extras(a) == extras(b)

############################################################################################

# Shared utility functions for all scan series display methods
module ScanSeriesDisplay

# Helper functions for value formatting and wrapping
function format_value_with_wrap(value_str::String, max_width::Int)
    if length(value_str) <= max_width
        return value_str
    else
        # Split long values at reasonable breakpoints
        return join([value_str[i:min(i+max_width-1, end)] for i in 1:max_width:length(value_str)], "\n")
    end
end

function print_wrapped_value(io::IO, prefix::String, value_str::String, continuation_indent::String)
    lines = split(value_str, '\n')
    if length(lines) == 1
        # No actual wrapping needed, just print normally
        println(io, "$prefix $value_str")
    else
        println(io, "$prefix $(lines[1])")
        for line in lines[2:end]
            println(io, "$continuation_indent$line")
        end
    end
end

# Check if a value is "complex" (needs nested display)
function is_complex_value(value)
    return value isa NamedTuple && length(value) > 0 ||
           value isa Dict && length(value) > 0 ||
           value isa AbstractVector && length(value) > 3 ||  # Show arrays with >3 elements as nested
           value isa AbstractMatrix
end

# Print a complex value with proper tree structure
function print_complex_value(io::IO, value, base_prefix::String, continuation_prefix::String, is_last::Bool)
    if value isa NamedTuple
        pairs_list = collect(pairs(value))
        for (i, (key, subvalue)) in enumerate(pairs_list)
            is_last_item = i == length(pairs_list)
            item_prefix = is_last_item ? "$(base_prefix)└─" : "$(base_prefix)├─"
            item_continuation = is_last_item ? "$(base_prefix) " : "$(base_prefix)│"
            
            if is_complex_value(subvalue)
                println(io, "$item_prefix $key:")
                print_complex_value(io, subvalue, "$(item_continuation)  ", "$(item_continuation)  ", true)
            else
                value_str = format_value_with_wrap(string(subvalue), 50)
                if '\n' in value_str
                    print_wrapped_value(io, "$item_prefix $key =", value_str, "$(item_continuation)     ")
                else
                    println(io, "$item_prefix $key = $value_str")
                end
            end
        end
    elseif value isa Dict
        pairs_list = collect(value)
        for (i, (key, subvalue)) in enumerate(pairs_list)
            is_last_item = i == length(pairs_list)
            item_prefix = is_last_item ? "$(base_prefix)└─" : "$(base_prefix)├─"
            item_continuation = is_last_item ? "$(base_prefix) " : "$(base_prefix)│"
            
            if is_complex_value(subvalue)
                println(io, "$item_prefix $key:")
                print_complex_value(io, subvalue, "$(item_continuation)  ", "$(item_continuation)  ", true)
            else
                value_str = format_value_with_wrap(string(subvalue), 50)
                if '\n' in value_str
                    print_wrapped_value(io, "$item_prefix $key =", value_str, "$(item_continuation)     ")
                else
                    println(io, "$item_prefix $key = $value_str")
                end
            end
        end
    elseif value isa AbstractVector
        if length(value) <= 3
            # Short arrays: display inline
            println(io, "$(base_prefix)└─ [$(join(value, ", "))]")
        else
            # Long arrays: show first few elements + length
            preview = join(value[1:min(3, end)], ", ")
            println(io, "$(base_prefix)├─ [$preview, ...]")
            println(io, "$(base_prefix)└─ length: $(length(value))")
        end
    elseif value isa AbstractMatrix
        println(io, "$(base_prefix)├─ $(size(value, 1))×$(size(value, 2)) $(typeof(value).name.name)")
        if size(value, 1) <= 2 && size(value, 2) <= 3
            # Small matrices: show content
            nrows = size(value, 1)
            for (i, row_idx) in enumerate(axes(value, 1))
                is_last_row = i == nrows
                row_prefix = is_last_row ? "$(base_prefix)└─" : "$(base_prefix)├─"
                println(io, "$row_prefix row $i: [$(join(value[row_idx, :], ", "))]")
            end
        else
            println(io, "$(base_prefix)└─ (content omitted)")
        end
    else
        # Fallback for other complex types
        println(io, "$(base_prefix)└─ $(typeof(value)): $(string(value))")
    end
end

# Print annotations section (shared between ChromScanSeries and MassScanSeries)
function print_annotations(io::IO, metadata_sections, additional_metadata)
    non_empty_sections = filter(section -> !isempty(section[2]), metadata_sections)
    has_metadata = !isempty(non_empty_sections) || !isempty(additional_metadata)
    
    if !has_metadata
        return false
    end
    
    println(io, "└─ Annotations:")
    
    # Print non-empty metadata sections
    for (i, (name, data)) in enumerate(non_empty_sections)
        is_last_section = (i == length(non_empty_sections)) && isempty(additional_metadata)
        section_prefix = is_last_section ? "   └─" : "   ├─"
        section_continuation = is_last_section ? "   " : "   │"
        
        # Format the metadata content
        if length(data) == 1
            key, value = first(pairs(data))
            if is_complex_value(value)
                println(io, "$section_prefix $name:")
                println(io, "$section_continuation  └─ $key:")
                print_complex_value(io, value, "$(section_continuation)     ", "$(section_continuation)     ", true)
            else
                value_str = format_value_with_wrap(string(value), 60)
                if '\n' in value_str
                    println(io, "$section_prefix $name:")
                    print_wrapped_value(io, "$section_continuation  └─ $key =", value_str, "$section_continuation     ")
                else
                    println(io, "$section_prefix $name: $key = $value_str")
                end
            end
        else
            println(io, "$section_prefix $name:")
            data_pairs = collect(pairs(data))
            for (j, (key, value)) in enumerate(data_pairs)
                is_last_item = j == length(data_pairs)
                item_prefix = is_last_item ? "$section_continuation  └─" : "$section_continuation  ├─"
                item_continuation = is_last_item ? "$section_continuation   " : "$section_continuation  │"
                
                if is_complex_value(value)
                    println(io, "$item_prefix $key:")
                    print_complex_value(io, value, "$(item_continuation)  ", "$(item_continuation)  ", true)
                else
                    value_str = format_value_with_wrap(string(value), 50)
                    if '\n' in value_str
                        print_wrapped_value(io, "$item_prefix $key =", value_str, "$(item_continuation)     ")
                    else
                        println(io, "$item_prefix $key = $value_str")
                    end
                end
            end
        end
    end
    
    # Print additional metadata if present
    if !isempty(additional_metadata)
        println(io, "   └─ Additional metadata:")
        metadata_pairs = collect(additional_metadata)
        for (i, (key, value)) in enumerate(metadata_pairs)
            is_last = i == length(metadata_pairs)
            prefix = is_last ? "      └─" : "      ├─"
            continuation = is_last ? "       " : "      │"
            
            if is_complex_value(value)
                println(io, "$prefix $key:")
                print_complex_value(io, value, "$(continuation)  ", "$(continuation)  ", true)
            else
                value_str = format_value_with_wrap(string(value), 50)
                if '\n' in value_str
                    print_wrapped_value(io, "$prefix $key =", value_str, "$(continuation)     ")
                else
                    if is_last
                        print(io, "$prefix $key = $value_str")
                    else
                        println(io, "$prefix $key = $value_str")
                    end
                end
            end
        end
    end
    
    return true
end

end # module ScanSeriesDisplay

# ChromScanSeries show method
function Base.show(io::IO, css::ChromScanSeries)
    # Get basic information
    n_scans = length(css.scans)
    scan_type = eltype(css.scans)
    
    # Extract retention information
    retentions = [rawretention(scan) for scan in css.scans]
    retention_unit = retentionunit(first(css.scans))
    retention_range = (minimum(retentions), maximum(retentions))
    
    # Extract intensity information
    intensities = [rawintensity(scan) for scan in css.scans]
    intensity_unit = intensityunit(first(css.scans))
    intensity_range = (minimum(intensities), maximum(intensities))
    
    # Format units for display
    ret_unit_str = retention_unit === nothing ? "unitless" : string(retention_unit)
    int_unit_str = intensity_unit === nothing ? "unitless" : string(intensity_unit)
    
    # Check what metadata sections have content
    metadata_sections = [
        ("Instrument", css.instrument),
        ("Acquisition", css.acquisition),
        ("User", css.user),
        ("Sample", css.sample)
    ]
    
    non_empty_sections = filter(section -> !isempty(section[2]), metadata_sections)
    has_metadata = !isempty(non_empty_sections) || !isempty(css.extras)
    
    # Print header
    scan_word = n_scans == 1 ? "scan" : "scans"
    println(io, "ChromScanSeries with $n_scans $scan_word")
    
    # Handle potentially long scan type
    scan_type_str = string(scan_type)
    if length(scan_type_str) > 80
        println(io, "├─ Scan type:")
        ScanSeriesDisplay.print_wrapped_value(io, "│  └─", scan_type_str, "     ")
    else
        println(io, "├─ Scan type: $scan_type_str")
    end
    
    # Print retention information
    println(io, "├─ Retention:")
    println(io, "│  ├─ Range: $(retention_range[1]) to $(retention_range[2]) ($ret_unit_str)")
    println(io, "│  └─ Type: $(eltype(retentions))")
    
    # Print intensity information (adjust connector based on whether metadata follows)
    intensity_prefix = has_metadata ? "├─" : "└─"
    println(io, "$intensity_prefix Intensity:")
    intensity_sub_prefix = has_metadata ? "│" : " "
    println(io, "$intensity_sub_prefix  ├─ Range: $(intensity_range[1]) to $(intensity_range[2]) ($int_unit_str)")
    println(io, "$intensity_sub_prefix  └─ Type: $(eltype(intensities))")
    
    # Print annotations using shared function
    ScanSeriesDisplay.print_annotations(io, metadata_sections, css.extras)
end

# MassScanSeries show method
function Base.show(io::IO, mss::MassScanSeries)
    # Get basic information
    n_scans = length(mss.scans)
    scan_type = eltype(mss.scans)
    
    # Extract retention information
    retentions = [rawretention(scan) for scan in mss.scans]
    retention_unit = retentionunit(first(mss.scans))
    retention_range = (minimum(retentions), maximum(retentions))
    
    # Extract m/z information
    all_mz = vcat([rawmzvalues(scan) for scan in mss.scans]...)
    mz_unit = mzunit(first(mss.scans))
    mz_range = (minimum(all_mz), maximum(all_mz))
    
    # Extract intensity information
    all_intensities = vcat([rawintensities(scan) for scan in mss.scans]...)
    intensity_unit = intensityunit(first(mss.scans))
    intensity_range = (minimum(all_intensities), maximum(all_intensities))
    
    # Extract MS levels
    levels = [level(scan) for scan in mss.scans]
    unique_levels = sort(unique(levels))
    
    # Format units for display
    ret_unit_str = retention_unit === nothing ? "unitless" : string(retention_unit)
    mz_unit_str = mz_unit === nothing ? "unitless" : string(mz_unit)
    int_unit_str = intensity_unit === nothing ? "unitless" : string(intensity_unit)
    
    # Check what metadata sections have content
    metadata_sections = [
        ("Instrument", mss.instrument),
        ("Acquisition", mss.acquisition),
        ("User", mss.user),
        ("Sample", mss.sample)
    ]
    
    non_empty_sections = filter(section -> !isempty(section[2]), metadata_sections)
    has_metadata = !isempty(non_empty_sections) || !isempty(mss.extras)
    
    # Print header
    scan_word = n_scans == 1 ? "scan" : "scans"
    println(io, "MassScanSeries with $n_scans $scan_word")
    
    # Handle potentially long scan type
    scan_type_str = string(scan_type)
    if length(scan_type_str) > 80
        println(io, "├─ Scan type:")
        ScanSeriesDisplay.print_wrapped_value(io, "│  └─", scan_type_str, "     ")
    else
        println(io, "├─ Scan type: $scan_type_str")
    end
    
    # Print retention information
    println(io, "├─ Retention:")
    println(io, "│  ├─ Range: $(retention_range[1]) to $(retention_range[2]) ($ret_unit_str)")
    println(io, "│  └─ Type: $(eltype(retentions))")
    
    # Print m/z information
    println(io, "├─ M/Z values:")
    println(io, "│  ├─ Range: $(mz_range[1]) to $(mz_range[2]) ($mz_unit_str)")
    println(io, "│  ├─ Total data points: $(length(all_mz))")
    println(io, "│  └─ Type: $(eltype(all_mz))")
    
    # Print intensity information (adjust connector based on whether metadata follows)
    intensity_prefix = has_metadata ? "├─" : "└─"
    println(io, "$intensity_prefix Intensity:")
    intensity_sub_prefix = has_metadata ? "│" : " "
    println(io, "$intensity_sub_prefix  ├─ Range: $(intensity_range[1]) to $(intensity_range[2]) ($int_unit_str)")
    println(io, "$intensity_sub_prefix  ├─ Type: $(eltype(all_intensities))")
    
    # Print MS levels
    if length(unique_levels) == 1
        println(io, "$intensity_sub_prefix  └─ MS Level: $(unique_levels[1])")
    else
        println(io, "$intensity_sub_prefix  └─ MS Levels: $(join(unique_levels, ", "))")
    end
    
    # Print annotations using shared function
    ScanSeriesDisplay.print_annotations(io, metadata_sections, mss.extras)
end

# Shared MIME show method for both types
function Base.show(io::IO, ::MIME"text/plain", series::Union{ChromScanSeries, MassScanSeries})
    show(io, series)
end
