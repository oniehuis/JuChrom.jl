# ── indextrim ─────────────────────────────────────────────────────────────────────────────

"""
    indextrim(
        series::AbstractScanSeries;
        start=firstindex(scans(series)), 
        stop=lastindex(scans(series))
    ) -> typeof(series)

Return a new scan series of the same concrete type as `series` containing only the
scans whose indices are within the inclusive range `[start, stop]`. The function
returns a copy of `series` with only the requested scan slice. It throws `BoundsError`
if the indices are out of bounds, and `ArgumentError` if `start > stop` or the
resulting scan list is empty. 
    
See also
[`AbstractScanSeries`](@ref JuChrom.AbstractScanSeries), 
[`indextrim!`](@ref), 
[`retentiontrim`](@ref),
[`retentiontrim!`](@ref), 
[`scancount`](@ref)
[`scans`](@ref).

# Examples
```jldoctest
julia> scan1 = ChromScan(1.0u"minute", 5.0)
       scan2 = ChromScan(2.0u"minute", 6.0)
       scan3 = ChromScan(3.0u"minute", 7.0)
       series = ChromScanSeries([scan1, scan2, scan3]);

julia> trimmed = indextrim(series; start=2, stop=2);

julia> scancount(trimmed)
1

julia> retentions(trimmed) == [2.0]u"minute"
true
```
"""
function indextrim(
    series::AbstractScanSeries; 
    start::T1=firstindex(scans(series)), 
    stop::T2=lastindex(scans(series))
    ) where {T1<:Integer, T2<:Integer}

    # Call the non-mutating helper to extract scans by index range
    filtered_scans = _indextrim(series, start, stop)

    # Return a new copy of the series with only the filtered scans included
    copy_with(series, (; scans=filtered_scans))
end


function _indextrim(series::AbstractScanSeries, start::Integer, stop::Integer)
    # Ensure indices are within valid bounds
    checkbounds(scans(series), start)
    checkbounds(scans(series), stop)

    # Validate index order
    start ≤ stop || throw(ArgumentError("`start` must be less than or equal to `stop`."))

    # Slice the scans in the inclusive index range
    filtered_scans = scans(series)[start:stop]

    # Ensure result is not empty
    isempty(filtered_scans) && throw(ArgumentError("No scans left after trimming."))

    filtered_scans
end

# ── indextrim! ────────────────────────────────────────────────────────────────────────────

"""
    indextrim!(
        series::AbstractScanSeries;
        start=firstindex(scans(series)), 
        stop=lastindex(scans(series))
    ) -> typeof(series)

Mutate `series` by retaining only scans whose indices fall within the inclusive range
`[start, stop]`. The underlying scan vector is trimmed in place, so the original
series is modified and no new series is allocated. It throws `BoundsError` if the
indices are out of bounds, and `ArgumentError` if `start > stop` or the resulting scan
list is empty. 

See also 
[`AbstractScanSeries`](@ref), 
[`indextrim`](@ref),
[`retentiontrim`](@ref), 
[`retentiontrim!`](@ref), 
[`scancount`](@ref)
[`scans`](@ref).

# Examples
```jldoctest
julia> scan1 = ChromScan(1.0u"minute", 5.0)
       scan2 = ChromScan(2.0u"minute", 6.0)
       scan3 = ChromScan(3.0u"minute", 7.0)
       series = ChromScanSeries([scan1, scan2, scan3]);

julia> indextrim!(series; start=2, stop=2);

julia> scancount(series)
1

julia> retentions(series) == [2.0]u"minute"
true
```
"""
function indextrim!(
    series::AbstractScanSeries; 
    start::Integer=firstindex(scans(series)), 
    stop::Integer=lastindex(scans(series))
)
    scanvec = scans(series)

    start ≤ stop || throw(ArgumentError("`start` must be less than or equal to `stop`."))
    checkbounds(scanvec, start:stop)

    # Remove elements after the stop index
    if stop < length(scanvec)
        deleteat!(scanvec, (stop + 1):length(scanvec))
    end

    # Remove elements before the start index
    if start > 1
        deleteat!(scanvec, 1:(start - 1))
    end

    # Ensure result is not empty
    isempty(scanvec) && throw(ArgumentError("No scans left after trimming."))

    series
end

# ── levelscans ────────────────────────────────────────────────────────────────────────────

"""
    levelscans(series::AbstractMassScanSeries, target_level::Integer=1) -> typeof(series)

Return a new `AbstractMassScanSeries` containing only scans with the specified m/z 
`target_level` (e.g., MS¹ only). The returned series preserves metadata and leaves the 
original unchanged. It throws `ArgumentError` if no scans exist at the requested level. 
    
See also
[`AbstractMassScanSeries`](@ref),
[`levels`](@ref levels(::AbstractMassScanSeries)),
[`scancount`](@ref), 
[`scans`](@ref).

# Examples
```jldoctest
julia> scan1 = MassScan(1.0u"minute", [100.0], [10.0]; level=1)
       scan2 = MassScan(2.0u"minute", [150.0], [15.0]; level=2)
       series = MassScanSeries([scan1, scan2]);

julia> filtered = levelscans(series, 2);

julia> scancount(filtered)
1

julia> level(scan(filtered, 1))
2
```
"""
function levelscans(series::AbstractMassScanSeries, target_level::Integer=1)
    # Filter to scans with the specified MS level
    level_scans = filter(scan -> level(scan) == target_level, scans(series))

    # Ensure at least one scan remains
    isempty(level_scans) && throw_noscans_error(target_level)

    # Return new series with filtered scans
    copy_with(series, (; scans=level_scans))
end

@noinline function throw_noscans_error(level::Integer)
    throw(ArgumentError("No scans found at level $level."))
end

# ── retentiontrim ─────────────────────────────────────────────────────────────────────────

"""
    retentiontrim(
        series::AbstractScanSeries; 
        start=first(retentions(series)), 
        stop=last(retentions(series))
    ) -> typeof(series)

Return a new scan series of the same concrete type as `series` containing only those
scans whose retention values fall within the inclusive range `[start, stop]`. The
method supports series with unitful or unitless retention values and defaults to the
range covered by the series. It throws `ArgumentError` if `start > stop` or no scans
remain after trimming. 

See also
[`AbstractScanSeries`](@ref), 
[`indextrim`](@ref),
[`indextrim!`](@ref), 
[`retentiontrim!`](@ref), 
[`scancount`](@ref),
[`scans`](@ref).

# Examples
```jldoctest
julia> scan1 = ChromScan(1.0u"minute", 5.0)
       scan2 = ChromScan(2.0u"minute", 6.0)
       scan3 = ChromScan(3.0u"minute", 7.0)
       series = ChromScanSeries([scan1, scan2, scan3]);

julia> trimmed = retentiontrim(series; start=1.5u"minute", stop=2.5u"minute");

julia> scancount(trimmed)
1

julia> retentions(trimmed) == [2.0]u"minute"
true
```
"""
function retentiontrim(
    series::AbstractScanSeries{<:AbstractScan, <:Unitful.Units}; 
    start::T1=first(retentions(series)), 
    stop::T2=last(retentions(series))
    ) where {T1<:Unitful.AbstractQuantity, T2<:Unitful.AbstractQuantity}

    # Call the non-mutating helper to get the filtered scans
    filtered_scans = _retentiontrim(series, start, stop)

    # Return a new copy of the series with only the filtered scans included
    copy_with(series, (; scans=filtered_scans))
end

function retentiontrim(
    series::AbstractScanSeries{<:AbstractScan, <:Nothing}; 
    start::T1=first(retentions(series)), 
    stop::T2=last(retentions(series))
    ) where {T1<:Real, T2<:Real}

    # Call the non-mutating helper to get the filtered scans
    filtered_scans = _retentiontrim(series, start, stop)

    # Return a new copy of the series with only the filtered scans included
    copy_with(series, (; scans=filtered_scans))
end

function _retentiontrim(series::AbstractScanSeries, start, stop)

    # Validate that the time range is in correct order
    start ≤ stop || throw(ArgumentError("`start` must be less than or equal to `stop`."))

    # Filter scans to include only those within the [start, stop] range
    filtered_scans = filter(scan -> start ≤ retention(scan) ≤ stop, scans(series))

    # Ensure that trimming did not remove all scans
    isempty(filtered_scans) && throw(ArgumentError("No scans left after trimming."))

    # Return the filtered collection of scans
    filtered_scans
end

# ── retentiontrim! ────────────────────────────────────────────────────────────────────────

"""
    retentiontrim!(
        series::AbstractScanSeries; 
        start=first(retentions(series)), 
        stop=last(retentions(series))
    ) -> typeof(series)

Mutate the `series` in place by removing scans whose retention values fall outside the
inclusive range `[start, stop]`. The method supports series with unitful or unitless
retention values and trims the underlying scan vector. It throws `ArgumentError` if
`start > stop` or no scans remain after filtering. 

See also 
[`AbstractScanSeries`](@ref),
[`indextrim`](@ref), 
[`indextrim!`](@ref), 
[`retentiontrim`](@ref), 
[`scancount`](@ref),
[`scans`](@ref).

# Examples
```jldoctest
julia> scan1 = ChromScan(1.0u"minute", 5.0)
       scan2 = ChromScan(2.0u"minute", 6.0)
       scan3 = ChromScan(3.0u"minute", 7.0)
       series = ChromScanSeries([scan1, scan2, scan3]);

julia> retentiontrim!(series; start=1.5u"minute", stop=2.5u"minute");

julia> scancount(series)
1

julia> retention(scan(series, 1)) == 2.0u"minute"
true
```
"""
function retentiontrim!(
    series::AbstractScanSeries{<:AbstractScan, <:Unitful.Units, <:Any}; 
    start::T1=first(retentions(series)), 
    stop::T2=last(retentions(series))
    ) where {T1<:Unitful.AbstractQuantity, T2<:Unitful.AbstractQuantity}

    start ≤ stop || throw(ArgumentError("`start` must be less than or equal to `stop`."))

    scanvec = scans(series)
    filter!(scan -> start ≤ retention(scan) ≤ stop, scanvec)

    # Ensure result is not empty
    isempty(scanvec) && throw(ArgumentError("No scans left after trimming."))

    series
end

function retentiontrim!(
    series::AbstractScanSeries{<:AbstractScan, <:Nothing, <:Any}; 
    start::T1=first(retentions(series)), 
    stop::T2=last(retentions(series))
    ) where {T1<:Real, T2<:Real}

    start ≤ stop || throw(ArgumentError("`start` must be less than or equal to `stop`."))

    scanvec = scans(series)
    filter!(scan -> start ≤ retention(scan) ≤ stop, scanvec)

    # Ensure result is not empty
    isempty(scanvec) && throw(ArgumentError("No scans left after trimming."))

    series
end
