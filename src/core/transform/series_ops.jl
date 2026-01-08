# ── indextrim ─────────────────────────────────────────────────────────────────────────────

"""
    indextrim(series::AbstractScanSeries; start=firstindex(scans(series)), 
              stop=lastindex(scans(series))) -> AbstractScanSeries

Return a new `AbstractScanSeries` containing only the scans whose indices are within the 
inclusive range `[start, stop]`.

# Arguments
- `series::AbstractScanSeries`: The scan series to trim.
- `start::Integer`: Starting index of the scan range to retain. Defaults to the first index.
- `stop::Integer`: Ending index of the scan range to retain. Defaults to the last index.

# Returns
A new copy of `series` with only the scans in the specified index range.

# Throws
- `BoundsError` if `start` or `stop` is outside the bounds of `scans(series)`.
- `ArgumentError` if `start > stop`.
- `ArgumentError` if the resulting scan list is empty.

See also [`AbstractScanSeries`](@ref), [`indextrim!`](@ref), [`retentiontrim`](@ref),
[`retentiontrim!`](@ref), [`scans`](@ref), [`scancount`](@ref).

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
````
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
    indextrim!(series::AbstractScanSeries; start=firstindex(scans(series)), 
               stop=lastindex(scans(series))) -> AbstractScanSeries

Mutate `series` by retaining only scans whose indices fall within the inclusive range
`[start, stop]`.

This operation trims the underlying scan vector in-place. It is non-allocating and 
modifies the original data structure.

# Arguments
- `series::AbstractScanSeries`: The scan series to modify.
- `start::Integer`: Starting index of the scan range to retain. Defaults to the first index.
- `stop::Integer`: Ending index of the scan range to retain. Defaults to the last index.

# Returns
- The mutated `series` with only scans from index `start` to `stop` remaining.

# Throws
- `BoundsError` if `start` or `stop` is outside the bounds of `scans(series)`.
- `ArgumentError` if `start > stop`.
- `ArgumentError` if the resulting scan list is empty.

See also [`AbstractScanSeries`](@ref), [`indextrim`](@ref), [`retentiontrim`](@ref),
[`retentiontrim!`](@ref), [`scans`](@ref), [`scancount`](@ref).

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
    levelscans(series::AbstractMassScanSeries, target_level::Integer=1) 
        -> AbstractMassScanSeries

Return a new `AbstractMassScanSeries` containing only scans with the specified 
`target_level`.

# Arguments
- `series::AbstractMassScanSeries`: A series of mass spectrometry scans.
- `target_level::Integer=1`: The MS level to filter for (e.g. 1 for MS1, 2 for MS2).

# Returns
A new mass scan series containing only the scans at the specified level. All metadata is 
preserved.

# Throws
- `ArgumentError` if no scans exist at the specified `target_level`.

See also [`AbstractMassScanSeries`](@ref), [`scans`](@ref), [`scan`](@ref), [`level`](@ref).

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
    retentiontrim(series::AbstractScanSeries; start=first(retentions(series)), 
                  stop=last(retentions(series))) -> AbstractScanSeries

Return a new `AbstractScanSeries` containing only those scans whose retention values fall 
within the inclusive range `[start, stop]`.

This method supports series with unitful or unitless retention values. It automatically 
uses the correct default range based on the retention values present in the series.

# Arguments
- `series::AbstractScanSeries`: The scan series to trim.
- `start`: The lower bound of the retention range. Defaults to the first retention value in 
  `series`.
- `stop`: The upper bound of the retention range. Defaults to the last retention value in 
  `series`.

# Returns
A new `AbstractScanSeries` with scans filtered by the specified retention range.

# Throws
- `ArgumentError` if `start > stop`.
- `ArgumentError` if no scans remain after trimming.

See also [`AbstractScanSeries`](@ref), [`retentiontrim!`](@ref), [`indextrim`](@ref), 
[`indextrim!`](@ref), [`scans`](@ref), [`scancount`](@ref).

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
    ) where {T1<:Unitful.Quantity, T2<:Unitful.Quantity}

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

function retentiontrim(
    msmatrix::MassScanMatrix{<:Any, <:Unitful.Units}; 
    start::T1=first(retentions(msmatrix)), 
    stop::T2=last(retentions(msmatrix))
    ) where {T1<:Unitful.Quantity, T2<:Unitful.Quantity}

    _retentiontrim(msmatrix::MassScanMatrix, start, stop)
end

function retentiontrim(
    msmatrix::MassScanMatrix{<:Any, <:Nothing}; 
    start::T1=first(retentions(msmatrix)), 
    stop::T2=last(retentions(msmatrix))
    ) where {T1<:Real, T2<:Real}

    _retentiontrim(msmatrix::MassScanMatrix, start, stop)
end

function _retentiontrim(msmatrix::MassScanMatrix, start, stop)

        # Validate that the time range is in correct order
    start ≤ stop || throw(ArgumentError("`start` must be less than or equal to `stop`."))

    # Filter scans to include only those within the [start, stop] range
    i_start = searchsortedfirst(retentions(msmatrix), start)
    i_stop = searchsortedlast(retentions(msmatrix), stop)
    indices = i_start:i_stop

    # Ensure that trimming did not remove all scans
    isempty(indices) && throw(ArgumentError("No scans left after trimming."))

    # Create and return a new MassScanMatrix with the filtered scans 
    MassScanMatrix(
        retentions(msmatrix)[indices],
        retentionunit(msmatrix),
        mzvalues(msmatrix),
        mzunit(msmatrix),
        rawintensities(msmatrix)[indices, :],
        intensityunit(msmatrix),
        level=level(msmatrix),
        instrument=deepcopy(instrument(msmatrix)),
        acquisition=deepcopy(acquisition(msmatrix)),
        user=deepcopy(user(msmatrix)),
        sample=deepcopy(sample(msmatrix)),
        extras=deepcopy(extras(msmatrix))
    )
end

# ── retentiontrim! ────────────────────────────────────────────────────────────────────────

"""
    retentiontrim!(series::AbstractScanSeries; start=first(retentions(series)), 
                   stop=last(retentions(series))) -> AbstractScanSeries

Mutate the `series` in-place by removing scans whose retention values fall outside the 
inclusive range `[start, stop]`.

This method supports series with either unitful or unitless retention values.

# Arguments
- `series::AbstractScanSeries`: The scan series to modify.
- `start`: Lower bound for retention filtering. Defaults to the first retention in the 
  series.
- `stop`: Upper bound for retention filtering. Defaults to the last retention in the 
  series.

# Returns
- The same `series` object, with its scans filtered in-place.

# Throws
- `ArgumentError` if `start > stop`.
- `ArgumentError` if no scans remain after filtering.

See also [`AbstractScanSeries`](@ref), [`retentiontrim`](@ref), [`indextrim`](@ref), 
[`indextrim!`](@ref), [`scans`](@ref), [`scancount`](@ref).

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
    ) where {T1<:Unitful.Quantity, T2<:Unitful.Quantity}

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
