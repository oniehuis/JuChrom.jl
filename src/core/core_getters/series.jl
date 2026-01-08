# ── acquisition ───────────────────────────────────────────────────────────────────────────

"""
    acquisition(series::AbstractScanSeries) -> NamedTuple

Returns the acquisition metadata associated with the scan series.

The acquisition metadata typically includes details such as scan mode, method parameters, 
instrument configuration, or other context relevant to how the data was acquired.

# Arguments
- `series`: A concrete subtype of `AbstractScanSeries`.

# Returns
A `NamedTuple` containing acquisition-related metadata.

See also [`AbstractScanSeries`](@ref), [`instrument`](@ref), [`user`](@ref), 
[`sample`](@ref).

# Examples
```jldoctest
julia> s = [ChromScan(0.0u"s", 12.3), ChromScan(1.0u"s", 45.6)];

julia> series = ChromScanSeries(s; acquisition=(mode="gradient", method="HPLC-QC-42"));

julia> acquisition(series)
(mode = "gradient", method = "HPLC-QC-42")

julia> acquisition(series).method
"HPLC-QC-42"
```
"""
acquisition(series::AbstractScanSeries) = series.acquisition

# ── extras ────────────────────────────────────────────────────────────────────────────────

"""
    extras(series::AbstractScanSeries) -> Dict{String, Any}

Returns the unstructured extras associated with the scan series.

This dictionary may contain arbitrary key-value pairs that are not captured by the 
structured fields such as `instrument`, `acquisition`, `user`, or `sample`.

# Arguments
- `series`: A concrete subtype of `AbstractScanSeries`.

# Returns
A `Dict{String, Any}` containing unstructured metadata.

See also [`instrument`](@ref), [`acquisition`](@ref), [`user`](@ref), [`sample`](@ref).

# Examples
```jldoctest
julia> s = [ChromScan(1.0u"s", 10.0u"V"), ChromScan(2.0u"s", 200.0u"V")];

julia> series = ChromScanSeries(s, extras=Dict("run"=>42, "comment"=>"Baseline drift"));

julia> extras(series)
Dict{String, Any} with 2 entries:
  "run"     => 42
  "comment" => "Baseline drift"

julia> extras(series)["qc_passed"] = true;

julia> extras(series)
Dict{String, Any} with 3 entries:
  "qc_passed" => true
  "run"       => 42
  "comment"   => "Baseline drift"
```
"""
extras(series::AbstractScanSeries) = series.extras

# ── intensities ───────────────────────────────────────────────────────────────────────────

"""
    intensities(series::AbstractChromScanSeries;
    unit::Union{Nothing, Unitful.Units}=nothing) -> Vector{<:Real} or Vector{<:Quantity}

Return the intensity values from each scan in a chromatographic scan series.  
If the scan series stores unitful intensities and a `unit` is specified, the values are 
converted accordingly.

# Arguments
- `series`: A concrete subtype of `AbstractChromScanSeries` containing chromatographic 
  scans.
- `unit`: *(optional)* A `Unitful.Units` object to convert intensities to. If omitted or 
  `nothing`, intensities are returned as stored.

# Returns
- A vector of intensity values, either as plain numbers or `Quantity`s, depending on the 
  stored data.

# Throws
- `ArgumentError` if the scan series stores unitless intensities but a unit conversion is 
  requested.
- `AssertionError` if the series contains no scans.

See also [`AbstractChromScanSeries`](@ref), [`scans`](@ref), [`scancount`](@ref).

# Examples
```jldoctest
julia> s1 = ChromScan(1.0u"s", 1.0u"pA");

julia> s2 = ChromScan(2.0u"s", 2.0u"pA");

julia> series = ChromScanSeries([s1, s2]);

julia> intensities(series) == [1.0, 2.0]u"pA"
true

julia> intensities(series; unit=u"fA") == [1e3, 2e3]u"fA"
true
```
"""
@inline function intensities(
    series::AbstractChromScanSeries; 
    unit::Union{Nothing, Unitful.Units}=nothing)

    intensity.(scans(series); unit=unit)
end

"""
    intensities(series::AbstractMassScanSeries, scanindex::Integer; 
                unit::Union{Nothing, Unitful.Units}=nothing) -> AbstractVector

Returns the intensity values for the scan at the specified index within the given mass scan 
series.

# Arguments
- `series::AbstractMassScanSeries`: The series of mass scans.
- `scanindex::Integer`: The index of the scan to extract intensities from.
- `unit::Union{Nothing, Unitful.Units}`: Optional unit to convert the returned 
  intensities into. If `nothing`, the scan's native unit is used.

# Returns
- An `AbstractVector` of intensity values (possibly unitful) corresponding to the scan at 
  the specified index.

# Throws
- `BoundsError` if `scanindex` is outside the valid range of the scan series.

# Examples
```jldoctest
julia> ms1 = MassScan(1.0u"s", [100.0, 150.0], [10.0, 20.0]u"pA");

julia> ms2 = MassScan(2.0u"s", [100.0, 150.0], [12.0, 22.0]u"pA");

julia> mss = MassScanSeries([ms1, ms2]; acquisition=(mode="FullScan",));

julia> intensities(mss, 1) == [10.0, 20.0]u"pA"
true

julia> intensities(mss, 2; unit=u"fA") == [12e3, 22e3]u"fA"
true
```
"""
@inline function intensities(
    series::AbstractMassScanSeries,
    scanindex::Integer; 
    unit::Union{Nothing, Unitful.Units}=nothing)

    intensities(scan(series, scanindex); unit=unit)
end

# ── intensity ─────────────────────────────────────────────────────────────────────────────

"""
    intensity(series::AbstractChromScanSeries, scanindex::Integer; 
              unit::Union{Nothing, Unitful.Units}=nothing) -> Real

Returns the intensity value for the chromatographic scan at the specified index within the 
given chromatogram scan series.

# Arguments
- `series::AbstractChromScanSeries`: The series of chromatographic scans.
- `scanindex::Integer`: The index of the scan to extract the intensity from.
- `unit::Union{Nothing, Unitful.Units}`: Optional unit to convert the returned intensity 
  into. If `nothing`, the scan's native unit is used.

# Returns
- A scalar `Real` (possibly unitful) representing the intensity of the scan at the 
  specified index.

# Throws
- `BoundsError` if `scanindex` is outside the valid range of the scan series.

# Examples
```jldoctest
julia> cs1 = ChromScan(1.0u"s", 1.5u"nA");

julia> cs2 = ChromScan(2.0u"s", 2.5u"nA");

julia> css = ChromScanSeries([cs1, cs2]; acquisition=(mode="TIC",));

julia> intensity(css, 1)
1.5 nA

julia> intensity(css, 2; unit=u"pA")
2500.0 pA
```
"""
@inline function intensity(
    series::AbstractChromScanSeries,
    scanindex::Integer; 
    unit::Union{Nothing, Unitful.Units}=nothing)

    intensity(scan(series, scanindex); unit=unit)
end

# ── intensityunit ─────────────────────────────────────────────────────────────────────────

"""
    intensityunit(series::AbstractScanSeries) -> Union{Unitful.Units, Nothing}

Return the unit associated with the intensity values of a scan series.

This function retrieves the intensity unit from the first scan in the series, assuming
all scans have consistent units.

# Arguments
- `series`: A subtype of `AbstractScanSeries`.

# Returns
- A `Unitful.Units` subtype representing the intensity unit, or `nothing` if unspecified.

# Examples
```julia
julia> ms = MassScan(1.0u"s", [100.0], [10.0]u"pA")
       series = MassScanSeries([ms]);

julia> intensityunit(series)
pA
```
"""
intensityunit(series::AbstractScanSeries) = intensityunit(scan(series, 1))

# ── instrument ────────────────────────────────────────────────────────────────────────────

"""
    instrument(series::AbstractScanSeries) -> NamedTuple

Returns the instrument metadata associated with the scan series.

This typically includes information about the hardware used for data acquisition, such as 
detector type, manufacturer, model, or configuration parameters.

# Arguments
- `series`: A concrete subtype of `AbstractScanSeries`.

# Returns
A `NamedTuple` containing instrument-specific metadata.

See also [`AbstractScanSeries`](@ref), [`acquisition`](@ref), [`user`](@ref), 
[`sample`](@ref).

# Examples
```jldoctest
julia> s1 = ChromScan(1.0u"s", 10.0);

julia> s2 = ChromScan(2.0u"s", 200.0);

julia> series = ChromScanSeries([s1, s2], instrument=(detector="UV-DAD",));

julia> instrument(series)
(detector = "UV-DAD",)

julia> instrument(series).detector
"UV-DAD"
```
"""
instrument(series::AbstractScanSeries) = series.instrument

# ── levels ────────────────────────────────────────────────────────────────────────────────

"""
    levels(series::AbstractMassScanSeries) -> Vector{<:Integer}

Returns a sorted list of unique scan levels present in the given mass scan series.

Scan levels typically indicate the stage of mass spectrometry (e.g. MS¹, MS²). This 
function extracts the level from each scan, filters duplicates, and returns them in
ascending order.

# Arguments
- `series::AbstractMassScanSeries`: The series of mass scans to analyze.

# Returns
- A sorted vector of unique scan levels.

# Examples
```jldoctest
julia> ms1 = MassScan(1.0u"s", [100.0], [10.0]u"pA"; level=1);

julia> ms2 = MassScan(2.0u"s", [150.0], [20.0]u"pA"; level=2);

julia> ms3 = MassScan(3.0u"s", [200.0], [30.0]u"pA"; level=1);

julia> mss = MassScanSeries([ms1, ms2, ms3]; acquisition=(mode="DDA",));

julia> levels(mss) == [1, 2]
true
```
"""
levels(series::AbstractMassScanSeries) = sort(unique(level.(scans(series))))

# ── mzunit ────────────────────────────────────────────────────────────────────────────────

"""
    mzunit(series::AbstractMassScanSeries) -> Union{Unitful.Units, Nothing}

Return the unit associated with the m/z values of a mass scan series.

This function retrieves the m/z unit from the first scan in the series, assuming
all scans have consistent units.

# Arguments
- `series`: A subtype of `AbstractMassScanSeries`.

# Returns
- A `Unitful.Units` subtype representing the m/z unit, or `nothing` if unspecified.

# Examples
```jldoctest
julia> ms = MassScan(1.0u"s", [100.0]u"Th", [10.0])
       series = MassScanSeries([ms]);

julia> mzunit(series)
Th
```
"""
mzunit(series::AbstractMassScanSeries) = mzunit(scan(series, 1))

# ── mzvalues ──────────────────────────────────────────────────────────────────────────────

"""
    mzvalues(series::AbstractMassScanSeries, scanindex::Integer; 
             unit::Union{Nothing, Unitful.Units}=nothing) -> AbstractVector{<:Real}

Returns the m/z (mass-to-charge ratio) values for the scan at the specified index in the 
given mass scan series.

If the scan stores unitful m/z values, they are optionally converted to the specified unit.
If no unit is specified, the scan's native unit is used. If the m/z values are unitless, 
they are returned as-is.

# Arguments
- `series`: A mass scan series (`AbstractMassScanSeries`) to index into.
- `scanindex`: The index of the scan to retrieve m/z values from.
- `unit`: *(optional)* Unit to which m/z values should be converted. Must be `nothing` for 
  scans with unitless m/z values.

# Returns
- An `AbstractVector{<:Real}` or `AbstractVector{<:Quantity}` containing the m/z values of 
  the specified scan.

# Throws
- `BoundsError` if the index is out of range.
- `ArgumentError` if a unit is requested for a scan with unitless m/z values.

# Examples
```jldoctest
julia> s1 = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]u"pA");

julia> series = MassScanSeries([s1]);

julia> mzvalues(series, 1) == [100.0, 200.0]
true
```
"""
@inline function mzvalues(series::AbstractMassScanSeries, scanindex::Integer; 
    unit::Union{Nothing, Unitful.Units}=nothing)

    mzvalues(scan(series, scanindex); unit=unit)
end

# ── rawintensities ────────────────────────────────────────────────────────────────────────

"""
    rawintensities(series::AbstractChromScanSeries;
    unit::Union{Nothing, Unitful.Units}=nothing) -> Vector

Return the numeric (unitless) intensity values from all scans in a chromatographic scan 
series.

If the scan’s `intensity_unit` is not `nothing`, the intensity values are optionally 
converted to the user-specified `unit` using `uconvert`, then stripped of units. If no 
`unit` is specified, the stored unit is used for stripping.

# Arguments
- `series`: A concrete subtype of `AbstractChromScanSeries`.
- `unit`: *(optional)* A `Unitful.Units` object to which raw intensities should be 
  converted. If the scan series does not define intensity units, this must remain `nothing`.

# Returns
A vector of raw intensity values, one per scan, with optional unit conversion.

# Throws
- `ArgumentError` if a unit is requested but the scan series has no defined intensity unit.

See also [`AbstractChromScanSeries`](@ref), [`intensities`](@ref), [`rawintensity`](@ref), 
[`scans`](@ref).

# Examples
```jldoctest
julia> s1 = ChromScan(1.0u"s", 1.0u"pA");

julia> s2 = ChromScan(2.0u"s", 2.0u"pA");

julia> series = ChromScanSeries([s1, s2]);

julia> rawintensities(series) == [1.0, 2.0]
true
```
"""
@inline function rawintensities(
    series::AbstractChromScanSeries;
    unit::Union{Nothing, Unitful.Units}=nothing)

    rawintensity.(scans(series); unit=unit)
end

"""
    rawintensities(series::AbstractMassScanSeries, scanindex::Integer; 
                   unit::Union{Nothing, Unitful.Units}=nothing) -> AbstractVector{<:Real}

Returns the numeric (unitless) intensity values for the scan at the specified index within
the given mass scan series.

If the scan stores unitful intensities, the values are optionally converted to the 
specified unit and then stripped of units. If the intensities are already unitless, 
no conversion is performed.

# Arguments
- `series::AbstractMassScanSeries`: The series of mass scans.
- `scanindex::Integer`: The index of the scan to extract raw intensities from.
- `unit::Union{Nothing, Unitful.Units}`: Optional unit to convert the intensities into 
  before stripping. If `nothing`, the stored unit (if any) is used. If intensities are 
  already unitless, `unit` must also be `nothing`.

# Returns
- An `AbstractVector{<:Real}` of unitless intensity values from the specified scan.

# Throws
- `BoundsError` if `scanindex` is outside the valid range of the scan series.
- `ArgumentError` if a unit is requested for a scan that stores unitless intensities.
- `AssertionError` if the scan claims to have a unit but `intensity_unit` is `nothing`.

# Examples
```jldoctest
julia> ms1 = MassScan(1.0u"s", [100.0, 150.0], [10.0, 20.0]u"pA");

julia> ms2 = MassScan(2.0u"s", [100.0, 150.0], [12.0, 22.0]u"pA");

julia> mss = MassScanSeries([ms1, ms2]; acquisition=(mode="FullScan",));

julia> rawintensities(mss, 1) == [10.0, 20.0]
true

julia> rawintensities(mss, 2; unit=u"nA") == [0.012, 0.022]
true
```
"""
@inline function rawintensities(
    series::AbstractMassScanSeries,
    scanindex::Integer; 
    unit::Union{Nothing, Unitful.Units}=nothing)

    rawintensities(scan(series, scanindex); unit=unit)
end

# ── rawretentions ─────────────────────────────────────────────────────────────────────────

"""
    rawretentions(series::AbstractScanSeries; 
                  unit::Union{Nothing, Unitful.Units}=nothing) -> Vector{<:Real}

Returns the raw (unitless) retention values for all scans in the given scan series.

If the scan stores unitful retention values, they are optionally converted to the specified 
unit and then stripped of units. If retention values are already unitless, no conversion is 
performed.

# Arguments
- `series`: A subtype of `AbstractScanSeries` whose scans define retention values.
- `unit`: *(optional)* A `Unitful.Units` object to which each scan's retention value should 
  be converted before stripping. Must be `nothing` for scans with unitless retention.

# Returns
- A `Vector{<:Real}` of raw, unitless retention values — one per scan in the series.

# Throws
- `ArgumentError` if a unit is requested for a scan that does not define a retention unit.

# Examples
```jldoctest
julia> s1 = ChromScan(1.0u"s", 1.0);

julia> s2 = ChromScan(2.0u"s", 2.0);

julia> series = ChromScanSeries([s1, s2]);

julia> rawretentions(series) == [1.0, 2.0]
true

julia> rawretentions(series; unit=u"ms") == [1e3, 2e3]
true
```
"""
@inline function rawretentions(
    series::AbstractScanSeries;
    unit::Union{Nothing, Unitful.Units}=nothing)

    rawretention.(scans(series); unit=unit)
end

# ── retentions ────────────────────────────────────────────────────────────────────────────

"""
    retentions(series::AbstractScanSeries; 
               unit::Union{Nothing, Unitful.Units}=nothing) -> Vector{<:Real}

Returns the retention values for all scans in the given scan series, preserving units.

If the scan stores unitful retention values, the values are optionally converted to the 
specified unit. If no unit is specified, the scan's stored unit is used. If the retention 
values are unitless, no conversion is performed and the values are returned as-is.

# Arguments
- `series`: A subtype of `AbstractScanSeries` whose scans define retention values.
- `unit`: *(optional)* A `Unitful.Units` object to which each scan's retention value should 
  be converted. Must be `nothing` for scans with unitless retention.

# Returns
- A `Vector{<:Real}` or `Vector{<:Quantity}` of retention values, one per scan.

# Throws
- `ArgumentError` if a unit is requested for a scan that does not define a retention unit.

# Examples
```jldoctest
julia> s1 = ChromScan(1.0u"s", 100.0);

julia> s2 = ChromScan(2.0u"s", 200.0);

julia> series = ChromScanSeries([s1, s2]);

julia> retentions(series) == [1.0, 2.0]u"s"
true

julia> retentions(series; unit=u"ms") == [1000.0, 2000.0]u"ms"
true
```
"""
@inline function retentions(
    series::AbstractScanSeries;
    unit::Union{Nothing, Unitful.Units}=nothing)

    retention.(scans(series); unit=unit)
end

# ── retentionunit ─────────────────────────────────────────────────────────────────────────

"""
    retentionunit(series::AbstractScanSeries) -> Union{Unitful.Units, Nothing}

Return the unit associated with the retention values of a scan series.

This function retrieves the retention unit from the first scan in the series, assuming
all scans have consistent units (which should be enforced by the series constructor).

# Arguments
- `series`: A subtype of `AbstractScanSeries`, such as `MassScanSeries` or `ChromScanSeries`.

# Returns
- A `Unitful.Units` subtype representing the retention unit, or `nothing` if unspecified.

# Examples
```jldoctest
julia> ms = MassScan(1.0u"s", [100.0], [10.0])
       series = MassScanSeries([ms]);

julia> retentionunit(series)
s
```
"""
retentionunit(series::AbstractScanSeries) = retentionunit(scan(series, 1))

# ── sample ────────────────────────────────────────────────────────────────────────────────

"""
    sample(series::AbstractScanSeries) -> NamedTuple

Returns the sample metadata associated with the scan series.

This typically includes information about the analyzed sample, such as sample ID, origin, 
preparation details, or treatment conditions.

# Arguments
- `series`: A concrete subtype of `AbstractScanSeries`.

# Returns
A `NamedTuple` containing sample-related metadata.

See also [`AbstractScanSeries`](@ref), [`instrument`](@ref), [`acquisition`](@ref), 
[`user`](@ref).

# Examples
```jldoctest
julia> s = [ChromScan(1.0u"s", 10.0), ChromScan(2.0u"s", 200.0)];

julia> series = ChromScanSeries(s, sample=(ID="Polistes dominula", locality="Germany"));

julia> sample(series)
(ID = "Polistes dominula", locality = "Germany")

julia> sample(series).ID
"Polistes dominula"
```
"""
sample(series::AbstractScanSeries) = series.sample

# ── scan ──────────────────────────────────────────────────────────────────────────────────

"""
    scan(series::AbstractScanSeries, scanindex::Integer) -> AbstractScan

Returns the scan at the specified index from the given scan series.

This is a direct accessor for individual scans, supporting bounds-checked indexed access.

# Arguments
- `series::AbstractScanSeries`: The scan series to index into.
- `scanindex::Integer`: The index of the scan to retrieve.

# Returns
- An `AbstractScan` object corresponding to the requested index.

# Throws
- `BoundsError` if `scanindex` is outside the valid range of the scan series.

# Examples
```jldoctest
julia> s1 = ChromScan(1.0u"s", 2.0u"nA");

julia> s2 = ChromScan(2.0u"s", 2.5u"nA");

julia> css = ChromScanSeries([s1, s2]; acquisition=(mode="TIC",));

julia> scan(css, 2) == s2
true
```
"""
@inline function scan(series::AbstractScanSeries, scanindex::Integer)
    vec = scans(series)
    checkbounds(vec, scanindex)
    vec[scanindex]
end

# ── scancount ─────────────────────────────────────────────────────────────────────────────

"""
    scancount(series::AbstractScanSeries) -> Int

Returns the number of individual scan entries in the scan series.

This is equivalent to `length(series)`, and reflects the total number of scan measurements 
stored in the series.

# Arguments
- `series`: A concrete subtype of `AbstractScanSeries`.

# Returns
An `Int` representing the number of scan elements.

See also [`AbstractScanSeries`](@ref), [`scans`](@ref).

# Examples
```jldoctest
julia> series = ChromScanSeries([ChromScan(1.0u"s", 10.0), ChromScan(2.0u"s", 200.0)]);

julia> scancount(series)
2

julia> scancount(series) == length(series)
true
```
"""
scancount(series::AbstractScanSeries) = length(scans(series))

# ── scans ─────────────────────────────────────────────────────────────────────────────────

"""
    scans(series::AbstractScanSeries) -> AbstractVector{<:AbstractScan}

Return the collection of scans contained in the given scan series.

# Arguments
- `series`: A subtype of `AbstractScanSeries`, such as `ChromScanSeries` or 
  `MassScanSeries`.

# Returns
- A vector of scan objects in acquisition order.

# Examples
```jldoctest
julia> s1 = ChromScan(1.0u"s", 5.0u"pA");

julia> s2 = ChromScan(2.0u"s", 7.5u"pA");

julia> series = ChromScanSeries([s1, s2]);

julia> length(scans(series)) == scancount(series)
true

julia> all(retention.(scans(series)) .== [1.0u"s", 2.0u"s"])
true

julia> all(intensity.(scans(series)) .== [5.0u"pA", 7.5u"pA"])
true
```
"""
scans(series::AbstractScanSeries) = series.scans

# ── uniquemzvalues ────────────────────────────────────────────────────────────────────────

"""
    uniquemzvalues(series::AbstractMassScanSeries, target_level::Integer=1; 
                   unit::Union{Nothing, Unitful.Units}=nothing) -> Vector{<:Real}

Returns a sorted vector of unique m/z values (ions) found in scans at the specified MS 
level within the given mass scan series.

If the scans store unitful m/z values, they are optionally converted to the specified 
`unit`. If no unit is provided, the scan's native unit is used. If m/z values are unitless, 
no  conversion is performed and values are returned as-is.

# Arguments
- `series`: The mass scan series to extract m/z values from.
- `target_level`: *(optional)* The MS level to filter scans by (default: 1).
- `unit`: *(optional)* Unit to which m/z values should be converted before deduplication. 
  Must be `nothing` for unitless scans.

# Returns
- A sorted vector of unique m/z values (`Real` or `Quantity`), depending on unit usage.

# Throws
- `ArgumentError` if no scans exist at the specified level.
- `ArgumentError` if a unit is requested for a scan with unitless m/z values.

# Examples
```jldoctest
julia> s1 = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]; level=1);

julia> s2 = MassScan(2.0u"s", [150.0, 200.0], [15.0, 25.0]; level=2);

julia> mss = MassScanSeries([s1, s2]);

julia> uniquemzvalues(mss) == [100.0, 200.0]
true
```
"""
@inline function uniquemzvalues(series::AbstractMassScanSeries, target_level::Integer=1;
    unit::Union{Nothing, Unitful.Units}=nothing)

    # Filter scans by target level
    scans_at_level = filter(scan -> level(scan) == target_level, scans(series))

    isempty(scans_at_level) && throw(
        ArgumentError("Mass scan series has no scan data for level $target_level"))

    # Extract, concatenate, and deduplicate ion values
    sort(unique(vcat(mzvalues.(scans_at_level, unit=unit)...)))
end

# ── user ──────────────────────────────────────────────────────────────────────────────────

"""
    user(series::AbstractScanSeries) -> NamedTuple

Returns the user metadata associated with the scan series.

This typically includes information about the operator, analyst, or laboratory personnel 
involved in the acquisition or processing of the scan data.

# Arguments
- `series`: A concrete subtype of `AbstractScanSeries`.

# Returns
A `NamedTuple` containing user-related metadata.

See also [`AbstractScanSeries`](@ref), [`acquisition`](@ref), [`instrument`](@ref), 
[`sample`](@ref).

# Examples
```jldoctest
julia> vec = [ChromScan(1.0u"s", 10.0), ChromScan(2.0u"s", 200.0)];

julia> series = ChromScanSeries(vec, user=(name="Alice", email=""));

julia> user(series).name
"Alice"
```
"""
user(series::AbstractScanSeries) = series.user
