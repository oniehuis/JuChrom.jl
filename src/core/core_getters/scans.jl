# ── attrs ──────────────────────────────────────────────────────────────────────────────

"""
    attrs(scan::AbstractScan)

Return the scan-level attrs associated with the given `scan` object.

The attrs is stored as a `NamedTuple` and may contain additional contextual or 
acquisition-specific information beyond the primary measurement fields (e.g., scan ID, 
instrument settings, annotations).

# Arguments
- `scan`: A subtype of `AbstractScan`

# Returns
A `NamedTuple` containing the attrs fields defined for the scan.

See also [`AbstractScan`](@ref), [`AbstractChromScan`](@ref), [`AbstractMassScan`](@ref).

# Examples
```jldoctest
julia> csc = ChromScan(5.0, 1200.0, attrs=(scan_id=42, comment="baseline"));

julia> attrs(csc)
(scan_id = 42, comment = "baseline")

julia> attrs(csc).scan_id
42
```
"""
attrs(scan::AbstractScan) = scan.attrs

# ── intensity ─────────────────────────────────────────────────────────────────────────────

"""
    intensity(scan::AbstractChromScan; unit::Union{Nothing, Unitful.Units}=nothing)

Return the intensity value of a chromatographic scan, optionally as a `Unitful.Quantity`.

If the scan’s `intensity_unit` is not `nothing`, the intensity is returned as a quantity 
in that unit. If a `unit` is specified, the value is converted using `uconvert`.

If `intensity_unit` is `nothing`:
- and a `unit` is provided, an error is thrown
- and no `unit` is provided, the raw numeric value is returned unchanged

This function does **not** assign units to raw values — it only performs conversions when a 
unit is already stored.

# Arguments
- `scan`: A subtype of `AbstractChromScan`
- `unit`: (Optional) target unit (e.g. `u"pA"`) — used only if a unit is stored

# Returns
- A `Quantity` in the scan’s stored or requested unit, if applicable
- A plain `Real` value if no unit is stored and no `unit` is requested

# Throws
- `ArgumentError` if `unit` is specified for a scan with no stored `intensity_unit`

See also [`AbstractChromScan`](@ref), [`AbstractScan`](@ref), [`rawintensity`](@ref), 
[`intensityunit`](@ref).

# Examples
```jldoctest
julia> csc = ChromScan(3.0u"minute", 100u"pA");

julia> intensity(csc)
100 pA

julia> intensity(csc; unit=u"fA")
100000 fA

julia> csc2 = ChromScan(5000, 100);  # intensity_unit == nothing

julia> intensity(csc2)
100

julia> intensity(csc2; unit=u"pA")
ERROR: ArgumentError: Cannot convert unitless intensity to a unit
[...]
```
"""
@inline function intensity(
    scan::AbstractChromScan{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(scan.intensity, unit, "intensity")
end

@inline function intensity(
    scan::AbstractChromScan{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert_scalar(scan.intensity, scan.intensity_unit, unit)
end

# ── intensityunit ─────────────────────────────────────────────────────────────────────────

"""
    intensityunit(scan::AbstractScan) -> Union{Unitful.Units, Nothing}

Return the unit associated with the intensity value(s) of a scan.

This function retrieves the `intensity_unit` field from any scan subtype, which indicates 
the physical unit used for the signal intensity or intensities (e.g. `u"pA"`, `u"counts"`). 
If no unit was specified at construction, returns `nothing`.

# Arguments
- `scan`: A subtype of `AbstractScan`, such as `ChromScan` or `MassScan`

# Returns
- A `Unitful.Units` subtype representing the intensity unit, or `nothing` if unspecified


See also [`AbstractScan`](@ref), [`AbstractChromScan`](@ref), [`AbstractMassScan`](@ref), 
[`intensity`](@ref), [`rawintensity`](@ref), [`intensities`](@ref), 
[`rawintensities`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]u"pA");

julia> intensityunit(msc)
pA

julia> csc = ChromScan(1.0u"minute", 42.0);

julia> intensityunit(csc) === nothing
true
```
"""
intensityunit(scan::AbstractScan) = scan.intensity_unit

# ── intensities ───────────────────────────────────────────────────────────────────────────

"""
    intensities(scan::AbstractMassScan; unit=nothing)

Return the intensity values from a mass spectrometric scan.

This method supports both unitless and unitful scans. If `unit` is specified, the 
intensities are converted to the desired unit using `Unitful.uconvert`. If the scan is 
unitless and a unit is requested, an error is thrown.

# Arguments
- `scan`: A subtype of `AbstractMassScan`
- `unit`: Desired unit to convert intensities into (`Unitful.Units` subtype or `nothing`, 
  default)

# Returns
A vector of intensity values. If a unit is provided, the values are returned in that unit.

# Throws
- `ArgumentError` if:
  - The scan has no intensity unit and a unit conversion is requested.
- `AssertionError` if:
  - The scan claims to have a unit but `intensity_unit` is `nothing` (internal 
    inconsistency)


See also [`AbstractMassScan`](@ref), [`AbstractScan`](@ref), [`rawintensities`](@ref), 
[`intensityunit`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]u"pA");

julia> intensities(msc) == [10.0, 20.0]u"pA"
true

julia> intensities(msc, unit=u"nA") == [0.01, 0.02]u"nA"
true
```
"""
@inline function intensities(
    scan::AbstractMassScan{<:Any, <:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(scan.intensities, unit, "intensities")
end

@inline function intensities(
    scan::AbstractMassScan{<:Any, <:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(scan.intensities, scan.intensity_unit, unit)
end

# ── level ────────────────────────────────────────────────────────────────────

"""
    level(scan::AbstractMassScan) -> Integer

Return the MS level of a mass spectrometric scan.

The MS level indicates the stage of mass spectrometry at which the scan was acquired:
- `1` for full MS1 scans (e.g. precursor spectra),
- `2` for MS/MS scans (e.g. product ion spectra),
- Higher values for advanced fragmentation strategies (e.g. MS^n).

# Arguments
- `scan`: A subtype of `AbstractMassScan`, such as `MassScan`

# Returns
- An `Integer` value representing the MS level (guaranteed to be ≥ 1)


See also [`AbstractMassScan`](@ref), [`AbstractScan`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]);

julia> level(msc)
1
```
"""
level(scan::AbstractMassScan) = scan.level

# ── mzcount ───────────────────────────────────────────────────────────────────────────────

"""
    mzcount(scan::AbstractMassScan) -> Int

Return the number of m/z values in the scan.

This is a simple utility function that returns the length of the `mz_values` vector stored 
in a mass spectrometric scan. It reflects the number of data points (peaks) measured at a 
given separation coordinate.

# Arguments
- `scan`: A subtype of `AbstractMassScan`, such as `MassScan`

# Returns
- An `Int` representing the number of m/z values in the scan


See also [`AbstractMassScan`](@ref), [`AbstractScan`](@ref), [`mzvalues`](@ref), 
[`rawmzvalues`](@ref), [`intensities`](@ref), [`rawintensities`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(0.5u"s", [100.0, 150.0, 200.0], [10.0, 20.0, 30.0]);

julia> mzcount(msc)
3
```
"""
mzcount(scan::AbstractMassScan) = length(scan.mz_values)

# ── mzunit ────────────────────────────────────────────────────────────────────────────────

"""
    mzunit(scan::AbstractMassScan) -> Union{Unitful.Units, Nothing}

Return the unit associated with the mass-to-charge (m/z) values of a scan.

Typically, this is `nothing`, since m/z values are conventionally unitless. If a unit was
explicitly provided during scan construction, it will be returned.

# Arguments
- `scan`: A subtype of `AbstractMassScan`, such as `MassScan`

# Returns
- A `Unitful.Units` subtype representing the m/z unit, or `nothing` if unspecified


See also [`AbstractMassScan`](@ref), [`AbstractScan`](@ref), [`mzvalues`](@ref), 
[`rawmzvalues`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]);

julia> mzunit(msc) === nothing
true
```
"""
mzunit(scan::AbstractMassScan) = scan.mz_unit

# ── mzvalues───────────────────────────────────────────────────────────────────────────────

"""
    mzvalues(scan::AbstractMassScan; unit=nothing)

Return the mass-to-charge (m/z) values from a mass spectrometric scan.

This method supports both unitless and unitful scans. If the scan has an associated m/z 
unit, values can be converted to a different unit using the optional `unit` argument. 
If the scan is unitless and a unit is requested, an error is thrown.

# Arguments
- `scan`: A subtype of `AbstractMassScan`
- `unit`: Desired unit to convert m/z values into (`Unitful.Units` subtype or `nothing`, 
  default)

# Returns
A vector of m/z values. If a unit is specified, values are converted to that unit; 
otherwise, values are returned in their stored form.

# Throws
- `ArgumentError` if:
  - The scan is unitless and a unit is requested.
- `AssertionError` if:
  - The scan claims to have a unit but `mz_unit` is `nothing` (internal inconsistency)

See also [`AbstractMassScan`](@ref), [`AbstractScan`](@ref), [`rawmzvalues`](@ref), 
[`mzunit`](@ref), [`mzcount`](@ref), [`intensities`](@ref), [`rawintensities`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]);

julia> mzvalues(msc) == [100.0, 200.0]
true
```
"""
@inline function mzvalues(
    scan::AbstractMassScan{<:Any, Nothing, <:Any};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(scan.mz_values, unit, "m/z values")
end

@inline function mzvalues(
    scan::AbstractMassScan{<:Any, <:Unitful.Units, <:Any};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(scan.mz_values, scan.mz_unit, unit)
end

# ── rawintensities ────────────────────────────────────────────────────────────────────────

"""
    rawintensities(scan::AbstractMassScan; unit=nothing)

Return the raw (unitless) numeric intensity values from a mass spectrometric scan.

If the scan stores unitful intensities, the unit is stripped using `Unitful.ustrip`. If a
specific unit is requested via the `unit` keyword, the values are first converted to that
unit before stripping.

# Arguments
- `scan`: A subtype of `AbstractMassScan`
- `unit`: Desired unit to convert intensities into before stripping 
  (`Unitful.Units` subtype or `nothing`, default)

# Returns
A vector of numeric (unitless) intensity values.

# Throws
- `ArgumentError` if:
  - The scan is unitless but a unit conversion is requested.
- `AssertionError` if:
  - The scan claims to have a unit but `intensity_unit` is `nothing` (internal 
    inconsistency)


See also [`AbstractMassScan`](@ref), [`AbstractScan`](@ref), [`intensities`](@ref), 
[`intensityunit`](@ref), [`mzvalues`](@ref), [`rawmzvalues`](@ref), [`mzcount`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]u"pA");

julia> rawintensities(msc) == [10.0, 20.0]
true

julia> rawintensities(msc, unit=u"nA") == [0.01, 0.02]
true
```
"""
@inline function rawintensities(
    scan::AbstractMassScan{<:Any, <:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(scan.intensities, unit, "intensities")
end

@inline function rawintensities(
    scan::AbstractMassScan{<:Any, <:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip(scan.intensities, scan.intensity_unit, unit)
end

# ── rawintensity ──────────────────────────────────────────────────────────────────────────

"""
    rawintensity(scan::AbstractChromScan; unit::Union{Nothing, Unitful.Units}=nothing)

Return the numeric intensity value of a chromatographic scan as a unitless `Real`.

If the scan’s `intensity_unit` is not `nothing`, the intensity value is optionally 
converted to the user-specified `unit` using `uconvert`, then stripped of units. If no 
`unit` is specified, the stored unit is used for stripping.

If `intensity_unit` is `nothing`:
- and a `unit` is requested, an error is thrown
- and no `unit` is requested, the raw numeric value is returned unchanged

# Arguments
- `scan`: A subtype of `AbstractChromScan`
- `unit`: (Optional) target unit for conversion before stripping — must not be used for 
  unitless scans

# Returns
A `Real` value with all units stripped.

# Throws
- `ArgumentError` if a `unit` is requested for a scan with no stored `intensity_unit`

# Behavior Summary
- If `intensity_unit ≠ nothing` and `unit` is provided → convert, then strip
- If `intensity_unit ≠ nothing` and `unit == nothing` → strip using stored unit
- If `intensity_unit == nothing` and `unit == nothing` → return raw value
- If `intensity_unit == nothing` and `unit ≠ nothing` → throw `ArgumentError`

See also [`AbstractChromScan`](@ref), [`AbstractScan`](@ref), [`intensity`](@ref), 
[`intensityunit`](@ref).

# Examples
```jldoctest
julia> csc = ChromScan(5.0u"minute", 100u"pA");

julia> rawintensity(csc)
100

julia> rawintensity(csc; unit=u"fA")
100000

julia> csc2 = ChromScan(5000, 200);  # unitless intensity

julia> rawintensity(csc2)
200

julia> rawintensity(csc2; unit=u"pA")
ERROR: ArgumentError: Cannot convert unitless intensity to a unit
[...]
```
"""
@inline function rawintensity(
    scan::AbstractChromScan{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(scan.intensity, unit, "intensity")
end

@inline function rawintensity(
    scan::AbstractChromScan{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip_scalar(scan.intensity, scan.intensity_unit, unit)
end

# ── rawmzvalues ───────────────────────────────────────────────────────────────────────────

"""
    rawmzvalues(scan::AbstractMassScan; unit=nothing)

Return the raw (unit-stripped) m/z values from a mass spectrometric scan.

This method returns the numeric representation of m/z values, optionally converted 
to a different unit. If the scan is unitless and a unit conversion is requested, 
an error is thrown.

# Arguments
- `scan`: A subtype of `AbstractMassScan`
- `unit`: Desired unit to convert values into before stripping units 
  (`Unitful.Units` subtype or `nothing`, default)

# Returns
A vector of plain numeric m/z values, optionally converted to the requested unit 
and stripped of units.

# Throws
- `ArgumentError` if:
  - The scan is unitless and a unit is requested.
- `AssertionError` if:
  - The scan claims to have a unit but `mz_unit` is `nothing` (internal inconsistency)


See also [`AbstractMassScan`](@ref), [`AbstractScan`](@ref), [`mzvalues`](@ref), 
[`mzunit`](@ref), [`mzcount`](@ref), [`intensities`](@ref), [`rawintensities`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]);

julia> rawmzvalues(msc) == [100.0, 200.0]
true
```
"""
@inline function rawmzvalues(
    scan::AbstractMassScan{<:Any, <:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(scan.mz_values, unit, "m/z values")
end

@inline function rawmzvalues(
    scan::AbstractMassScan{<:Any, <:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip(scan.mz_values, scan.mz_unit, unit)
end

# ── rawretention ──────────────────────────────────────────────────────────────────────────

"""
    rawretention(scan::AbstractScan; unit::Union{Nothing, Unitful.Units}=nothing)

Return the numeric retention value of a scan, always as a unitless `Real`.

If the scan’s `retention_unit` is not `nothing`, the value is optionally converted to the 
user-specified `unit` using `uconvert`, then stripped of units. If no `unit` is specified, 
the stored unit is used for stripping.

If `retention_unit` is `nothing`:
- and a `unit` is requested, an error is thrown
- and no `unit` is requested, the raw numeric value is returned unchanged

# Arguments
- `scan`: A subtype of `AbstractScan`
- `unit`: (Optional) target unit to convert to before unit stripping — must not be used for 
  unitless scans

# Returns
A `Real` value with all units stripped.

# Throws
- `ArgumentError` if a `unit` is requested for a scan with no stored `retention_unit`

# Behavior Summary
- If `retention_unit ≠ nothing` and `unit` is provided → convert, then strip
- If `retention_unit ≠ nothing` and `unit == nothing` → strip using stored unit
- If `retention_unit == nothing` and `unit == nothing` → return raw value
- If `retention_unit == nothing` and `unit ≠ nothing` → throw `ArgumentError`

See also [`AbstractScan`](@ref), [`AbstractChromScan`](@ref), [`AbstractMassScan`](@ref), 
[`retention`](@ref), [`retentionunit`](@ref).

# Examples
```jldoctest
julia> csc = ChromScan(5.0u"minute", 1000);

julia> rawretention(csc)
5.0

julia> rawretention(csc; unit=u"s")
300.0

julia> csc2 = ChromScan(5000, 200);  # unitless retention

julia> rawretention(csc2)
5000

julia> rawretention(csc2; unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention to a unit
[...]
```
"""
@inline function rawretention(
    scan::AbstractScan{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(scan.retention, unit, "retention")
end

@inline function rawretention(
    scan::AbstractScan{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip_scalar(scan.retention, scan.retention_unit, unit)
end

# ── retention ─────────────────────────────────────────────────────────────────────────────

"""
    retention(scan::AbstractScan; unit::Union{Nothing, Unitful.Units}=nothing)

Return the retention value of a scan, optionally converted to a specified unit.

If the scan’s `retention_unit` is not `nothing`, the return value is a `Unitful.Quantity`. 
If a `unit` is provided, the value is converted to that unit using `uconvert`.

If `retention_unit` is `nothing`:
- and a `unit` is requested, an error is thrown
- and no `unit` is requested, the raw numeric value is returned unchanged

This function does **not** assign units to unitless scans — it only performs conversions 
when a unit is already stored in `retention_unit`.

# Arguments
- `scan`: A subtype of `AbstractScan`
- `unit`: (Optional) target unit (e.g. `u"s"` or `u"minute"`) — used only if a unit is 
  stored

# Returns
- A `Quantity` in the scan’s stored or requested unit, if applicable
- A plain `Real` value if the scan is unitless and no `unit` is requested

# Throws
- `ArgumentError` if a `unit` is requested for a unitless scan

See also [`AbstractScan`](@ref), [`AbstractChromScan`](@ref), [`AbstractMassScan`](@ref), 
[`rawretention`](@ref), [`retentionunit`](@ref).

# Examples
```jldoctest
julia> csc = ChromScan(3.0u"minute", 1000);

julia> retention(csc)
3.0 minute

julia> retention(csc; unit=u"s")
180.0 s

julia> csc2 = ChromScan(5000, 100);  # unitless retention

julia> retention(csc2)
5000

julia> retention(csc2; unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention to a unit
[...]
```
"""
@inline function retention(
    scan::AbstractScan{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(scan.retention, unit, "retention")
end

@inline function retention(
    scan::AbstractScan{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert_scalar(scan.retention, scan.retention_unit, unit)
end

# ── retentionunit ─────────────────────────────────────────────────────────────────────────

"""
    retentionunit(scan::AbstractScan) -> Union{Unitful.Units, Nothing}

Return the unit associated with the retention value of a scan.

This function retrieves the `retention_unit` field from any scan subtype, which indicates 
the physical unit used along the separation axis (e.g. `u"s"`, `u"minute"`). If no unit 
was specified at construction, returns `nothing`.

# Arguments
- `scan`: A subtype of `AbstractScan`, such as `ChromScan` or `MassScan`

# Returns
- A `Unitful.Units` subtype representing the retention unit, or `nothing` if unspecified

See also [`AbstractScan`](@ref), [`AbstractChromScan`](@ref), [`AbstractMassScan`](@ref), 
[`retention`](@ref), [`rawretention`](@ref).

# Examples
```jldoctest
julia> msc = MassScan(1.0u"s", [100.0, 200.0], [10.0, 20.0]);

julia> retentionunit(msc)
s

julia> csc = ChromScan(1.0u"minute", 42.0);

julia> retentionunit(csc)
minute
```
"""
retentionunit(scan::AbstractScan) = scan.retention_unit
