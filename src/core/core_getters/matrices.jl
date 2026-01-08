# ── acquisition ───────────────────────────────────────────────────────────────────────────

"""
    acquisition(msm::AbstractMassScanMatrix) -> NamedTuple

Returns the acquisition metadata associated with the mass scan matrix.

The acquisition metadata typically includes details such as scan mode, method parameters,
instrument configuration, or other context relevant to how the mass spectrometry data was 
acquired.

# Arguments
- `msm`: A concrete subtype of `AbstractMassScanMatrix`.

# Returns
A `NamedTuple` containing acquisition-related metadata.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), [`instrument`](@ref), 
[`sample`](@ref), [`user`](@ref), [`extras`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0], 
                            acquisition=(mode = "DIA", method = "Orbitrap-MS"));

julia> acquisition(msm)
(mode = "DIA", method = "Orbitrap-MS")

julia> acquisition(msm).method
"Orbitrap-MS"
```
"""
acquisition(msm::AbstractMassScanMatrix) = msm.acquisition

# ── extras ────────────────────────────────────────────────────────────────────────────────

"""
    extras(msm::AbstractMassScanMatrix) -> Dict{String, Any}

Returns the unstructured extras associated with the mass scan matrix.

This dictionary may contain arbitrary key-value pairs that are not captured by the 
structured fields such as `instrument`, `acquisition`, `user`, or `sample`.

# Arguments
- `msm`: A concrete subtype of `AbstractMassScanMatrix`.

# Returns
A `Dict{String, Any}` containing unstructured metadata.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), [`acquisition`](@ref), 
[`instrument`](@ref), [`sample`](@ref), [`user`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]; 
                             extras=Dict("run" => 42, "comment" => "Baseline drift"));

julia> extras(msm)
Dict{String, Any} with 2 entries:
  "run"     => 42
  "comment" => "Baseline drift"

julia> extras(msm)["qc_passed"] = true;

julia> extras(msm)
Dict{String, Any} with 3 entries:
  "qc_passed" => true
  "run"       => 42
  "comment"   => "Baseline drift"
```
"""
extras(msm::AbstractMassScanMatrix) = msm.extras

# ── instrument ────────────────────────────────────────────────────────────────────────────

"""
    instrument(msm::AbstractMassScanMatrix) -> NamedTuple

Returns the instrument metadata associated with the mass scan matrix.

This typically includes information about the hardware used for data acquisition, such as 
detector type, manufacturer, model, or configuration parameters.

# Arguments
- `msm`: A concrete subtype of `AbstractMassScanMatrix`.

# Returns
A `NamedTuple` containing instrument-specific metadata.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), [`acquisition`](@ref), 
[`sample`](@ref), [`user`](@ref), [`extras`](@ref).

# Examples
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]; 
                            instrument=(detector = "Orbitrap", manufacturer = "Thermo"));
julia> instrument(msm)
(detector = "Orbitrap", manufacturer = "Thermo")

julia> instrument(msm).detector
"Orbitrap"
```
"""
instrument(msm::AbstractMassScanMatrix) = msm.instrument

# ── intensities ───────────────────────────────────────────────────────────────────────────

"""
    intensities(msm::AbstractMassScanMatrix;
                unit::Union{Nothing, Unitful.Units}=nothing) -> Matrix{<:Real}

Return the intensity values from all scans in a mass scan matrix.

This function retrieves the intensity values from each scan in the mass scan matrix.  
If the matrix stores unitful intensities and a `unit` is specified, the values are 
converted accordingly. If the matrix stores unitless intensities, a unit conversion 
request will throw an `ArgumentError`.

# Arguments
- `msm`: A concrete subtype of `AbstractMassScanMatrix` containing mass scans.
- `unit`: *(optional)* A `Unitful.Units` object to convert intensities to. If omitted or 
  `nothing`, intensities are returned as stored.

# Returns
- A matrix of intensity values, either as plain numbers or `Quantity`s, depending on the 
  stored data.

# Throws
- `ArgumentError` if the scan matrix stores unitless intensities but a unit conversion is 
  requested.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), 
[`intensityunit`](@ref), [`rawintensities`](@ref), [`mzcount`](@ref), [`scancount`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1, 2]u"s", [100.0, 150.0], [10.0 20.0; 15.0 25.0]u"pA");

julia> intensities(msm) ≈ [10.0 20.0; 15.0 25.0]u"pA"
true

julia> intensities(msm, unit=u"nA") ≈ [0.01 0.02; 0.015 0.025]u"nA"
true

julia> msm2 = MassScanMatrix([1, 2]u"s", [100.0, 150.0], [10.0 20.0; 15.0 25.0]);

julia> intensities(msm2) ≈ [10.0 20.0; 15.0 25.0]
true

julia> intensities(msm2, unit=u"nA")
ERROR: ArgumentError: Cannot convert unitless intensities to a unit
[...]
```
"""
@inline function intensities(
    msm::AbstractMassScanMatrix{<:Any, <:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(msm.intensities, unit, "intensities")
end

@inline function intensities(
    msm::AbstractMassScanMatrix{<:Any, <:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(msm.intensities, msm.intensity_unit, unit)
end

# ── intensityunit ─────────────────────────────────────────────────────────────────────────

"""
    intensityunit(msm::AbstractMassScanMatrix) -> Union{Unitful.Units, Nothing}

Return the unit associated with the intensity values of a mass scan matrix.

This function retrieves the `intensity_unit` field from any mass scan matrix subtype, 
which indicates the physical unit used for the signal intensities (e.g. `u"pA"`). If 
no unit was specified at construction, returns `nothing`.

# Arguments
- `msm`: A subtype of `AbstractMassScanMatrix`, such as `MassScanMatrix`.

# Returns
- A `Unitful.Units` subtype representing the intensity unit, or `nothing` if unspecified.

See also [`AbstractMassScanMatrix`](@ref), [`AbstractMassScan`](@ref), 
[`intensities`](@ref), [`rawintensities`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]u"pA");

julia> intensityunit(msm)
pA

julia> msm2 = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> intensityunit(msm2) === nothing
true
```
"""
intensityunit(msm::AbstractMassScanMatrix) = msm.intensity_unit

# ── level ─────────────────────────────────────────────────────────────────────────────────

"""
    level(msm::AbstractMassScanMatrix) -> Integer

Return the MS level of a mass scan matrix.

The MS level indicates the stage of mass spectrometry at which the scans were acquired:
- `1` for full MS1 scans (e.g. precursor spectra),
- `2` for MS/MS scans (e.g. product ion spectra),
- Higher values for advanced fragmentation strategies (e.g. MSⁿ).

# Arguments
- `msm`: A concrete subtype of `AbstractMassScanMatrix`, such as `MassScanMatrix`.

# Returns
An `Integer` value representing the MS level (guaranteed to be ≥ 1).

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), [`levels`](@ref), 
[`levelscans`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]; level=2);

julia> level(msm)
2
```
"""
level(msm::AbstractMassScanMatrix) = msm.level

# ── mzcount ───────────────────────────────────────────────────────────────────────────────

"""
    mzcount(msm::AbstractMassScanMatrix) -> Int

Return the number of m/z values in the mass scan matrix.

This utility function returns the length of the `mz_values` vector stored in a mass 
spectrometric scan matrix. It reflects the number of unique m/z data points (peaks) 
measured across all scans in the matrix, typically corresponding to the columns of the 
intensity matrix.

# Arguments
- `msm::AbstractMassScanMatrix`: A mass scan matrix object, such as `MassScanMatrix`, 
  containing multiple scans and their associated m/z values.

# Returns
- `Int`: The number of m/z values in the scan matrix.

See also [`AbstractMassScanMatrix`](@ref), [`AbstractMassScan`](@ref), [`mzvalues`](@ref), 
[`rawmzvalues`](@ref).

# Examples
julia> msm = MassScanMatrix([1.0]u"s", [100.0, 200.0], [1.0 2.0]);

julia> mzcount(msm)
2

julia> mzcount(msm) == size(intensities(msm), 2)
true
```
"""
mzcount(msm::AbstractMassScanMatrix) = length(msm.mz_values)

# ── mzunit ────────────────────────────────────────────────────────────────────────────────

"""
    mzunit(msm::AbstractMassScanMatrix) -> Union{Unitful.Units, Nothing}

Return the unit associated with the mass-to-charge (m/z) values of a mass scan matrix.

Typically, this is `nothing`, since m/z values are conventionally unitless. If a unit was
explicitly provided during matrix construction, it will be returned.

# Arguments
- `msm`: A subtype of `AbstractMassScanMatrix`, such as `MassScanMatrix`.

# Returns
- A `Unitful.Units` subtype representing the m/z unit, or `nothing` if unspecified.

See also [`AbstractMassScanMatrix`](@ref), [`AbstractMassScan`](@ref), [`mzvalues`](@ref), 
[`rawmzvalues`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0]u"Th", [1.0 2.0; 3.0 4.0]);

julia> mzunit(msm)
Th

julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> mzunit(msm) == nothing
true
```
"""
mzunit(msm::AbstractMassScanMatrix) = msm.mz_unit

# ── mzvalues ──────────────────────────────────────────────────────────────────────────────

"""
    mzvalues(msm::AbstractMassScanMatrix; unit::Union{Nothing, Unitful.Units}=nothing)

Return the mass-to-charge (m/z) values from a mass scan matrix.

This method supports both unitless and unitful matrices. If the matrix has an associated m/z 
unit, values can be converted to a different unit using the optional `unit` argument. 
If the matrix is unitless and a unit is requested, an error is thrown.

# Arguments
- `msm`: A subtype of `AbstractMassScanMatrix`.
- `unit`: Desired unit to convert m/z values into (`Unitful.Units` subtype or `nothing`, default).

# Returns
A vector of m/z values. If a unit is specified, values are converted to that unit; 
otherwise, values are returned in their stored form.

# Throws
- `ArgumentError` if:
    - The matrix is unitless and a unit is requested.
- `AssertionError` if:
    - The matrix claims to have a unit but `mz_unit` is `nothing` (internal inconsistency).

See also [`AbstractMassScanMatrix`](@ref), [`AbstractMassScan`](@ref), [`rawmzvalues`](@ref), 
[`mzunit`](@ref), [`mzcount`](@ref), [`intensities`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0]u"Th", [1.0 2.0; 3.0 4.0]);

julia> mzvalues(msm) == [100.0, 200.0]u"Th"
true

julia> mzvalues(msm, unit=u"kTh") == [0.1, 0.2]u"kTh"
true

julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> mzvalues(msm) == [100.0, 200.0]
true

julia> mzvalues(msm, unit=u"Th")
ERROR: ArgumentError: Cannot convert unitless m/z values to a unit
[...]
```
"""
@inline function mzvalues(
    msm::AbstractMassScanMatrix{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(msm.mz_values, unit, "m/z values")
end

@inline function mzvalues(
    msm::AbstractMassScanMatrix{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(msm.mz_values, msm.mz_unit, unit)
end

# ── rawintensities ────────────────────────────────────────────────────────────────────────

"""
    rawintensities(msm::AbstractMassScanMatrix;
                   unit::Union{Nothing, Unitful.Units}=nothing) -> Matrix

Return the numeric (unitless) intensity values from all scans in a mass scan matrix.

If the scan matrix’s `intensity_unit` is not `nothing`, the intensity values are optionally 
converted to the user-specified `unit` using `uconvert`, then stripped of units. If no 
`unit` is specified, the stored unit is used for stripping.

If the matrix does not define intensity units, requesting a unit conversion will throw an 
`ArgumentError`.

# Arguments
- `msm`: A concrete subtype of `AbstractMassScanMatrix`, such as `MassScanMatrix`.
- `unit`: *(optional)* A `Unitful.Units` object to which the raw intensities should be 
  converted. If the scan matrix does not define intensity units, this must remain `nothing`.

# Returns
A matrix of (unitless) intensity values, with optional unit conversion. Each row typically 
corresponds to a scan, and each column to an m/z value.

# Throws
- `ArgumentError` if a unit is requested but the scan matrix has no defined intensity unit.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), [`intensities`](@ref), 
[`intensityunit`](@ref), [`scancount`](@ref), [`mzcount`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1, 2]u"s", [100.0, 150.0], [10.0 20.0; 15.0 25.0]u"pA");

julia> rawintensities(msm, unit=u"nA") ≈ [0.01 0.02; 0.015 0.025]
true

julia> msm2 = MassScanMatrix([1, 2]u"s", [100.0, 150.0], [10.0 20.0; 15.0 25.0]);

julia> rawintensities(msm2) ≈ [10.0 20.0; 15.0 25.0]
true

julia> rawintensities(msm2, unit=u"pA")
ERROR: ArgumentError: Cannot convert unitless intensities to a unit
[...]
```
"""
@inline function rawintensities(
    msm::AbstractMassScanMatrix{<:Any, <:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(msm.intensities, unit, "intensities")
end

@inline function rawintensities(
    msm::AbstractMassScanMatrix{<:Any, <:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip(msm.intensities, msm.intensity_unit, unit)
end

# ── rawmzvalues ───────────────────────────────────────────────────────────────────────────

"""
    rawmzvalues(msm::AbstractMassScanMatrix; unit::Union{Nothing, Unitful.Units}=nothing)

Return the raw (unit-stripped) m/z values from a mass scan matrix.

This method returns the numeric representation of m/z values, optionally converted 
to a different unit. If the matrix is unitless and a unit conversion is requested, 
an error is thrown.

# Arguments
- `msm`: A subtype of `AbstractMassScanMatrix`.
- `unit`: Desired unit to convert values into before stripping units 
  (`Unitful.Units` subtype or `nothing`, default).

# Returns
A vector of plain numeric m/z values, optionally converted to the requested unit 
and stripped of units.

# Throws
- `ArgumentError` if:
    - The matrix is unitless and a unit is requested.
- `AssertionError` if:
    - The matrix claims to have a unit but `mz_unit` is `nothing` (internal inconsistency).

See also [`AbstractMassScanMatrix`](@ref), [`AbstractMassScan`](@ref), [`mzvalues`](@ref), 
[`mzunit`](@ref), [`mzcount`](@ref), [`intensities`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0]u"Th", [1.0 2.0; 3.0 4.0]);

julia> rawmzvalues(msm) == [100.0, 200.0]
true

julia> msm2 = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> rawmzvalues(msm2) == [100.0, 200.0]
true

julia> rawmzvalues(msm2, unit=u"kTh")
ERROR: ArgumentError: Cannot convert unitless m/z values to a unit
[...]
```
"""
@inline function rawmzvalues(
    msm::AbstractMassScanMatrix{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(msm.mz_values, unit, "m/z values")
end

@inline function rawmzvalues(
    msm::AbstractMassScanMatrix{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip(msm.mz_values, msm.mz_unit, unit)
end

# ── rawretentions ─────────────────────────────────────────────────────────────────────────

"""
    rawretentions(msm::AbstractMassScanMatrix{Nothing}; 
                  unit::Union{Nothing, Unitful.Units}=nothing) -> Vector{<:Real}

Returns the raw (unitless) retention values for all scans in the given mass scan matrix.

If the scan matrix stores unitless retention values, they are returned as-is. If a unit is 
requested via the `unit` keyword, an error is thrown since conversion is not possible.

# Arguments
- `msm`: A subtype of `AbstractMassScanMatrix` whose retentions are unitless.
- `unit`: *(optional)* Must be `nothing` for unitless retention values; otherwise, 
  an error is thrown.

# Returns
A `Vector{<:Real}` of raw, unitless retention values for the scan matrix.

# Throws
- `ArgumentError` if a unit is requested for a scan matrix that does not define a retention 
  unit.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), [`retentions`](@ref), 
[`retentionunit`](@ref), [`scancount`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> rawretentions(msm) == [1.0, 2.0]
true

julia> msm2 = MassScanMatrix([1.0, 2.0], [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> rawretentions(msm2, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retentions to a unit
[...]
```
"""
@inline function rawretentions(
    msm::AbstractMassScanMatrix{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(msm.retentions, unit, "retentions")
end

@inline function rawretentions(
    msm::AbstractMassScanMatrix{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip(msm.retentions, msm.retention_unit, unit)
end

# ── retentions ────────────────────────────────────────────────────────────────────────────

"""
    retentions(msm::AbstractMassScanMatrix; unit::Union{Nothing, Unitful.Units}=nothing)

Return the vector of retention values for a scan object, optionally converted to a 
specified unit.

If the scan’s `retention_unit` is not `nothing`, the returned vector consists of 
`Unitful.Quantity` values. If a `unit` is provided, all values are converted to that unit 
using `uconvert.`.

If `retention_unit` is `nothing`:
- and a `unit` is requested, an error is thrown
- and no `unit` is requested, the raw numeric vector is returned unchanged

This function does **not** assign units to unitless scans — it only performs conversions 
when a unit is already stored in `retention_unit`.

# Arguments
- `msm`: A subtype of `AbstractScan` whose `retention` field is a vector of values.
- `unit`: (Optional) target unit (e.g. `u"s"` or `u"minute"`) — used only if a unit is 
  stored.

# Returns
- A vector of `Quantity` values in the scan’s stored or requested unit, if applicable.
- A plain vector of `Real` values if the scan is unitless and no `unit` is requested.

# Throws
- `ArgumentError` if a `unit` is requested for a unitless scan.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), 
[`rawretentions`](@ref), [`retentionunit`](@ref), [`scancount`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> retentions(msm) == [1.0, 2.0]u"s"
true

julia> retentions(msm; unit=u"ms") == [1000.0, 2000.0]u"ms"
true

julia> msm2 = MassScanMatrix([1.0, 2.0], [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> retentions(msm2) == [1.0, 2.0]  # unitless
true

julia> retentions(msm2; unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retentions to a unit
[...]
```
"""
@inline function retentions(
    msm::AbstractMassScanMatrix{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(msm.retentions, unit, "retentions")
end

@inline function retentions(
    msm::AbstractMassScanMatrix{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(msm.retentions, msm.retention_unit, unit)
end

# ── retentionunit ─────────────────────────────────────────────────────────────────────────

"""
    retentionunit(scan::AbstractMassScanMatrix) -> Union{Unitful.Units, Nothing}

Return the unit associated with the retention value of a scan.

This function retrieves the `retention_unit` field from any scan subtype, which indicates 
the physical unit used along the separation axis (e.g. `u"s"`, `u"minute"`). If no unit 
was specified at construction, returns `nothing`.

# Arguments
- `scan`: A subtype of `AbstractScan`, such as `ChromScan` or `MassScan`.

# Returns
- A `Unitful.Units` subtype representing the retention unit, or `nothing` if unspecified.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), 
[`retentions`](@ref), [`rawretentions`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> retentionunit(msm)
s

julia> msm2 = MassScanMatrix([1.0, 2.0]u"minute", [100.0, 200.0], [1.0 2.0; 3.0 4.0]);

julia> retentionunit(msm2)
minute
```
"""
retentionunit(msm::AbstractMassScanMatrix) = msm.retention_unit

# ── sample ────────────────────────────────────────────────────────────────────────────────

"""
    sample(msm::AbstractMassScanMatrix) -> NamedTuple

Returns the sample metadata associated with the mass scan matrix.

This typically includes information about the analyzed sample, such as sample ID, origin, 
preparation details, or treatment conditions.

# Arguments
- `msm`: A concrete subtype of `AbstractMassScanMatrix`.

# Returns
A `NamedTuple` containing sample-related metadata.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), [`acquisition`](@ref), 
[`instrument`](@ref), [`user`](@ref), [`extras`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0],
                            sample=(ID = "Polistes dominula", locality = "Germany"));

julia> sample(msm)
(ID = "Polistes dominula", locality = "Germany")

julia> sample(msm).ID
"Polistes dominula"
```
"""
sample(msm::AbstractMassScanMatrix) = msm.sample

"""
    scancount(msm::AbstractMassScanMatrix) -> Int

Return the number of individual scans in the mass scan matrix.

This function returns the number of scan entries (e.g. time points or separation 
coordinates) in the mass scan matrix. It reflects the total number of scans stored in the 
intensity matrix.

# Arguments
- `msm::AbstractMassScanMatrix`: A mass scan matrix object, such as `MassScanMatrix`, 
  containing multiple scans and their associated retention times or separation coordinates.

# Returns
- `Int`: The number of scan elements in the matrix.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), [`retentions`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([0.5, 1.0]u"s", [100.0, 150.0, 200.0], [10 20 30; 15 25 35]);

julia> scancount(msm)
2

julia> scancount(msm) == size(intensities(msm), 1)
true
```
"""
scancount(msm::AbstractMassScanMatrix) = length(msm.retentions)

# ── user ──────────────────────────────────────────────────────────────────────────────────

"""
    user(msm::AbstractMassScanMatrix) -> NamedTuple

Returns the user metadata associated with the mass scan matrix.

This typically includes information about the operator, data analyst, or other
user-specific annotations relevant to the acquisition or processing of the data.

# Arguments
- `msm`: A concrete subtype of `AbstractMassScanMatrix`.

# Returns
A `NamedTuple` containing user-related metadata.

See also [`AbstractMassScanMatrix`](@ref), [`MassScanMatrix`](@ref), [`acquisition`](@ref), 
[`instrument`](@ref), [`sample`](@ref), [`extras`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0],
                            user=(operator = "Alice", project = "QC-2025"));

julia> user(msm)
(operator = "Alice", project = "QC-2025")

julia> user(msm).operator
"Alice"
```
"""
user(msm::AbstractMassScanMatrix) = msm.user
