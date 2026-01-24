# ── acquisition ───────────────────────────────────────────────────────────────────────────

"""
    acquisition(msm::AbstractMassScanMatrix) -> NamedTuple

Returns the acquisition metadata associated with the mass scan matrix.

The acquisition metadata typically includes details such as scan mode, method parameters,
instrument configuration, or other context relevant to how the mass spectrometry data was 
acquired.

`msm` is a concrete subtype of `AbstractMassScanMatrix`. Returns a `NamedTuple` containing
acquisition-related metadata.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`extras`](@ref JuChrom.extras(::AbstractMassScanMatrix)), 
[`instrument`](@ref JuChrom.instrument(::AbstractMassScanMatrix)), 
[`level`](@ref JuChrom.level(::AbstractMassScanMatrix)), 
[`sample`](@ref JuChrom.sample(::AbstractMassScanMatrix)), 
[`user`](@ref JuChrom.user(::AbstractMassScanMatrix)).

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

`msm` is a concrete subtype of `AbstractMassScanMatrix`. Returns a `Dict{String, Any}`
containing unstructured metadata.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`acquisition`](@ref JuChrom.acquisition(::AbstractMassScanMatrix)), 
[`instrument`](@ref JuChrom.instrument(::AbstractMassScanMatrix)), 
[`level`](@ref JuChrom.level(::AbstractMassScanMatrix)), 
[`sample`](@ref JuChrom.sample(::AbstractMassScanMatrix)), 
[`user`](@ref JuChrom.user(::AbstractMassScanMatrix)).

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

`msm` is a concrete subtype of `AbstractMassScanMatrix`. Returns a `NamedTuple` containing
instrument-specific metadata.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`acquisition`](@ref JuChrom.acquisition(::AbstractMassScanMatrix)), 
[`extras`](@ref JuChrom.extras(::AbstractMassScanMatrix)), 
[`level`](@ref JuChrom.level(::AbstractMassScanMatrix)), 
[`sample`](@ref JuChrom.sample(::AbstractMassScanMatrix)), 
[`user`](@ref JuChrom.user(::AbstractMassScanMatrix)).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0];
                            instrument=(detector="Orbitrap", manufacturer="Thermo"));

julia> instrument(msm)
(detector = "Orbitrap", manufacturer = "Thermo")

julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [1.0 2.0; 3.0 4.0];
                            instrument=(detector="Orbitrap", manufacturer="Thermo"));


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

`msm` is a concrete subtype of `AbstractMassScanMatrix` containing mass scans. `unit` is
an optional `Unitful.Units` object to convert intensities to; if omitted or `nothing`,
intensities are returned as stored. Returns a matrix of intensity values, either plain
numbers or `AbstractQuantity`s depending on stored data. Throws `ArgumentError` if the scan matrix
stores unitless intensities but a unit conversion is requested.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassScanMatrix)), 
[`mzcount`](@ref JuChrom.mzcount(::AbstractMassScanMatrix)), 
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanMatrix)), 
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})), 
[`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})), 
[`rawretentions`](@ref JuChrom.rawretentions(::AbstractMassScanMatrix{Nothing})), 
[`retentions`](@ref JuChrom.retentions(::AbstractMassScanMatrix{Nothing})), 
[`retentionunit`](@ref JuChrom.retentionunit(::AbstractMassScanMatrix)), 
[`scancount`](@ref JuChrom.scancount(::AbstractMassScanMatrix)).

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
    _handle_unitful_convert(msm.intensities, msm.intensityunit, unit)
end

# ── intensityunit ─────────────────────────────────────────────────────────────────────────

"""
    intensityunit(msm::AbstractMassScanMatrix) -> Union{Unitful.Units, Nothing}

Return the unit associated with the intensity values of a mass scan matrix.

This function retrieves the `intensityunit` field from any mass scan matrix subtype, 
which indicates the physical unit used for the signal intensities (e.g. `u"pA"`). If 
no unit was specified at construction, returns `nothing`.

`msm` is a subtype of `AbstractMassScanMatrix` such as `MassScanMatrix`. Returns a
`Unitful.Units` subtype representing the intensity unit, or `nothing` if unspecified.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`intensities`](@ref JuChrom.intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})).

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
intensityunit(msm::AbstractMassScanMatrix) = msm.intensityunit

# ── level ─────────────────────────────────────────────────────────────────────────────────

"""
    level(msm::AbstractMassScanMatrix) -> Integer

Return the MS level of a mass scan matrix.

The MS level indicates the stage of mass spectrometry at which the scans were acquired
(for example `1` for MS¹, `2` for MS², higher values for MSⁿ). `msm` is a concrete
subtype of `AbstractMassScanMatrix` such as `MassScanMatrix`. Returns an `Integer` value
representing the MS level (guaranteed to be ≥ 1).

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`acquisition`](@ref JuChrom.acquisition(::AbstractMassScanMatrix)), 
[`extras`](@ref JuChrom.extras(::AbstractMassScanMatrix)), 
[`instrument`](@ref JuChrom.instrument(::AbstractMassScanMatrix)), 
[`level`](@ref JuChrom.level(::AbstractMassScanMatrix)), 
[`sample`](@ref JuChrom.sample(::AbstractMassScanMatrix)), 
[`user`](@ref JuChrom.user(::AbstractMassScanMatrix)).

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

This utility function returns the length of the `mzvalues` vector stored in a mass 
spectrometric scan matrix. It reflects the number of unique m/z data points (peaks) 
measured across all scans in the matrix, typically corresponding to the columns of the 
intensity matrix.

`msm` is a mass scan matrix object such as `MassScanMatrix` containing multiple scans and 
their associated m/z values. Returns an `Int` number of m/z values in the scan matrix.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`intensities`](@ref JuChrom.intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanMatrix)), 
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})), 
[`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0]u"s", [100.0, 200.0], [1.0 2.0]);

julia> mzcount(msm)
2

julia> mzcount(msm) == size(intensities(msm), 2)
true
```
"""
mzcount(msm::AbstractMassScanMatrix) = length(msm.mzvalues)

# ── mzindex ───────────────────────────────────────────────────────────────────────────────

"""
    mzindex(msm::AbstractMassScanMatrix, mzvalue::Union{Real, Unitful.AbstractQuantity}) -> Int

Return the index of the m/z value in the mass scan matrix.

This function searches the `mzvalues` vector stored in the mass scan matrix and returns the
1-based index of the requested `mzvalue`. The return value is `nothing` if the value is not
present. For unitless m/z vectors, `mzvalue` must be a `Real`. For unitful m/z vectors,
`mzvalue` must be a `Unitful.AbstractQuantity` with compatible units.

`msm` is a mass scan matrix object such as `MassScanMatrix`. The index corresponds to the
column position in the intensity matrix.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`intensities`](@ref JuChrom.intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanMatrix)), 
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})), 
[`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0]u"s", [100.0, 200.0], [1.0 2.0]);

julia> mzindex(msm, 100.0) == 1
true

julia> mzindex.(msm, [100.0, 200.0]) == [1, 2]
true
```
"""
mzindex(msm::AbstractMassScanMatrix{<:Any, Nothing}, mzvalue::Real) = 
    findfirst(==(mzvalue), mzvalues(msm))    

mzindex(msm::AbstractMassScanMatrix{<:Any, <:Unitful.Units}, mzvalue::Unitful.AbstractQuantity) = 
    findfirst(==(mzvalue), mzvalues(msm)) 

# ── mzunit ────────────────────────────────────────────────────────────────────────────────

"""
    mzunit(msm::AbstractMassScanMatrix) -> Union{Unitful.Units, Nothing}

Return the unit associated with the mass-to-charge (m/z) values of a mass scan matrix.

Typically, this is `nothing`, since m/z values are conventionally unitless. If a unit was
explicitly provided during matrix construction, it will be returned.

`msm` is a subtype of `AbstractMassScanMatrix` such as `MassScanMatrix`. Returns a
`Unitful.Units` subtype representing the m/z unit, or `nothing` if unspecified.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`mzcount`](@ref JuChrom.mzcount(::AbstractMassScanMatrix)), 
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})), 
[`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})).

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
mzunit(msm::AbstractMassScanMatrix) = msm.mzunit

# ── mzvalues ──────────────────────────────────────────────────────────────────────────────

"""
    mzvalues(msm::AbstractMassScanMatrix; unit::Union{Nothing, Unitful.Units}=nothing)

Return the mass-to-charge (m/z) values from a mass scan matrix.

This method supports both unitless and unitful matrices. If the matrix has an associated m/z 
unit, values can be converted to a different unit using the optional `unit` argument. 
If the matrix is unitless and a unit is requested, an error is thrown.

`msm` is a subtype of `AbstractMassScanMatrix`. `unit` is an optional target unit.
Returns a vector of m/z values, converted to the requested unit when specified, otherwise
returned in stored form. Throws `ArgumentError` if the matrix is unitless and a unit is
requested, and `AssertionError` if the matrix claims to have a unit but `mzunit` is
`nothing`.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`intensities`](@ref JuChrom.intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`mzcount`](@ref JuChrom.mzcount(::AbstractMassScanMatrix)), 
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanMatrix)), 
[`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})).

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
    _handle_unitless(msm.mzvalues, unit, "m/z values")
end

@inline function mzvalues(
    msm::AbstractMassScanMatrix{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert(msm.mzvalues, msm.mzunit, unit)
end

# ── rawintensities ────────────────────────────────────────────────────────────────────────

"""
    rawintensities(msm::AbstractMassScanMatrix;
                   unit::Union{Nothing, Unitful.Units}=nothing) -> Matrix

Return the numeric (unitless) intensity values from all scans in a mass scan matrix.

If the scan matrix’s `intensityunit` is not `nothing`, the intensity values are optionally 
converted to the user-specified `unit` using `uconvert`, then stripped of units. If no 
`unit` is specified, the stored unit is used for stripping.

If the matrix does not define intensity units, requesting a unit conversion will throw an 
`ArgumentError`.

`msm` is a concrete subtype of `AbstractMassScanMatrix` such as `MassScanMatrix`. `unit`
is an optional `Unitful.Units` object to convert raw intensities before stripping; it must
be `nothing` for unitless matrices. Returns a matrix of unitless intensity values with
optional conversion (rows are scans, columns are m/z values). Throws `ArgumentError` if a
unit is requested but the scan matrix has no defined intensity unit.

[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`intensities`](@ref JuChrom.intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassScanMatrix)), 
[`mzcount`](@ref JuChrom.mzcount(::AbstractMassScanMatrix)), 
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanMatrix)), 
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})), 
[`rawmzvalues`](@ref JuChrom.rawmzvalues(::AbstractMassScanMatrix{<:Any, Nothing})), 
[`rawretentions`](@ref JuChrom.rawretentions(::AbstractMassScanMatrix{Nothing})), 
[`retentions`](@ref JuChrom.retentions(::AbstractMassScanMatrix{Nothing})), 
[`retentionunit`](@ref JuChrom.retentionunit(::AbstractMassScanMatrix)), 
[`scancount`](@ref JuChrom.scancount(::AbstractMassScanMatrix)).

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
    _handle_unitful_strip(msm.intensities, msm.intensityunit, unit)
end

# ── rawmzvalues ───────────────────────────────────────────────────────────────────────────

"""
    rawmzvalues(msm::AbstractMassScanMatrix; unit::Union{Nothing, Unitful.Units}=nothing)

Return the raw (unit-stripped) m/z values from a mass scan matrix.

This method returns the numeric representation of m/z values, optionally converted 
to a different unit. If the matrix is unitless and a unit conversion is requested, 
an error is thrown.

`msm` is a subtype of `AbstractMassScanMatrix`. `unit` is an optional target unit.
Returns a vector of numeric m/z values, optionally converted and stripped of units. Throws
`ArgumentError` if the matrix is unitless and a unit is requested, and `AssertionError` if
the matrix claims to have a unit but `mzunit` is `nothing`.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`intensities`](@ref JuChrom.intensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})), 
[`mzcount`](@ref JuChrom.mzcount(::AbstractMassScanMatrix)), 
[`mzunit`](@ref JuChrom.mzunit(::AbstractMassScanMatrix)), 
[`mzvalues`](@ref JuChrom.mzvalues(::AbstractMassScanMatrix{<:Any, Nothing})), 
[`rawintensities`](@ref JuChrom.rawintensities(::AbstractMassScanMatrix{<:Any, <:Any, Nothing})).

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
    _handle_unitless(msm.mzvalues, unit, "m/z values")
end

@inline function rawmzvalues(
    msm::AbstractMassScanMatrix{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip(msm.mzvalues, msm.mzunit, unit)
end

# ── rawretentions ─────────────────────────────────────────────────────────────────────────

"""
    rawretentions(msm::AbstractMassScanMatrix{Nothing}; 
                  unit::Union{Nothing, Unitful.Units}=nothing) -> Vector{<:Real}

Returns the raw (unitless) retention values for all scans in the given mass scan matrix.

If the scan matrix stores unitless retention values, they are returned as-is. If a unit is 
requested via the `unit` keyword, an error is thrown since conversion is not possible.

`msm` is a subtype of `AbstractMassScanMatrix` with unitless retentions. `unit` must be
`nothing` for unitless retention values. Returns a vector of raw, unitless retention
values. Throws `ArgumentError` if a unit is requested for a scan matrix that does not
define a retention unit.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`retentions`](@ref JuChrom.retentions(::AbstractMassScanMatrix{Nothing})), 
[`retentionunit`](@ref JuChrom.retentionunit(::AbstractMassScanMatrix)), 
[`scancount`](@ref JuChrom.scancount(::AbstractMassScanMatrix)).

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
    _handle_unitful_strip(msm.retentions, msm.retentionunit, unit)
end

# ── retentions ────────────────────────────────────────────────────────────────────────────

"""
    retentions(msm::AbstractMassScanMatrix; unit::Union{Nothing, Unitful.Units}=nothing)

Return the vector of retention values from a mass scan matrix, optionally converted to a
specified unit.

If the matrix stores a `retentionunit`, values are returned as `Unitful.AbstractQuantity`s and
any requested `unit` is applied via `uconvert`. If the matrix is unitless, the raw numeric
vector is returned unless a unit conversion is requested, which throws `ArgumentError`.
This function does not assign units to unitless matrices; it only converts when a unit is
stored.

`msm` is a subtype of `AbstractMassScanMatrix` whose retentions are stored as a vector.
`unit` is an optional target unit (for example `u"s"` or `u"minute"`) used only when a
unit is stored.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`rawretentions`](@ref JuChrom.rawretentions(::AbstractMassScanMatrix{Nothing})), 
[`retentionunit`](@ref JuChrom.retentionunit(::AbstractMassScanMatrix)), 
[`scancount`](@ref JuChrom.scancount(::AbstractMassScanMatrix)).

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
    _handle_unitful_convert(msm.retentions, msm.retentionunit, unit)
end

# ── retentionunit ─────────────────────────────────────────────────────────────────────────

"""
    retentionunit(scan::AbstractMassScanMatrix) -> Union{Unitful.Units, Nothing}

Return the unit associated with the retention value of a scan.

This function retrieves the `retentionunit` field from any scan subtype, which indicates 
the physical unit used along the separation axis (e.g. `u"s"`, `u"minute"`). If no unit 
was specified at construction, returns `nothing`.

`msm` is a subtype of `AbstractMassScanMatrix`. Returns a `Unitful.Units` subtype
representing the retention unit, or `nothing` if unspecified.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`rawretentions`](@ref JuChrom.rawretentions(::AbstractMassScanMatrix{Nothing})), 
[`retentions`](@ref JuChrom.retentions(::AbstractMassScanMatrix{Nothing})).

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
retentionunit(msm::AbstractMassScanMatrix) = msm.retentionunit

# ── sample ────────────────────────────────────────────────────────────────────────────────

"""
    sample(msm::AbstractMassScanMatrix) -> NamedTuple

Returns the sample metadata associated with the mass scan matrix.

This typically includes information about the analyzed sample, such as sample ID, origin, 
preparation details, or treatment conditions.

`msm` is a concrete subtype of `AbstractMassScanMatrix`. Returns a `NamedTuple` containing
sample-related metadata.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`acquisition`](@ref JuChrom.acquisition(::AbstractMassScanMatrix)), 
[`extras`](@ref JuChrom.extras(::AbstractMassScanMatrix)), 
[`instrument`](@ref JuChrom.instrument(::AbstractMassScanMatrix)), 
[`level`](@ref JuChrom.level(::AbstractMassScanMatrix)), 
[`user`](@ref JuChrom.user(::AbstractMassScanMatrix)).

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

`msm` is a mass scan matrix object such as `MassScanMatrix` containing multiple scans and
their associated retention times or separation coordinates. Returns an `Int` number of
scan elements in the matrix.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`rawretentions`](@ref JuChrom.rawretentions(::AbstractMassScanMatrix{Nothing})), 
[`retentions`](@ref JuChrom.retentions(::AbstractMassScanMatrix{Nothing})).

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

`msm` is a concrete subtype of `AbstractMassScanMatrix`. Returns a `NamedTuple` containing
user-related metadata.

See also
[`AbstractMassScanMatrix`](@ref JuChrom.AbstractMassScanMatrix), 
[`MassScanMatrix`](@ref JuChrom.MassScanMatrix), 
[`acquisition`](@ref JuChrom.acquisition(::AbstractMassScanMatrix)), 
[`extras`](@ref JuChrom.extras(::AbstractMassScanMatrix)), 
[`instrument`](@ref JuChrom.instrument(::AbstractMassScanMatrix)), 
[`level`](@ref JuChrom.level(::AbstractMassScanMatrix)), 
[`sample`](@ref JuChrom.sample(::AbstractMassScanMatrix)).

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
