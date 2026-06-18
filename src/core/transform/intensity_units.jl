# ── dwellnormalize ────────────────────────────────────────────────────────────────────────

"""
    dwellnormalize(msm::MassScanMatrix, dwell::AbstractVector{<:Real}, unit::Unitful.Units)

Return a copy of `msm` whose intensity values are divided by m/z-specific dwell intervals.

The input matrix is interpreted as containing signal accumulated over m/z dwell intervals.
The values may be raw, vendor-scaled, or otherwise transformed intensities. The relevant
assumption is that the stored values scale with the duration of the dwell interval. The
input intensity unit must be absent. `dwell` provides one dwell interval per m/z value, in
the unit specified by `unit`, and must therefore have length `mzcount(msm)`. Each intensity
column is divided by the corresponding dwell value. The sum of the dwell intervals must
not exceed the shortest scan interval implied by the retention axis, allowing for floating
point roundoff. The returned matrix uses `inverse(unit)` as its intensity unit. Retention
coordinates, m/z values, MS level, and metadata are preserved.

Throws `ArgumentError` if `intensityunit(msm)` is not `nothing`, if the retention axis has
no unit or fewer than two scans, if any dwell value is non-finite or non-positive, or if
the total dwell duration exceeds the shortest scan interval. Throws `DimensionMismatch` if
`dwell` does not have length `mzcount(msm)`.

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [10.0 20.0; 30.0 40.0]);

julia> normalized = dwellnormalize(msm, [0.4, 0.6], u"s");

julia> rawintensities(normalized)
2×2 Matrix{Float64}:
 25.0  33.3333
 75.0  66.6667

julia> intensityunit(normalized)
s^-1
```
"""
function dwellnormalize(
    msm::MassScanMatrix,
    dwell::AbstractVector{<:Real},
    unit::Unitful.Units
)
    inputunit = intensityunit(msm)
    isnothing(inputunit) ||
        throw(ArgumentError("MassScanMatrix intensityunit must be nothing."))

    length(dwell) == mzcount(msm) ||
        throw(DimensionMismatch("dwell must have length mzcount(msm)."))

    all(isfinite, dwell) ||
        throw(ArgumentError("all dwell values must be finite."))

    all(t -> t > zero(t), dwell) ||
        throw(ArgumentError("all dwell values must be positive."))

    validate_dwell_fits_scan_interval(msm, dwell, unit)

    normalized = rawintensities(msm) * Diagonal(inv.(dwell))
    normalizedunit = inverse(unit)

    MassScanMatrix(
        copy(rawretentions(msm)),
        retentionunit(msm),
        copy(rawmzvalues(msm)),
        mzunit(msm),
        normalized,
        normalizedunit;
        level=level(msm),
        instrument=deepcopy(instrument(msm)),
        acquisition=deepcopy(acquisition(msm)),
        user=deepcopy(user(msm)),
        sample=deepcopy(sample(msm)),
        extras=deepcopy(extras(msm)),
    )
end

function validate_dwell_fits_scan_interval(
    msm::MassScanMatrix,
    dwell::AbstractVector{<:Real},
    unit::Unitful.Units
)
    rtunit = retentionunit(msm)
    !isnothing(rtunit) ||
        throw(ArgumentError("retentionunit is required to validate dwell intervals."))
    scancount(msm) > 1 ||
        throw(ArgumentError("at least two scans are required to validate dwell intervals."))

    scanintervals = diff(rawretentions(msm))
    all(isfinite, scanintervals) ||
        throw(ArgumentError("all scan intervals must be finite."))
    all(interval -> interval > zero(interval), scanintervals) ||
        throw(ArgumentError("all scan intervals must be positive."))

    totaldwell = try
        ustrip(rtunit, sum(dwell) * unit)
    catch
        throw(ArgumentError(
            "dwell unit must be compatible with retentionunit(msm)."))
    end
    shortestinterval = minimum(scanintervals)
    tolerance = sqrt(eps(Float64)) * max(abs(totaldwell), abs(shortestinterval), 1.0)
    totaldwell ≤ shortestinterval + tolerance ||
        throw(ArgumentError(
            "sum(dwell) must not exceed the shortest scan interval."))

    nothing
end

"""
    dwellnormalize(msm::MassScanMatrix, dwell::AbstractVector{<:AbstractQuantity})

Return a copy of `msm` whose intensity values are divided by unitful m/z-specific dwell
intervals.

This method accepts dwell intervals with Unitful units and delegates to
`dwellnormalize(msm, dwellvalues, dwellunit)` after checking and stripping the units. The
unitful `dwell` vector must contain one dwell interval per m/z value and use a consistent
unit. The returned matrix stores intensity values normalized by those dwell intervals and
uses the reciprocal dwell unit as its intensity unit.

Throws `ArgumentError` if units in `dwell` are inconsistent, if `intensityunit(msm)` is
not `nothing`, if any dwell value is non-finite or non-positive, or if the total dwell
duration exceeds the shortest scan interval. Throws `DimensionMismatch` if `dwell` does
not have length `mzcount(msm)`.

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [10.0 20.0; 30.0 40.0]);

julia> normalized = dwellnormalize(msm, [0.4, 0.6]u"s");

julia> rawintensities(normalized)
2×2 Matrix{Float64}:
 25.0  33.3333
 75.0  66.6667

julia> intensityunit(normalized)
s^-1
```
"""
function dwellnormalize(
    msm::MassScanMatrix,
    dwell::AbstractVector{<:AbstractQuantity{<:Real}},
)
    dwellvalues, dwellunit = strip_units_checked(dwell, "dwell")
    dwellnormalize(msm, dwellvalues, dwellunit)
end

"""
    dwellnormalize(msm::MassScanMatrix)

Return a copy of `msm` whose intensity values are normalized by an inferred uniform dwell
interval.

This convenience method assumes that each scan interval is divided evenly among all m/z
values in the matrix. The scan interval is estimated as the mean spacing between adjacent
retention coordinates, and the dwell interval is computed as that mean scan interval
divided by `mzcount(msm)`. This method is intended for data where the stored intensities
represent signal accumulated over dwell intervals and the scan schedule is known to be
uniform across m/z values. If the input intensity unit is absent, the returned matrix uses
the reciprocal retention unit as its intensity unit.

Throws `ArgumentError` if `msm` has no retention unit, if it has fewer than two scans, if
`intensityunit(msm)` is not `nothing`, or if the inferred scan intervals are non-finite or
non-positive.

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [10.0 20.0; 30.0 40.0]);

julia> normalized = dwellnormalize(msm);

julia> rawintensities(normalized)
2×2 Matrix{Float64}:
 20.0  40.0
 60.0  80.0

julia> intensityunit(normalized)
s^-1
```
"""
function dwellnormalize(msm::MassScanMatrix)
    !isnothing(retentionunit(msm)) ||
        throw(ArgumentError("retentionunit is required to infer dwell intervals."))
    scancount(msm) > 1 ||
        throw(ArgumentError("at least two scans are required to infer dwell intervals."))

    scanintervals = diff(rawretentions(msm))
    all(isfinite, scanintervals) ||
        throw(ArgumentError("all inferred scan intervals must be finite."))
    all(t -> t > zero(t), scanintervals) ||
        throw(ArgumentError("all inferred scan intervals must be positive."))

    mean_scaninterval = mean(scanintervals)
    dwell = fill(mean_scaninterval / mzcount(msm), mzcount(msm))
    dwellnormalize(msm, dwell, retentionunit(msm))
end

"""
    dwellnormalize(vmsm::VarianceMassScanMatrix, dwell::AbstractVector{<:Real}, unit::Unitful.Units)

Return a copy of `vmsm` whose intensity values and variances are normalized by
m/z-specific dwell intervals.

The parent mass-scan matrix is normalized with
`dwellnormalize(parent(vmsm), dwell, unit)`, which also copies the retention axis, m/z
axis, units, and metadata. Variances are divided by the squared dwell interval of the
corresponding m/z column, so they remain compatible with the normalized intensities. The
returned object is a new `VarianceMassScanMatrix`.
"""
function dwellnormalize(
    vmsm::VarianceMassScanMatrix,
    dwell::AbstractVector{<:Real},
    unit::Unitful.Units
)
    normalizedmsm = dwellnormalize(parent(vmsm), dwell, unit)
    normalizedvariances = rawvariances(vmsm) * Diagonal(abs2.(inv.(dwell)))
    VarianceMassScanMatrix(normalizedmsm, normalizedvariances)
end

"""
    dwellnormalize(vmsm::VarianceMassScanMatrix, dwell::AbstractVector{<:AbstractQuantity})

Return a copy of `vmsm` whose intensity values and variances are normalized by unitful
m/z-specific dwell intervals.
"""
function dwellnormalize(
    vmsm::VarianceMassScanMatrix,
    dwell::AbstractVector{<:AbstractQuantity{<:Real}}
)
    dwellvalues, dwellunit = strip_units_checked(dwell, "dwell")
    dwellnormalize(vmsm, dwellvalues, dwellunit)
end

"""
    dwellnormalize(vmsm::VarianceMassScanMatrix)

Return a copy of `vmsm` whose intensity values and variances are normalized by an inferred
uniform dwell interval.
"""
function dwellnormalize(vmsm::VarianceMassScanMatrix)
    !isnothing(retentionunit(vmsm)) ||
        throw(ArgumentError("retentionunit is required to infer dwell intervals."))
    scancount(vmsm) > 1 ||
        throw(ArgumentError("at least two scans are required to infer dwell intervals."))

    scanintervals = diff(rawretentions(vmsm))
    all(isfinite, scanintervals) ||
        throw(ArgumentError("all inferred scan intervals must be finite."))
    all(t -> t > zero(t), scanintervals) ||
        throw(ArgumentError("all inferred scan intervals must be positive."))

    mean_scaninterval = mean(scanintervals)
    dwell = fill(mean_scaninterval / mzcount(vmsm), mzcount(vmsm))
    dwellnormalize(vmsm, dwell, retentionunit(vmsm))
end

# ── withintensityunit ─────────────────────────────────────────────────────────────────────

"""
    withintensityunit(msm::MassScanMatrix, unit::Unitful.Units)

Return a copy of `msm` with the supplied intensity unit attached without changing the
stored intensity values.

This function is intended for cases where the numeric intensity values already have a known
physical unit, but the matrix currently stores `intensityunit(msm) ≡ nothing`. It only
changes the intensity unit metadata. Retention coordinates, m/z values, intensity values,
MS level, and metadata are preserved.

Throws `ArgumentError` if `msm` already has an intensity unit.

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0], [10.0 20.0; 30.0 40.0]);

julia> annotated = withintensityunit(msm, u"pA");

julia> rawintensities(annotated) == rawintensities(msm)
true

julia> intensityunit(annotated)
pA
```
"""
function withintensityunit(msm::MassScanMatrix, unit::Unitful.Units)
    isnothing(intensityunit(msm)) ||
        throw(ArgumentError("MassScanMatrix intensityunit must be nothing."))

    MassScanMatrix(
        copy(rawretentions(msm)),
        retentionunit(msm),
        copy(rawmzvalues(msm)),
        mzunit(msm),
        copy(rawintensities(msm)),
        unit;
        level=level(msm),
        instrument=deepcopy(instrument(msm)),
        acquisition=deepcopy(acquisition(msm)),
        user=deepcopy(user(msm)),
        sample=deepcopy(sample(msm)),
        extras=deepcopy(extras(msm)),
    )
end
