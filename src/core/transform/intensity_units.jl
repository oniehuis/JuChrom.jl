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

    dwell_normalize_unchecked(msm, dwell, unit)
end

function dwell_normalize_unchecked(
    msm::MassScanMatrix,
    dwell::AbstractVector{<:Real},
    unit::Unitful.Units
)
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
        ustrip(rtunit, sum(Float64, dwell) * unit)
    catch
        throw(ArgumentError(
            "dwell unit must be compatible with retentionunit(msm)."))
    end
    shortestinterval = Float64(minimum(scanintervals))
    tolerance = max(
        length(dwell) * eps(Float64) * max(abs(totaldwell), 
        abs(shortestinterval)),
        eps(Float64)
    )
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
    dwellnormalize(msm::MassScanMatrix, dwell::Real, unit::Unitful.Units)

Return a copy of `msm` whose intensity values are divided by a single dwell interval
shared by all m/z channels.

This scalar form represents simultaneous acquisition with a known effective dwell time.
The dwell interval must be finite, positive, and no larger than the shortest scan interval
on the retention axis. The returned matrix uses `inverse(unit)` as its intensity unit.
"""
function dwellnormalize(
    msm::MassScanMatrix,
    dwell::Real,
    unit::Unitful.Units
)
    inputunit = intensityunit(msm)
    isnothing(inputunit) ||
        throw(ArgumentError("MassScanMatrix intensityunit must be nothing."))
    validate_scalar_dwell(msm, dwell, unit)

    dwell_normalize_unchecked(msm, fill(Float64(dwell), mzcount(msm)), unit)
end

"""
    dwellnormalize(msm::MassScanMatrix, dwell::Unitful.Quantity)

Return a copy of `msm` whose intensity values are divided by a single unitful dwell
interval shared by all m/z channels.
"""
function dwellnormalize(
    msm::MassScanMatrix,
    dwell::AbstractQuantity{<:Real}
)
    dwellnormalize(msm, Unitful.ustrip(dwell), Unitful.unit(dwell))
end

"""
    dwellnormalize(msm::MassScanMatrix; acquisition=:sequential)

Return a copy of `msm` whose intensity values are normalized by an inferred uniform dwell
interval.

This convenience method infers a uniform dwell interval from the shortest adjacent
retention spacing. For `acquisition=:sequential`, that shortest scan interval is divided
by `mzcount(msm)`. For `acquisition=:simultaneous`, the full shortest scan interval is
used for every m/z value. The returned matrix uses the reciprocal retention unit as its
intensity unit.

Throws `ArgumentError` if `msm` has no retention unit, if it has fewer than two scans, if
`intensityunit(msm)` is not `nothing`, or if the inferred scan intervals are non-finite or
non-positive. `acquisition` must be `:sequential` or `:simultaneous`.

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
function dwellnormalize(msm::MassScanMatrix; acquisition::Symbol=:sequential)
    dwell = infer_uniform_dwell_intervals(msm, acquisition)
    dwell_normalize_unchecked(msm, dwell, retentionunit(msm))
end

function infer_uniform_dwell_intervals(
    msm::Union{MassScanMatrix, VarianceMassScanMatrix},
    acquisition::Symbol
)
    validate_dwell_acquisition(acquisition)
    !isnothing(retentionunit(msm)) ||
        throw(ArgumentError("retentionunit is required to infer dwell intervals."))
    scancount(msm) > 1 ||
        throw(ArgumentError("at least two scans are required to infer dwell intervals."))

    scanintervals = diff(rawretentions(msm))
    all(isfinite, scanintervals) ||
        throw(ArgumentError("all inferred scan intervals must be finite."))
    all(t -> t > zero(t), scanintervals) ||
        throw(ArgumentError("all inferred scan intervals must be positive."))

    shortestscaninterval = Float64(minimum(scanintervals))
    dwell = acquisition ≡ :sequential ?
        shortestscaninterval / mzcount(msm) :
        shortestscaninterval
    fill(dwell, mzcount(msm))
end

function validate_dwell_acquisition(acquisition::Symbol)
    acquisition in (:sequential, :simultaneous) || throw(ArgumentError(
        "acquisition must be :sequential or :simultaneous"))

    nothing
end

function validate_scalar_dwell(
    msm::MassScanMatrix,
    dwell::Real,
    unit::Unitful.Units
)
    isfinite(dwell) || throw(ArgumentError("dwell must be finite."))
    dwell > zero(dwell) || throw(ArgumentError("dwell must be positive."))

    rtunit = retentionunit(msm)
    !isnothing(rtunit) ||
        throw(ArgumentError("retentionunit is required to validate dwell interval."))
    scancount(msm) > 1 ||
        throw(ArgumentError("at least two scans are required to validate dwell interval."))

    scanintervals = diff(rawretentions(msm))
    all(isfinite, scanintervals) ||
        throw(ArgumentError("all scan intervals must be finite."))
    all(interval -> interval > zero(interval), scanintervals) ||
        throw(ArgumentError("all scan intervals must be positive."))

    dwellretention = try
        ustrip(rtunit, Float64(dwell) * unit)
    catch
        throw(ArgumentError(
            "dwell unit must be compatible with retentionunit(msm)."))
    end
    shortestinterval = Float64(minimum(scanintervals))
    tolerance = max(
        eps(Float64) * max(abs(dwellretention), 
        abs(shortestinterval)),
        eps(Float64)
    )
    dwellretention ≤ shortestinterval + tolerance ||
        throw(ArgumentError("dwell must not exceed the shortest scan interval."))

    nothing
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
    dwell_normalize_variances(vmsm, normalizedmsm, dwell)
end

function dwell_normalize_variances(
    vmsm::VarianceMassScanMatrix,
    normalizedmsm::MassScanMatrix,
    dwell::AbstractVector{<:Real}
)
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
    dwellnormalize(vmsm::VarianceMassScanMatrix, dwell::Real, unit::Unitful.Units)

Return a copy of `vmsm` whose intensity values and variances are normalized by a single
dwell interval shared by all m/z channels.
"""
function dwellnormalize(
    vmsm::VarianceMassScanMatrix,
    dwell::Real,
    unit::Unitful.Units
)
    normalizedmsm = dwellnormalize(parent(vmsm), dwell, unit)
    dwell_normalize_variances(
        vmsm,
        normalizedmsm,
        fill(Float64(dwell), mzcount(vmsm))
    )
end

"""
    dwellnormalize(vmsm::VarianceMassScanMatrix, dwell::Unitful.Quantity)

Return a copy of `vmsm` whose intensity values and variances are normalized by a single
unitful dwell interval shared by all m/z channels.
"""
function dwellnormalize(
    vmsm::VarianceMassScanMatrix,
    dwell::AbstractQuantity{<:Real}
)
    dwellnormalize(vmsm, Unitful.ustrip(dwell), Unitful.unit(dwell))
end

"""
    dwellnormalize(vmsm::VarianceMassScanMatrix; acquisition=:sequential)

Return a copy of `vmsm` whose intensity values and variances are normalized by an inferred
uniform dwell interval.
"""
function dwellnormalize(
    vmsm::VarianceMassScanMatrix;
    acquisition::Symbol=:sequential
)
    dwell = infer_uniform_dwell_intervals(vmsm, acquisition)
    normalizedmsm = dwell_normalize_unchecked(parent(vmsm), dwell, retentionunit(vmsm))
    dwell_normalize_variances(vmsm, normalizedmsm, dwell)
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

"""
    withintensityunit(vmsm::VarianceMassScanMatrix, unit::Unitful.Units)

Return a copy of `vmsm` with the supplied intensity unit attached without changing the
stored intensity or variance values.

The parent mass-scan matrix is annotated with `unit`. The variance matrix is copied
unchanged and its unit is set to `unit^2`, matching the annotated intensity unit.

Throws `ArgumentError` if `vmsm` already has an intensity unit or variance unit.
"""
function withintensityunit(vmsm::VarianceMassScanMatrix, unit::Unitful.Units)
    isnothing(intensityunit(vmsm)) ||
        throw(ArgumentError("VarianceMassScanMatrix intensityunit must be nothing."))
    isnothing(varianceunit(vmsm)) ||
        throw(ArgumentError("VarianceMassScanMatrix varianceunit must be nothing."))

    annotatedmsm = withintensityunit(parent(vmsm), unit)
    VarianceMassScanMatrix(annotatedmsm, copy(rawvariances(vmsm)), unit^2)
end
