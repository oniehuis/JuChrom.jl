# ── dwellnormalize ────────────────────────────────────────────────────────────────────────

"""
    dwellnormalize(msm::Union{MassScanMatrix, VarianceMassScanMatrix}, dwell)

Return a copy with intensities divided by explicitly supplied dwell durations.

`dwell` is either a single unitful duration shared by every m/z column, or a tuple/vector
of `mz => duration` pairs or a dictionary associating m/z values with unitful durations.
Every duration must be real-valued, finite, positive, and have units of time. Durations
are converted to seconds; output intensities have unit `s^-1`, independently of the
units or ordering of the specification. A shared duration does not imply simultaneous
acquisition.

For m/z-specific durations, every matrix column, including zero-filled columns, needs an
entry. Additional ions are allowed. Matching is exact, without rounding or tolerance.
Keys must be finite and positive, and either all plain real numbers or all real-valued
quantities. Plain keys use the matrix's m/z unit, if any. Unitful keys are converted to
that unit before matching and require a unitful matrix m/z axis. Compatible key units
may differ. Duplicate keys are rejected after unit conversion, including duplicates for
ions absent from the matrix. All entries, including additional ions, are validated.

The input intensity unit must be absent, and its values must represent signal proportional
to dwell duration. If the stored values already represent a rate, use
[`withintensityunit`](@ref) to attach its known unit without rescaling. This function does
not infer dwell from recorded ions or retention spacing, and does not check durations
against scan intervals. Single scans and transformed or unitless retention axes work.

For `VarianceMassScanMatrix`, variances are divided by squared dwell durations and have
unit `s^-2`. Dense/sparse storage is preserved. Axes, MS level, and metadata are copied.

Throws `ArgumentError` for invalid durations or keys, duplicate or missing keys, malformed
pair collections, or an input with an intensity unit. Positional duration vectors,
unitless durations, automatic inference, and the `acquisition` keyword are not supported.

See also
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassScanMatrix)),
[`varianceunit`](@ref JuChrom.varianceunit(::AbstractVarianceMassScanMatrix)),
[`withintensityunit`](@ref).

# Examples
```jldoctest
julia> msm = MassScanMatrix([1.0, 2.0]u"s", [57, 71], [10.0 20.0; 30.0 40.0]);

julia> shared = dwellnormalize(msm, 10u"ms");

julia> rawintensities(shared) == [1000 2000; 3000 4000]
true

julia> normalized = dwellnormalize(msm, (71 => 20u"ms", 57 => 0.01u"s"));

julia> rawintensities(normalized) == [1000 1000; 3000 2000]
true

julia> intensityunit(normalized) == u"s^-1"
true

julia> rates = dwellnormalize(msm, Dict(57 => 10u"ms", 71 => 20u"ms", 85 => 5u"ms"));

julia> rawintensities(rates) == rawintensities(normalized)
true

julia> vmsm = VarianceMassScanMatrix(msm, [1.0 4.0; 9.0 16.0]);

julia> normalizedvar = dwellnormalize(vmsm, [57 => 10u"ms", 71 => 20u"ms"]);

julia> rawvariances(normalizedvar) == [10000 10000; 90000 40000]
true

julia> varianceunit(normalizedvar) == u"s^-2"
true
```
"""
function dwellnormalize(
    msm::Union{MassScanMatrix, VarianceMassScanMatrix},
    dwell::Union{AbstractQuantity, AbstractVector, Tuple, AbstractDict}
)
    isnothing(intensityunit(msm)) || throw(ArgumentError(
        "dwellnormalize requires intensities without an attached unit."))
    durations = resolve_dwelltimes(msm, dwell)
    dwell_normalize_unchecked(msm, durations)
end

# Validate a duration and return its numeric value in the internal time unit (seconds).
function validate_dwelltime(duration)
    duration isa AbstractQuantity{<:Real} || throw(ArgumentError(
        "Each dwell duration must be a real-valued unitful time quantity."))
    Unitful.dimension(duration) == Unitful.dimension(u"s") || throw(ArgumentError(
        "Dwell durations must have units of time."))
    seconds = ustrip(u"s", duration)
    isfinite(seconds) && seconds > zero(seconds) || throw(ArgumentError(
        "Dwell durations must be finite and positive in seconds."))
    # Integer/rational inverse durations can overflow when squared for variance scaling.
    float(seconds)
end

function resolve_dwelltimes(msm, dwell::AbstractQuantity)
    fill(validate_dwelltime(dwell), mzcount(msm))
end

function resolve_dwelltimes(msm, dwell::Union{AbstractVector, Tuple, AbstractDict})
    isempty(dwell) && throw(ArgumentError("The m/z-to-dwell specification must not be empty."))
    entries = collect(dwell)
    all(entry -> entry isa Pair, entries) || throw(ArgumentError(
        "Supply mz => duration pairs, not a positional duration vector."))
    keys = dwell_mz_keys(entries, mzunit(msm))
    by_mz = Dict{Real, Real}()
    for (mz, entry) in zip(keys, entries)
        haskey(by_mz, mz) && throw(ArgumentError("Duplicate dwell specification for m/z $mz."))
        by_mz[mz] = validate_dwelltime(last(entry))
    end
    durations = map(rawmzvalues(msm)) do mz
        haskey(by_mz, mz) || throw(ArgumentError("Missing dwell duration for m/z $mz."))
        by_mz[mz]
    end
    # A concrete common duration type keeps mixed-precision matrix arithmetic numeric.
    durationtype = mapreduce(typeof, promote_type, durations)
    convert(Vector{durationtype}, durations)
end

function dwell_mz_keys(entries, mz_unit)
    keys = map(first, entries)
    if all(x -> x isa Real, keys)
        values = keys
    elseif all(x -> x isa AbstractQuantity{<:Real}, keys)
        isnothing(mz_unit) && throw(ArgumentError(
            "Unitful dwell m/z keys require a matrix with an m/z unit."))
        all(x -> Unitful.dimension(x) == Unitful.dimension(mz_unit), keys) ||
            throw(ArgumentError("Dwell m/z key units must be compatible with the matrix's m/z unit."))
        values = map(x -> ustrip(mz_unit, x), keys)
    else
        throw(ArgumentError(
            "Dwell m/z keys must be either all real numbers or all real-valued quantities."))
    end
    all(x -> isfinite(x) && x > zero(x), values) || throw(ArgumentError(
        "Dwell m/z keys must be finite and positive."))
    values
end

function dwell_normalize_unchecked(msm::MassScanMatrix, dwell::AbstractVector{<:Real})
    normalized = rawintensities(msm) * Diagonal(inv.(dwell))
    MassScanMatrix(
        copy(rawretentions(msm)), retentionunit(msm),
        copy(rawmzvalues(msm)), mzunit(msm), normalized, u"s^-1";
        level=level(msm), instrument=deepcopy(instrument(msm)),
        acquisition=deepcopy(acquisition(msm)), user=deepcopy(user(msm)),
        sample=deepcopy(sample(msm)), extras=deepcopy(extras(msm)))
end

function dwell_normalize_unchecked(vmsm::VarianceMassScanMatrix, dwell::AbstractVector{<:Real})
    normalizedmsm = dwell_normalize_unchecked(parent(vmsm), dwell)
    normalizedvariances = rawvariances(vmsm) * Diagonal(abs2.(inv.(dwell)))
    VarianceMassScanMatrix(normalizedmsm, normalizedvariances)
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

See also
[`dwellnormalize`](@ref),
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassScanMatrix)).

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

See also
[`dwellnormalize`](@ref),
[`intensityunit`](@ref JuChrom.intensityunit(::AbstractMassScanMatrix)),
[`varianceunit`](@ref JuChrom.varianceunit(::AbstractVarianceMassScanMatrix)).
"""
function withintensityunit(vmsm::VarianceMassScanMatrix, unit::Unitful.Units)
    isnothing(intensityunit(vmsm)) ||
        throw(ArgumentError("VarianceMassScanMatrix intensityunit must be nothing."))
    isnothing(varianceunit(vmsm)) ||
        throw(ArgumentError("VarianceMassScanMatrix varianceunit must be nothing."))

    annotatedmsm = withintensityunit(parent(vmsm), unit)
    VarianceMassScanMatrix(annotatedmsm, copy(rawvariances(vmsm)), unit^2)
end
