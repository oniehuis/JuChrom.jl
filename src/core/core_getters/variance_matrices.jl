# Getter forwarding for VarianceMassScanMatrix keeps the wrapper compatible with the
# AbstractMassScanMatrix interface while exposing variance-specific accessors.

acquisition(vmsm::VarianceMassScanMatrix) = acquisition(parent(vmsm))
extras(vmsm::VarianceMassScanMatrix) = extras(parent(vmsm))
instrument(vmsm::VarianceMassScanMatrix) = instrument(parent(vmsm))
intensityunit(vmsm::VarianceMassScanMatrix) = intensityunit(parent(vmsm))
level(vmsm::VarianceMassScanMatrix) = level(parent(vmsm))
mzcount(vmsm::VarianceMassScanMatrix) = mzcount(parent(vmsm))
mzunit(vmsm::VarianceMassScanMatrix) = mzunit(parent(vmsm))
retentionunit(vmsm::VarianceMassScanMatrix) = retentionunit(parent(vmsm))
sample(vmsm::VarianceMassScanMatrix) = sample(parent(vmsm))
scancount(vmsm::VarianceMassScanMatrix) = scancount(parent(vmsm))
user(vmsm::VarianceMassScanMatrix) = user(parent(vmsm))

@inline function intensities(
    vmsm::VarianceMassScanMatrix;
    unit::Union{Nothing, Unitful.Units}=nothing)
    intensities(parent(vmsm); unit=unit)
end

@inline function rawintensities(
    vmsm::VarianceMassScanMatrix;
    unit::Union{Nothing, Unitful.Units}=nothing)
    rawintensities(parent(vmsm); unit=unit)
end

@inline function mzvalues(
    vmsm::VarianceMassScanMatrix;
    unit::Union{Nothing, Unitful.Units}=nothing)
    mzvalues(parent(vmsm); unit=unit)
end

@inline function rawmzvalues(
    vmsm::VarianceMassScanMatrix;
    unit::Union{Nothing, Unitful.Units}=nothing)
    rawmzvalues(parent(vmsm); unit=unit)
end

@inline function retentions(
    vmsm::VarianceMassScanMatrix;
    unit::Union{Nothing, Unitful.Units}=nothing)
    retentions(parent(vmsm); unit=unit)
end

@inline function rawretentions(
    vmsm::VarianceMassScanMatrix;
    unit::Union{Nothing, Unitful.Units}=nothing)
    rawretentions(parent(vmsm); unit=unit)
end

"""
    varianceunit(vmsm::AbstractVarianceMassScanMatrix) -> Union{Unitful.Units, Nothing}

Return the unit associated with `variances(vmsm)`, or `nothing` for unitless variances.
"""
varianceunit(vmsm::AbstractVarianceMassScanMatrix) = vmsm.varianceunit

"""
    variances(vmsm::AbstractVarianceMassScanMatrix; unit=nothing)

Return the variance matrix, attaching or converting units when `varianceunit(vmsm)` is set.
"""
@inline function variances(
    vmsm::AbstractVarianceMassScanMatrix;
    unit::Union{Nothing, Unitful.Units}=nothing)

    vunit = varianceunit(vmsm)
    isnothing(vunit) && return _handle_unitless(vmsm.variances, unit, "variances")
    _handle_unitful_convert(vmsm.variances, vunit, unit)
end

"""
    rawvariances(vmsm::AbstractVarianceMassScanMatrix; unit=nothing)

Return the numeric variance matrix, optionally converting unitful variances before stripping.
"""
@inline function rawvariances(
    vmsm::AbstractVarianceMassScanMatrix;
    unit::Union{Nothing, Unitful.Units}=nothing)

    vunit = varianceunit(vmsm)
    isnothing(vunit) && return _handle_unitless(vmsm.variances, unit, "variances")
    _handle_unitful_strip(vmsm.variances, vunit, unit)
end
