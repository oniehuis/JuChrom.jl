"""
    varpred(intensity, model::LinearObservedIntensityVarianceModel;
        varfloor=0, extrapolation=:allow)

Evaluate a linear observed-intensity variance model:

    variance = model.intercept + model.slope * max(intensity - model.intensity_offset, 0)

`intensity` may be a scalar or array. Intensities below `model.intensity_offset` predict
the baseline variance `model.intercept`; only intensity above the offset contributes to
the linear term. `varfloor` is applied elementwise and can be unitful when the model is
unitful.

`model.intensity_min` and `model.intensity_max` store the calibrated intensity range.
Predictions outside that range are linear extrapolations. `extrapolation` controls the
reaction to finite values outside the calibrated range:

- `:allow` ignores the range and returns the extrapolated prediction.
- `:warn` emits a warning and returns the extrapolated prediction.
- `:throw` throws an `ArgumentError`.
"""
@inline function varpred(
    intensity,
    model::LinearObservedIntensityVarianceModel;
    varfloor::Number=0,
    extrapolation::Symbol=:allow,
)
    alkane_variance_check_extrapolation(intensity, model, extrapolation)
    floor = LinearObservedIntensityVarianceFloor(varfloor)
    @. apply_variance_floor(
        linear_observed_intensity_variance(
            linear_variance_intensity_delta(intensity, model.intensity_offset),
            model.intercept,
            model.slope
        ),
        floor
    )
end

function alkane_variance_check_extrapolation(
    intensity,
    model::LinearObservedIntensityVarianceModel,
    extrapolation::Symbol,
)
    extrapolation in (:allow, :warn, :throw) || throw(ArgumentError(
        "extrapolation must be :allow, :warn, or :throw"))
    extrapolation === :allow && return nothing

    values = alkane_variance_calibration_range_values(intensity, model)
    finite = Float64[value for value in values if isfinite(value)]
    isempty(finite) && return nothing

    lower = model.intensity_min
    upper = model.intensity_max
    below = count(<(lower), finite)
    above = count(>(upper), finite)
    below == 0 && above == 0 && return nothing

    observed_min = minimum(finite)
    observed_max = maximum(finite)
    message = string(
        "variance prediction uses intensities outside the calibrated range ",
        "[", lower, ", ", upper, "]; observed range is [",
        observed_min, ", ", observed_max, "]",
    )
    if extrapolation === :warn
        @warn message below above
        return nothing
    end

    throw(ArgumentError(message))
end

function alkane_variance_calibration_range_values(
    intensity::Number,
    model::LinearObservedIntensityVarianceModel,
)
    value = alkane_variance_calibration_range_value(intensity, model)
    isfinite(value) ? Float64[value] : Float64[]
end

function alkane_variance_calibration_range_values(
    intensities,
    model::LinearObservedIntensityVarianceModel,
)
    values = Float64[]
    for intensity in intensities
        value = alkane_variance_calibration_range_value(intensity, model)
        isfinite(value) && push!(values, value)
    end

    values
end

function alkane_variance_calibration_range_value(
    intensity,
    model::LinearObservedIntensityVarianceModel,
)
    if isunitful(intensity)
        if isunitful(model.intensity_offset)
            return Float64(ustrip(Unitful.unit(model.intensity_offset), intensity))
        else
            return NaN
        end
    end

    Float64(intensity)
end

@inline linear_variance_intensity_delta(intensity::Integer, intensity_offset) =
    linear_variance_intensity_delta(float(intensity), intensity_offset)
@inline function linear_variance_intensity_delta(intensity, intensity_offset)
    delta = intensity - intensity_offset
    max(delta, zero(delta))
end

@inline linear_observed_intensity_variance(intensity, intercept, slope) =
    muladd(slope, intensity, intercept)

struct LinearObservedIntensityVarianceFloor{T<:Number}
    value::T
end

Base.broadcastable(floor::LinearObservedIntensityVarianceFloor) = Ref(floor)

@inline function apply_variance_floor(val, floor::LinearObservedIntensityVarianceFloor)
    floor_adj = linear_observed_intensity_variance_coerce_floor(
        floor.value,
        val,
    )
    linear_observed_intensity_variance_floorfinite(val, floor_adj)
end

@inline function linear_observed_intensity_variance_coerce_floor(value, reference)
    if isunitful(reference)
        if isunitful(value)
            Unitful.dimension(value) == Unitful.dimension(reference) ||
                throw(Unitful.DimensionError(value, reference))
            value
        else
            value * Unitful.oneunit(reference)
        end
    else
        isunitful(value) && throw(ArgumentError(
            "varfloor has units but predicted variances are unitless"))
        value
    end
end

@inline function linear_observed_intensity_variance_floorfinite(
    val::T1,
    floor::T2
) where {T1<:Number, T2<:Number}
    T3 = promote_type(T1, T2)
    val_p = T3(val)
    floor_p = T3(floor)
    (isfinite(val) && val_p ≥ floor_p) ? val_p : floor_p
end

"""
    varpred(msm, model::LinearObservedIntensityVarianceModel;
        varfloor=0, extrapolation=:allow)

Return a [`VarianceMassScanMatrix`](@ref JuChrom.VarianceMassScanMatrix) with variances predicted from
`rawintensities(msm)` using `model`.

This is the scan-matrix counterpart to `varpred(intensity, model)`. The model must be on
the same intensity scale as `msm`. Unitful models are evaluated after attaching
`intensityunit(msm)` to the numeric raw intensity matrix. `extrapolation` has the same
meaning as for scalar or array inputs.
"""
function varpred(
    msm::AbstractMassScanMatrix,
    model::LinearObservedIntensityVarianceModel;
    varfloor::Number=0,
    extrapolation::Symbol=:allow,
)
    model_has_units = isunitful(model.intercept) ||
        isunitful(model.slope) ||
        isunitful(model.intensity_offset)
    intensity_unit = intensityunit(msm)
    intensities_for_model = if model_has_units && !isnothing(intensity_unit)
        rawintensities(msm) .* intensity_unit
    else
        rawintensities(msm)
    end

    predicted = varpred(
        intensities_for_model,
        model;
        varfloor=varfloor,
        extrapolation=extrapolation,
    )
    VarianceMassScanMatrix(msm, predicted)
end
