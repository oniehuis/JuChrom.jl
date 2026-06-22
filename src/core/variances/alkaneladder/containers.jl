"""
    LinearObservedIntensityVarianceModel(
        intercept,
        slope,
        intensity_offset,
        intensity_min,
        intensity_max,
        rho_lag1
    )

Empirical variance model calibrated from alkane ladder peak residuals:

    variance(I) = intercept + slope * max(I - intensity_offset, 0)

`I` is the baseline-inclusive observed intensity on the same preprocessing scale used
during calibration. `intercept` is the variance at `intensity_offset`, typically the
flat-region baseline intensity. `intensity_min` and `intensity_max` record the calibrated
range used for extrapolation checks; for this model the lower bound is normally `0`
because all intensities below `intensity_offset` predict the baseline variance. `rho_lag1`
stores the estimated lag-1 residual autocorrelation for variance inflation during scan
pooling.

The model can be unitless or unitful. For unitful models, `intercept` must have
intensity² units, and `slope` and `intensity_offset` must have intensity units.
`rho_lag1` is always unitless.
"""
struct LinearObservedIntensityVarianceModel{
    T1<:Number,
    T2<:Number,
    T3<:Number,
    T4<:Real,
    T5<:Real,
    T6<:Real
}
    intercept::T1
    slope::T2
    intensity_offset::T3
    intensity_min::T4
    intensity_max::T5
    rho_lag1::T6

    function LinearObservedIntensityVarianceModel(
        intercept::T1,
        slope::T2,
        intensity_min::T3,
        intensity_max::T4,
        rho_lag1::T5
    ) where {
        T1<:Number,
        T2<:Number,
        T3<:Number,
        T4<:Real,
        T5<:Real
    }
        LinearObservedIntensityVarianceModel(
            intercept,
            slope,
            zero(slope),
            intensity_min,
            intensity_max,
            rho_lag1,
        )
    end

    function LinearObservedIntensityVarianceModel(
        intercept::T1,
        slope::T2,
        intensity_offset::T3,
        intensity_min::T4,
        intensity_max::T5,
        rho_lag1::T6
    ) where {
        T1<:Number,
        T2<:Number,
        T3<:Number,
        T4<:Real,
        T5<:Real,
        T6<:Real
    }
        validate_linear_observed_intensity_variance_units(
            intercept,
            slope,
            intensity_offset,
        )
        isfinite(intensity_offset) || throw(ArgumentError(
            "intensity_offset must be finite"))
        isfinite(intensity_min) || throw(ArgumentError("intensity_min must be finite"))
        isfinite(intensity_max) || throw(ArgumentError("intensity_max must be finite"))
        intensity_min ≤ intensity_max || throw(ArgumentError(
            "intensity_min must be less than or equal to intensity_max"))
        isfinite(rho_lag1) || throw(ArgumentError("rho_lag1 must be finite"))
        abs(rho_lag1) < 1 || throw(ArgumentError("rho_lag1 must satisfy abs(rho_lag1) < 1"))

        zero_intercept = zero(intercept)
        zero_slope = zero(slope)
        intercept ≥ zero_intercept || throw(ArgumentError("intercept must be nonnegative"))
        slope ≥ zero_slope || throw(ArgumentError("slope must be nonnegative"))

        new{T1, T2, T3, T4, T5, T6}(
            intercept,
            slope,
            intensity_offset,
            intensity_min,
            intensity_max,
            rho_lag1
        )
    end
end

Base.broadcastable(model::LinearObservedIntensityVarianceModel) = Base.RefValue(model)

"""
    AlkaneVarianceFit

Result returned by [`fitalkanevariancemodel`](@ref). The fields expose the fitted model,
quality control, and peak-level calibration details directly; display is compact for REPL
use.
"""
struct AlkaneVarianceFit{
    TModel,
    TQC,
    TMSM,
    TResult,
    TSignal,
    TExcluded,
    TIncluded,
    TSettings,
    TSpectra,
    TFailures,
    TPeakInputs,
    TPeakInputFailures,
    TPeakFits,
    TPeakFitFailures,
    TPeakQC,
    TResidualRecords,
    TFlatRecords,
    TFlatQC,
    TVarianceFit,
}
    success::Bool
    status::Symbol
    model::TModel
    qc::TQC
    msm::TMSM
    result::TResult
    signal::TSignal
    excludeladdersteps::TExcluded
    includedladdersteps::TIncluded
    settings::TSettings
    spectra::TSpectra
    failures::TFailures
    peakinputs::TPeakInputs
    peakinputfailures::TPeakInputFailures
    peakfits::TPeakFits
    peakfitfailures::TPeakFitFailures
    peakqc::TPeakQC
    residualrecords::TResidualRecords
    flatrecords::TFlatRecords
    flatqc::TFlatQC
    variancefit::TVarianceFit
end

function validate_linear_observed_intensity_variance_units(intercept, slope, intensity_offset)
    intercept_has_units = isunitful(intercept)
    slope_has_units = isunitful(slope)
    offset_has_units = isunitful(intensity_offset)

    if intercept_has_units
        intensity_ref = sqrt(Unitful.oneunit(intercept))
        slope_has_units || throw(Unitful.DimensionError(slope, intensity_ref))
        Unitful.dimension(slope) == Unitful.dimension(intensity_ref) ||
            throw(Unitful.DimensionError(slope, intensity_ref))
        offset_has_units || throw(Unitful.DimensionError(intensity_offset, intensity_ref))
        Unitful.dimension(intensity_offset) == Unitful.dimension(intensity_ref) ||
            throw(Unitful.DimensionError(intensity_offset, intensity_ref))
    elseif slope_has_units
        throw(Unitful.DimensionError(slope, one(Float64)))
    elseif offset_has_units
        throw(Unitful.DimensionError(intensity_offset, one(Float64)))
    end

    nothing
end

function Base.show(io::IO, model::LinearObservedIntensityVarianceModel)
    print(io, "LinearObservedIntensityVarianceModel(")
    print(io, "intercept=", model.intercept)
    print(io, ", slope=", model.slope)
    print(io, ", intensity_offset=")
    print(io, model.intensity_offset)
    print(io, ", intensity_range=(")
    print(io, model.intensity_min)
    print(io, ", ")
    print(io, model.intensity_max)
    print(io, "), rho_lag1=")
    print(io, model.rho_lag1)
    print(io, ")")
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    model::LinearObservedIntensityVarianceModel,
)
    println(io, "LinearObservedIntensityVarianceModel")
    print(io, "  variance(I) = ")
    print(io, model.intercept)
    print(io, " + ")
    print(io, model.slope)
    print(io, " * max(I - ")
    print(io, model.intensity_offset)
    println(io, ", 0)")
    print(io, "  intensity range: [")
    print(io, model.intensity_min)
    print(io, ", ")
    print(io, model.intensity_max)
    println(io, "]")
    print(io, "  lag-1 residual autocorrelation: ")
    print(io, model.rho_lag1)
end
