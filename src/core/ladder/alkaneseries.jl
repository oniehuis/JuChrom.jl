
defaultalkanestandard() = nothing

struct AlkaneSeriesResult{S,V,I}
    standard::S
    variances::V
    varianceinfo::I
end

function findalkaneseries(
    msm::MassScanMatrix;
    standard=defaultalkanestandard(),
    variances=nothing,

    carbonrange=8:40,

    variancewindowsize=13,
    variancemintransitioncount=7,
    variancepositivecountquantile=0.01,
    variancezerothresholdquantile=0.99,
    varianceintensityfloor=nothing,

    baselineλ=1e7,
    baselinenonnegative=true,
    baselinepeakthreshold=6.0,
    baselinepeakslope=0.5,
    baselinezerothreshold=0.5,

    mzpeakvariancefloor=1.0,
    mzpeakmaxwidth=101,
    mzpeakzfloor=1.0,
    mzpeakzrise=3.0,
    mzpeakzmin=4.0,

    smoothing=3,
    minrelativeintensity=0.05,
    molecularionwindow=2,
    spectralpower=0.5,
    includedistance=true,

    minsteps=5,
    stepreward=0.05,
    minrelativepeakheight=1e-4,
    minpeakheight=1e-4,
    maxcandidatespertrace=100,
    maxsnapscans=3,
    spacingweight=1.0,
    gapincreaseweight=0.25,

    apexscanwindow=2,
    apexweighting=:variance,
    maxapexoffsetscans=1.0,

    extractstepspectra=false,
    stepspectrascanwindow=2,
    stepspectraweighting=:variance,
    stepspectraallownegative=true,
)

    # Validate inputs
    isnothing(standard) && throw(ArgumentError(
        "findalkaneseries requires an alkane standard; defaultalkanestandard is not implemented yet"))

    # Validate variances if provided
    variances isa AbstractMatrix{<:Real} || throw(ArgumentError(
        "variances must be a matrix matching rawintensities(msm)"))

    expectedsize = size(rawintensities(msm))
    size(variances) == expectedsize || throw(DimensionMismatch(
        "variances must have size $(expectedsize), matching rawintensities(msm)"))

    for variance in variances
        isfinite(variance) || throw(ArgumentError("variances must be finite"))
        variance ≥ 0 || throw(ArgumentError("variances must be nonnegative"))
    end

    # Estimate variances if not provided
    if isnothing(variances)
        countvars = countvariances(
            msm;
            windowsize=variancewindowsize,
            mintransitioncount=variancemintransitioncount,
            positivecountquantile=variancepositivecountquantile,
            zerothresholdquantile=variancezerothresholdquantile,
            intensityfloor=varianceintensityfloor,
        )
        σ² = countvars.variances
        varianceinfo = countvars
    else
        validatealkaneseriesvariances(msm, variances)
        σ² = variances
        varianceinfo = nothing
    end


    # Subtract baseline

    AlkaneSeriesResult(standard, σ², varianceinfo)
end
