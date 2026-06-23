"""
    AlkaneReferenceChannels

Matched reference-spectrum channels for one alkane ladder step.
"""
struct AlkaneReferenceChannels{T1<:NamedTuple, T2<:Real}
    carbon::Int
    referenceattrs::T1
    mzindices::Vector{Int}
    mzvalues::Vector{T2}
    referenceintensities::Vector{Float64}
end

"""
    AlkaneChannelInfo

Matched alkane reference channels on the measured m/z grid.
"""
struct AlkaneChannelInfo{
    T1<:Real,
    T2<:AbstractVector{<:AlkaneReferenceChannels},
    T3<:Union{Nothing, Unitful.Units}
}
    mzindices::Vector{Int}
    mzvalues::Vector{T1}
    references::T2
    carbonrange::Vector{Int}
    minrelativeintensity::Float64
    mzunit::T3
end

struct AlkaneReferenceChannelMatch{
    T1<:AbstractVector{<:Integer},
    T2<:AbstractVector{<:Real},
    T3<:AbstractVector{<:Real}
}
    mzindices::T1
    mzvalues::T2
    referenceintensities::T3
end

struct AlkaneBaselineInfo{
    T1<:MassScanMatrix,
    T2<:Real,
    T3<:Real,
    T4<:Real,
    T5<:Real,
    T6<:Real
}
    baselines::T1
    estimator::Symbol
    λ::T2
    nonnegative::Bool
    peakthreshold::T3
    peakslope::T4
    zerothreshold::T5
    zerofractionthreshold::T6
end

abstract type AbstractAlkaneAbundanceInfo end
abstract type AbstractAlkaneMolecularIonInfo end
abstract type AbstractAlkaneLadderPathInfo end
abstract type AbstractAlkaneLadderApexInfo end
abstract type AbstractAlkaneLadderAdditionInfo end
abstract type AbstractAlkaneSeriesDataInfo end
abstract type AbstractAlkaneLadderApex end

"""
    AlkaneSeriesResult

Result container returned by [`findalkanes`](@ref) and [`findalkaneseries`](@ref).

The ladder apices in `apexinfo` and `additioninfo` use raw numeric retention values. The
corresponding unit of those values is stored in `retentionunit`, matching the analyzed
`MassScanMatrix`. Reference-spectrum RI values remain part of the `standard` via
`attrs(spectrum).ri`.

`success` and `status` summarize whether the result is ready for default RT -> RI mapper
construction. `success` is `true` and `status === :ok` only when the default calibration
point selection yields at least three mapper anchors.
"""
struct AlkaneSeriesResult{
    T1<:Union{Nothing, AlkaneStandard},
    T2<:AbstractMatrix{<:Real},
    T3<:Union{Nothing, NamedTuple, CountVarianceEstimate, AbstractVarianceMassScanMatrix},
    T4<:Union{Nothing, AlkaneBaselineInfo},
    T5<:Union{Nothing, AlkaneChannelInfo},
    T6<:Union{Nothing, AbstractAlkaneAbundanceInfo},
    T7<:Union{Nothing, AbstractAlkaneMolecularIonInfo},
    T8<:Union{Nothing, AbstractAlkaneLadderPathInfo},
    T9<:Union{Nothing, AbstractAlkaneLadderApexInfo},
    T10<:Union{Nothing, AbstractAlkaneLadderAdditionInfo},
    T11<:Union{Nothing, AbstractAlkaneSeriesDataInfo},
    T12<:Union{Nothing, Unitful.Units}
}
    success::Bool
    status::Symbol
    standard::T1
    variances::T2
    varianceinfo::T3
    baselineinfo::T4
    channelinfo::T5
    abundanceinfo::T6
    molecularioninfo::T7
    pathinfo::T8
    apexinfo::T9
    additioninfo::T10
    datainfo::T11
    retentionunit::T12
end

retentionunit(result::AlkaneSeriesResult) = result.retentionunit

const MIN_ALKANE_MAPPER_CALIBRATION_POINTS = 3

function AlkaneSeriesResult(
    success::Bool,
    status::Symbol,
    standard::Union{Nothing, AlkaneStandard},
    variances::AbstractMatrix{<:Real},
    varianceinfo::Union{
        Nothing,
        NamedTuple,
        CountVarianceEstimate,
        AbstractVarianceMassScanMatrix
    },
    baselineinfo::Union{Nothing, AlkaneBaselineInfo},
    channelinfo::Union{Nothing, AlkaneChannelInfo},
    abundanceinfo::Union{Nothing, AbstractAlkaneAbundanceInfo},
    molecularioninfo::Union{Nothing, AbstractAlkaneMolecularIonInfo},
    pathinfo::Union{Nothing, AbstractAlkaneLadderPathInfo},
    apexinfo::Union{Nothing, AbstractAlkaneLadderApexInfo},
    additioninfo::Union{Nothing, AbstractAlkaneLadderAdditionInfo}
)
    AlkaneSeriesResult(
        success,
        status,
        standard,
        variances,
        varianceinfo,
        baselineinfo,
        channelinfo,
        abundanceinfo,
        molecularioninfo,
        pathinfo,
        apexinfo,
        additioninfo,
        nothing,
        nothing
    )
end

abstract type AbstractAlkaneLadderCandidate end

"""
    AlkaneLadderStep

Stable public view of one refined alkane ladder step returned by
[`alkaneladdersteps`](@ref).

`apexretention` is unitless and has the unit stored in `result.retentionunit`.
"""
struct AlkaneLadderStep{T<:AbstractAlkaneLadderApex}
    ladderstep::Int
    apexscanindex::Float64
    apexretention::Float64
    source::Symbol
    massspectrumcosine::Float64
    requiredcosine::Float64
    goodforcalibration::Bool
    apex::T
end

"""
    AlkaneLadderCalibrationPoint

Stable RT -> RI calibration point derived from an identified alkane ladder step or added
manually before fitting a retention mapper.

`retention` is the raw numeric retention value. `retentionunit` is the unit associated
with that value, or `nothing` for unitless retentions. `retentionindex` is the Kováts
retention index used as the mapper target.
"""
struct AlkaneLadderCalibrationPoint
    ladderstep::Int
    retention::Float64
    retentionunit::Union{Nothing, Unitful.Units}
    retentionindex::Float64
    source::Symbol
    goodforcalibration::Bool
end

"""
    alkaneladdersteps(result; molecularion=true, gapfilled=true, edgeextended=true)

Return a sorted stable view of already-refined alkane ladder steps.

The returned values are `AlkaneLadderStep` objects. No apex refinement is recomputed:
molecular-ion-supported steps are taken from `result.apexinfo.apexes`, and accepted
gap/edge additions are taken from `addition.apex`.
"""
function alkaneladdersteps(
    result::AlkaneSeriesResult;
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true
)
    alkane_ladder_steps(
        result.apexinfo,
        result.additioninfo;
        molecularion=molecularion,
        gapfilled=gapfilled,
        edgeextended=edgeextended
    )
end

function alkane_ladder_steps(
    apexinfo::AbstractAlkaneLadderApexInfo,
    additioninfo::AbstractAlkaneLadderAdditionInfo;
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true
)
    steps = AlkaneLadderStep[]
    if molecularion
        append!(steps, alkane_ladder_steps_from_apexes(
            apexinfo.apexes,
            :molecularion
        ))
    end

    if gapfilled
        append!(steps, alkane_ladder_steps_from_additions(
            additioninfo.gapfilled,
            :gapfilled
        ))
    end

    if edgeextended
        append!(steps, alkane_ladder_steps_from_additions(
            additioninfo.leftextended,
            :leftextended
        ))
        append!(steps, alkane_ladder_steps_from_additions(
            additioninfo.rightextended,
            :rightextended
        ))
    end

    sort!(steps; by=step -> step.ladderstep)

    steps
end

function alkane_ladder_steps_from_apexes(
    apexes::AbstractVector{<:AbstractAlkaneLadderApex},
    source::Symbol
)
    steps = AlkaneLadderStep[]
    for apex in apexes
        Bool(apex.success) || continue
        push!(steps, alkane_ladder_step_from_apex(apex, source))
    end

    steps
end

function alkane_ladder_steps_from_additions(additions::AbstractVector, source::Symbol)
    steps = AlkaneLadderStep[]
    for addition in additions
        apex = addition.apex
        Bool(apex.success) || continue
        push!(steps, alkane_ladder_step_from_apex(apex, source))
    end

    steps
end

function alkane_ladder_step_from_apex(apex::AbstractAlkaneLadderApex, source::Symbol)
    AlkaneLadderStep(
        apex.ladderstep,
        apex.apexscanindex,
        apex.apexretention,
        source,
        apex.mass_spectrum_cosine,
        apex.required_cosine,
        apex.good_for_calibration,
        apex
    )
end

"""
    alkaneladdercalibrationpoints(result;
        molecularion=true,
        gapfilled=true,
        edgeextended=true,
        goodforcalibration=true,
        exclude=Int[],
        include=Int[],
        extra=AlkaneLadderCalibrationPoint[])

Return RT -> RI calibration points from an [`AlkaneSeriesResult`](@ref).

By default, only ladder steps that passed the calibration-quality gate are returned.
`exclude` removes algorithm-derived ladder steps by carbon number. `include` keeps
algorithm-derived steps with those carbon numbers even when `goodforcalibration=true` and
the step failed the gate. It does not override `molecularion`, `gapfilled`, or
`edgeextended` source filters.

`extra` is appended after `exclude` is applied, so manual points can replace algorithmic
points by excluding the algorithmic step and adding a manual point with the same
`ladderstep`. The returned points are sorted by `ladderstep` and validated to contain no
duplicate ladder steps.
"""
function alkaneladdercalibrationpoints(
    result::AlkaneSeriesResult;
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true,
    goodforcalibration::Bool=true,
    exclude::AbstractVector{<:Integer}=Int[],
    include::AbstractVector{<:Integer}=Int[],
    extra::AbstractVector{<:AlkaneLadderCalibrationPoint}=AlkaneLadderCalibrationPoint[]
)
    alkane_ladder_calibration_points(
        result.standard,
        result.apexinfo,
        result.additioninfo,
        result.retentionunit;
        molecularion=molecularion,
        gapfilled=gapfilled,
        edgeextended=edgeextended,
        goodforcalibration=goodforcalibration,
        exclude=exclude,
        include=include,
        extra=extra
    )
end

function alkane_ladder_calibration_points(
    standard::Union{Nothing, AlkaneStandard},
    apexinfo::AbstractAlkaneLadderApexInfo,
    additioninfo::AbstractAlkaneLadderAdditionInfo,
    retentionunit::Union{Nothing, Unitful.Units};
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true,
    goodforcalibration::Bool=true,
    exclude::AbstractVector{<:Integer}=Int[],
    include::AbstractVector{<:Integer}=Int[],
    extra::AbstractVector{<:AlkaneLadderCalibrationPoint}=AlkaneLadderCalibrationPoint[]
)
    retentionindices = alkane_ladder_retention_indices(standard)
    excludedsteps = Set(exclude)
    includedsteps = Set(include)
    points = AlkaneLadderCalibrationPoint[]
    for step in alkane_ladder_steps(
            apexinfo,
            additioninfo;
            molecularion=molecularion,
            gapfilled=gapfilled,
            edgeextended=edgeextended
        )
        step.ladderstep in excludedsteps && continue
        if goodforcalibration &&
                !step.goodforcalibration &&
                !(step.ladderstep in includedsteps)
            continue
        end
        haskey(retentionindices, step.ladderstep) || throw(ArgumentError(
            "alkane standard has no retention index for C$(step.ladderstep)"))
        push!(points, AlkaneLadderCalibrationPoint(
            step.ladderstep,
            step.apexretention,
            retentionunit,
            retentionindices[step.ladderstep],
            step.source,
            step.goodforcalibration
        ))
    end

    append!(points, extra)
    sort!(points; by=point -> point.ladderstep)
    validate_alkane_ladder_calibration_points(points)

    points
end

function alkane_series_mapper_status(
    standard::Union{Nothing, AlkaneStandard},
    apexinfo::Union{Nothing, AbstractAlkaneLadderApexInfo},
    additioninfo::Union{Nothing, AbstractAlkaneLadderAdditionInfo},
    retentionunit::Union{Nothing, Unitful.Units}
)
    isnothing(standard) && return false, :missing_standard
    (isnothing(apexinfo) || isnothing(additioninfo)) && return false, :too_few_mapper_steps

    points = alkane_ladder_calibration_points(standard, apexinfo, additioninfo, retentionunit)
    length(points) >= MIN_ALKANE_MAPPER_CALIBRATION_POINTS ||
        return false, :too_few_mapper_steps

    true, :ok
end

function alkane_ladder_retention_indices(standard::AlkaneStandard)
    spectra = alkane_standard_spectra(standard)
    retentionindices = Dict{Int, Float64}()
    for (carbon, spectrum) in alkane_spectra_by_carbon(spectra)
        retentionindices[carbon] = attrs(spectrum).ri
    end

    retentionindices
end

function alkane_ladder_retention_indices(::Nothing)
    throw(ArgumentError(
        "alkane ladder calibration points require result.standard"))
end

function validate_alkane_ladder_calibration_points(
    points::AbstractVector{<:AlkaneLadderCalibrationPoint}
)
    seen = Set{Int}()
    for point in points
        isfinite(point.retention) || throw(ArgumentError(
            "alkane ladder calibration point C$(point.ladderstep) has nonfinite retention"))
        isfinite(point.retentionindex) || throw(ArgumentError(
            "alkane ladder calibration point C$(point.ladderstep) has nonfinite RI"))
        point.ladderstep in seen && throw(ArgumentError(
            "duplicate alkane ladder calibration point for C$(point.ladderstep); " *
            "use exclude to replace an algorithmic point with an extra point"))
        push!(seen, point.ladderstep)
    end

    nothing
end

function alkane_ladder_calibration_vectors(
    points::AbstractVector{<:AlkaneLadderCalibrationPoint}
)
    validate_alkane_ladder_calibration_points(points)
    sorted = sort(collect(points); by=point -> point.retention)
    isempty(sorted) && return Float64[], Float64[]
    retention_unit = first(sorted).retentionunit
    for point in sorted
        point.retentionunit == retention_unit || throw(ArgumentError(
            "all alkane ladder calibration points must have the same retention unit"))
    end

    retentionvalues = [point.retention for point in sorted]
    retentions = isnothing(retention_unit) ?
        retentionvalues :
        retentionvalues .* retention_unit
    retentionindices = [point.retentionindex for point in sorted]

    retentions, retentionindices
end

function fitmap(
    points::AbstractVector{<:AlkaneLadderCalibrationPoint};
    λ::Real=3e-9,
    kwargs...
)
    retentions, retentionindices = alkane_ladder_calibration_vectors(points)
    fitmap(retentions, retentionindices; λ=λ, kwargs...)
end

function fitmap(
    result::AlkaneSeriesResult;
    molecularion::Bool=true,
    gapfilled::Bool=true,
    edgeextended::Bool=true,
    goodforcalibration::Bool=true,
    exclude::AbstractVector{<:Integer}=Int[],
    include::AbstractVector{<:Integer}=Int[],
    extra::AbstractVector{<:AlkaneLadderCalibrationPoint}=AlkaneLadderCalibrationPoint[],
    λ::Real=3e-9,
    kwargs...
)
    points = alkaneladdercalibrationpoints(
        result;
        molecularion=molecularion,
        gapfilled=gapfilled,
        edgeextended=edgeextended,
        goodforcalibration=goodforcalibration,
        exclude=exclude,
        include=include,
        extra=extra
    )

    fitmap(points; λ=λ, kwargs...)
end

"""
    findalkaneseries(msm::MassScanMatrix, varianceestimate::CountVarianceEstimate; kwargs...)
    findalkaneseries(vmsm::AbstractVarianceMassScanMatrix; kwargs...)

Identify an n-alkane ladder in a preprocessed GC-MS variance mass-scan matrix.

`findalkaneseries` is the low-level entry point. It assumes that `msm` is already the
signal to analyze and that the variance estimate carries one variance per intensity value.
Use [`findalkanes`](@ref) when count variances should be estimated or an ARPLS baseline
should be subtracted first. Bare variance matrices are not accepted; pass a
[`CountVarianceEstimate`](@ref JuChrom.CountVarianceEstimate) or
[`VarianceMassScanMatrix`](@ref JuChrom.VarianceMassScanMatrix) so the variance unit is
checked against the signal intensity unit.

The algorithm matches reference spectra to the integer-binned m/z grid, fits scan-wise
alkane abundance tracks, detects local abundance windows, scores molecular-ion contrast,
selects a monotone ladder path by dynamic programming, refines selected apices with
ion-level m/z-retention timing correction, and then proposes/refines single-step gap fills
and edge extensions.

Regular tuning keywords:

  * `standard`: alkane reference spectra. The default is `defaultalkanestandard()`. Each
    spectrum must carry `attrs(spectrum).order` and a finite Kováts retention index in
    `attrs(spectrum).ri`.
  * `carbonrange`: carbon numbers considered for the ladder.
  * `minrelativeintensity`: minimum relative reference intensity retained during m/z
    channel matching; `0.0` keeps all reference ions.
  * `variancefloor`: minimum variance used in variance-weighted abundance, molecular-ion,
    path-spectrum, apex, and addition calculations.
  * `nonnegative`: if `true`, abundance estimates are floored at zero.
  * `thresholdfraction`: abundance fraction of the local apex at which abundance-window
    extension stops.
  * `minrisez`: minimum abundance rise, in abundance standard-error units, required for
    accepting a local abundance window.
  * `molecularioncenterzmin`: minimum z-score for the molecular-ion center signal.
  * `molecularionisolationzmin`: minimum z-score for molecular-ion contrast over flanking
    control ions.
  * `pathminsteps`: minimum number of steps required for a successful ladder path.
  * `pathmaxmissingsteps`: maximum number of internal carbon steps that a path may skip.
  * `pathmassspectrummatch`: enables full-spectrum reference matching in the path score.
  * `apexmzscanorder`: m/z scan timing model. Use `:ascending`, `:descending`, or
    `:simultaneous` when known. The default `:inferdirection` compares `:ascending` and
    `:descending`; it does not test `:simultaneous`.
  * `apexmzretentionkwargs`: explicit keyword arguments for [`mzretention`](@ref). If it
    contains `order`, it must be consistent with `apexmzscanorder`.
  * `apexionmzvalues`: explicit m/z values used for apex fitting. When `nothing`, ions are
    selected from the reference spectrum.
  * `apexmaxshiftfromguess`: maximum allowed apex shift, in scans, away from the path
    candidate during ion-level apex refinement.
  * `additionmaxextensionsteps`: maximum number of iterative edge-extension steps before
    the end of `carbonrange` is reached.
  * `additionmassspectrummatch`: enables mass-spectrum gates for proposed gap fills and
    edge extensions.

Advanced or usually-stable keywords:

  * `varianceinfo`, `baselineinfo`, `datainfo`: provenance objects stored in the result.
    They do not change detection, except that `datainfo` is used instead of recomputing
    checksums when supplied.
  * `molecularionwindow`: m/z half-window around molecular-ion and control-ion locations.
  * `molecularionstepmass`: nominal alkane molecular-ion mass increment.
  * `pathstepreward`: reward added per selected path step.
  * `pathmaxcandidatesperstep`: candidate-window cap per carbon number before dynamic
    programming.
  * `pathspacingweight`: penalty for nonsmooth retention spacing.
  * `pathgapincreaseweight`: penalty for increasing spacing between neighboring path gaps.
  * `pathmaxgapratio`: largest allowed adjacent scan-spacing ratio.
  * `pathmissingsteppenalty`: penalty for each skipped internal carbon step.
  * `pathmassspectrummatchdistanceweight`: strength of the path full-spectrum mismatch
    penalty.
  * `pathmassspectrumvariancefloor`: variance floor used only for path-stage local
    full-spectrum matching.
  * `apexscanwindow`: scans on each side of the candidate used for ion-level apex fitting.
  * `apexlogfloorfraction`: local maximum fraction added before log-quadratic apex fits.
  * `apexionexcludemzvalues`: m/z values ignored during automatic apex-ion selection.
  * `apexionminrelativeintensity`: minimum relative reference intensity for automatic
    apex-ion selection.
  * `apexminioncount`: minimum number of ions required for a refined apex.
  * `apexcenteredscantolerance`: scan-offset tolerance for treating an apex fit as
    centered during recentering.
  * `apexfitqualityoutlierz`: robust z-score threshold for excluding poor apex fits from
    calibration; `nothing` disables this exclusion.
  * `apexfitqualityminsteps`: minimum number of refined steps required before robust
    apex-fit-quality z-scores are applied.
  * `apexmzscanordermaxpeaks`: number of spread-out path peaks used in the first scan-order
    inference attempt; `nothing` means all peaks.
  * `apexmzscanorderminpeaks`: minimum usable peaks needed to accept scan-order inference.
  * `apexmzscanorderminapexvarianceratio`: evidence ratio required between the better and
    worse scan-order hypotheses.
  * `apexmzscanordershapeioncount`: number of regularly spaced ions used to fit each
    provisional peak shape during scan-order inference.
  * `apexmzscanordershapemzspacing`: nominal m/z spacing of those provisional shape ions.
  * `apexmzscanorderextremeioncount`: number of low/high m/z ion pairs used for fast
    scan-order contrast.
  * `apexmzscanorderminioncount`: minimum successful ion apex fits required per peak for
    scan-order evidence.
  * `additionminradius`, `additionradiusfraction`, `additionpositionsigmafraction`: scan
    search radius and position penalty for single-step gap fills.
  * `additiongapmincosinefloor`, `additiongapcosinetolerance`: mass-spectrum gate for gap
    fills.
  * `additionedgemaxanchors`, `additionedgeminradius`, `additionedgeradiusfraction`,
    `additionedgemincosinefloor`, `additionedgecosinetolerance`,
    `additionedgecosineanchorcount`, `additionedgepositionsigmafraction`: edge-extension
    anchoring, search-radius, cosine, and position-penalty controls.
  * `additionmassspectrumvariancefloor`: variance floor used for addition-stage
    mass-spectrum matching.

Expected failure modes:

  * The function throws if the m/z axis is not integer-binned or if too few requested
    carbon numbers have observable molecular-ion channels.
  * A path may fail when molecular-ion contrast is too weak, when too few steps pass the
    gates, or when retention spacing is too irregular for the path constraints.
  * Scan-order inference can remain ambiguous when signals are very low, when only early
    ladder steps are present and selected ions have little acquisition-time separation, or
    when peaks are overloaded/saturated so ion apices differ for reasons other than m/z
    scan timing.

The returned `AlkaneSeriesResult` stores reference channel matching in `channelinfo`,
abundance tracks and windows in `abundanceinfo`, molecular-ion evidence in
`molecularioninfo`, the selected path in `pathinfo`, refined apices in `apexinfo`,
gap/edge additions in `additioninfo`, and preprocessing provenance in `varianceinfo`,
`baselineinfo`, and `datainfo`. The raw retention unit of the analyzed matrix is stored in
`retentionunit`. The top-level `success`/`status` fields summarize mapper readiness:
`success` is `true` and `status === :ok` only when the default calibration-point selection
yields at least three RT -> RI mapper anchors.
"""
function findalkaneseries(
    msm::MassScanMatrix,
    variance_source::Union{CountVarianceEstimate, AbstractVarianceMassScanMatrix};
    standard=defaultalkanestandard(),
    varianceinfo=nothing,
    baselineinfo=nothing,
    datainfo=nothing,

    carbonrange=8:40,

    minrelativeintensity=0.0,

    variancefloor=1.0,
    nonnegative::Bool=false,
    thresholdfraction=0.05,
    minrisez=10.0,
    molecularionwindow=1,
    molecularionstepmass=14,
    molecularioncenterzmin=1.645,
    molecularionisolationzmin=1.645,
    pathminsteps=5,
    pathstepreward=0.05,
    pathmaxcandidatesperstep=100,
    pathspacingweight=25.0,
    pathgapincreaseweight=5.0,
    pathmaxgapratio=2.5,
    pathmaxmissingsteps=1,
    pathmissingsteppenalty=2.0,
    pathmassspectrummatch=true,
    pathmassspectrummatchdistanceweight=5.0,
    pathmassspectrumvariancefloor=1.0,
    apexscanwindow=2,
    apexlogfloorfraction=1e-3,
    apexmzscanorder=:inferdirection,
    apexmzretentionkwargs=nothing,
    apexionexcludemzvalues=DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
    apexionmzvalues=nothing,
    apexionminrelativeintensity=0.1,
    apexminioncount=3,
    apexmaxshiftfromguess=3.0,
    apexcenteredscantolerance=0.25,
    apexfitqualityoutlierz=3.0,
    apexfitqualityminsteps=6,
    apexmzscanordermaxpeaks=3,
    apexmzscanorderminpeaks=5,
    apexmzscanorderminapexvarianceratio=2.0,
    apexmzscanordershapeioncount=3,
    apexmzscanordershapemzspacing=14,
    apexmzscanorderextremeioncount=5,
    apexmzscanorderminioncount=8,
    additionminradius=5.0,
    additionradiusfraction=0.15,
    additionpositionsigmafraction=0.05,
    additionmaxextensionsteps=100,
    additionmassspectrummatch=true,
    additiongapmincosinefloor=0.85,
    additiongapcosinetolerance=0.03,
    additionedgemaxanchors=6,
    additionedgeminradius=5.0,
    additionedgeradiusfraction=0.2,
    additionedgemincosinefloor=0.9,
    additionedgecosinetolerance=0.03,
    additionedgecosineanchorcount=3,
    additionedgepositionsigmafraction=0.1,
    additionmassspectrumvariancefloor=1.0

)
    validate_alkane_series_variance_source(msm, variance_source)
    variances, stored_varianceinfo = extract_alkane_series_variances(msm, variance_source)
    result_varianceinfo = isnothing(varianceinfo) ? stored_varianceinfo : varianceinfo

    # Validate inputs
    isnothing(standard) && throw(ArgumentError(
        "findalkaneseries requires an alkane standard"))

    validate_alkane_series_variances(msm, variances)
    channelinfo = alkane_mz_channels(
        msm,
        standard,
        carbonrange,
        minrelativeintensity
    )
    alkane_ladder_require_molecular_ion_channels(
        msm;
        carbonrange=channelinfo.carbonrange,
        ionwindow=molecularionwindow,
        stepmass=molecularionstepmass,
        minsteps=pathminsteps
    )

    # Infer the abundance of each ladder step for each scan, delineate local peak windows,
    # and estimate the variance of each abundance estimate
    abundanceinfo = alkane_abundance_info(
        msm,
        variances,
        channelinfo,
        variancefloor,
        nonnegative,
        thresholdfraction,
        minrisez
    )
    molecularioninfo = alkane_molecular_ion_info(
        msm,
        variances,
        abundanceinfo,
        molecularionwindow,
        molecularionstepmass,
        variancefloor,
        molecularioncenterzmin,
        molecularionisolationzmin
    )
    pathinfo = alkaneladderpath(
        molecularioninfo;
        centerzmin=molecularioncenterzmin,
        isolationzmin=molecularionisolationzmin,
        minsteps=pathminsteps,
        stepreward=pathstepreward,
        maxcandidatesperstep=pathmaxcandidatesperstep,
        spacingweight=pathspacingweight,
        gapincreaseweight=pathgapincreaseweight,
        maxgapratio=pathmaxgapratio,
        maxmissingsteps=pathmaxmissingsteps,
        missingsteppenalty=pathmissingsteppenalty,
        msm=msm,
        variances=variances,
        standard=standard,
        massspectrummatch=pathmassspectrummatch,
        massspectrummatchdistanceweight=pathmassspectrummatchdistanceweight,
        massspectrumvariancefloor=pathmassspectrumvariancefloor
    )
    apexsettings = AlkaneLadderApexSettings(
        standard,
        apexscanwindow,
        variancefloor,
        apexlogfloorfraction,
        apexmzretentionkwargs,
        apexmzscanorder,
        apexionexcludemzvalues,
        apexionmzvalues,
        apexionminrelativeintensity,
        apexminioncount,
        apexmaxshiftfromguess,
        apexcenteredscantolerance,
        apexfitqualityoutlierz,
        apexfitqualityminsteps,
        apexmzscanordermaxpeaks,
        apexmzscanorderminpeaks,
        apexmzscanorderminapexvarianceratio,
        apexmzscanordershapeioncount,
        apexmzscanordershapemzspacing,
        apexmzscanorderextremeioncount,
        apexmzscanorderminioncount
    )
    apexinfo = alkaneladderapexes(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexsettings
    )
    additionsettings = alkane_ladder_addition_settings(
        msm,
        variances,
        abundanceinfo,
        apexinfo,
        additionminradius,
        additionradiusfraction,
        additionpositionsigmafraction,
        additionmaxextensionsteps,
        additionmassspectrummatch,
        additiongapmincosinefloor,
        additiongapcosinetolerance,
        additionedgemaxanchors,
        additionedgeminradius,
        additionedgeradiusfraction,
        additionedgemincosinefloor,
        additionedgecosinetolerance,
        additionedgecosineanchorcount,
        additionedgepositionsigmafraction,
        additionmassspectrumvariancefloor,
        apexscanwindow,
        variancefloor,
        apexlogfloorfraction,
        apexionexcludemzvalues,
        apexionmzvalues,
        apexionminrelativeintensity,
        apexminioncount,
        apexinfo.settings.mzretentionkwargs,
        apexmaxshiftfromguess,
        apexcenteredscantolerance,
        apexfitqualityoutlierz,
        apexfitqualityminsteps,
        channelinfo.carbonrange
    )
    additioninfo = alkaneladderadditions(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexinfo,
        additionsettings,
        standard
    )

    result_datainfo = isnothing(datainfo) ?
        alkane_series_datainfo(msm, msm, variances) :
        datainfo
    result_success, result_status = alkane_series_mapper_status(
        standard,
        apexinfo,
        additioninfo,
        retentionunit(msm)
    )

    AlkaneSeriesResult(
        result_success,
        result_status,
        standard,
        variances,
        result_varianceinfo,
        baselineinfo,
        channelinfo,
        abundanceinfo,
        molecularioninfo,
        pathinfo,
        apexinfo,
        additioninfo,
        result_datainfo,
        retentionunit(msm)
    )
end

function findalkaneseries(vmsm::AbstractVarianceMassScanMatrix; kwargs...)
    findalkaneseries(parent(vmsm), vmsm; kwargs...)
end

function findalkaneseries(msm::MassScanMatrix, variances; kwargs...)
    throw(ArgumentError(
        "findalkaneseries no longer accepts variances as a separate matrix; " *
        "use findalkaneseries(msm, CountVarianceEstimate(...)) or " *
        "findalkaneseries(VarianceMassScanMatrix(msm, variances))"))
end

function alkane_mz_channels(
    msm::MassScanMatrix,
    standard::AlkaneStandard,
    carbonrange,
    minrelativeintensity::Real
)
    isnothing(standard) && throw(ArgumentError(
        "alkane m/z channel matching requires an alkane standard"))
    validate_alkane_channel_settings(minrelativeintensity)

    carbon_numbers = collect(carbonrange)
    isempty(carbon_numbers) && throw(ArgumentError("carbonrange must not be empty"))
    all(carbon -> carbon isa Integer, carbon_numbers) || throw(ArgumentError(
        "carbonrange must contain integer carbon numbers"))

    spectra = alkane_standard_spectra(standard)
    spectra_by_carbon = alkane_spectra_by_carbon(spectra)
    grid_bins = alkane_mz_bins(msm)
    grid_mzs = mzvalues(msm)
    grid_index_by_mz = Dict(mz => index for (index, mz) in pairs(grid_bins))

    referenceinfo = map(carbon_numbers) do carbon
        spectrum = get(spectra_by_carbon, Int(carbon), nothing)
        isnothing(spectrum) && throw(ArgumentError(
            "standard does not contain a reference spectrum for C$(carbon)"))

        match = alkane_reference_channel_match(
            spectrum,
            grid_index_by_mz,
            grid_mzs,
            minrelativeintensity
        )
        isempty(match.mzindices) && throw(ArgumentError(
            "reference spectrum for C$(carbon) has no m/z channels in common with msm"))

        AlkaneReferenceChannels(
            Int(carbon),
            attrs(spectrum),
            match.mzindices,
            match.mzvalues,
            match.referenceintensities
        )
    end

    common_mz_indices = sort!(unique!(reduce(
        vcat,
        (info.mzindices for info in referenceinfo);
        init=Int[]
    )))

    AlkaneChannelInfo(
        common_mz_indices,
        collect(grid_mzs[common_mz_indices]),
        referenceinfo,
        Int.(carbon_numbers),
        Float64(minrelativeintensity),
        mzunit(msm)
    )
end

"""
    findalkanes(msm::MassScanMatrix; kwargs...)
    findalkanes(vmsm::AbstractVarianceMassScanMatrix; kwargs...)

Identify an n-alkane ladder from ordinary raw-count GC-MS data.

`findalkanes` is the high-level entry point. It can estimate count variances, subtract an
ARPLS baseline, store preprocessing provenance, and then call [`findalkaneseries`](@ref)
on the prepared signal. If variances were estimated externally, pass them as a
[`CountVarianceEstimate`](@ref JuChrom.CountVarianceEstimate) or
[`VarianceMassScanMatrix`](@ref JuChrom.VarianceMassScanMatrix). Use `findalkaneseries`
directly when the signal has already been preprocessed outside this function.

Preprocessing keywords:

  * `standard`: alkane reference spectra passed to ladder detection.
  * `variances`: caller-supplied [`CountVarianceEstimate`](@ref JuChrom.CountVarianceEstimate)
    or [`AbstractVarianceMassScanMatrix`](@ref JuChrom.AbstractVarianceMassScanMatrix). If
    `nothing`, variances are estimated from raw counts with [`countvariances`](@ref). Bare
    variance matrices are rejected here because the variance unit must be checked against
    the intensity unit of the analyzed matrix.
  * `subtractbaseline`: if `true`, fit and subtract an ARPLS baseline before ladder
    detection. If `false`, `variances` must be supplied.
  * `variancewindowsize`: scan-window size used by `countvariances`.
  * `variancemintransitioncount`: minimum transition count required for a usable variance
    window.
  * `variancepositivecountquantile`: positive-count quantile used for the count floor.
  * `variancezerothresholdquantile`: quantile used to define near-zero windows during
    variance estimation.
  * `varianceintensityfloor`: explicit count-variance intensity floor, or `nothing` for
    automatic selection.
  * `baselineλ`: ARPLS smoothness parameter.
  * `baselinenonnegative`: constrain fitted baselines to nonnegative intensities.
  * `baselinepeakthreshold`: ARPLS peak-mask threshold.
  * `baselinepeakslope`: ARPLS peak-mask slope parameter.
  * `baselinezerothreshold`: intensity threshold used to identify zeros.
  * `baselinezerofractionthreshold`: if at least this fraction of a channel is zero, zeros
    are treated as real observations rather than sparse dropouts.

Ladder detection keywords are the same named controls documented for
[`findalkaneseries`](@ref): `carbonrange`, `minrelativeintensity`, abundance thresholds,
molecular-ion gates, path scoring, apex m/z scan-order handling, apex fit quality, gap
filling, and edge extension.

The returned `AlkaneSeriesResult` stores estimated variances in `result.varianceinfo`,
baseline provenance in `result.baselineinfo`, and SHA-256 checksums for the raw matrix,
analyzed signal, and variance matrix in `result.datainfo`. Use `result.success` as the
simple pipeline gate before calling `fitmap(result)`; it is `true` only when the default
calibration-point selection yields at least three RT -> RI mapper anchors.
"""
function findalkanes(
    msm::MassScanMatrix;
    standard=defaultalkanestandard(),
    variances=nothing,

    carbonrange=8:40,
    minrelativeintensity=0.0,

    variancewindowsize=13,
    variancemintransitioncount=7,
    variancepositivecountquantile=0.01,
    variancezerothresholdquantile=0.99,
    varianceintensityfloor=nothing,

    baselineλ=1e6,
    subtractbaseline::Bool=true,
    baselinenonnegative=true,
    baselinepeakthreshold=10.0,
    baselinepeakslope=0.5,
    baselinezerothreshold=1.0,
    baselinezerofractionthreshold=0.2,

    variancefloor=1.0,
    nonnegative::Bool=false,
    thresholdfraction=0.05,
    minrisez=10.0,

    molecularionwindow=1,
    molecularionstepmass=14,
    molecularioncenterzmin=1.645,
    molecularionisolationzmin=1.645,

    pathminsteps=5,
    pathstepreward=0.05,
    pathmaxcandidatesperstep=100,
    pathspacingweight=25.0,
    pathgapincreaseweight=5.0,
    pathmaxgapratio=2.5,
    pathmaxmissingsteps=1,
    pathmissingsteppenalty=2.0,
    pathmassspectrummatch=true,
    pathmassspectrummatchdistanceweight=5.0,
    pathmassspectrumvariancefloor=1.0,

    apexscanwindow=2,
    apexlogfloorfraction=1e-3,
    apexmzscanorder=:inferdirection,
    apexmzretentionkwargs=nothing,
    apexionexcludemzvalues=DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
    apexionmzvalues=nothing,
    apexionminrelativeintensity=0.1,
    apexminioncount=3,
    apexmaxshiftfromguess=3.0,
    apexcenteredscantolerance=0.25,
    apexfitqualityoutlierz=3.0,
    apexfitqualityminsteps=6,
    apexmzscanordermaxpeaks=3,
    apexmzscanorderminpeaks=5,
    apexmzscanorderminapexvarianceratio=2.0,
    apexmzscanordershapeioncount=3,
    apexmzscanordershapemzspacing=14,
    apexmzscanorderextremeioncount=5,
    apexmzscanorderminioncount=8,

    additionminradius=5.0,
    additionradiusfraction=0.15,
    additionpositionsigmafraction=0.05,
    additionmaxextensionsteps=100,
    additionmassspectrummatch=true,
    additiongapmincosinefloor=0.85,
    additiongapcosinetolerance=0.03,
    additionedgemaxanchors=6,
    additionedgeminradius=5.0,
    additionedgeradiusfraction=0.2,
    additionedgemincosinefloor=0.9,
    additionedgecosinetolerance=0.03,
    additionedgecosineanchorcount=3,
    additionedgepositionsigmafraction=0.1,
    additionmassspectrumvariancefloor=1.0
)

    isnothing(standard) && throw(ArgumentError("findalkanes requires an alkane standard"))

    if !subtractbaseline && isnothing(variances)
        throw(ArgumentError("variances must be provided when subtractbaseline=false"))
    end

    alkane_ladder_require_molecular_ion_channels(
        msm;
        carbonrange=carbonrange,
        ionwindow=molecularionwindow,
        stepmass=molecularionstepmass,
        minsteps=pathminsteps
    )

    if isnothing(variances)
        varianceestimate = countvariances(
            msm;
            windowsize=variancewindowsize,
            mintransitioncount=variancemintransitioncount,
            positivecountquantile=variancepositivecountquantile,
            zerothresholdquantile=variancezerothresholdquantile,
            intensityfloor=varianceintensityfloor
        )
        σ², varianceinfo = extract_alkane_series_variances(msm, varianceestimate)
        validate_alkane_series_variances(msm, σ²)
    else
        validate_alkane_series_variance_source(msm, variances)
        σ², varianceinfo = extract_alkane_series_variances(msm, variances)
    end

    baselineinfo = if subtractbaseline
        baselines = arpls(
            msm;
            λ=baselineλ,
            nonnegative=baselinenonnegative,
            variances=σ²,
            peakthreshold=baselinepeakthreshold,
            peakslope=baselinepeakslope,
            zerothreshold=baselinezerothreshold,
            zerofractionthreshold=baselinezerofractionthreshold
        )
        validate_alkane_series_baselines(msm, baselines)
        AlkaneBaselineInfo(
            baselines,
            :arpls,
            baselineλ,
            baselinenonnegative,
            baselinepeakthreshold,
            baselinepeakslope,
            baselinezerothreshold,
            baselinezerofractionthreshold
        )
    else
        nothing
    end

    signal = isnothing(baselineinfo) ? msm : msm - baselineinfo.baselines

    datainfo = alkane_series_datainfo(msm, signal, σ²)

    findalkaneseries(
        signal,
        varianceinfo;
        standard=standard,
        varianceinfo=varianceinfo,
        baselineinfo=baselineinfo,
        datainfo=datainfo,
        carbonrange=carbonrange,
        minrelativeintensity=minrelativeintensity,
        variancefloor=variancefloor,
        nonnegative=nonnegative,
        thresholdfraction=thresholdfraction,
        minrisez=minrisez,
        molecularionwindow=molecularionwindow,
        molecularionstepmass=molecularionstepmass,
        molecularioncenterzmin=molecularioncenterzmin,
        molecularionisolationzmin=molecularionisolationzmin,
        pathminsteps=pathminsteps,
        pathstepreward=pathstepreward,
        pathmaxcandidatesperstep=pathmaxcandidatesperstep,
        pathspacingweight=pathspacingweight,
        pathgapincreaseweight=pathgapincreaseweight,
        pathmaxgapratio=pathmaxgapratio,
        pathmaxmissingsteps=pathmaxmissingsteps,
        pathmissingsteppenalty=pathmissingsteppenalty,
        pathmassspectrummatch=pathmassspectrummatch,
        pathmassspectrummatchdistanceweight=pathmassspectrummatchdistanceweight,
        pathmassspectrumvariancefloor=pathmassspectrumvariancefloor,
        apexscanwindow=apexscanwindow,
        apexlogfloorfraction=apexlogfloorfraction,
        apexmzscanorder=apexmzscanorder,
        apexmzretentionkwargs=apexmzretentionkwargs,
        apexionexcludemzvalues=apexionexcludemzvalues,
        apexionmzvalues=apexionmzvalues,
        apexionminrelativeintensity=apexionminrelativeintensity,
        apexminioncount=apexminioncount,
        apexmaxshiftfromguess=apexmaxshiftfromguess,
        apexcenteredscantolerance=apexcenteredscantolerance,
        apexfitqualityoutlierz=apexfitqualityoutlierz,
        apexfitqualityminsteps=apexfitqualityminsteps,
        apexmzscanordermaxpeaks=apexmzscanordermaxpeaks,
        apexmzscanorderminpeaks=apexmzscanorderminpeaks,
        apexmzscanorderminapexvarianceratio=apexmzscanorderminapexvarianceratio,
        apexmzscanordershapeioncount=apexmzscanordershapeioncount,
        apexmzscanordershapemzspacing=apexmzscanordershapemzspacing,
        apexmzscanorderextremeioncount=apexmzscanorderextremeioncount,
        apexmzscanorderminioncount=apexmzscanorderminioncount,
        additionminradius=additionminradius,
        additionradiusfraction=additionradiusfraction,
        additionpositionsigmafraction=additionpositionsigmafraction,
        additionmaxextensionsteps=additionmaxextensionsteps,
        additionmassspectrummatch=additionmassspectrummatch,
        additiongapmincosinefloor=additiongapmincosinefloor,
        additiongapcosinetolerance=additiongapcosinetolerance,
        additionedgemaxanchors=additionedgemaxanchors,
        additionedgeminradius=additionedgeminradius,
        additionedgeradiusfraction=additionedgeradiusfraction,
        additionedgemincosinefloor=additionedgemincosinefloor,
        additionedgecosinetolerance=additionedgecosinetolerance,
        additionedgecosineanchorcount=additionedgecosineanchorcount,
        additionedgepositionsigmafraction=additionedgepositionsigmafraction,
        additionmassspectrumvariancefloor=additionmassspectrumvariancefloor
    )
end

function findalkanes(vmsm::AbstractVarianceMassScanMatrix; kwargs...)
    forwarded = (; kwargs...)
    haskey(forwarded, :variances) && throw(ArgumentError(
        "do not pass variances when calling findalkanes on an AbstractVarianceMassScanMatrix"))

    findalkanes(parent(vmsm); variances=vmsm, kwargs...)
end

function extract_alkane_series_variances(
    msm::MassScanMatrix,
    estimate::CountVarianceEstimate
)
    rawvariances(estimate; unit=default_varianceunit(msm)), estimate
end

function extract_alkane_series_variances(
    msm::MassScanMatrix,
    vmsm::AbstractVarianceMassScanMatrix
)
    rawvariances(vmsm; unit=default_varianceunit(msm)), vmsm
end

function validate_alkane_series_variance_source(
    msm::MassScanMatrix,
    estimate::CountVarianceEstimate
)
    validate_varianceunit(msm, varianceunit(estimate))
    size(rawvariances(estimate; unit=default_varianceunit(msm))) ==
        size(rawintensities(msm)) || throw(DimensionMismatch(
            "variance estimate matrix must match rawintensities(msm)"))

    nothing
end

function validate_alkane_series_variance_source(
    msm::MassScanMatrix,
    variance_source::AbstractVarianceMassScanMatrix
)
    retentionunit(variance_source) == retentionunit(msm) || throw(ArgumentError(
        "variance source retention unit must match msm retention unit"))
    mzunit(variance_source) == mzunit(msm) || throw(ArgumentError(
        "variance source m/z unit must match msm m/z unit"))
    intensityunit(variance_source) == intensityunit(msm) || throw(ArgumentError(
        "variance source intensity unit must match msm intensity unit"))
    rawretentions(variance_source) == rawretentions(msm) || throw(DimensionMismatch(
        "variance source retentions must match msm retentions"))
    rawmzvalues(variance_source) == rawmzvalues(msm) || throw(DimensionMismatch(
        "variance source m/z values must match msm m/z values"))
    size(rawvariances(variance_source)) == size(rawintensities(msm)) ||
        throw(DimensionMismatch(
            "variance source variance matrix must match rawintensities(msm)"))

    nothing
end

function validate_alkane_series_variance_source(msm::MassScanMatrix, variance_source)
    throw(ArgumentError(
        "variances must be nothing, CountVarianceEstimate, or AbstractVarianceMassScanMatrix"))
end

function alkane_standard_spectra(standard::AlkaneStandard)
    validate_alkane_standard_spectra(standard.spectra)
end

alkane_standard_spectra(spectra::AbstractVector{<:AbstractMassSpectrum}) =
    validate_alkane_standard_spectra(spectra)

function alkane_standard_spectra(standard)
    throw(ArgumentError(
        "alkane standard must be an AlkaneStandard or a vector of mass spectra"))
end

function validate_alkane_standard_spectra(spectra::AbstractVector{<:AbstractMassSpectrum})
    isempty(spectra) && throw(ArgumentError("alkane standard spectra must not be empty"))

    spectra
end

function alkane_spectra_by_carbon(spectra::AbstractVector{<:AbstractMassSpectrum})
    spectra_by_carbon = Dict{Int, AbstractMassSpectrum}()
    for spectrum in spectra
        spectrum_attrs = attrs(spectrum)
        hasproperty(spectrum_attrs, :order) || throw(ArgumentError(
            "each alkane reference spectrum must have attrs(spectrum).order"))
        carbon = getproperty(spectrum_attrs, :order)
        carbon isa Integer || throw(ArgumentError(
            "attrs(spectrum).order must be an integer carbon number"))
        :ri in keys(spectrum_attrs) || throw(ArgumentError(
            "each alkane reference spectrum must have attrs(spectrum).ri"))
        ri = spectrum_attrs.ri
        ri isa Real && isfinite(ri) || throw(ArgumentError(
            "attrs(spectrum).ri must be a finite Kováts retention index"))
        haskey(spectra_by_carbon, Int(carbon)) && throw(ArgumentError(
            "standard contains more than one reference spectrum for C$(carbon)"))
        spectra_by_carbon[Int(carbon)] = spectrum
    end

    spectra_by_carbon
end

function validate_alkane_channel_settings(minrelativeintensity::Real)
    isfinite(minrelativeintensity) && 0 ≤ minrelativeintensity < 1 || throw(
        ArgumentError("minrelativeintensity must be finite and in [0, 1)"))

    nothing
end

function alkane_reference_channel_match(
    spectrum::AbstractMassSpectrum,
    grid_index_by_mz::AbstractDict{Int, <:Integer},
    grid_mzs::AbstractVector{<:Real},
    minrelativeintensity::Real
)
    spectrum_mzs = integer_mz_values(mzvalues(spectrum), "reference spectrum m/z values")
    spectrum_intensities = Float64.(intensities(spectrum))
    all(isfinite, spectrum_intensities) || throw(ArgumentError(
        "reference spectrum contains nonfinite intensities"))
    all(≥(0), spectrum_intensities) || throw(ArgumentError(
        "reference spectrum contains negative intensities"))
    max_intensity = maximum(spectrum_intensities)
    max_intensity > 0 || throw(ArgumentError(
        "reference spectrum must contain at least one positive intensity"))

    reference_by_index = Dict{Int, Float64}()
    for (mz, intensity) in zip(spectrum_mzs, spectrum_intensities)
        intensity > 0 || continue
        intensity / max_intensity > minrelativeintensity || continue

        mz_index = get(grid_index_by_mz, mz, nothing)
        isnothing(mz_index) && continue
        reference_by_index[mz_index] = get(reference_by_index, mz_index, 0.0) + intensity
    end

    mz_indices = sort!(collect(keys(reference_by_index)))

    AlkaneReferenceChannelMatch(
        mz_indices,
        collect(grid_mzs[mz_indices]),
        [reference_by_index[index] for index in mz_indices]
    )
end

function alkane_mz_bins(msm::MassScanMatrix)
    validate_alkane_mz_unit(mzunit(msm))
    mzs = isnothing(mzunit(msm)) ? rawmzvalues(msm) : rawmzvalues(msm; unit=Th)

    integer_mz_values(mzs, "msm m/z channels")
end

function validate_alkane_mz_unit(unit::Union{Nothing, Unitful.Units})
    isnothing(unit) && return nothing

    try
        uconvert(Th, 1.0 * unit)
    catch
        throw(ArgumentError(
            "alkane ladder m/z values must be unitless or compatible with u\"Th\""))
    end

    nothing
end

function integer_mz_values(mzs::AbstractVector{<:Real}, label::AbstractString)
    integer_mzs = Vector{Int}(undef, length(mzs))
    for (index, mz) in pairs(mzs)
        isfinite(mz) || throw(ArgumentError("$(label) must be finite"))
        isinteger(mz) || throw(ArgumentError("$(label) must be integer binned"))
        integer_mzs[index] = Int(mz)
    end

    integer_mzs
end

function validate_alkane_series_variances(
    msm::AbstractMassScanMatrix,
    variances::AbstractMatrix{<:Real}
)
    expectedsize = size(rawintensities(msm))
    size(variances) == expectedsize || throw(DimensionMismatch(
        "variances must have size $(expectedsize), matching rawintensities(msm)"))

    for variance in variances
        isfinite(variance) || throw(ArgumentError("variances must be finite"))
        variance ≥ 0 || throw(ArgumentError("variances must be nonnegative"))
    end

    nothing
end

function validate_alkane_series_baselines(
    msm::MassScanMatrix,
    baselines::Union{Nothing, MassScanMatrix}
)
    isnothing(baselines) && return nothing
    size(rawintensities(baselines)) == size(rawintensities(msm)) || throw(
        DimensionMismatch("baselines must match the size of rawintensities(msm)"))
    retentions(baselines) == retentions(msm) || throw(DimensionMismatch(
        "baselines must match msm retention coordinates"))
    mzvalues(baselines) == mzvalues(msm) || throw(DimensionMismatch(
        "baselines must match msm m/z values"))
    level(baselines) == level(msm) || throw(DimensionMismatch(
        "baselines must match msm MS level"))

    nothing
end
