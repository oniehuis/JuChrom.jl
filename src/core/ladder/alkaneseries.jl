struct AlkaneSeriesResult{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11}
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
end

function AlkaneSeriesResult(
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
)
    AlkaneSeriesResult(
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
    )
end

struct AlkaneLadderStep{T}
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
    steps = AlkaneLadderStep[]
    if molecularion
        append!(steps, alkane_ladder_steps_from_apexes(
            result.apexinfo.apexes,
            :molecularion,
        ))
    end

    if gapfilled
        append!(steps, alkane_ladder_steps_from_additions(
            result.additioninfo.gapfilled,
            :gapfilled,
        ))
    end

    if edgeextended
        append!(steps, alkane_ladder_steps_from_additions(
            result.additioninfo.leftextended,
            :leftextended,
        ))
        append!(steps, alkane_ladder_steps_from_additions(
            result.additioninfo.rightextended,
            :rightextended,
        ))
    end

    sort!(steps; by=step -> step.ladderstep)

    steps
end

function alkane_ladder_steps_from_apexes(apexes::AbstractVector, source::Symbol)
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

function alkane_ladder_step_from_apex(apex, source::Symbol)
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
    findalkaneseries(msm::MassScanMatrix, variances;
        standard=defaultalkanestandard(), varianceinfo=nothing,
        baselineinfo=nothing, datainfo=nothing, ...)

Identify an n-alkane series in a preprocessed GC-MS mass-scan matrix.

This is the low-level analysis entry point. It assumes that `msm` is already the signal
to analyze and therefore requires per-cell `variances`. It does not estimate variances
or fit/subtract a baseline. If upstream preprocessing should be retained with the result,
pass it as `baselineinfo`.

The function first matches the alkane reference spectra to the integer-binned m/z grid of
`msm`. For every requested carbon number it then estimates a scan-wise alkane abundance
track by fitting the reference ion pattern to the measured scan intensities with the
provided variances. Local abundance windows are detected on those tracks. Inside each
window, molecular-ion evidence is computed by comparing the expected molecular ion with
neighboring control ions. A dynamic-programming path then selects a monotone ladder through
the carbon series, balancing molecular-ion score, retention spacing, missing steps, and
optional local mass-spectrum match evidence. Finally, selected path steps are refined with
ion-level apex fitting and m/z-retention timing correction. Single missing internal steps
and edge extensions are proposed separately and are apex-refined before being returned.

`standard` supplies the reference spectra and defaults to `defaultalkanestandard()`.
`varianceinfo`, `baselineinfo`, and `datainfo` are stored in the result so that callers can
keep preprocessing provenance together with the ladder result; they do not change the
algorithm except that `datainfo` is used as supplied instead of being generated from `msm`
and `variances`. `carbonrange` selects the carbon numbers to consider. `minrelativeintensity`
filters reference ions before m/z channel matching, with `0.0` retaining all reference ions.

`variancefloor` is the minimum variance used in variance-weighted calculations, preventing
near-zero variances from dominating fits. `nonnegative` controls whether abundance fits are
floored at zero. `thresholdfraction` controls how far an abundance track may fall from its
local apex before a peak window stops, expressed as a fraction of the apex abundance.
`minrisez` is the minimum abundance rise, in abundance-standard-error units, required for
accepting a local abundance window.

`molecularionwindow` controls the m/z half-window around the molecular-ion and control-ion
locations. `molecularionstepmass` is the nominal m/z spacing between neighboring alkane
fragment or molecular-ion positions. `molecularioncenterzmin` requires the molecular-ion
center signal to be above noise, and `molecularionisolationzmin` requires it to be stronger
than the flanking control ions. The molecular-ion score used downstream is based on the
contrast after these gates.

`pathminsteps` is the minimum number of steps required for a successful path.
`pathstepreward` rewards longer paths. `pathmaxcandidatesperstep` limits how many candidate
windows per carbon number are retained for dynamic programming. `pathspacingweight`
penalizes deviations from smooth retention spacing. `pathgapincreaseweight` penalizes
changes in spacing between neighboring ladder gaps. `pathmaxgapratio` rejects paths with
implausibly large adjacent spacing ratios. `pathmaxmissingsteps` controls how many internal
carbon steps may be skipped, and `pathmissingsteppenalty` penalizes such skips.
`pathmassspectrummatch` enables local full-spectrum matching in the path objective, and
`pathmassspectrummatchdistanceweight` controls how strongly poor mass-spectrum matches are
penalized.

`apexscanwindow` sets the number of scans on each side of the candidate scan used for
ion-level apex fitting. `apexlogfloorfraction` adds a small fraction of the local maximum
before log fitting to stabilize low intensities. `apexmzretentionkwargs` may provide
explicit keyword arguments for [`mzretention`](@ref). `apexmzscanorder` controls the m/z
timing model used during apex refinement. Use `:ascending` or `:descending` for known
sequential quadrupole scan directions, and `:simultaneous` for TOF-like or other
full-spectrum acquisition where all m/z values in a scan are treated as observed at the
scan-level retention. The default `:inferdirection` infers only the sequential quadrupole
direction by comparing `:ascending` and `:descending`; it does not test `:simultaneous`.
`apexmzscanordermaxpeaks` is the initial number of spread-out ladder peaks used for this
direction inference. If the subset is insufficient or ambiguous, more peaks may be used as
defined by the scan-order inference logic. `apexmzscanorderminpeaks` is the minimum number
of usable peaks needed to judge direction. `apexmzscanorderminapexvarianceratio` is the
required evidence ratio between the better and worse scan-order hypotheses.
`apexmzscanordershapeioncount` and `apexmzscanordershapemzspacing` control the small
regularly spaced ion set used to fit each provisional peak shape during scan-order
inference. `apexmzscanorderextremeioncount` controls how many low- and high-m/z contrast
ions are selected for fixed-shape apex alignment, and `apexmzscanorderminioncount` caps the
minimum number of successful ion apex fits required per candidate peak.

`additionminradius`, `additionradiusfraction`, and `additionpositionsigmafraction` control
the search radius and position penalty for proposed missing internal steps. The effective
radius is the maximum of `additionminradius` and `additionradiusfraction` times the local
step spacing, while `additionpositionsigmafraction` sets the Gaussian position scale used
when ranking candidates. `additionmaxextensionsteps` limits iterative edge extension at
both ends of the ladder. The effective per-edge value is additionally capped by the
remaining carbon numbers between the current refined edge and the corresponding end of
`carbonrange`.

The returned `AlkaneSeriesResult` stores reference channel matching in `channelinfo`,
abundance tracks and windows in `abundanceinfo`, molecular-ion evidence in
`molecularioninfo`, the selected dynamic-programming path in `pathinfo`, refined path
apexes in `apexinfo`, gap and edge additions in `additioninfo`, and preprocessing metadata
in `varianceinfo`, `baselineinfo`, and `datainfo`.
"""
function findalkaneseries(
    msm::MassScanMatrix,
    variances;
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
    apexscanwindow=2,
    apexlogfloorfraction=1e-3,
    apexmzscanorder=:inferdirection,
    apexmzretentionkwargs=nothing,
    apexmzscanordermaxpeaks=3,
    apexmzscanorderminpeaks=3,
    apexmzscanorderminapexvarianceratio=1.25,
    apexmzscanordershapeioncount=3,
    apexmzscanordershapemzspacing=14,
    apexmzscanorderextremeioncount=5,
    apexmzscanorderminioncount=8,
    additionminradius=5.0,
    additionradiusfraction=0.15,
    additionpositionsigmafraction=0.05,
    additionmaxextensionsteps=100,

)

    # Validate inputs
    isnothing(standard) && throw(ArgumentError(
        "findalkaneseries requires an alkane standard"))

    validate_alkane_series_variances(msm, variances)
    channelinfo = alkane_mz_channels(
        msm;
        standard=standard,
        carbonrange=carbonrange,
        minrelativeintensity=minrelativeintensity,
    )
    alkane_ladder_require_molecular_ion_channels(
        msm;
        carbonrange=channelinfo.carbonrange,
        ionwindow=molecularionwindow,
        stepmass=molecularionstepmass,
        minsteps=pathminsteps,
    )

    # Infer the abundance of each ladder step for each scan, delineate local peak windows,
    # and estimate the variance of each abundance estimate
    abundanceinfo = alkaneabundanceinfo(
        msm,
        variances,
        channelinfo;
        variancefloor=variancefloor,
        nonnegative=nonnegative,
        thresholdfraction=thresholdfraction,
        minrisez=minrisez,
    )
    molecularioninfo = alkanemolecularioninfo(
        msm,
        variances,
        abundanceinfo;
        ionwindow=molecularionwindow,
        stepmass=molecularionstepmass,
        variancefloor=variancefloor,
        centerzmin=molecularioncenterzmin,
        isolationzmin=molecularionisolationzmin,
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
    )
    apexinfo = alkaneladderapexes(
        msm,
        variances,
        abundanceinfo,
        pathinfo;
        standard=standard,
        scanwindow=apexscanwindow,
        variancefloor=variancefloor,
        logfloorfraction=apexlogfloorfraction,
        mzretentionkwargs=apexmzretentionkwargs,
        mzscanorder=apexmzscanorder,
        mzscanordermaxpeaks=apexmzscanordermaxpeaks,
        mzscanorderminpeaks=apexmzscanorderminpeaks,
        mzscanorderminapexvarianceratio=apexmzscanorderminapexvarianceratio,
        mzscanordershapeioncount=apexmzscanordershapeioncount,
        mzscanordershapemzspacing=apexmzscanordershapemzspacing,
        mzscanorderextremeioncount=apexmzscanorderextremeioncount,
        mzscanorderminioncount=apexmzscanorderminioncount,
    )
    additioninfo = alkaneladderadditions(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexinfo;
        minradius=additionminradius,
        radiusfraction=additionradiusfraction,
        positionsigmafraction=additionpositionsigmafraction,
        maxextensionsteps=additionmaxextensionsteps,
        standard=standard,
        apexscanwindow=apexscanwindow,
        apexvariancefloor=variancefloor,
        apexlogfloorfraction=apexlogfloorfraction,
        apexmzretentionkwargs=apexinfo.settings.mzretentionkwargs,
        carbonrange=channelinfo.carbonrange,
    )

    result_datainfo = isnothing(datainfo) ?
        alkane_series_datainfo(msm, msm, variances) :
        datainfo

    AlkaneSeriesResult(
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
        result_datainfo,
    )
end

function alkane_mz_channels(
    msm::MassScanMatrix;
    standard=defaultalkanestandard(),
    carbonrange=8:40,
    minrelativeintensity=0.0,
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
            grid_mzs;
            minrelativeintensity=minrelativeintensity,
        )
        isempty(match.mzindices) && throw(ArgumentError(
            "reference spectrum for C$(carbon) has no m/z channels in common with msm"))

        (
            carbon=Int(carbon),
            referenceattrs=attrs(spectrum),
            mzindices=match.mzindices,
            mzvalues=match.mzvalues,
            referenceintensities=match.referenceintensities,
        )
    end

    common_mz_indices = sort!(unique!(reduce(
        vcat,
        (info.mzindices for info in referenceinfo);
        init=Int[],
    )))

    (
        mzindices=common_mz_indices,
        mzvalues=collect(grid_mzs[common_mz_indices]),
        references=referenceinfo,
        carbonrange=Int.(carbon_numbers),
        minrelativeintensity=Float64(minrelativeintensity),
        mzunit=mzunit(msm),
    )
end

"""
    findalkanes(msm::MassScanMatrix; standard=defaultalkanestandard(),
        variances=nothing, subtractbaseline=true, kwargs...)

Convenience wrapper for raw-count GC-MS data.

This is the high-level entry point for ordinary uncorrected GC-MS count data. It estimates
or accepts variances, optionally fits and subtracts a baseline, records preprocessing
metadata, and then calls [`findalkaneseries`](@ref) on the prepared signal. Use
`findalkaneseries` directly when you already have a baseline-corrected signal, a custom
variance model, or nonstandard preprocessing that should not be repeated by this wrapper.

`standard` supplies the reference alkane spectra and is passed through to
`findalkaneseries`. `variances` may be an already prepared variance matrix with the same
shape as `rawintensities(msm)`. When `variances === nothing`, variances are estimated from
the raw counts with [`countvariances`](@ref). `subtractbaseline` controls whether an ARPLS
baseline is fitted and subtracted before ladder detection. If `subtractbaseline=false`,
`variances` must be provided, because count-based variances are intended for uncorrected
raw counts and should not be inferred from caller-preprocessed signals.

`variancewindowsize`, `variancemintransitioncount`, `variancepositivecountquantile`, and
`variancezerothresholdquantile` are passed to [`countvariances`](@ref) when `variances` is
not supplied. They control the scan-window size, the minimum number of level crossings
needed for a usable noise window, the positive-count quantile used for the count floor, and
the quantile used to define near-zero windows. `varianceintensityfloor` overrides the
automatically chosen intensity floor when it is not `nothing`.

`baselineλ`, `baselinenonnegative`, `baselinepeakthreshold`, `baselinepeakslope`, and
`baselinezerothreshold` are passed to [`arpls`](@ref) when `subtractbaseline=true`.
`baselineλ` controls baseline smoothness. `baselinenonnegative` constrains the fitted
baseline to nonnegative intensities. `baselinepeakthreshold`, `baselinepeakslope`, and
`baselinezerothreshold` control peak masking and zero-level handling during baseline
estimation.

All remaining keyword arguments are forwarded unchanged to [`findalkaneseries`](@ref).
This includes the carbon range, abundance-window thresholds, molecular-ion gates, path
selection weights, apex m/z scan-order handling, and gap or edge extension settings.

The returned result is the `AlkaneSeriesResult` from `findalkaneseries`. If variances were
estimated, the full variance estimate is stored in `result.varianceinfo`. If a baseline was
subtracted, the fitted baseline matrix and ARPLS settings are stored in
`result.baselineinfo`, and the checksums in `result.datainfo` identify the raw matrix,
baseline-corrected signal, and variance matrix used for detection.
"""
function findalkanes(
    msm::MassScanMatrix;
    standard=defaultalkanestandard(),
    variances=nothing,

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

    kwargs...,
)

    isnothing(standard) && throw(ArgumentError("findalkanes requires an alkane standard"))

    if !subtractbaseline && isnothing(variances)
        throw(ArgumentError("variances must be provided when subtractbaseline=false"))
    end

    forwarded = (; kwargs...)
    alkane_ladder_require_molecular_ion_channels(
        msm;
        carbonrange=get(forwarded, :carbonrange, 8:40),
        ionwindow=get(forwarded, :molecularionwindow, 1),
        stepmass=get(forwarded, :molecularionstepmass, 14),
        minsteps=get(forwarded, :pathminsteps, 5),
    )

    if isnothing(variances)
        varianceestimate = countvariances(
            msm;
            windowsize=variancewindowsize,
            mintransitioncount=variancemintransitioncount,
            positivecountquantile=variancepositivecountquantile,
            zerothresholdquantile=variancezerothresholdquantile,
            intensityfloor=varianceintensityfloor,
        )
        σ², varianceinfo = extract_alkane_series_variances(varianceestimate)
        validate_alkane_series_variances(msm, σ²)
    else
        validate_alkane_series_variances(msm, variances)
        σ² = variances
        varianceinfo = nothing
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
        )
        validate_alkane_series_baselines(msm, baselines)
        (
            baselines=baselines,
            estimator=:arpls,
            λ=baselineλ,
            nonnegative=baselinenonnegative,
            peakthreshold=baselinepeakthreshold,
            peakslope=baselinepeakslope,
            zerothreshold=baselinezerothreshold,
        )
    else
        nothing
    end

    signal = isnothing(baselineinfo) ? msm : msm - baselineinfo.baselines

    datainfo = alkane_series_datainfo(msm, signal, σ²)

    findalkaneseries(
        signal,
        σ²;
        standard=standard,
        varianceinfo=varianceinfo,
        baselineinfo=baselineinfo,
        datainfo=datainfo,
        kwargs...,
    )
end

extract_alkane_series_variances(variances::AbstractMatrix{<:Real}) = variances, nothing

extract_alkane_series_variances(varianceestimate) =
    varianceestimate.variances, varianceestimate

function alkane_standard_spectra(standard::AlkaneStandard)
    spectra = standard.spectra
    spectra isa AbstractVector{<:AbstractMassSpectrum} || throw(ArgumentError(
        "standard.spectra must be a vector of mass spectra"))

    validate_alkane_standard_spectra(spectra)
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
        haskey(spectra_by_carbon, Int(carbon)) && throw(ArgumentError(
            "standard contains more than one reference spectrum for C$(carbon)"))
        spectra_by_carbon[Int(carbon)] = spectrum
    end

    spectra_by_carbon
end

function validate_alkane_channel_settings(minrelativeintensity)
    minrelativeintensity isa Real || throw(ArgumentError(
        "minrelativeintensity must be real"))
    isfinite(minrelativeintensity) && 0 ≤ minrelativeintensity < 1 || throw(
        ArgumentError("minrelativeintensity must be finite and in [0, 1)"))

    nothing
end

function alkane_reference_channel_match(
    spectrum::AbstractMassSpectrum,
    grid_index_by_mz::AbstractDict{Int, <:Integer},
    grid_mzs::AbstractVector;
    minrelativeintensity::Real,
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

    (
        mzindices = mz_indices,
        mzvalues = collect(grid_mzs[mz_indices]),
        referenceintensities = [reference_by_index[index] for index in mz_indices],
    )
end

function alkane_mz_bins(msm::MassScanMatrix)
    validate_alkane_mz_unit(mzunit(msm))
    mzs = isnothing(mzunit(msm)) ? rawmzvalues(msm) : rawmzvalues(msm; unit=Th)

    integer_mz_values(mzs, "msm m/z channels")
end

function validate_alkane_mz_unit(unit)
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

function validate_alkane_series_variances(msm::MassScanMatrix, variances)
    variances isa AbstractMatrix{<:Real} || throw(ArgumentError(
        "variances must be a matrix matching rawintensities(msm)"))

    expectedsize = size(rawintensities(msm))
    size(variances) == expectedsize || throw(DimensionMismatch(
        "variances must have size $(expectedsize), matching rawintensities(msm)"))

    for variance in variances
        isfinite(variance) || throw(ArgumentError("variances must be finite"))
        variance ≥ 0 || throw(ArgumentError("variances must be nonnegative"))
    end

    nothing
end

function validate_alkane_series_baselines(msm::MassScanMatrix, baselines)
    isnothing(baselines) && return nothing
    baselines isa MassScanMatrix || throw(ArgumentError(
        "baselines must be nothing or a MassScanMatrix matching msm"))
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
