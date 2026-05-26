struct AlkaneSeriesResult{S,V,I,B,C,T,P,A,Q}
    standard::S
    variances::V
    varianceinfo::I
    baselineinfo::B
    channelinfo::C
    traces::T
    seriespath::P
    apexes::A
    stepspectra::Q
end

AlkaneSeriesResult(standard, variances, varianceinfo) =
    AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )

AlkaneSeriesResult(standard, variances, varianceinfo, baselineinfo) =
    AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        baselineinfo,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )

AlkaneSeriesResult(standard, variances, varianceinfo, baselineinfo, channelinfo) =
    AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        baselineinfo,
        channelinfo,
        nothing,
        nothing,
        nothing,
        nothing,
    )

AlkaneSeriesResult(standard, variances, varianceinfo, baselineinfo, channelinfo, traces) =
    AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        baselineinfo,
        channelinfo,
        traces,
        nothing,
        nothing,
        nothing,
    )

AlkaneSeriesResult(
    standard,
    variances,
    varianceinfo,
    baselineinfo,
    channelinfo,
    traces,
    seriespath,
) =
    AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        baselineinfo,
        channelinfo,
        traces,
        seriespath,
        nothing,
        nothing,
    )

AlkaneSeriesResult(
    standard,
    variances,
    varianceinfo,
    baselineinfo,
    channelinfo,
    traces,
    seriespath,
    apexes,
) =
    AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        baselineinfo,
        channelinfo,
        traces,
        seriespath,
        apexes,
        nothing,
    )

"""
    findalkaneseries(msm::MassScanMatrix, variances;
        standard=defaultalkanestandard(), varianceinfo=nothing,
        baselineinfo=nothing, kwargs...) -> AlkaneSeriesResult

Identify an n-alkane series in a preprocessed GC-MS mass-scan matrix.

This is the low-level analysis entry point. It assumes that `msm` is already the signal
to analyze and therefore requires per-cell `variances`. It does not estimate variances
or fit/subtract a baseline. If upstream preprocessing should be retained with the result,
pass it as `baselineinfo`.

Reference-spectrum m/z channels shared with `msm` are identified before downstream
analysis and stored in `result.channelinfo`. Match, molecular-ion, m/z-peak-distance, and
combined evidence traces are stored in `result.traces`. The dynamic-programming series
path and its candidates are stored in `result.seriespath`, and variance-weighted apex
refinements are stored in `result.apexes`. If requested, full-grid ladder step mass
spectra are stored in `result.stepspectra`. Set `stepspectracarbons` to `:all`, one
carbon number, or a collection such as `10:20` to choose which ladder step spectra are
extracted. The `msm` m/z axis must be integer binned.
"""
function findalkaneseries(
    msm::MassScanMatrix,
    variances;
    standard=defaultalkanestandard(),
    varianceinfo=nothing,
    baselineinfo=nothing,

    carbonrange=8:40,

    mzpeakvariancefloor=1.0,
    mzpeakzmin=4.0,
    mzpeakeps=1e-12,
    mzretentionkwargs=nothing,

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
    spacingweight=1.0,
    gapincreaseweight=0.25,

    apexscanwindow=2,
    apexweighting=:variance,
    maxapexoffsetscans=1.0,
    apexlogfloorfraction=1e-3,
    apexvariancefloor=1.0,

    extractstepspectra=false,
    stepspectracarbons=:all,
    stepspectrascanwindow=2,
    stepspectraweighting=:variance,
    stepspectraallownegative=true,
)

    # Validate inputs
    isnothing(standard) && throw(ArgumentError(
        "findalkaneseries requires an alkane standard"))
    apexweighting == :variance || throw(ArgumentError(
        "apexweighting is fixed to :variance"))
    stepspectraweighting == :variance || throw(ArgumentError(
        "stepspectraweighting is fixed to :variance"))

    validate_alkane_series_variances(msm, variances)
    channelinfo = alkane_mz_channels(
        msm;
        standard=standard,
        carbonrange=carbonrange,
        minrelativeintensity=minrelativeintensity,
    )
    traces = alkane_series_traces(
        msm,
        variances,
        channelinfo.carbonrange;
        standard=standard,
        channelinfo=channelinfo,
        minrelativeintensity=minrelativeintensity,
        molecularionwindow=molecularionwindow,
        spectralpower=spectralpower,
        smoothing=smoothing,
        includedistance=includedistance,
        mzpeakzmin=mzpeakzmin,
        mzpeakeps=mzpeakeps,
        mzpeakvariancefloor=mzpeakvariancefloor,
        mzretentionkwargs=mzretentionkwargs,
    )
    seriespath = alkaneseriespath(
        traces.evidence;
        minsteps=minsteps,
        stepreward=stepreward,
        minrelativepeakheight=minrelativepeakheight,
        minpeakheight=minpeakheight,
        maxcandidatespertrace=maxcandidatespertrace,
        spacingweight=spacingweight,
        gapincreaseweight=gapincreaseweight,
    )
    apexes = refine_alkane_series_apexes(
        msm,
        variances,
        seriespath;
        channelinfo=channelinfo,
        scanwindow=apexscanwindow,
        maxapexoffsetscans=maxapexoffsetscans,
        mzretentionkwargs=mzretentionkwargs,
        logfloorfraction=apexlogfloorfraction,
        variancefloor=apexvariancefloor,
    )
    stepspectra = extractstepspectra ?
        alkane_series_step_spectra(
            msm,
            variances,
            apexes;
            scanwindow=stepspectrascanwindow,
            carbons=stepspectracarbons,
            mzretentionkwargs=mzretentionkwargs,
            variancefloor=apexvariancefloor,
            allownegative=stepspectraallownegative,
        ) :
        nothing

    AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        baselineinfo,
        channelinfo,
        traces,
        seriespath,
        apexes,
        stepspectra,
    )
end

function alkane_mz_channels(
    msm::MassScanMatrix;
    standard=defaultalkanestandard(),
    carbonrange=8:40,
    minrelativeintensity=0.05,
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
        variances=nothing, subtractbaseline=true, kwargs...) -> AlkaneSeriesResult

Convenience wrapper for raw-count GC-MS data.

This wrapper handles count-based variance estimation and optional baseline subtraction,
then delegates to [`findalkaneseries`](@ref) for the downstream analysis. By default it
uses [`countvariances`](@ref) for variances and [`arpls`](@ref) for the baseline. If a
different variance model or baseline method is needed, run that preprocessing explicitly
and call [`findalkaneseries`](@ref) with the prepared signal and variances.

When a baseline is subtracted, the fitted baseline matrix and baseline settings are stored
in `result.baselineinfo`. The subtracted matrix is available as
`result.baselineinfo.baselines`.

If `subtractbaseline=false`, `variances` must be supplied because count-based variances
cannot be reliably inferred from caller-preprocessed signals.
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

    baselineλ=1e7,
    subtractbaseline::Bool=true,
    baselinenonnegative=true,
    baselinepeakthreshold=6.0,
    baselinepeakslope=0.5,
    baselinezerothreshold=0.5,

    kwargs...,
)

    isnothing(standard) && throw(ArgumentError(
        "findalkanes requires an alkane standard"))

    if !subtractbaseline && isnothing(variances)
        throw(ArgumentError(
            "variances must be provided when subtractbaseline=false"))
    end

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

    findalkaneseries(
        signal,
        σ²;
        standard=standard,
        varianceinfo=varianceinfo,
        baselineinfo=baselineinfo,
        kwargs...,
    )
end

function extract_alkane_series_variances(varianceestimate)
    if hasproperty(varianceestimate, :variances)
        return getproperty(varianceestimate, :variances), varianceestimate
    end

    varianceestimate, nothing
end

function alkane_standard_spectra(standard)
    spectra = if standard isa AbstractVector{<:AbstractMassSpectrum}
        standard
    elseif hasproperty(standard, :spectra)
        getproperty(standard, :spectra)
    else
        throw(ArgumentError(
            "alkane standard must provide reference spectra in a `spectra` field"))
    end

    spectra isa AbstractVector{<:AbstractMassSpectrum} || throw(ArgumentError(
        "standard.spectra must be a vector of mass spectra"))
    isempty(spectra) && throw(ArgumentError("standard.spectra must not be empty"))

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
    isfinite(minrelativeintensity) && 0 <= minrelativeintensity < 1 || throw(
        ArgumentError("minrelativeintensity must be finite and in [0, 1)"))

    nothing
end

function alkane_reference_channel_match(
    spectrum::AbstractMassSpectrum,
    grid_index_by_mz::AbstractDict{Int,<:Integer},
    grid_mzs::AbstractVector;
    minrelativeintensity::Real,
)
    spectrum_mzs = integer_mz_values(mzvalues(spectrum), "reference spectrum m/z values")
    spectrum_intensities = Float64.(intensities(spectrum))
    all(isfinite, spectrum_intensities) || throw(ArgumentError(
        "reference spectrum contains nonfinite intensities"))
    all(>=(0), spectrum_intensities) || throw(ArgumentError(
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
        mzindices=mz_indices,
        mzvalues=collect(grid_mzs[mz_indices]),
        referenceintensities=[reference_by_index[index] for index in mz_indices],
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
