"""
    AlkaneMolecularIonContrast

Molecular-ion contrast result for one alkane abundance window.
"""
struct AlkaneMolecularIonContrast{T<:Real}
    ladderstep::Int
    leftindex::Int
    apexindex::Int
    rightindex::Int
    scanindices::Vector{Int}
    peakmodel::Vector{Float64}
    apexabundance::Float64
    molecularion::Int
    centerindices::Vector{Int}
    lowerindices::Vector{Int}
    upperindices::Vector{Int}
    centerions::Vector{T}
    lowerions::Vector{T}
    upperions::Vector{T}
    lowercontrolused::Bool
    uppercontrolused::Bool
    centerabundance::Float64
    lowerabundance::Float64
    upperabundance::Float64
    lowercontrolpenalty::Float64
    uppercontrolpenalty::Float64
    controlpenalty::Float64
    controlvariance::Float64
    controlstderr::Float64
    controlcount::Int
    signedcontrastabundance::Float64
    contrastabundance::Float64
    contrastvariance::Float64
    contraststderr::Float64
    centerstderr::Float64
    centerz::Float64
    lowercontrolvetoz::Float64
    uppercontrolvetoz::Float64
    centervslowerz::Float64
    centervsupperz::Float64
    isolationz::Float64
    molecularionscore::Float64
    z::Float64
end

"""
    AlkaneMolecularIonInfo

Result of the molecular-ion evidence stage.
"""
struct AlkaneMolecularIonInfo{
    T1<:AbstractDict{Int,<:AbstractVector{<:AlkaneMolecularIonContrast}},
    T2<:NamedTuple
} <: AbstractAlkaneMolecularIonInfo
    contrasts::T1
    zscorevectors::Dict{Int, Vector{Float64}}
    settings::T2
end

struct AlkaneFittedIonAbundance
    abundance::Float64
    variance::Float64
end

struct AlkaneMolecularIonContrastEntry{T<:Real}
    carbon::Int
    contrasts::Vector{AlkaneMolecularIonContrast{T}}
end

struct AlkaneMolecularIonZScoreEntry
    carbon::Int
    zscores::Vector{Float64}
end

function alkane_molecular_ion_info(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    ionwindow::Integer,
    stepmass::Integer,
    variancefloor::Real,
    centerzmin::Real,
    isolationzmin::Real
)
    contrasts = alkane_molecular_ion_window_contrasts(
        msm,
        variances,
        abundanceinfo.abundances,
        abundanceinfo.windows,
        ionwindow,
        stepmass,
        variancefloor,
        centerzmin,
        isolationzmin
    )

    zscorevectors = alkane_molecular_ion_zscore_vectors(
        abundanceinfo.abundances,
        contrasts,
        0.0
    )

    AlkaneMolecularIonInfo(
        contrasts,
        zscorevectors,
        (
            ionwindow=Int(ionwindow),
            stepmass=Int(stepmass),
            variancefloor=Float64(variancefloor),
            centerzmin=Float64(centerzmin),
            isolationzmin=Float64(isolationzmin)
        )
    )
end

function alkane_molecular_ion_window_contrasts(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundances::AbstractDict{Int, <:AbstractVector{<:Real}},
    windows::AbstractDict{Int, <:AbstractVector{AlkaneAbundanceWindow}},
    ionwindow::Integer,
    stepmass::Integer,
    variancefloor::Real,
    centerzmin::Real,
    isolationzmin::Real
)
    validate_alkane_molecular_ion_settings(
        ionwindow,
        stepmass,
        variancefloor,
        centerzmin,
        isolationzmin
    )
    validate_alkane_series_variances(msm, variances)

    X = rawintensities(msm)
    mzbins = alkane_mz_bins(msm)
    grid_mzs = rawmzvalues(msm)
    contrasttype = AlkaneMolecularIonContrast{eltype(grid_mzs)}
    windowpairs = collect(pairs(windows))
    for (carbon, stepwindows) in windowpairs
        step = carbon
        haskey(abundances, step) || throw(ArgumentError(
            "abundances does not contain a vector for C$(step)"))
        molecularion = alkane_molecular_ion(step, stepmass)
        centerindices = alkane_molecular_ion_group_indices(
            mzbins,
            molecularion,
            ionwindow
        )
        alkane_molecular_ion_group_is_annotatable(
            mzbins,
            molecularion,
            centerindices
        ) || continue

        abundance = alkane_abundance_values(abundances[step], step)
        length(abundance) == size(X, 1) || throw(DimensionMismatch(
            "abundance vector for C$(step) must have length $(size(X, 1))"))
        for window in stepwindows
            1 ≤ window.leftindex ≤ window.rightindex ≤ length(abundance) ||
                throw(DimensionMismatch(
                    "window indices for C$(step) must fit abundance length"))
        end
    end

    results = Vector{AlkaneMolecularIonContrastEntry{eltype(grid_mzs)}}(
        undef,
        length(windowpairs)
    )

    Base.Threads.@threads for index in eachindex(windowpairs)
        carbon = first(windowpairs[index])
        stepwindows = last(windowpairs[index])
        step = carbon
        haskey(abundances, step) || throw(ArgumentError(
            "abundances does not contain a vector for C$(step)"))

        molecularion = alkane_molecular_ion(step, stepmass)
        centerindices = alkane_molecular_ion_group_indices(mzbins, molecularion, ionwindow)
        lowerindices = alkane_mz_window_indices(mzbins, molecularion - stepmass, ionwindow)
        upperindices = alkane_mz_window_indices(mzbins, molecularion + stepmass, ionwindow)

        if !alkane_molecular_ion_group_is_annotatable(mzbins, molecularion, centerindices)
            results[index] = AlkaneMolecularIonContrastEntry(step, contrasttype[])
            continue
        end

        abundance = alkane_abundance_values(abundances[step], step)
        length(abundance) == size(X, 1) || throw(DimensionMismatch(
            "abundance vector for C$(step) must have length $(size(X, 1))"))

        stepcontrasts = contrasttype[]
        for window in stepwindows
            alkane_abundance_window_has_positive_peak_model(abundance, window) ||
                continue
            push!(stepcontrasts, alkane_molecular_ion_window_contrast(
                    X,
                    variances,
                    grid_mzs,
                    abundance,
                    window,
                    molecularion,
                    centerindices,
                    lowerindices,
                    upperindices,
                    variancefloor,
                    centerzmin,
                    isolationzmin
                )
            )
        end
        results[index] = AlkaneMolecularIonContrastEntry(step, stepcontrasts)
    end

    contrasts = Dict{Int, Vector{AlkaneMolecularIonContrast{eltype(grid_mzs)}}}()
    for result in results
        contrasts[result.carbon] = result.contrasts
    end

    contrasts
end

function alkane_molecular_ion_zscore_vectors(
    abundances::AbstractDict{Int, <:AbstractVector{<:Real}},
    contrasts::AbstractDict{Int, <:AbstractVector{<:AlkaneMolecularIonContrast}},
    fillvalue::Real
)
    isfinite(fillvalue) || throw(ArgumentError("fillvalue must be finite"))
    fill = Float64(fillvalue)

    abundancepairs = collect(pairs(abundances))
    for (carbon, abundance) in abundancepairs
        step = carbon
        abundancevalues = alkane_abundance_values(abundance, step)
        stepcontrasts = get(contrasts, step, AlkaneMolecularIonContrast[])
        for contrast in stepcontrasts
            1 ≤ contrast.leftindex ≤ contrast.rightindex ≤ length(abundancevalues) ||
                throw(DimensionMismatch(
                    "window indices for C$(step) must fit abundance length"))
        end
    end

    results = Vector{AlkaneMolecularIonZScoreEntry}(undef, length(abundancepairs))

    Base.Threads.@threads for pairindex in eachindex(abundancepairs)
        carbon = first(abundancepairs[pairindex])
        abundance = last(abundancepairs[pairindex])
        step = carbon
        abundancevalues = alkane_abundance_values(abundance, step)
        zscores = Base.fill(fill, length(abundancevalues))
        stepcontrasts = get(contrasts, step, AlkaneMolecularIonContrast[])

        for contrast in stepcontrasts
            1 ≤ contrast.leftindex ≤ contrast.rightindex ≤ length(abundancevalues) ||
                throw(DimensionMismatch(
                    "window indices for C$(step) must fit abundance length"))
            score = contrast.molecularionscore
            z = isfinite(score) ? max(score, fill) : fill
            for index in contrast.leftindex:contrast.rightindex
                zscores[index] = max(zscores[index], z)
            end
        end

        results[pairindex] = AlkaneMolecularIonZScoreEntry(step, zscores)
    end

    zscorevectors = Dict{Int, Vector{Float64}}()
    for result in results
        zscorevectors[result.carbon] = result.zscores
    end

    zscorevectors
end

function validate_alkane_molecular_ion_settings(
    ionwindow::Integer,
    stepmass::Integer,
    variancefloor::Real,
    centerzmin::Real,
    isolationzmin::Real
)
    ionwindow ≥ 0 || throw(ArgumentError("ionwindow must be nonnegative"))
    stepmass > 0 || throw(ArgumentError("stepmass must be positive"))
    validate_alkane_abundance_variancefloor(variancefloor)
    isfinite(centerzmin) && centerzmin ≥ 0 || throw(ArgumentError(
        "centerzmin must be finite and nonnegative"))
    isfinite(isolationzmin) && isolationzmin ≥ 0 || throw(ArgumentError(
        "isolationzmin must be finite and nonnegative"))

    nothing
end

function alkane_ladder_require_molecular_ion_channels(
    msm::MassScanMatrix;
    carbonrange=8:40,
    ionwindow::Integer=1,
    stepmass::Integer=14,
    minsteps::Integer=5
)
    validate_alkane_molecular_ion_channel_precheck_settings(
        carbonrange,
        ionwindow,
        stepmass,
        minsteps
    )

    annotatable = alkane_ladder_annotatable_molecular_ion_steps(
        alkane_mz_bins(msm),
        carbonrange,
        ionwindow,
        stepmass
    )
    length(annotatable) ≥ minsteps && return annotatable

    throw(ArgumentError(
        "alkane ladder identification requires molecular-ion channels for at least " *
        "$(minsteps) carbon steps, but only $(length(annotatable)) are " *
        "annotatable on the measured m/z grid. The data may be SIM data or may not " *
        "cover the alkane molecular-ion m/z values. Annotatable steps: " *
        string(annotatable)))
end

function alkane_ladder_annotatable_molecular_ion_steps(
    mzbins::AbstractVector{<:Integer},
    carbonrange::AbstractVector{Int},
    ionwindow::Integer,
    stepmass::Integer
)
    annotatable = Int[]
    for carbon in carbonrange
        molecularion = alkane_molecular_ion(carbon, stepmass)
        centerindices = alkane_molecular_ion_group_indices(
            mzbins,
            molecularion,
            ionwindow
        )
        if alkane_molecular_ion_group_is_annotatable(
            mzbins,
            molecularion,
            centerindices
        )
            push!(annotatable, carbon)
        end
    end

    annotatable
end

function validate_alkane_molecular_ion_channel_precheck_settings(
    carbonrange,
    ionwindow::Integer,
    stepmass::Integer,
    minsteps::Integer
)
    carbons = collect(carbonrange)
    isempty(carbons) && throw(ArgumentError("carbonrange must not be empty"))
    all(carbon -> carbon isa Integer, carbons) || throw(ArgumentError(
        "carbonrange must contain integer carbon numbers"))
    ionwindow ≥ 0 || throw(ArgumentError("ionwindow must be nonnegative"))
    stepmass > 0 || throw(ArgumentError("stepmass must be positive"))
    minsteps ≥ 1 || throw(ArgumentError("minsteps must be at least 1"))

    nothing
end

alkane_molecular_ion(carbon::Integer) = alkane_molecular_ion(carbon, 14)

alkane_molecular_ion(carbon::Integer, stepmass::Integer)::Int = stepmass * carbon + 2

function alkane_mz_window_indices(
    mzs::AbstractVector{<:Integer},
    center::Integer,
    halfwidth::Integer
)
    halfwidth ≥ 0 || throw(ArgumentError("halfwidth must be nonnegative"))

    findall(mz -> abs(mz - center) ≤ halfwidth, mzs)
end

function alkane_molecular_ion_group_indices(
    mzs::AbstractVector{<:Integer},
    molecularion::Integer,
    forwardwindow::Integer
)
    forwardwindow ≥ 0 || throw(ArgumentError("forwardwindow must be nonnegative"))
    lower = molecularion
    upper = molecularion + forwardwindow

    findall(mz -> lower ≤ mz ≤ upper, mzs)
end

function alkane_molecular_ion_group_is_annotatable(
    mzs::AbstractVector{<:Integer},
    molecularion::Integer,
    centerindices::AbstractVector{<:Integer}
)
    mi = molecularion
    hasmi = any(index -> mzs[index] == mi, centerindices)
    hasmi || return false

    gridhasmiplusone = any(mz -> mz == mi + 1, mzs)
    hasmiplusone = any(index -> mzs[index] == mi + 1, centerindices)

    hasmiplusone || !gridhasmiplusone
end

function alkane_abundance_window_has_positive_peak_model(
    abundance::AbstractVector{<:Real},
    window::AlkaneAbundanceWindow
)
    1 ≤ window.leftindex ≤ window.rightindex ≤ length(abundance) ||
        throw(DimensionMismatch(
            "window indices for C$(window.ladderstep) must fit abundance length"))

    any(index -> abundance[index] > 0, window.leftindex:window.rightindex)
end

function alkane_molecular_ion_window_contrast(
    X::AbstractMatrix{<:Real},
    variances::AbstractMatrix{<:Real},
    grid_mzs::AbstractVector{<:Real},
    abundance::AbstractVector{<:Real},
    window::AlkaneAbundanceWindow,
    molecularion::Integer,
    centerindices::AbstractVector{<:Integer},
    lowerindices::AbstractVector{<:Integer},
    upperindices::AbstractVector{<:Integer},
    variancefloor::Real,
    centerzmin::Real,
    isolationzmin::Real
)
    scanindices = window.leftindex:window.rightindex
    peakmodel = abundance[scanindices]
    peakmodel .= max.(peakmodel, 0.0)
    peakheight, peakoffset = findmax(peakmodel)
    peakheight > 0 || throw(ArgumentError(
        "abundance window for C$(window.ladderstep) has no positive peak model"))
    modelapexindex = first(scanindices) + peakoffset - 1
    peakmodel ./= peakheight

    centerfit = alkane_fitted_ion_group_abundance(
        X,
        variances,
        scanindices,
        centerindices,
        peakmodel,
        variancefloor
    )
    lowerused = !isempty(lowerindices)
    upperused = !isempty(upperindices)
    lowerfit = lowerused ?
        alkane_fitted_ion_group_abundance(
            X,
            variances,
            scanindices,
            lowerindices,
            peakmodel,
            variancefloor
        ) :
        AlkaneFittedIonAbundance(0.0, 0.0)
    upperfit = upperused ?
        alkane_fitted_ion_group_abundance(
            X,
            variances,
            scanindices,
            upperindices,
            peakmodel,
            variancefloor
        ) :
        AlkaneFittedIonAbundance(0.0, 0.0)

    lowerpenalty = max(lowerfit.abundance, 0.0)
    upperpenalty = max(upperfit.abundance, 0.0)
    controlpenalty = lowerpenalty + upperpenalty
    signedcontrast = centerfit.abundance - lowerfit.abundance - upperfit.abundance
    contrast = centerfit.abundance - controlpenalty
    contrastvariance = centerfit.variance + lowerfit.variance + upperfit.variance
    contraststderr = sqrt(contrastvariance)
    z = contraststderr > 0 ? contrast / contraststderr : NaN
    centerstderr = sqrt(centerfit.variance)
    centerz = centerstderr > 0 ? centerfit.abundance / centerstderr : NaN
    lowercontrolvetoz = alkane_molecular_ion_control_veto_z(lowerfit, centerfit)
    uppercontrolvetoz = alkane_molecular_ion_control_veto_z(upperfit, centerfit)
    centervslowerz = lowerused ?
        alkane_molecular_ion_center_vs_control_z(centerfit, lowerfit) :
        NaN
    centervsupperz = upperused ?
        alkane_molecular_ion_center_vs_control_z(centerfit, upperfit) :
        NaN
    isolationz = alkane_molecular_ion_isolation_z(centervslowerz, centervsupperz)
    molecularionscore = isfinite(centerz) &&
        centerz ≥ centerzmin &&
        isfinite(isolationz) &&
        isolationz ≥ isolationzmin &&
        isfinite(z) ?
        max(z, 0.0) :
        0.0

    AlkaneMolecularIonContrast(
        window.ladderstep,
        window.leftindex,
        modelapexindex,
        window.rightindex,
        collect(scanindices),
        peakmodel,
        window.apexabundance,
        molecularion,
        centerindices,
        lowerindices,
        upperindices,
        grid_mzs[centerindices],
        grid_mzs[lowerindices],
        grid_mzs[upperindices],
        lowerused,
        upperused,
        centerfit.abundance,
        lowerfit.abundance,
        upperfit.abundance,
        lowerpenalty,
        upperpenalty,
        controlpenalty,
        lowerfit.variance + upperfit.variance,
        sqrt(lowerfit.variance + upperfit.variance),
        count(identity, (lowerused, upperused)),
        signedcontrast,
        contrast,
        contrastvariance,
        contraststderr,
        centerstderr,
        centerz,
        lowercontrolvetoz,
        uppercontrolvetoz,
        centervslowerz,
        centervsupperz,
        isolationz,
        molecularionscore,
        z
    )
end

function alkane_fitted_ion_group_abundance(
    X::AbstractMatrix{<:Real},
    variances::AbstractMatrix{<:Real},
    scanindices::AbstractVector{<:Integer},
    mzindices::AbstractVector{<:Integer},
    peakmodel::AbstractVector{<:Real},
    variancefloor::Real
)
    abundance = 0.0
    abundancevariance = 0.0

    for mzindex in mzindices
        fit = alkane_fitted_ion_abundance(
            X,
            variances,
            scanindices,
            mzindex,
            peakmodel,
            variancefloor
        )
        abundance += fit.abundance
        abundancevariance += fit.variance
    end

    AlkaneFittedIonAbundance(abundance, abundancevariance)
end

function alkane_fitted_ion_abundance(
    X::AbstractMatrix{<:Real},
    variances::AbstractMatrix{<:Real},
    scanindices::AbstractVector{<:Integer},
    mzindex::Integer,
    peakmodel::AbstractVector{<:Real},
    variancefloor::Real
)
    size(variances) == size(X) || throw(DimensionMismatch(
        "variances must have size $(size(X)), matching X"))
    length(scanindices) == length(peakmodel) || throw(DimensionMismatch(
        "scanindices and peakmodel must have the same length"))

    numerator = denominator = 0.0
    for (modelindex, scanindex) in pairs(scanindices)
        weight = inv(max(variances[scanindex, mzindex], variancefloor))
        modelvalue = peakmodel[modelindex]
        numerator += weight * modelvalue * X[scanindex, mzindex]
        denominator += weight * abs2(modelvalue)
    end

    denominator > 0 || return AlkaneFittedIonAbundance(0.0, 0.0)

    AlkaneFittedIonAbundance(numerator / denominator, inv(denominator))
end

function alkane_molecular_ion_control_veto_z(
    controlfit::AlkaneFittedIonAbundance,
    centerfit::AlkaneFittedIonAbundance
)
    variance = controlfit.variance + centerfit.variance
    variance > 0 || return NaN

    (controlfit.abundance - centerfit.abundance) / sqrt(variance)
end

function alkane_molecular_ion_center_vs_control_z(
    centerfit::AlkaneFittedIonAbundance,
    controlfit::AlkaneFittedIonAbundance
)
    variance = centerfit.variance + controlfit.variance
    variance > 0 || return NaN

    (centerfit.abundance - controlfit.abundance) / sqrt(variance)
end

function alkane_molecular_ion_isolation_z(
    centervslowerz::Real,
    centervsupperz::Real
)
    values = Float64[]
    isfinite(centervslowerz) && push!(values, centervslowerz)
    isfinite(centervsupperz) && push!(values, centervsupperz)
    isempty(values) && return NaN

    minimum(values)
end
