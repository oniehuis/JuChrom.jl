function alkanemolecularioninfo(
    msm::MassScanMatrix,
    variances,
    abundanceinfo;
    ionwindow::Integer=1,
    stepmass::Integer=14,
    variancefloor::Real=1.0,
    centerzmin::Real=1.645,
    isolationzmin::Real=1.645
)
    contrasts = alkanemolecularionwindowcontrasts(
        msm,
        variances,
        abundanceinfo.abundances,
        abundanceinfo.windows;
        ionwindow=ionwindow,
        stepmass=stepmass,
        variancefloor=variancefloor,
        centerzmin=centerzmin,
        isolationzmin=isolationzmin
    )
    zscorevectors = alkanemolecularionzscorevectors(
        abundanceinfo.abundances,
        contrasts
    )

    (
        contrasts=contrasts,
        zscorevectors=zscorevectors,
        settings=(
            ionwindow=Int(ionwindow),
            stepmass=Int(stepmass),
            variancefloor=Float64(variancefloor),
            centerzmin=Float64(centerzmin),
            isolationzmin=Float64(isolationzmin)
        )
    )
end

function alkanemolecularionwindowcontrasts(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundances::AbstractDict,
    windows::AbstractDict;
    ionwindow::Integer=1,
    stepmass::Integer=14,
    variancefloor::Real=1.0,
    centerzmin::Real=1.645,
    isolationzmin::Real=1.645
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
    grid_mzs = mzvalues(msm)
    contrasts = Dict{Int, Vector{NamedTuple}}()

    for (carbon, stepwindows) in pairs(windows)
        step = Int(carbon)
        haskey(abundances, step) || throw(ArgumentError(
            "abundances does not contain a vector for C$(step)"))

        molecularion = alkane_molecular_ion(step; stepmass=Int(stepmass))
        centerindices = alkane_molecular_ion_group_indices(
            mzbins,
            molecularion,
            Int(ionwindow)
        )
        lowerindices = alkane_mz_window_indices(
            mzbins,
            molecularion - Int(stepmass),
            Int(ionwindow)
        )
        upperindices = alkane_mz_window_indices(
            mzbins,
            molecularion + Int(stepmass),
            Int(ionwindow)
        )
        if !alkane_molecular_ion_group_is_annotatable(
                mzbins,
                molecularion,
                centerindices
            )
            contrasts[step] = NamedTuple[]
            continue
        end

        abundance = alkane_abundance_values(abundances[step], step)
        length(abundance) == size(X, 1) || throw(DimensionMismatch(
            "abundance vector for C$(step) must have length $(size(X, 1))"))

        stepcontrasts = NamedTuple[]
        for window in stepwindows
            alkane_abundance_window_has_positive_peak_model(abundance, window) ||
                continue
            push!(
                stepcontrasts,
                alkane_molecular_ion_window_contrast(
                    X,
                    variances,
                    grid_mzs,
                    abundance,
                    window,
                    molecularion,
                    centerindices,
                    lowerindices,
                    upperindices;
                    variancefloor=Float64(variancefloor),
                    centerzmin=Float64(centerzmin),
                    isolationzmin=Float64(isolationzmin)
                )
            )
        end
        contrasts[step] = stepcontrasts
    end

    contrasts
end

function alkanemolecularionzscorevectors(
    abundances::AbstractDict,
    contrasts::AbstractDict;
    fillvalue::Real=0.0
)
    isfinite(fillvalue) || throw(ArgumentError("fillvalue must be finite"))

    zscorevectors = Dict{Int, Vector{Float64}}()
    for (carbon, abundance) in pairs(abundances)
        step = Int(carbon)
        abundancevalues = alkane_abundance_values(abundance, step)
        zscores = fill(Float64(fillvalue), length(abundancevalues))
        stepcontrasts = get(contrasts, step, NamedTuple[])

        for contrast in stepcontrasts
            1 <= contrast.leftindex <= contrast.rightindex <= length(abundancevalues) ||
                throw(DimensionMismatch(
                    "window indices for C$(step) must fit abundance length"))
            score = hasproperty(contrast, :molecularionscore) ?
                Float64(contrast.molecularionscore) :
                Float64(contrast.z)
            z = isfinite(score) ? max(score, Float64(fillvalue)) : Float64(fillvalue)
            for index in contrast.leftindex:contrast.rightindex
                zscores[index] = max(zscores[index], z)
            end
        end

        zscorevectors[step] = zscores
    end

    zscorevectors
end

function validate_alkane_molecular_ion_settings(
    ionwindow,
    stepmass,
    variancefloor,
    centerzmin=1.645,
    isolationzmin=1.645
)
    ionwindow isa Integer || throw(ArgumentError("ionwindow must be an integer"))
    ionwindow >= 0 || throw(ArgumentError("ionwindow must be nonnegative"))
    stepmass isa Integer || throw(ArgumentError("stepmass must be an integer"))
    stepmass > 0 || throw(ArgumentError("stepmass must be positive"))
    validate_alkane_abundance_variancefloor(variancefloor)
    isfinite(centerzmin) && centerzmin >= 0 || throw(ArgumentError(
        "centerzmin must be finite and nonnegative"))
    isfinite(isolationzmin) && isolationzmin >= 0 || throw(ArgumentError(
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

    mzbins = alkane_mz_bins(msm)
    annotatable = alkane_ladder_annotatable_molecular_ion_steps(
        mzbins;
        carbonrange=carbonrange,
        ionwindow=ionwindow,
        stepmass=stepmass
    )
    length(annotatable) >= Int(minsteps) && return annotatable

    throw(ArgumentError(
        "alkane ladder identification requires molecular-ion channels for at least " *
        "$(Int(minsteps)) carbon steps, but only $(length(annotatable)) are " *
        "annotatable on the measured m/z grid. The data may be SIM data or may not " *
        "cover the alkane molecular-ion m/z values. Annotatable steps: " *
        string(annotatable)))
end

function alkane_ladder_annotatable_molecular_ion_steps(
    mzbins::AbstractVector{<:Integer};
    carbonrange=8:40,
    ionwindow::Integer=1,
    stepmass::Integer=14
)
    [
        Int(carbon) for carbon in carbonrange
        if alkane_molecular_ion_group_is_annotatable(
            mzbins,
            alkane_molecular_ion(Int(carbon); stepmass=Int(stepmass)),
            alkane_molecular_ion_group_indices(
                mzbins,
                alkane_molecular_ion(Int(carbon); stepmass=Int(stepmass)),
                Int(ionwindow)
            )
        )
    ]
end

function validate_alkane_molecular_ion_channel_precheck_settings(
    carbonrange,
    ionwindow,
    stepmass,
    minsteps
)
    carbons = collect(carbonrange)
    isempty(carbons) && throw(ArgumentError("carbonrange must not be empty"))
    all(carbon -> carbon isa Integer, carbons) || throw(ArgumentError(
        "carbonrange must contain integer carbon numbers"))
    ionwindow isa Integer || throw(ArgumentError("ionwindow must be an integer"))
    ionwindow >= 0 || throw(ArgumentError("ionwindow must be nonnegative"))
    stepmass isa Integer || throw(ArgumentError("stepmass must be an integer"))
    stepmass > 0 || throw(ArgumentError("stepmass must be positive"))
    minsteps isa Integer || throw(ArgumentError("minsteps must be an integer"))
    minsteps >= 1 || throw(ArgumentError("minsteps must be at least 1"))

    nothing
end

alkane_molecular_ion(carbon::Integer; stepmass::Integer=14) =
    Int(stepmass) * Int(carbon) + 2

function alkane_mz_window_indices(
    mzs::AbstractVector{<:Integer},
    center::Integer,
    halfwidth::Integer
)
    halfwidth >= 0 || throw(ArgumentError("halfwidth must be nonnegative"))

    findall(mz -> abs(Int(mz) - Int(center)) <= halfwidth, mzs)
end

function alkane_molecular_ion_group_indices(
    mzs::AbstractVector{<:Integer},
    molecularion::Integer,
    forwardwindow::Integer
)
    forwardwindow >= 0 || throw(ArgumentError("forwardwindow must be nonnegative"))
    lower = Int(molecularion)
    upper = Int(molecularion) + Int(forwardwindow)

    findall(mz -> lower <= Int(mz) <= upper, mzs)
end

function alkane_molecular_ion_group_is_annotatable(
    mzs::AbstractVector{<:Integer},
    molecularion::Integer,
    centerindices::AbstractVector{<:Integer}
)
    mi = Int(molecularion)
    hasmi = any(index -> Int(mzs[index]) == mi, centerindices)
    hasmi || return false

    gridhasmiplusone = any(mz -> Int(mz) == mi + 1, mzs)
    hasmiplusone = any(index -> Int(mzs[index]) == mi + 1, centerindices)

    hasmiplusone || !gridhasmiplusone
end

function alkane_abundance_window_has_positive_peak_model(
    abundance::AbstractVector{<:Real},
    window
)
    1 <= window.leftindex <= window.rightindex <= length(abundance) ||
        throw(DimensionMismatch(
            "window indices for C$(window.ladderstep) must fit abundance length"))

    maximum(Float64.(abundance[window.leftindex:window.rightindex])) > 0
end

function alkane_molecular_ion_window_contrast(
    X::AbstractMatrix{<:Real},
    variances::AbstractMatrix{<:Real},
    grid_mzs::AbstractVector,
    abundance::AbstractVector{<:Real},
    window,
    molecularion::Integer,
    centerindices::AbstractVector{<:Integer},
    lowerindices::AbstractVector{<:Integer},
    upperindices::AbstractVector{<:Integer};
    variancefloor::Float64,
    centerzmin::Float64,
    isolationzmin::Float64
)
    scanindices = window.leftindex:window.rightindex
    peakmodel = Float64.(abundance[scanindices])
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
        peakmodel;
        variancefloor=variancefloor
    )
    lowerused = !isempty(lowerindices)
    upperused = !isempty(upperindices)
    lowerfit = lowerused ?
        alkane_fitted_ion_group_abundance(
            X,
            variances,
            scanindices,
            lowerindices,
            peakmodel;
            variancefloor=variancefloor
        ) :
        (abundance=0.0, variance=0.0)
    upperfit = upperused ?
        alkane_fitted_ion_group_abundance(
            X,
            variances,
            scanindices,
            upperindices,
            peakmodel;
            variancefloor=variancefloor
        ) :
        (abundance=0.0, variance=0.0)

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
        centerz >= centerzmin &&
        isfinite(isolationz) &&
        isolationz >= isolationzmin &&
        isfinite(z) ?
        max(z, 0.0) :
        0.0

    (
        ladderstep=Int(window.ladderstep),
        leftindex=Int(window.leftindex),
        apexindex=modelapexindex,
        rightindex=Int(window.rightindex),
        scanindices=collect(scanindices),
        peakmodel=collect(peakmodel),
        apexabundance=hasproperty(window, :apexabundance) ?
            Float64(window.apexabundance) :
            NaN,
        molecularion=Int(molecularion),
        centerindices=collect(centerindices),
        lowerindices=collect(lowerindices),
        upperindices=collect(upperindices),
        centerions=collect(grid_mzs[centerindices]),
        lowerions=collect(grid_mzs[lowerindices]),
        upperions=collect(grid_mzs[upperindices]),
        lowercontrolused=lowerused,
        uppercontrolused=upperused,
        centerabundance=centerfit.abundance,
        lowerabundance=lowerfit.abundance,
        upperabundance=upperfit.abundance,
        lowercontrolpenalty=lowerpenalty,
        uppercontrolpenalty=upperpenalty,
        controlpenalty=controlpenalty,
        controlvariance=lowerfit.variance + upperfit.variance,
        controlstderr=sqrt(lowerfit.variance + upperfit.variance),
        controlcount=Int(lowerused) + Int(upperused),
        signedcontrastabundance=signedcontrast,
        contrastabundance=contrast,
        contrastvariance=contrastvariance,
        contraststderr=contraststderr,
        centerstderr=centerstderr,
        centerz=centerz,
        lowercontrolvetoz=lowercontrolvetoz,
        uppercontrolvetoz=uppercontrolvetoz,
        centervslowerz=centervslowerz,
        centervsupperz=centervsupperz,
        isolationz=isolationz,
        molecularionscore=molecularionscore,
        z=z
    )
end

function alkane_fitted_ion_group_abundance(
    X::AbstractMatrix{<:Real},
    variances::AbstractMatrix{<:Real},
    scanindices,
    mzindices::AbstractVector{<:Integer},
    peakmodel::AbstractVector{<:Real};
    variancefloor::Float64
)
    abundance = 0.0
    abundancevariance = 0.0

    for mzindex in mzindices
        fit = alkane_fitted_ion_abundance(
            X,
            variances,
            scanindices,
            mzindex,
            peakmodel;
            variancefloor=variancefloor
        )
        abundance += fit.abundance
        abundancevariance += fit.variance
    end

    (abundance=abundance, variance=abundancevariance)
end

function alkane_fitted_ion_abundance(
    X::AbstractMatrix{<:Real},
    variances::AbstractMatrix{<:Real},
    scanindices,
    mzindex::Integer,
    peakmodel::AbstractVector{<:Real};
    variancefloor::Float64
)
    size(variances) == size(X) || throw(DimensionMismatch(
        "variances must have size $(size(X)), matching X"))
    length(scanindices) == length(peakmodel) || throw(DimensionMismatch(
        "scanindices and peakmodel must have the same length"))

    numerator = 0.0
    denominator = 0.0
    for (modelindex, scanindex) in pairs(scanindices)
        weight = inv(max(Float64(variances[scanindex, mzindex]), variancefloor))
        modelvalue = Float64(peakmodel[modelindex])
        numerator += weight * modelvalue * Float64(X[scanindex, mzindex])
        denominator += weight * abs2(modelvalue)
    end

    denominator > 0 || return (abundance=0.0, variance=0.0)

    (abundance=numerator / denominator, variance=inv(denominator))
end

function alkane_molecular_ion_control_veto_z(controlfit, centerfit)
    variance = controlfit.variance + centerfit.variance
    variance > 0 || return NaN

    (controlfit.abundance - centerfit.abundance) / sqrt(variance)
end

function alkane_molecular_ion_center_vs_control_z(centerfit, controlfit)
    variance = centerfit.variance + controlfit.variance
    variance > 0 || return NaN

    (centerfit.abundance - controlfit.abundance) / sqrt(variance)
end

function alkane_molecular_ion_isolation_z(
    centervslowerz::Real,
    centervsupperz::Real
)
    values = Float64[]
    isfinite(centervslowerz) && push!(values, Float64(centervslowerz))
    isfinite(centervsupperz) && push!(values, Float64(centervsupperz))
    isempty(values) && return NaN

    minimum(values)
end
