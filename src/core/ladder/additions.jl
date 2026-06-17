function alkaneladderadditions(
    msm::MassScanMatrix,
    variances,
    abundanceinfo,
    pathinfo,
    apexinfo;
    minradius::Real=5.0,
    radiusfraction::Real=0.15,
    positionsigmafraction::Real=0.05,
    maxextensionsteps::Integer=100,
    standard=defaultalkanestandard(),
    massspectrummatch::Bool=true,
    gapmincosinefloor::Real=0.85,
    gapcosinetolerance::Real=0.03,
    edgemaxanchors::Integer=6,
    edgeminradius::Real=5.0,
    edgeradiusfraction::Real=0.2,
    edgemincosinefloor::Real=0.9,
    edgecosinetolerance::Real=0.03,
    edgecosineanchorcount::Integer=3,
    edgepositionsigmafraction::Real=0.1,
    massspectrumvariancefloor::Real=1.0,
    apexscanwindow::Integer=2,
    apexvariancefloor::Real=1.0,
    apexlogfloorfraction::Real=1e-3,
    apexionexcludemzvalues=DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
    apexionmzvalues=nothing,
    apexionminrelativeintensity::Real=DEFAULT_ALKANE_APEX_ION_MIN_RELATIVE_INTENSITY,
    apexminioncount::Integer=3,
    apexmzretentionkwargs=nothing,
    apexmaxshiftfromguess::Real=3.0,
    carbonrange=nothing
)
    validate_alkane_series_variances(msm, variances)
    validate_alkane_ladder_addition_settings(
        minradius,
        radiusfraction,
        positionsigmafraction,
        maxextensionsteps,
        gapmincosinefloor,
        gapcosinetolerance,
        edgemaxanchors,
        edgeminradius,
        edgeradiusfraction,
        edgemincosinefloor,
        edgecosinetolerance,
        edgecosineanchorcount,
        edgepositionsigmafraction,
        massspectrumvariancefloor
    )
    validate_alkane_ladder_apex_settings(
        msm,
        variances,
        apexscanwindow,
        apexminioncount,
        apexionminrelativeintensity,
        apexvariancefloor,
        apexlogfloorfraction,
        apexmaxshiftfromguess
    )
    hasproperty(abundanceinfo, :abundances) || throw(ArgumentError(
        "abundanceinfo must contain abundances"))
    hasproperty(abundanceinfo, :windows) || throw(ArgumentError(
        "abundanceinfo must contain windows"))
    hasproperty(pathinfo, :path) || throw(ArgumentError("pathinfo must contain path"))
    hasproperty(apexinfo, :apexes) || throw(ArgumentError("apexinfo must contain apexes"))

    resolved_apex_mzretentionkwargs = alkane_ladder_addition_apex_mzretentionkwargs(
        apexinfo,
        apexmzretentionkwargs
    )
    resolved_carbonrange = alkane_ladder_addition_carbonrange(abundanceinfo, carbonrange)
    settings = (
        minradius=Float64(minradius),
        radiusfraction=Float64(radiusfraction),
        positionsigmafraction=Float64(positionsigmafraction),
        maxextensionsteps=Int(maxextensionsteps),
        massspectrummatch=Bool(massspectrummatch),
        gapmincosinefloor=Float64(gapmincosinefloor),
        gapcosinetolerance=Float64(gapcosinetolerance),
        edgemaxanchors=Int(edgemaxanchors),
        edgeminradius=Float64(edgeminradius),
        edgeradiusfraction=Float64(edgeradiusfraction),
        edgemincosinefloor=Float64(edgemincosinefloor),
        edgecosinetolerance=Float64(edgecosinetolerance),
        edgecosineanchorcount=Int(edgecosineanchorcount),
        edgepositionsigmafraction=Float64(edgepositionsigmafraction),
        massspectrumvariancefloor=Float64(massspectrumvariancefloor),
        apexscanwindow=Int(apexscanwindow),
        apexvariancefloor=Float64(apexvariancefloor),
        apexlogfloorfraction=Float64(apexlogfloorfraction),
        apexionexcludemzvalues=apexionexcludemzvalues,
        apexionmzvalues=apexionmzvalues,
        apexionminrelativeintensity=Float64(apexionminrelativeintensity),
        apexminioncount=Int(apexminioncount),
        apexmzretentionkwargs=resolved_apex_mzretentionkwargs,
        apexmaxshiftfromguess=Float64(apexmaxshiftfromguess),
        carbonrange=resolved_carbonrange
    )
    massspectrumcontext = alkane_ladder_mass_spectrum_match_context(
        msm,
        variances,
        standard;
        enabled=settings.massspectrummatch,
        variancefloor=settings.massspectrumvariancefloor
    )
    anchors = alkane_successful_ladder_apexes(apexinfo)
    gapfilled, gapdiagnostics = alkane_single_gap_additions(
        msm,
        variances,
        abundanceinfo,
        anchors,
        scancount(msm);
        settings=settings,
        massspectrumcontext=massspectrumcontext,
        standard=standard
    )
    edgeanchors = NamedTuple[]
    append!(edgeanchors, anchors)
    append!(edgeanchors, (alkane_ladder_addition_anchor(addition) for addition in gapfilled))
    sort!(edgeanchors; by=anchor -> Int(anchor.ladderstep))
    leftextended, leftdiagnostics = alkane_iterative_edge_additions(
        msm,
        variances,
        abundanceinfo,
        edgeanchors,
        :left,
        scancount(msm);
        settings=settings,
        massspectrumcontext=massspectrumcontext,
        standard=standard
    )
    rightextended, rightdiagnostics = alkane_iterative_edge_additions(
        msm,
        variances,
        abundanceinfo,
        edgeanchors,
        :right,
        scancount(msm);
        settings=settings,
        massspectrumcontext=massspectrumcontext,
        standard=standard
    )

    additions = NamedTuple[]
    append!(additions, gapfilled)
    append!(additions, leftextended)
    append!(additions, rightextended)
    sort!(additions; by=addition -> addition.ladderstep)

    (
        status=isempty(additions) ? :empty : :success,
        additions=additions,
        gapfilled=gapfilled,
        leftextended=leftextended,
        rightextended=rightextended,
        diagnostics=(
            gapfilled=gapdiagnostics,
            leftextended=leftdiagnostics,
            rightextended=rightdiagnostics
        ),
        settings=settings
    )
end

function alkane_ladder_addition_apex_mzretentionkwargs(apexinfo, apexmzretentionkwargs)
    !isnothing(apexmzretentionkwargs) && return apexmzretentionkwargs
    if hasproperty(apexinfo, :settings) &&
            hasproperty(apexinfo.settings, :mzretentionkwargs)
        return apexinfo.settings.mzretentionkwargs
    end

    throw(ArgumentError(
        "apexmzretentionkwargs must be provided when apexinfo does not contain " *
        "resolved mzretentionkwargs"))
end

function alkane_ladder_addition_carbonrange(abundanceinfo, carbonrange)
    carbons = if isnothing(carbonrange)
        carbonkeys = Int[]
        if hasproperty(abundanceinfo, :abundances)
            append!(carbonkeys, Int.(collect(keys(abundanceinfo.abundances))))
        end
        if hasproperty(abundanceinfo, :windows)
            append!(carbonkeys, Int.(collect(keys(abundanceinfo.windows))))
        end
        sort!(unique!(carbonkeys))
    else
        collect(carbonrange)
    end

    isempty(carbons) && throw(ArgumentError("carbonrange must not be empty"))
    all(carbon -> carbon isa Integer, carbons) || throw(ArgumentError(
        "carbonrange must contain integer carbon numbers"))

    Int.(carbons)
end

function validate_alkane_ladder_addition_settings(
    minradius,
    radiusfraction,
    positionsigmafraction,
    maxextensionsteps=100,
    gapmincosinefloor=0.85,
    gapcosinetolerance=0.03,
    edgemaxanchors=6,
    edgeminradius=5.0,
    edgeradiusfraction=0.2,
    edgemincosinefloor=0.9,
    edgecosinetolerance=0.03,
    edgecosineanchorcount=3,
    edgepositionsigmafraction=0.1,
    massspectrumvariancefloor=1.0
)
    isfinite(minradius) && minradius >= 0 || throw(ArgumentError(
        "minradius must be finite and nonnegative"))
    isfinite(radiusfraction) && radiusfraction >= 0 || throw(ArgumentError(
        "radiusfraction must be finite and nonnegative"))
    isfinite(positionsigmafraction) && positionsigmafraction > 0 ||
        throw(ArgumentError(
            "positionsigmafraction must be finite and positive"))
    maxextensionsteps isa Integer || throw(ArgumentError(
        "maxextensionsteps must be an integer"))
    maxextensionsteps >= 0 || throw(ArgumentError(
        "maxextensionsteps must be nonnegative"))
    isfinite(gapmincosinefloor) && 0 <= gapmincosinefloor <= 1 ||
        throw(ArgumentError("gapmincosinefloor must be finite and in [0, 1]"))
    isfinite(gapcosinetolerance) && gapcosinetolerance >= 0 ||
        throw(ArgumentError("gapcosinetolerance must be finite and nonnegative"))
    edgemaxanchors isa Integer || throw(ArgumentError(
        "edgemaxanchors must be an integer"))
    edgemaxanchors >= 2 || throw(ArgumentError(
        "edgemaxanchors must be at least 2"))
    isfinite(edgeminradius) && edgeminradius >= 0 || throw(ArgumentError(
        "edgeminradius must be finite and nonnegative"))
    isfinite(edgeradiusfraction) && edgeradiusfraction >= 0 || throw(ArgumentError(
        "edgeradiusfraction must be finite and nonnegative"))
    isfinite(edgemincosinefloor) && 0 <= edgemincosinefloor <= 1 ||
        throw(ArgumentError("edgemincosinefloor must be finite and in [0, 1]"))
    isfinite(edgecosinetolerance) && edgecosinetolerance >= 0 ||
        throw(ArgumentError("edgecosinetolerance must be finite and nonnegative"))
    edgecosineanchorcount isa Integer || throw(ArgumentError(
        "edgecosineanchorcount must be an integer"))
    edgecosineanchorcount >= 1 || throw(ArgumentError(
        "edgecosineanchorcount must be at least 1"))
    isfinite(edgepositionsigmafraction) && edgepositionsigmafraction > 0 ||
        throw(ArgumentError(
            "edgepositionsigmafraction must be finite and positive"))
    validate_alkane_abundance_variancefloor(massspectrumvariancefloor)

    nothing
end

function alkane_successful_ladder_apexes(apexinfo)
    anchors = [
        alkane_ladder_apex_anchor(apex) for apex in apexinfo.apexes
        if hasproperty(apex, :success) && Bool(apex.success)
    ]
    sort!(anchors; by=apex -> Int(apex.ladderstep))

    anchors
end

function alkane_ladder_apex_anchor(apex)
    (
        success=true,
        ladderstep=Int(apex.ladderstep),
        apexscanindex=Float64(apex.apexscanindex),
        apexretention=alkane_ladder_optional_float(apex, :apexretention, :apex_retention),
        source=hasproperty(apex, :source) ? Symbol(apex.source) : :molecularion,
        mass_spectrum_cosine=alkane_ladder_optional_float(
            apex,
            :mass_spectrum_cosine,
            :massspectrumcosine
        ),
        required_cosine=alkane_ladder_optional_float(
            apex,
            :required_cosine,
            :requiredcosine
        )
    )
end

function alkane_single_gap_additions(
    msm::MassScanMatrix,
    variances,
    abundanceinfo,
    anchors::AbstractVector,
    scancount::Integer;
    settings,
    massspectrumcontext,
    standard
)
    additions = NamedTuple[]
    diagnostics = NamedTuple[]
    length(anchors) >= 2 || return additions, diagnostics

    for index in 1:(length(anchors) - 1)
        left = anchors[index]
        right = anchors[index + 1]
        Int(right.ladderstep) - Int(left.ladderstep) == 2 || continue

        evaluation = alkane_ladder_addition_between_anchors(
            msm,
            variances,
            abundanceinfo,
            Int(left.ladderstep) + 1,
            :gapfilled,
            left,
            right,
            scancount;
            settings=settings,
            massspectrumcontext=massspectrumcontext,
            standard=standard
        )
        push!(diagnostics, evaluation.diagnostic)
        isnothing(evaluation.addition) || push!(additions, evaluation.addition)
    end

    additions, diagnostics
end

function alkane_edge_addition(
    msm::MassScanMatrix,
    variances,
    abundanceinfo,
    anchors::AbstractVector,
    direction::Symbol,
    scancount::Integer;
    settings,
    massspectrumcontext,
    standard
)
    direction in (:left, :right) || throw(ArgumentError(
        "direction must be :left or :right"))
    length(anchors) >= 2 || return NamedTuple[], nothing

    edgeanchors = alkane_ladder_edge_extension_anchor_results(
        anchors,
        direction;
        maxanchors=settings.edgemaxanchors
    )
    length(edgeanchors) >= 2 || return NamedTuple[], merge(
        alkane_ladder_addition_diagnostic_base(
            direction == :left ?
                Int(first(sort(collect(anchors); by=apex -> Int(apex.ladderstep))).ladderstep) - 1 :
                Int(last(sort(collect(anchors); by=apex -> Int(apex.ladderstep))).ladderstep) + 1,
            direction == :left ? :leftextended : :rightextended,
            NaN,
            NaN
        ),
        (reason=:insufficient_refined_anchors,)
    )
    edge = direction == :left ? first(edgeanchors) : last(edgeanchors)
    targetstep = direction == :left ?
        Int(edge.ladderstep) - 1 :
        Int(edge.ladderstep) + 1

    source = direction == :left ? :leftextended : :rightextended
    evaluation = alkane_ladder_addition_from_edge(
        msm,
        variances,
        abundanceinfo,
        targetstep,
        source,
        edge,
        edgeanchors,
        direction,
        scancount;
        settings=settings,
        massspectrumcontext=massspectrumcontext,
        standard=standard
    )
    additions = isnothing(evaluation.addition) ?
        NamedTuple[] :
        NamedTuple[evaluation.addition]

    additions, evaluation.diagnostic
end

function alkane_iterative_edge_additions(
    msm::MassScanMatrix,
    variances,
    abundanceinfo,
    anchors::AbstractVector,
    direction::Symbol,
    scancount::Integer;
    settings,
    massspectrumcontext,
    standard
)
    additions = NamedTuple[]
    diagnostics = NamedTuple[]
    workinganchors = NamedTuple[]
    append!(workinganchors, anchors)

    maxsteps = min(
        settings.maxextensionsteps,
        alkane_ladder_edge_extension_step_cap(
            workinganchors,
            direction,
            settings.carbonrange
        )
    )
    for _ in 1:maxsteps
        stepadditions, diagnostic = alkane_edge_addition(
            msm,
            variances,
            abundanceinfo,
            workinganchors,
            direction,
            scancount;
            settings=settings,
            massspectrumcontext=massspectrumcontext,
            standard=standard
        )
        isnothing(diagnostic) || push!(diagnostics, diagnostic)
        isempty(stepadditions) && break

        addition = only(stepadditions)
        push!(additions, addition)
        push!(workinganchors, alkane_ladder_addition_anchor(addition))
        sort!(workinganchors; by=anchor -> Int(anchor.ladderstep))
    end

    additions, diagnostics
end

function alkane_ladder_edge_extension_step_cap(
    anchors::AbstractVector,
    direction::Symbol,
    carbonrange
)
    isempty(anchors) && return 0
    direction in (:left, :right) || throw(ArgumentError(
        "direction must be :left or :right"))
    carbons = Int.(collect(carbonrange))
    isempty(carbons) && throw(ArgumentError("carbonrange must not be empty"))

    sorted = sort(collect(anchors); by=anchor -> Int(anchor.ladderstep))
    mincarbon = minimum(carbons)
    maxcarbon = maximum(carbons)
    edgecarbon = direction == :left ?
        Int(first(sorted).ladderstep) :
        Int(last(sorted).ladderstep)

    direction == :left ?
        max(0, edgecarbon - mincarbon) :
        max(0, maxcarbon - edgecarbon)
end

function alkane_ladder_addition_anchor(addition)
    (
        success=true,
        ladderstep=Int(addition.ladderstep),
        apexscanindex=Float64(addition.apexscanindex),
        apexretention=Float64(addition.apexretention),
        source=addition.source,
        mass_spectrum_cosine=alkane_ladder_optional_float(
            addition,
            :mass_spectrum_cosine,
            :massspectrumcosine
        ),
        required_cosine=alkane_ladder_optional_float(
            addition,
            :required_cosine,
            :requiredcosine
        )
    )
end

function alkane_ladder_optional_float(object, fields::Symbol...)
    for field in fields
        hasproperty(object, field) || continue
        value = getproperty(object, field)
        value isa Real || continue
        return Float64(value)
    end

    NaN
end

function alkane_ladder_addition_between_anchors(
    msm::MassScanMatrix,
    variances,
    abundanceinfo,
    targetstep::Integer,
    source::Symbol,
    leftanchor,
    rightanchor,
    scancount::Integer;
    settings,
    massspectrumcontext,
    standard
)
    leftscan = Float64(leftanchor.apexscanindex)
    rightscan = Float64(rightanchor.apexscanindex)
    localstepgap = (rightscan - leftscan) /
        (Int(rightanchor.ladderstep) - Int(leftanchor.ladderstep))
    expectedscan = (leftscan + rightscan) / 2
    diagnostic = alkane_ladder_addition_diagnostic_base(
        Int(targetstep),
        source,
        expectedscan,
        localstepgap
    )
    localstepgap > 0 && isfinite(localstepgap) || return (
        addition=nothing,
        diagnostic=merge(diagnostic, (reason=:invalid_anchor_spacing,))
    )

    evaluation = alkane_best_ladder_addition_window(
        msm,
        variances,
        abundanceinfo,
        Int(targetstep),
        source,
        expectedscan,
        localstepgap,
        scancount;
        settings=settings,
        massspectrumcontext=massspectrumcontext,
        standard=standard,
        direction=:between,
        leftanchor=leftanchor,
        rightanchor=rightanchor
    )

    evaluation
end

function alkane_ladder_addition_from_edge(
    msm::MassScanMatrix,
    variances,
    abundanceinfo,
    targetstep::Integer,
    source::Symbol,
    edgeanchor,
    edgeanchors,
    direction::Symbol,
    scancount::Integer;
    settings,
    massspectrumcontext,
    standard
)
    prediction = alkane_ladder_edge_extension_scan_prediction(edgeanchors, targetstep)
    isnothing(prediction) && return (
        addition=nothing,
        diagnostic=merge(
            alkane_ladder_addition_diagnostic_base(
                Int(targetstep),
                source,
                NaN,
                NaN
            ),
            (reason=:no_scan_prediction,)
        )
    )
    expectedscan = prediction.expectedscan
    localstepgap = prediction.localstepgap
    diagnostic = alkane_ladder_addition_diagnostic_base(
        Int(targetstep),
        source,
        expectedscan,
        localstepgap
    )
    localstepgap > 0 && isfinite(localstepgap) || return (
        addition=nothing,
        diagnostic=merge(diagnostic, (reason=:invalid_anchor_spacing,))
    )

    alkane_best_ladder_addition_window(
        msm,
        variances,
        abundanceinfo,
        Int(targetstep),
        source,
        expectedscan,
        localstepgap,
        scancount;
        settings=settings,
        massspectrumcontext=massspectrumcontext,
        standard=standard,
        direction=direction,
        edgeanchor=edgeanchor,
        edgeanchors=edgeanchors
    )
end

function alkane_best_ladder_addition_window(
    msm::MassScanMatrix,
    variances,
    abundanceinfo,
    targetstep::Integer,
    source::Symbol,
    expectedscan::Real,
    localstepgap::Real,
    scancount::Integer;
    settings,
    massspectrumcontext,
    standard,
    direction::Symbol,
    kwargs...
)
    diagnostic = alkane_ladder_addition_diagnostic_base(
        targetstep,
        source,
        expectedscan,
        localstepgap
    )
    isfinite(expectedscan) || return (
        addition=nothing,
        diagnostic=merge(diagnostic, (reason=:nonfinite_expected_scan,))
    )

    isedge = source in (:leftextended, :rightextended)
    minradius = isedge ? settings.edgeminradius : settings.minradius
    radiusfraction = isedge ? settings.edgeradiusfraction : settings.radiusfraction
    positionsigmafraction = isedge ?
        settings.edgepositionsigmafraction :
        settings.positionsigmafraction
    searchradius = max(minradius, radiusfraction * localstepgap)
    positionsigma = positionsigmafraction * localstepgap
    diagnostic = merge(
        diagnostic,
        (
            searchradius=Float64(searchradius),
            positionsigma=Float64(positionsigma)
        )
    )
    searchradius > 0 && isfinite(searchradius) || return (
        addition=nothing,
        diagnostic=merge(diagnostic, (reason=:invalid_search_radius,))
    )
    positionsigma > 0 && isfinite(positionsigma) || return (
        addition=nothing,
        diagnostic=merge(diagnostic, (reason=:invalid_position_sigma,))
    )

    windows = get(abundanceinfo.windows, Int(targetstep), NamedTuple[])
    abundance = get(abundanceinfo.abundances, Int(targetstep), nothing)
    diagnostic = merge(diagnostic, (candidatewindows=length(windows),))
    isempty(windows) && return (
        addition=nothing,
        diagnostic=merge(diagnostic, (reason=:no_abundance_windows,))
    )
    isnothing(abundance) && return (
        addition=nothing,
        diagnostic=merge(diagnostic, (reason=:no_abundance_vector,))
    )

    best = nothing
    bestrank = nothing
    inradius = 0
    directionok = 0
    finitecosine = 0
    passingcosine = 0
    completewindows = 0
    validpeakmodelwindows = 0
    apexconsistentwindows = 0
    requiredcosine = if source == :gapfilled
        alkane_single_gap_required_cosine(
            get(kwargs, :leftanchor, nothing),
            get(kwargs, :rightanchor, nothing);
            mincosinefloor=settings.gapmincosinefloor,
            cosinetolerance=settings.gapcosinetolerance
        )
    else
        alkane_edge_extension_required_cosine(
            get(kwargs, :edgeanchors, NamedTuple[]),
            direction;
            mincosinefloor=settings.edgemincosinefloor,
            cosinetolerance=settings.edgecosinetolerance,
            cosineanchorcount=settings.edgecosineanchorcount
        )
    end
    for window in windows
        apexindex = Int(window.apexindex)
        alkane_ladder_addition_direction_passes(
            apexindex,
            direction,
            expectedscan,
            scancount;
            kwargs...
        ) || continue
        directionok += 1
        scanerror = Float64(apexindex) - Float64(expectedscan)
        abs(scanerror) <= searchradius || continue
        inradius += 1
        if isedge && alkane_ladder_edge_window_apex_is_boundary_truncated(
                abundance,
                window,
                direction,
                scancount;
                maxboundaryapexdistance=settings.apexmaxshiftfromguess
            )
            continue
        end
        completewindows += 1

        apexabundance = hasproperty(window, :apexabundance) ?
            Float64(window.apexabundance) :
            Float64(abundance[apexindex])
        spectrumcandidate = alkane_ladder_addition_spectrum_candidate(
            abundance,
            window,
            Int(targetstep)
        )
        isnothing(spectrumcandidate) && continue
        validpeakmodelwindows += 1
        if isedge && Int(spectrumcandidate.apexindex) != apexindex
            continue
        end
        apexconsistentwindows += 1
        match = alkane_ladder_candidate_mass_spectrum_match(
            spectrumcandidate,
            massspectrumcontext
        )
        if isfinite(match.cosine)
            finitecosine += 1
        end
        if settings.massspectrummatch
            isfinite(match.cosine) || continue
            match.cosine >= requiredcosine || continue
            passingcosine += 1
            score = alkane_ladder_position_penalized_cosine_score(
                match.cosine,
                scanerror,
                positionsigma
            )
        else
            score = alkane_ladder_addition_score(
                apexabundance,
                scanerror,
                positionsigma
            )
        end
        isfinite(score) || continue
        candidate = merge(
            spectrumcandidate,
            (
                ladderstep=Int(targetstep),
                source=source,
                scanindex=apexindex,
                expectedscan=Float64(expectedscan),
                scanerror=Float64(scanerror),
                searchradius=Float64(searchradius),
                positionsigma=Float64(positionsigma),
                localstepgap=Float64(localstepgap),
                score=Float64(score),
                apexabundance=Float64(apexabundance),
                massspectrumcosine=Float64(match.cosine),
                massspectrumdistance=Float64(match.distance),
                massspectrumioncount=Int(match.ioncount),
                requiredcosine=Float64(requiredcosine),
                window=window
            )
        )
        apex = alkaneladderapex(
            msm,
            variances,
            abundanceinfo,
            candidate;
            scanwindow=settings.apexscanwindow,
            standard=standard,
            apexionexcludemzvalues=settings.apexionexcludemzvalues,
            apexionmzvalues=settings.apexionmzvalues,
            apexionminrelativeintensity=settings.apexionminrelativeintensity,
            minioncount=settings.apexminioncount,
            mzretentionkwargs=settings.apexmzretentionkwargs,
            variancefloor=settings.apexvariancefloor,
            logfloorfraction=settings.apexlogfloorfraction,
            maxapexshiftfromguess=settings.apexmaxshiftfromguess
        )
        apex.success || continue
        rank = (Float64(score), Float64(match.cosine), Float64(apexabundance))
        (isnothing(bestrank) || rank > bestrank) || continue
        bestrank = rank
        best = merge(
            candidate,
            (
                apex=apex,
                apexsuccess=true,
                apexscanindex=Float64(apex.apexscanindex),
                apexretention=Float64(apex.apexretention),
                apexrefinementreason=apex.reason
            )
        )
    end

    diagnostic = merge(
        diagnostic,
        (
            directionwindows=directionok,
            inradiuswindows=inradius,
            completewindows=completewindows,
            validpeakmodelwindows=validpeakmodelwindows,
            apexconsistentwindows=apexconsistentwindows,
            finitecosinewindows=finitecosine,
            passingcosinewindows=passingcosine,
            requiredcosine=Float64(requiredcosine)
        )
    )

    if !isnothing(best)
        return (
            addition=best,
            diagnostic=merge(
                diagnostic,
                (
                    status=:accepted,
                    reason=:passed,
                    bestscanindex=best.scanindex,
                    bestscanerror=best.scanerror,
                    bestscore=best.score,
                    bestapexabundance=best.apexabundance,
                    bestapexscanindex=best.apexscanindex,
                    bestcosine=best.massspectrumcosine
                )
            )
        )
    end

    reason = directionok == 0 ? :no_window_in_extension_direction :
        inradius == 0 ? :no_window_within_search_radius :
        completewindows == 0 && isedge ? :abundance_window_truncated_by_run_boundary :
        validpeakmodelwindows == 0 ? :no_valid_abundance_peak_model :
        apexconsistentwindows == 0 && isedge ? :abundance_window_apex_mismatch :
        settings.massspectrummatch && finitecosine == 0 ? :no_finite_mass_spectrum_cosine :
        settings.massspectrummatch && passingcosine == 0 ? :cosine_below_required_threshold :
        :no_candidate_selected

    (
        addition=nothing,
        diagnostic=merge(diagnostic, (reason=reason,))
    )
end

function alkane_ladder_addition_direction_passes(
    apexindex::Integer,
    direction::Symbol,
    expectedscan::Real,
    scancount::Integer;
    kwargs...
)
    1 <= apexindex <= scancount || return false
    if direction == :between
        leftanchor = get(kwargs, :leftanchor, nothing)
        rightanchor = get(kwargs, :rightanchor, nothing)
        isnothing(leftanchor) && return true
        isnothing(rightanchor) && return true
        return alkane_ladder_result_scanindex(leftanchor) <
            Float64(apexindex) <
            alkane_ladder_result_scanindex(rightanchor)
    end
    edgeanchor = get(kwargs, :edgeanchor, nothing)
    if direction == :left
        isnothing(edgeanchor) && return apexindex <= ceil(Int, expectedscan)
        return apexindex < Float64(alkane_ladder_result_scanindex(edgeanchor))
    elseif direction == :right
        isnothing(edgeanchor) && return apexindex >= floor(Int, expectedscan)
        return apexindex > Float64(alkane_ladder_result_scanindex(edgeanchor))
    end

    false
end

function alkane_ladder_result_scanindex(result)
    if hasproperty(result, :output_scan_index)
        return Float64(result.output_scan_index)
    elseif hasproperty(result, :outputscanindex)
        return Float64(result.outputscanindex)
    elseif hasproperty(result, :apexscanindex)
        return Float64(result.apexscanindex)
    elseif hasproperty(result, :scanindex)
        return Float64(result.scanindex)
    end

    NaN
end

function alkane_ladder_edge_extension_anchor_results(
    anchors::AbstractVector,
    direction::Symbol;
    maxanchors::Integer
)
    sorted = sort(collect(anchors); by=anchor -> Int(anchor.ladderstep))
    if direction == :left
        return first(sorted, min(length(sorted), Int(maxanchors)))
    elseif direction == :right
        return last(sorted, min(length(sorted), Int(maxanchors)))
    end

    throw(ArgumentError("direction must be :left or :right"))
end

function alkane_single_gap_required_cosine(
    left,
    right;
    mincosinefloor::Real,
    cosinetolerance::Real
)
    localcosines = Float64[]
    leftcosine = alkane_ladder_anchor_mass_spectrum_cosine(left)
    rightcosine = alkane_ladder_anchor_mass_spectrum_cosine(right)
    isfinite(leftcosine) && push!(localcosines, leftcosine)
    isfinite(rightcosine) && push!(localcosines, rightcosine)
    isempty(localcosines) && return Float64(mincosinefloor)

    max(Float64(mincosinefloor), minimum(localcosines) - Float64(cosinetolerance))
end

function alkane_edge_extension_required_cosine(
    anchors::AbstractVector,
    direction::Symbol;
    mincosinefloor::Real,
    cosinetolerance::Real,
    cosineanchorcount::Integer
)
    localanchors = direction == :left ?
        first(anchors, min(length(anchors), Int(cosineanchorcount))) :
        last(anchors, min(length(anchors), Int(cosineanchorcount)))
    localcosines = Float64[]
    for anchor in localanchors
        cosine = alkane_ladder_anchor_mass_spectrum_cosine(anchor)
        isfinite(cosine) && push!(localcosines, cosine)
    end
    isempty(localcosines) && return Float64(mincosinefloor)

    max(Float64(mincosinefloor), minimum(localcosines) - Float64(cosinetolerance))
end

function alkane_ladder_anchor_mass_spectrum_cosine(anchor)
    hasproperty(anchor, :mass_spectrum_cosine) &&
        return Float64(anchor.mass_spectrum_cosine)
    hasproperty(anchor, :massspectrumcosine) &&
        return Float64(anchor.massspectrumcosine)
    if hasproperty(anchor, :candidate)
        candidate = anchor.candidate
        hasproperty(candidate, :massspectrumcosine) &&
            return Float64(candidate.massspectrumcosine)
    end
    if hasproperty(anchor, :addition)
        addition = anchor.addition
        hasproperty(addition, :massspectrumcosine) &&
            return Float64(addition.massspectrumcosine)
    end

    NaN
end

function alkane_ladder_edge_extension_scan_prediction(
    anchors::AbstractVector,
    targetstep::Integer
)
    length(anchors) >= 2 || return nothing
    sorted = sort(collect(anchors); by=anchor -> Int(anchor.ladderstep))
    target = Int(targetstep)
    if target > Int(last(sorted).ladderstep)
        edge = last(sorted)
        neighbor = sorted[end - 1]
        edge_step = Int(edge.ladderstep)
        neighbor_step = Int(neighbor.ladderstep)
        step_delta = edge_step - neighbor_step
        step_delta > 0 || return nothing
        latest_gap =
            (alkane_ladder_result_scanindex(edge) -
             alkane_ladder_result_scanindex(neighbor)) / step_delta
        isfinite(latest_gap) && latest_gap > 0 || return nothing
        predicted_gap = latest_gap
        if length(sorted) >= 3
            previous = sorted[end - 2]
            previous_step = Int(previous.ladderstep)
            previous_step_delta = neighbor_step - previous_step
            if previous_step_delta > 0
                previous_gap =
                    (alkane_ladder_result_scanindex(neighbor) -
                     alkane_ladder_result_scanindex(previous)) / previous_step_delta
                gap_ratio = latest_gap / previous_gap
                if isfinite(gap_ratio) && gap_ratio > 0
                    predicted_gap = latest_gap * gap_ratio
                end
            end
        end
        isfinite(predicted_gap) && predicted_gap > 0 || return nothing
        expected_scan =
            alkane_ladder_result_scanindex(edge) +
            predicted_gap * Float64(target - edge_step)
    elseif target < Int(first(sorted).ladderstep)
        edge = first(sorted)
        neighbor = sorted[2]
        edge_step = Int(edge.ladderstep)
        neighbor_step = Int(neighbor.ladderstep)
        step_delta = neighbor_step - edge_step
        step_delta > 0 || return nothing
        latest_gap =
            (alkane_ladder_result_scanindex(neighbor) -
             alkane_ladder_result_scanindex(edge)) / step_delta
        isfinite(latest_gap) && latest_gap > 0 || return nothing
        predicted_gap = latest_gap
        if length(sorted) >= 3
            previous = sorted[3]
            previous_step = Int(previous.ladderstep)
            previous_step_delta = previous_step - neighbor_step
            if previous_step_delta > 0
                previous_gap =
                    (alkane_ladder_result_scanindex(previous) -
                     alkane_ladder_result_scanindex(neighbor)) / previous_step_delta
                gap_ratio = latest_gap / previous_gap
                if isfinite(gap_ratio) && gap_ratio > 0
                    predicted_gap = latest_gap * gap_ratio
                end
            end
        end
        isfinite(predicted_gap) && predicted_gap > 0 || return nothing
        expected_scan =
            alkane_ladder_result_scanindex(edge) -
            predicted_gap * Float64(edge_step - target)
    else
        return nothing
    end

    (
        expectedscan=Float64(expected_scan),
        rmse=0.0,
        localstepgap=Float64(predicted_gap),
        degree=1
    )
end

function alkane_ladder_edge_window_apex_is_boundary_truncated(
    abundance::AbstractVector{<:Real},
    window,
    direction::Symbol,
    scancount::Integer;
    maxboundaryapexdistance::Real
)
    isfinite(maxboundaryapexdistance) && maxboundaryapexdistance >= 0 ||
        throw(ArgumentError(
            "maxboundaryapexdistance must be finite and nonnegative"))
    1 <= window.leftindex <= window.rightindex <= length(abundance) ||
        return false

    scanindices = Int(window.leftindex):Int(window.rightindex)
    peakmodel = max.(Float64.(abundance[scanindices]), 0.0)
    maximum(peakmodel) > 0 || return false
    _, peakoffset = findmax(peakmodel)
    windowmaxscan = first(scanindices) + peakoffset - 1

    if direction == :right
        touchesboundary = Int(window.rightindex) >= scancount ||
            (hasproperty(window, :rightstop) && window.rightstop == :boundary)
        boundarydistance = scancount - windowmaxscan
        return touchesboundary &&
            boundarydistance <= Float64(maxboundaryapexdistance) + 1e-9
    elseif direction == :left
        touchesboundary = Int(window.leftindex) <= 1 ||
            (hasproperty(window, :leftstop) && window.leftstop == :boundary)
        boundarydistance = windowmaxscan - 1
        return touchesboundary &&
            boundarydistance <= Float64(maxboundaryapexdistance) + 1e-9
    end

    throw(ArgumentError("direction must be :left or :right"))
end

function alkane_ladder_addition_score(
    apexabundance::Real,
    scanerror::Real,
    positionsigma::Real
)
    isfinite(apexabundance) && apexabundance > 0 || return -Inf
    isfinite(scanerror) || return -Inf
    isfinite(positionsigma) && positionsigma > 0 || return -Inf

    Float64(apexabundance) * exp(-0.5 * abs2(Float64(scanerror) / positionsigma))
end

function alkane_ladder_position_penalized_cosine_score(
    cosine::Real,
    scanerror::Real,
    positionsigma::Real
)
    isfinite(cosine) && cosine > 0 || return -Inf
    isfinite(scanerror) || return -Inf
    isfinite(positionsigma) && positionsigma > 0 || return -Inf

    z = Float64(scanerror) / Float64(positionsigma)

    log(Float64(cosine)) - 0.5 * abs2(z)
end

function alkane_ladder_addition_spectrum_candidate(
    abundance::AbstractVector{<:Real},
    window,
    step::Integer
)
    1 <= window.leftindex <= window.rightindex <= length(abundance) ||
        return nothing

    scanindices = Int(window.leftindex):Int(window.rightindex)
    peakmodel = max.(Float64.(abundance[scanindices]), 0.0)
    peakheight, peakoffset = findmax(peakmodel)
    peakheight > 0 && isfinite(peakheight) || return nothing
    peakmodel ./= peakheight

    (
        ladderstep=Int(step),
        leftindex=Int(window.leftindex),
        apexindex=first(scanindices) + peakoffset - 1,
        rightindex=Int(window.rightindex),
        scanindices=collect(scanindices),
        peakmodel=collect(peakmodel)
    )
end

function alkane_ladder_mass_spectrum_match_context(
    msm::MassScanMatrix,
    variances,
    standard;
    enabled::Bool,
    variancefloor::Real
)
    enabled || return nothing
    validate_alkane_series_variances(msm, variances)
    validate_alkane_abundance_variancefloor(variancefloor)

    references = alkane_ladder_mass_spectrum_reference_vectors(msm, standard)
    isempty(references) && return nothing

    (
        X=rawintensities(msm),
        variances=variances,
        variancefloor=Float64(variancefloor),
        references=references
    )
end

function alkane_ladder_mass_spectrum_reference_vectors(
    msm::MassScanMatrix,
    standard
)
    spectra = alkane_standard_spectra(standard)
    mzbins = alkane_mz_bins(msm)
    indexbymz = Dict(mz => index for (index, mz) in pairs(mzbins))
    references = Dict{Int, Vector{Float64}}()

    for spectrum in spectra
        spectrumattrs = attrs(spectrum)
        hasproperty(spectrumattrs, :order) || continue
        carbon = getproperty(spectrumattrs, :order)
        carbon isa Integer || continue

        reference = zeros(Float64, mzcount(msm))
        spectrum_mzs = integer_mz_values(
            mzvalues(spectrum),
            "reference spectrum m/z values"
        )
        spectrum_intensities = Float64.(intensities(spectrum))
        for (mz, intensity) in zip(spectrum_mzs, spectrum_intensities)
            intensity > 0 || continue
            mzindex = get(indexbymz, mz, nothing)
            isnothing(mzindex) && continue
            reference[mzindex] += intensity
        end
        sum(abs2, reference) > 0 || continue
        references[Int(carbon)] = reference
    end

    references
end

function alkane_ladder_candidate_mass_spectrum_match(candidate, context)
    isnothing(context) && return (
        cosine=NaN,
        distance=NaN,
        ioncount=0
    )
    reference = get(context.references, Int(candidate.ladderstep), nothing)
    isnothing(reference) && return (
        cosine=NaN,
        distance=NaN,
        ioncount=0
    )

    observed = Vector{Float64}(undef, length(reference))
    for mzindex in eachindex(reference)
        fit = alkane_fitted_ion_abundance(
            context.X,
            context.variances,
            candidate.scanindices,
            mzindex,
            candidate.peakmodel;
            variancefloor=context.variancefloor
        )
        observed[mzindex] = max(Float64(fit.abundance), 0.0)
    end

    cosine = alkane_ladder_cosine_similarity(observed, reference)
    distance = isfinite(cosine) ? 1.0 - clamp(Float64(cosine), 0.0, 1.0) : NaN

    (
        cosine=cosine,
        distance=distance,
        ioncount=Base.count(>(0.0), reference)
    )
end

function alkane_ladder_cosine_similarity(
    observed::AbstractVector{<:Real},
    reference::AbstractVector{<:Real}
)
    length(observed) == length(reference) || throw(DimensionMismatch(
        "observed and reference vectors must have the same length"))
    isempty(observed) && return NaN

    numerator = 0.0
    observednorm2 = 0.0
    referencenorm2 = 0.0
    for (observedvalue, referencevalue) in zip(observed, reference)
        observedfloat = Float64(observedvalue)
        referencefloat = Float64(referencevalue)
        numerator += observedfloat * referencefloat
        observednorm2 += abs2(observedfloat)
        referencenorm2 += abs2(referencefloat)
    end
    referencenorm2 > 0 || return NaN
    observednorm2 > 0 || return 0.0

    clamp(numerator / sqrt(observednorm2 * referencenorm2), 0.0, 1.0)
end

function alkane_ladder_addition_diagnostic_base(
    targetstep::Integer,
    source::Symbol,
    expectedscan,
    localstepgap
)
    (
        status=:failed,
        reason=:not_evaluated,
        ladderstep=Int(targetstep),
        source=source,
        expectedscan=Float64(expectedscan),
        localstepgap=Float64(localstepgap),
        searchradius=NaN,
        positionsigma=NaN,
        candidatewindows=0,
        directionwindows=0,
        inradiuswindows=0,
        finitecosinewindows=0,
        passingcosinewindows=0,
        requiredcosine=NaN,
        bestscanindex=missing,
        bestscanerror=NaN,
        bestscore=NaN,
        bestapexabundance=NaN,
        bestapexscanindex=NaN,
        bestcosine=NaN
    )
end
