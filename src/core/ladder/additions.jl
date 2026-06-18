struct AlkaneLadderAdditionSettings{
    T1<:Real,
    T2<:Real,
    T3<:Real,
    T4<:Integer,
    T5<:Real,
    T6<:Real,
    T7<:Integer,
    T8<:Real,
    T9<:Real,
    T10<:Real,
    T11<:Real,
    T12<:Integer,
    T13<:Real,
    T14<:Real,
    T15<:Integer,
    T16<:Real,
    T17<:Real,
    T18<:Union{AbstractVector, Tuple},
    T19<:Union{Nothing, AbstractVector},
    T20<:Real,
    T21<:Integer,
    T22<:NamedTuple,
    T23<:Real,
    T24<:AbstractVector{<:Integer}
}
    minradius::T1
    radiusfraction::T2
    positionsigmafraction::T3
    maxextensionsteps::T4
    massspectrummatch::Bool
    gapmincosinefloor::T5
    gapcosinetolerance::T6
    edgemaxanchors::T7
    edgeminradius::T8
    edgeradiusfraction::T9
    edgemincosinefloor::T10
    edgecosinetolerance::T11
    edgecosineanchorcount::T12
    edgepositionsigmafraction::T13
    massspectrumvariancefloor::T14
    apexscanwindow::T15
    apexvariancefloor::T16
    apexlogfloorfraction::T17
    apexionexcludemzvalues::T18
    apexionmzvalues::T19
    apexionminrelativeintensity::T20
    apexminioncount::T21
    apexmzretentionkwargs::T22
    apexmaxshiftfromguess::T23
    carbonrange::T24
end

struct AlkaneLadderAdditionAnchor
    ladderstep::Int
    apexscanindex::Float64
    apexretention::Float64
    source::Symbol
    mass_spectrum_cosine::Float64
    required_cosine::Float64
end

struct AlkaneLadderAdditionCandidate{T<:AlkaneAbundanceWindow} <:
        AbstractAlkaneLadderCandidate
    ladderstep::Int
    leftindex::Int
    apexindex::Int
    rightindex::Int
    scanindices::Vector{Int}
    peakmodel::Vector{Float64}
    source::Symbol
    scanindex::Int
    expectedscan::Float64
    scanerror::Float64
    searchradius::Float64
    positionsigma::Float64
    localstepgap::Float64
    score::Float64
    apexabundance::Float64
    massspectrumcosine::Float64
    massspectrumdistance::Float64
    massspectrumioncount::Int
    requiredcosine::Float64
    window::T
end

struct AlkaneLadderAddition{T1<:AlkaneAbundanceWindow, T2<:AlkaneLadderApex}
    ladderstep::Int
    leftindex::Int
    apexindex::Int
    rightindex::Int
    scanindices::Vector{Int}
    peakmodel::Vector{Float64}
    source::Symbol
    scanindex::Int
    expectedscan::Float64
    scanerror::Float64
    searchradius::Float64
    positionsigma::Float64
    localstepgap::Float64
    score::Float64
    apexabundance::Float64
    massspectrumcosine::Float64
    massspectrumdistance::Float64
    massspectrumioncount::Int
    requiredcosine::Float64
    window::T1
    apex::T2
    apexsuccess::Bool
    apexscanindex::Float64
    apexretention::Float64
    apexrefinementreason::Symbol
end

function AlkaneLadderAddition(
    candidate::AlkaneLadderAdditionCandidate,
    apex::AlkaneLadderApex
)
    AlkaneLadderAddition(
        candidate.ladderstep,
        candidate.leftindex,
        candidate.apexindex,
        candidate.rightindex,
        candidate.scanindices,
        candidate.peakmodel,
        candidate.source,
        candidate.scanindex,
        candidate.expectedscan,
        candidate.scanerror,
        candidate.searchradius,
        candidate.positionsigma,
        candidate.localstepgap,
        candidate.score,
        candidate.apexabundance,
        candidate.massspectrumcosine,
        candidate.massspectrumdistance,
        candidate.massspectrumioncount,
        candidate.requiredcosine,
        candidate.window,
        apex,
        true,
        apex.apexscanindex,
        apex.apexretention,
        apex.reason
    )
end

alkane_ladder_candidate_step(candidate::AlkaneLadderAdditionCandidate) =
    candidate.ladderstep

alkane_ladder_candidate_source(candidate::AlkaneLadderAdditionCandidate) =
    candidate.source

alkane_ladder_candidate_is_gapfilled(candidate::AlkaneLadderAdditionCandidate) =
    candidate.source ≡ :gapfilled

alkane_ladder_candidate_is_edgeextended(candidate::AlkaneLadderAdditionCandidate) =
    candidate.source in (:leftextended, :rightextended)

alkane_ladder_candidate_mass_spectrum_cosine(candidate::AlkaneLadderAdditionCandidate) = 
    candidate.massspectrumcosine

alkane_ladder_candidate_required_cosine(candidate::AlkaneLadderAdditionCandidate) =
    candidate.requiredcosine

alkane_ladder_input_scan_index(candidate::AlkaneLadderAdditionCandidate) =
    candidate.scanindex

function alkane_ladder_fallback_scanindices(
    retentions::AbstractVector{<:Real},
    scanindex::Integer,
    scanwindow::Integer,
    candidate::AlkaneLadderAdditionCandidate
)
    leftindex = max(1, scanindex - scanwindow)
    rightindex = min(length(retentions), scanindex + scanwindow)
    leftindex = max(leftindex, candidate.window.leftindex)
    rightindex = min(rightindex, candidate.window.rightindex)

    collect(leftindex:rightindex)
end

struct AlkaneLadderAdditionDiagnostic
    status::Symbol
    reason::Symbol
    ladderstep::Int
    source::Symbol
    expectedscan::Float64
    localstepgap::Float64
    searchradius::Float64
    positionsigma::Float64
    candidatewindows::Int
    directionwindows::Int
    inradiuswindows::Int
    completewindows::Int
    validpeakmodelwindows::Int
    apexconsistentwindows::Int
    finitecosinewindows::Int
    passingcosinewindows::Int
    requiredcosine::Float64
    bestscanindex::Union{Missing, Int}
    bestscanerror::Float64
    bestscore::Float64
    bestapexabundance::Float64
    bestapexscanindex::Float64
    bestcosine::Float64
end

struct AlkaneLadderAdditionEvaluation{T<:Union{Nothing, AlkaneLadderAddition}}
    addition::T
    diagnostic::AlkaneLadderAdditionDiagnostic
end

struct AlkaneLadderAdditionDiagnostics
    gapfilled::Vector{AlkaneLadderAdditionDiagnostic}
    leftextended::Vector{AlkaneLadderAdditionDiagnostic}
    rightextended::Vector{AlkaneLadderAdditionDiagnostic}
end

struct AlkaneLadderAdditionInfo{T1<:AlkaneLadderAddition, T2<:AlkaneLadderAdditionSettings}
    status::Symbol
    additions::Vector{T1}
    gapfilled::Vector{T1}
    leftextended::Vector{T1}
    rightextended::Vector{T1}
    diagnostics::AlkaneLadderAdditionDiagnostics
    settings::T2
end

struct AlkaneLadderAdditionSpectrumCandidate{T<:AlkaneAbundanceWindow}
    ladderstep::Int
    leftindex::Int
    apexindex::Int
    rightindex::Int
    scanindices::Vector{Int}
    peakmodel::Vector{Float64}
    window::T
end

struct AlkaneLadderMassSpectrumMatchContext{
    T1<:AbstractMatrix{<:Real},
    T2<:AbstractMatrix{<:Real}
}
    X::T1
    variances::T2
    variancefloor::Float64
    references::Dict{Int, Vector{Float64}}
end

struct AlkaneLadderMassSpectrumMatch
    cosine::Float64
    distance::Float64
    ioncount::Int
end

struct AlkaneLadderEdgeScanPrediction
    expectedscan::Float64
    rmse::Float64
    localstepgap::Float64
    degree::Int
end

function alkane_ladder_addition_settings(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    apexinfo::Union{AlkaneLadderApexInfo, NamedTuple},
    minradius::Real,
    radiusfraction::Real,
    positionsigmafraction::Real,
    maxextensionsteps::Integer,
    massspectrummatch::Bool,
    gapmincosinefloor::Real,
    gapcosinetolerance::Real,
    edgemaxanchors::Integer,
    edgeminradius::Real,
    edgeradiusfraction::Real,
    edgemincosinefloor::Real,
    edgecosinetolerance::Real,
    edgecosineanchorcount::Integer,
    edgepositionsigmafraction::Real,
    massspectrumvariancefloor::Real,
    apexscanwindow::Integer,
    apexvariancefloor::Real,
    apexlogfloorfraction::Real,
    apexionexcludemzvalues::Union{AbstractVector, Tuple},
    apexionmzvalues::Union{Nothing, AbstractVector},
    apexionminrelativeintensity::Real,
    apexminioncount::Integer,
    apexmzretentionkwargs::Union{Nothing, NamedTuple},
    apexmaxshiftfromguess::Real,
    carbonrange::Union{Nothing, AbstractVector{<:Integer}, AbstractRange{<:Integer}}
)
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

    resolved_apex_mzretentionkwargs = alkane_ladder_addition_apex_mzretentionkwargs(
        apexinfo,
        apexmzretentionkwargs
    )
    resolved_carbonrange = alkane_ladder_addition_carbonrange(abundanceinfo, carbonrange)
    AlkaneLadderAdditionSettings(
        minradius,
        radiusfraction,
        positionsigmafraction,
        maxextensionsteps,
        massspectrummatch,
        gapmincosinefloor,
        gapcosinetolerance,
        edgemaxanchors,
        edgeminradius,
        edgeradiusfraction,
        edgemincosinefloor,
        edgecosinetolerance,
        edgecosineanchorcount,
        edgepositionsigmafraction,
        massspectrumvariancefloor,
        apexscanwindow,
        apexvariancefloor,
        apexlogfloorfraction,
        apexionexcludemzvalues,
        apexionmzvalues,
        apexionminrelativeintensity,
        apexminioncount,
        resolved_apex_mzretentionkwargs,
        apexmaxshiftfromguess,
        resolved_carbonrange
    )
end

function alkaneladderadditions(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    pathinfo::Union{AlkaneLadderPathInfo, NamedTuple},
    apexinfo::Union{AlkaneLadderApexInfo, NamedTuple},
    settings::AlkaneLadderAdditionSettings,
    standard::AlkaneStandard
)
    validate_alkane_series_variances(msm, variances)
    massspectrumcontext = alkane_ladder_mass_spectrum_match_context(
        msm,
        variances,
        standard,
        settings.massspectrummatch,
        settings.massspectrumvariancefloor
    )
    anchors = alkane_successful_ladder_apexes(apexinfo)
    gapfilled, gapdiagnostics = alkane_single_gap_additions(
        msm,
        variances,
        abundanceinfo,
        anchors,
        scancount(msm),
        settings,
        massspectrumcontext,
        standard
    )
    edgeanchors = AlkaneLadderAdditionAnchor[]
    append!(edgeanchors, anchors)
    append!(edgeanchors, (alkane_ladder_addition_anchor(addition) for addition in gapfilled))
    sort!(edgeanchors; by=anchor -> anchor.ladderstep)
    leftextended, leftdiagnostics = alkane_iterative_edge_additions(
        msm,
        variances,
        abundanceinfo,
        edgeanchors,
        :left,
        scancount(msm),
        settings,
        massspectrumcontext,
        standard
    )
    rightextended, rightdiagnostics = alkane_iterative_edge_additions(
        msm,
        variances,
        abundanceinfo,
        edgeanchors,
        :right,
        scancount(msm),
        settings,
        massspectrumcontext,
        standard
    )

    additions = AlkaneLadderAddition[]
    append!(additions, gapfilled)
    append!(additions, leftextended)
    append!(additions, rightextended)
    sort!(additions; by=addition -> addition.ladderstep)

    AlkaneLadderAdditionInfo(
        isempty(additions) ? :empty : :success,
        additions,
        gapfilled,
        leftextended,
        rightextended,
        AlkaneLadderAdditionDiagnostics(
            gapdiagnostics,
            leftdiagnostics,
            rightdiagnostics
        ),
        settings
    )
end

function alkane_ladder_addition_apex_mzretentionkwargs(
    apexinfo::Union{AlkaneLadderApexInfo, NamedTuple},
    apexmzretentionkwargs::Union{Nothing, NamedTuple}
)
    !isnothing(apexmzretentionkwargs) && return apexmzretentionkwargs
    apexinfo.settings.mzretentionkwargs
end

function alkane_ladder_addition_apex_settings(
    settings::AlkaneLadderAdditionSettings,
    standard::AlkaneStandard
)
    AlkaneLadderApexSettings(
        standard,
        settings.apexscanwindow,
        settings.apexvariancefloor,
        settings.apexlogfloorfraction,
        settings.apexmzretentionkwargs,
        get(settings.apexmzretentionkwargs, :order, :inferdirection),
        settings.apexionexcludemzvalues,
        settings.apexionmzvalues,
        settings.apexionminrelativeintensity,
        settings.apexminioncount,
        settings.apexmaxshiftfromguess,
        nothing,
        nothing,
        0,
        1.25,
        0,
        14,
        0,
        0
    )
end

function alkane_ladder_addition_carbonrange(
    abundanceinfo::AlkaneAbundanceInfo,
    carbonrange::Union{Nothing, AbstractVector{<:Integer}, AbstractRange{<:Integer}}
)
    carbons = if isnothing(carbonrange)
        carbonkeys = collect(keys(abundanceinfo.abundances))
        append!(carbonkeys, keys(abundanceinfo.windows))
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
    minradius::Real,
    radiusfraction::Real,
    positionsigmafraction::Real,
    maxextensionsteps::Integer,
    gapmincosinefloor::Real,
    gapcosinetolerance::Real,
    edgemaxanchors::Integer,
    edgeminradius::Real,
    edgeradiusfraction::Real,
    edgemincosinefloor::Real,
    edgecosinetolerance::Real,
    edgecosineanchorcount::Integer,
    edgepositionsigmafraction::Real,
    massspectrumvariancefloor::Real
)
    isfinite(minradius) && minradius ≥ 0 || throw(ArgumentError(
        "minradius must be finite and nonnegative"))
    isfinite(radiusfraction) && radiusfraction ≥ 0 || throw(ArgumentError(
        "radiusfraction must be finite and nonnegative"))
    isfinite(positionsigmafraction) && positionsigmafraction > 0 ||
        throw(ArgumentError("positionsigmafraction must be finite and positive"))
    maxextensionsteps isa Integer || throw(ArgumentError(
        "maxextensionsteps must be an integer"))
    maxextensionsteps ≥ 0 || throw(ArgumentError(
        "maxextensionsteps must be nonnegative"))
    isfinite(gapmincosinefloor) && 0 ≤ gapmincosinefloor ≤ 1 ||
        throw(ArgumentError("gapmincosinefloor must be finite and in [0, 1]"))
    isfinite(gapcosinetolerance) && gapcosinetolerance ≥ 0 ||
        throw(ArgumentError("gapcosinetolerance must be finite and nonnegative"))
    edgemaxanchors isa Integer || throw(ArgumentError(
        "edgemaxanchors must be an integer"))
    edgemaxanchors ≥ 2 || throw(ArgumentError(
        "edgemaxanchors must be at least 2"))
    isfinite(edgeminradius) && edgeminradius ≥ 0 || throw(ArgumentError(
        "edgeminradius must be finite and nonnegative"))
    isfinite(edgeradiusfraction) && edgeradiusfraction ≥ 0 || throw(ArgumentError(
        "edgeradiusfraction must be finite and nonnegative"))
    isfinite(edgemincosinefloor) && 0 ≤ edgemincosinefloor ≤ 1 ||
        throw(ArgumentError("edgemincosinefloor must be finite and in [0, 1]"))
    isfinite(edgecosinetolerance) && edgecosinetolerance ≥ 0 ||
        throw(ArgumentError("edgecosinetolerance must be finite and nonnegative"))
    edgecosineanchorcount isa Integer || throw(ArgumentError(
        "edgecosineanchorcount must be an integer"))
    edgecosineanchorcount ≥ 1 || throw(ArgumentError(
        "edgecosineanchorcount must be at least 1"))
    isfinite(edgepositionsigmafraction) && edgepositionsigmafraction > 0 ||
        throw(ArgumentError("edgepositionsigmafraction must be finite and positive"))
    validate_alkane_abundance_variancefloor(massspectrumvariancefloor)

    nothing
end

function alkane_successful_ladder_apexes(
    apexinfo::Union{AlkaneLadderApexInfo, NamedTuple}
)
    anchors = AlkaneLadderAdditionAnchor[]
    for apex in apexinfo.apexes
        apex.success || continue
        push!(anchors, alkane_ladder_apex_anchor(apex))
    end
    sort!(anchors; by=apex -> apex.ladderstep)

    anchors
end

function alkane_ladder_apex_anchor(apex::AlkaneLadderApex)
    AlkaneLadderAdditionAnchor(
        apex.ladderstep,
        apex.apexscanindex,
        apex.apexretention,
        apex.source,
        apex.mass_spectrum_cosine,
        apex.required_cosine
    )
end

function alkane_ladder_apex_anchor(apex::NamedTuple)
    AlkaneLadderAdditionAnchor(
        apex.ladderstep,
        apex.apexscanindex,
        get(apex, :apexretention, get(apex, :apex_retention, NaN)),
        get(apex, :source, :molecularion),
        get(apex, :mass_spectrum_cosine, get(apex, :massspectrumcosine, NaN)),
        get(apex, :required_cosine, get(apex, :requiredcosine, NaN))
    )
end

function alkane_single_gap_additions(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    anchors::AbstractVector{AlkaneLadderAdditionAnchor},
    scancount::Integer,
    settings::AlkaneLadderAdditionSettings,
    massspectrumcontext::Union{Nothing, AlkaneLadderMassSpectrumMatchContext},
    standard::AlkaneStandard
)
    additions = AlkaneLadderAddition[]
    diagnostics = AlkaneLadderAdditionDiagnostic[]
    length(anchors) ≥ 2 || return additions, diagnostics

    for index in 1:(length(anchors) - 1)
        left = anchors[index]
        right = anchors[index + 1]
        right.ladderstep - left.ladderstep == 2 || continue

        evaluation = alkane_ladder_addition_between_anchors(
            msm,
            variances,
            abundanceinfo,
            left.ladderstep + 1,
            :gapfilled,
            left,
            right,
            scancount,
            settings,
            massspectrumcontext,
            standard
        )
        push!(diagnostics, evaluation.diagnostic)
        isnothing(evaluation.addition) || push!(additions, evaluation.addition)
    end

    additions, diagnostics
end

function alkane_edge_addition(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    anchors::AbstractVector{AlkaneLadderAdditionAnchor},
    direction::Symbol,
    scancount::Integer,
    settings::AlkaneLadderAdditionSettings,
    massspectrumcontext::Union{Nothing, AlkaneLadderMassSpectrumMatchContext},
    standard::AlkaneStandard
)
    direction in (:left, :right) || throw(ArgumentError(
        "direction must be :left or :right"))
    length(anchors) ≥ 2 || return AlkaneLadderAddition[], nothing

    edgeanchors = alkane_ladder_edge_extension_anchor_results(
        anchors,
        direction,
        settings.edgemaxanchors
    )
    length(edgeanchors) ≥ 2 || return AlkaneLadderAddition[],
        alkane_ladder_addition_diagnostic_base(
            direction ≡ :left ?
                first(sort(collect(anchors); by=apex -> apex.ladderstep)).ladderstep - 1 :
                last(sort(collect(anchors); by=apex -> apex.ladderstep)).ladderstep + 1,
            direction ≡ :left ? :leftextended : :rightextended,
            NaN,
            NaN,
            :insufficient_refined_anchors
        )
    edge = direction ≡ :left ? first(edgeanchors) : last(edgeanchors)
    targetstep = direction ≡ :left ?
        edge.ladderstep - 1 :
        edge.ladderstep + 1

    source = direction ≡ :left ? :leftextended : :rightextended
    evaluation = alkane_ladder_addition_from_edge(
        msm,
        variances,
        abundanceinfo,
        targetstep,
        source,
        edge,
        edgeanchors,
        direction,
        scancount,
        settings,
        massspectrumcontext,
        standard
    )
    additions = isnothing(evaluation.addition) ?
        AlkaneLadderAddition[] :
        AlkaneLadderAddition[evaluation.addition]

    additions, evaluation.diagnostic
end

function alkane_iterative_edge_additions(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    anchors::AbstractVector{AlkaneLadderAdditionAnchor},
    direction::Symbol,
    scancount::Integer,
    settings::AlkaneLadderAdditionSettings,
    massspectrumcontext::Union{Nothing, AlkaneLadderMassSpectrumMatchContext},
    standard::AlkaneStandard
)
    additions = AlkaneLadderAddition[]
    diagnostics = AlkaneLadderAdditionDiagnostic[]
    workinganchors = AlkaneLadderAdditionAnchor[]
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
            scancount,
            settings,
            massspectrumcontext,
            standard
        )
        isnothing(diagnostic) || push!(diagnostics, diagnostic)
        isempty(stepadditions) && break

        addition = only(stepadditions)
        push!(additions, addition)
        push!(workinganchors, alkane_ladder_addition_anchor(addition))
        sort!(workinganchors; by=anchor -> anchor.ladderstep)
    end

    additions, diagnostics
end

function alkane_ladder_edge_extension_step_cap(
    anchors::AbstractVector{AlkaneLadderAdditionAnchor},
    direction::Symbol,
    carbonrange::Union{AbstractVector{<:Integer}, AbstractRange{<:Integer}}
)
    isempty(anchors) && return 0
    direction in (:left, :right) || throw(ArgumentError(
        "direction must be :left or :right"))
    carbons = collect(carbonrange)
    isempty(carbons) && throw(ArgumentError("carbonrange must not be empty"))

    sorted = sort(collect(anchors); by=anchor -> anchor.ladderstep)
    mincarbon = minimum(carbons)
    maxcarbon = maximum(carbons)
    edgecarbon = direction ≡ :left ? first(sorted).ladderstep : last(sorted).ladderstep

    direction ≡ :left ? max(0, edgecarbon - mincarbon) : max(0, maxcarbon - edgecarbon)
end

function alkane_ladder_addition_anchor(addition::AlkaneLadderAddition)
    AlkaneLadderAdditionAnchor(
        addition.ladderstep,
        addition.apexscanindex,
        addition.apexretention,
        addition.source,
        addition.massspectrumcosine,
        addition.requiredcosine
    )
end

function alkane_ladder_addition_between_anchors(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    targetstep::Integer,
    source::Symbol,
    leftanchor::AlkaneLadderAdditionAnchor,
    rightanchor::AlkaneLadderAdditionAnchor,
    scancount::Integer,
    settings::AlkaneLadderAdditionSettings,
    massspectrumcontext::Union{Nothing, AlkaneLadderMassSpectrumMatchContext},
    standard::AlkaneStandard
)
    leftscan = leftanchor.apexscanindex
    rightscan = rightanchor.apexscanindex
    localstepgap = (rightscan - leftscan) / (rightanchor.ladderstep - leftanchor.ladderstep)
    expectedscan = (leftscan + rightscan) / 2
    diagnostic = alkane_ladder_addition_diagnostic_base(
        targetstep,
        source,
        expectedscan,
        localstepgap,
        :not_evaluated
    )
    localstepgap > 0 && isfinite(localstepgap) || return AlkaneLadderAdditionEvaluation(
        nothing,
        alkane_ladder_addition_diagnostic_reason(diagnostic, :invalid_anchor_spacing)
    )

    evaluation = alkane_best_ladder_addition_window(
        msm,
        variances,
        abundanceinfo,
        targetstep,
        source,
        expectedscan,
        localstepgap,
        scancount,
        settings,
        massspectrumcontext,
        standard,
        :between,
        leftanchor,
        rightanchor,
        nothing,
        AlkaneLadderAdditionAnchor[]
    )

    evaluation
end

function alkane_ladder_addition_from_edge(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    targetstep::Integer,
    source::Symbol,
    edgeanchor::AlkaneLadderAdditionAnchor,
    edgeanchors::AbstractVector{AlkaneLadderAdditionAnchor},
    direction::Symbol,
    scancount::Integer,
    settings::AlkaneLadderAdditionSettings,
    massspectrumcontext::Union{Nothing, AlkaneLadderMassSpectrumMatchContext},
    standard::AlkaneStandard
)
    prediction = alkane_ladder_edge_extension_scan_prediction(edgeanchors, targetstep)
    isnothing(prediction) && return AlkaneLadderAdditionEvaluation(
        nothing,
        alkane_ladder_addition_diagnostic_base(
            targetstep,
            source,
            NaN,
            NaN,
            :no_scan_prediction
        )
    )
    expectedscan = prediction.expectedscan
    localstepgap = prediction.localstepgap
    diagnostic = alkane_ladder_addition_diagnostic_base(
        targetstep,
        source,
        expectedscan,
        localstepgap,
        :not_evaluated
    )
    localstepgap > 0 && isfinite(localstepgap) || return AlkaneLadderAdditionEvaluation(
        nothing,
        alkane_ladder_addition_diagnostic_reason(diagnostic, :invalid_anchor_spacing)
    )

    alkane_best_ladder_addition_window(
        msm,
        variances,
        abundanceinfo,
        targetstep,
        source,
        expectedscan,
        localstepgap,
        scancount,
        settings,
        massspectrumcontext,
        standard,
        direction,
        nothing,
        nothing,
        edgeanchor,
        edgeanchors
    )
end

function alkane_best_ladder_addition_window(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    abundanceinfo::AlkaneAbundanceInfo,
    targetstep::Integer,
    source::Symbol,
    expectedscan::Real,
    localstepgap::Real,
    scancount::Integer,
    settings::AlkaneLadderAdditionSettings,
    massspectrumcontext::Union{Nothing, AlkaneLadderMassSpectrumMatchContext},
    standard::AlkaneStandard,
    direction::Symbol,
    leftanchor::Union{Nothing, AlkaneLadderAdditionAnchor},
    rightanchor::Union{Nothing, AlkaneLadderAdditionAnchor},
    edgeanchor::Union{Nothing, AlkaneLadderAdditionAnchor},
    edgeanchors::AbstractVector{AlkaneLadderAdditionAnchor}
)
    diagnostic = alkane_ladder_addition_diagnostic_base(
        targetstep,
        source,
        expectedscan,
        localstepgap,
        :not_evaluated
    )
    isfinite(expectedscan) || return AlkaneLadderAdditionEvaluation(
        nothing,
        alkane_ladder_addition_diagnostic_reason(diagnostic, :nonfinite_expected_scan)
    )

    isedge = source in (:leftextended, :rightextended)
    minradius = isedge ? settings.edgeminradius : settings.minradius
    radiusfraction = isedge ? settings.edgeradiusfraction : settings.radiusfraction
    positionsigmafraction = isedge ?
        settings.edgepositionsigmafraction :
        settings.positionsigmafraction
    searchradius = max(minradius, radiusfraction * localstepgap)
    positionsigma = positionsigmafraction * localstepgap
    diagnostic = alkane_ladder_addition_diagnostic_search(
        diagnostic,
        searchradius,
        positionsigma
    )
    searchradius > 0 && isfinite(searchradius) || return AlkaneLadderAdditionEvaluation(
        nothing,
        alkane_ladder_addition_diagnostic_reason(diagnostic, :invalid_search_radius)
    )
    positionsigma > 0 && isfinite(positionsigma) || return AlkaneLadderAdditionEvaluation(
        nothing,
        alkane_ladder_addition_diagnostic_reason(diagnostic, :invalid_position_sigma)
    )

    windows = get(abundanceinfo.windows, targetstep, AlkaneAbundanceWindow[])
    abundance = get(abundanceinfo.abundances, targetstep, nothing)
    diagnostic = alkane_ladder_addition_diagnostic_candidate_windows(
        diagnostic,
        length(windows)
    )
    isempty(windows) && return AlkaneLadderAdditionEvaluation(
        nothing,
        alkane_ladder_addition_diagnostic_reason(diagnostic, :no_abundance_windows)
    )
    isnothing(abundance) && return AlkaneLadderAdditionEvaluation(
        nothing,
        alkane_ladder_addition_diagnostic_reason(diagnostic, :no_abundance_vector)
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
    requiredcosine = if source ≡ :gapfilled
        alkane_single_gap_required_cosine(
            leftanchor,
            rightanchor,
            settings.gapmincosinefloor,
            settings.gapcosinetolerance
        )
    else
        alkane_edge_extension_required_cosine(
            edgeanchors,
            direction,
            settings.edgemincosinefloor,
            settings.edgecosinetolerance,
            settings.edgecosineanchorcount
        )
    end
    for window in windows
        apexindex = window.apexindex
        alkane_ladder_addition_direction_passes(
            apexindex,
            direction,
            expectedscan,
            scancount,
            leftanchor,
            rightanchor,
            edgeanchor
        ) || continue
        directionok += 1
        scanerror = apexindex - expectedscan
        abs(scanerror) ≤ searchradius || continue
        inradius += 1
        if isedge && alkane_ladder_edge_window_apex_is_boundary_truncated(
                abundance,
                window,
                direction,
                scancount,
                settings.apexmaxshiftfromguess
            )
            continue
        end
        completewindows += 1

        apexabundance = window.apexabundance
        spectrumcandidate = alkane_ladder_addition_spectrum_candidate(
            abundance,
            window,
            targetstep
        )
        isnothing(spectrumcandidate) && continue
        validpeakmodelwindows += 1
        if isedge && spectrumcandidate.apexindex != apexindex
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
            match.cosine ≥ requiredcosine || continue
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
        candidate = AlkaneLadderAdditionCandidate(
            targetstep,
            spectrumcandidate.leftindex,
            spectrumcandidate.apexindex,
            spectrumcandidate.rightindex,
            spectrumcandidate.scanindices,
            spectrumcandidate.peakmodel,
            source,
            apexindex,
            expectedscan,
            scanerror,
            searchradius,
            positionsigma,
            localstepgap,
            score,
            apexabundance,
            match.cosine,
            match.distance,
            match.ioncount,
            requiredcosine,
            window
        )
        apexsettings = alkane_ladder_addition_apex_settings(settings, standard)
        apex = alkaneladderapex(
            msm,
            variances,
            abundanceinfo,
            candidate,
            apexsettings
        )
        apex.success || continue
        rank = (score, match.cosine, apexabundance)
        (isnothing(bestrank) || rank > bestrank) || continue
        bestrank = rank
        best = AlkaneLadderAddition(
            candidate,
            apex
        )
    end

    diagnostic = alkane_ladder_addition_diagnostic_counts(
        diagnostic,
        directionok,
        inradius,
        completewindows,
        validpeakmodelwindows,
        apexconsistentwindows,
        finitecosine,
        passingcosine,
        requiredcosine
    )

    if !isnothing(best)
        return AlkaneLadderAdditionEvaluation(
            best,
            alkane_ladder_addition_diagnostic_accepted(diagnostic, best)
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

    AlkaneLadderAdditionEvaluation(
        nothing,
        alkane_ladder_addition_diagnostic_reason(diagnostic, reason)
    )
end

function alkane_ladder_addition_direction_passes(
    apexindex::Integer,
    direction::Symbol,
    expectedscan::Real,
    scancount::Integer,
    leftanchor::Union{Nothing, AlkaneLadderAdditionAnchor},
    rightanchor::Union{Nothing, AlkaneLadderAdditionAnchor},
    edgeanchor::Union{Nothing, AlkaneLadderAdditionAnchor}
)
    1 ≤ apexindex ≤ scancount || return false
    if direction ≡ :between
        isnothing(leftanchor) && return true
        isnothing(rightanchor) && return true
        return alkane_ladder_result_scanindex(leftanchor) <
            apexindex <
            alkane_ladder_result_scanindex(rightanchor)
    end
    if direction ≡ :left
        isnothing(edgeanchor) && return apexindex ≤ ceil(Int, expectedscan)
        return apexindex < alkane_ladder_result_scanindex(edgeanchor)
    elseif direction ≡ :right
        isnothing(edgeanchor) && return apexindex ≥ floor(Int, expectedscan)
        return apexindex > alkane_ladder_result_scanindex(edgeanchor)
    end

    false
end

alkane_ladder_result_scanindex(anchor::AlkaneLadderAdditionAnchor) =
    anchor.apexscanindex

function alkane_ladder_edge_extension_anchor_results(
    anchors::AbstractVector{AlkaneLadderAdditionAnchor},
    direction::Symbol,
    maxanchors::Integer
)
    sorted = sort(collect(anchors); by=anchor -> anchor.ladderstep)
    if direction ≡ :left
        return first(sorted, min(length(sorted), maxanchors))
    elseif direction ≡ :right
        return last(sorted, min(length(sorted), maxanchors))
    end

    throw(ArgumentError("direction must be :left or :right"))
end

function alkane_single_gap_required_cosine(
    left::Union{Nothing, AlkaneLadderAdditionAnchor},
    right::Union{Nothing, AlkaneLadderAdditionAnchor},
    mincosinefloor::Real,
    cosinetolerance::Real
)
    localcosines = Float64[]
    leftcosine = alkane_ladder_anchor_mass_spectrum_cosine(left)
    rightcosine = alkane_ladder_anchor_mass_spectrum_cosine(right)
    isfinite(leftcosine) && push!(localcosines, leftcosine)
    isfinite(rightcosine) && push!(localcosines, rightcosine)
    isempty(localcosines) && return mincosinefloor

    max(mincosinefloor, minimum(localcosines) - cosinetolerance)
end

function alkane_edge_extension_required_cosine(
    anchors::AbstractVector{AlkaneLadderAdditionAnchor},
    direction::Symbol,
    mincosinefloor::Real,
    cosinetolerance::Real,
    cosineanchorcount::Integer
)
    localanchors = direction ≡ :left ?
        first(anchors, min(length(anchors), cosineanchorcount)) :
        last(anchors, min(length(anchors), cosineanchorcount))
    localcosines = Float64[]
    for anchor in localanchors
        cosine = alkane_ladder_anchor_mass_spectrum_cosine(anchor)
        isfinite(cosine) && push!(localcosines, cosine)
    end
    isempty(localcosines) && return mincosinefloor

    max(mincosinefloor, minimum(localcosines) - cosinetolerance)
end

alkane_ladder_anchor_mass_spectrum_cosine(::Nothing) = NaN

alkane_ladder_anchor_mass_spectrum_cosine(anchor::AlkaneLadderAdditionAnchor) =
    anchor.mass_spectrum_cosine

function alkane_ladder_edge_extension_scan_prediction(
    anchors::AbstractVector{AlkaneLadderAdditionAnchor},
    targetstep::Integer
)
    length(anchors) ≥ 2 || return nothing
    sorted = sort(collect(anchors); by=anchor -> anchor.ladderstep)
    if targetstep > last(sorted).ladderstep
        edge = last(sorted)
        neighbor = sorted[end - 1]
        edge_step = edge.ladderstep
        neighbor_step = neighbor.ladderstep
        step_delta = edge_step - neighbor_step
        step_delta > 0 || return nothing
        latest_gap =
            (alkane_ladder_result_scanindex(edge) -
             alkane_ladder_result_scanindex(neighbor)) / step_delta
        isfinite(latest_gap) && latest_gap > 0 || return nothing
        predicted_gap = latest_gap
        if length(sorted) ≥ 3
            previous = sorted[end - 2]
            previous_step = previous.ladderstep
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
            predicted_gap * (targetstep - edge_step)
    elseif targetstep < first(sorted).ladderstep
        edge = first(sorted)
        neighbor = sorted[2]
        edge_step = edge.ladderstep
        neighbor_step = neighbor.ladderstep
        step_delta = neighbor_step - edge_step
        step_delta > 0 || return nothing
        latest_gap =
            (alkane_ladder_result_scanindex(neighbor) -
             alkane_ladder_result_scanindex(edge)) / step_delta
        isfinite(latest_gap) && latest_gap > 0 || return nothing
        predicted_gap = latest_gap
        if length(sorted) ≥ 3
            previous = sorted[3]
            previous_step = previous.ladderstep
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
            predicted_gap * (edge_step - targetstep)
    else
        return nothing
    end

    AlkaneLadderEdgeScanPrediction(
        expected_scan,
        0.0,
        predicted_gap,
        1
    )
end

function alkane_ladder_edge_window_apex_is_boundary_truncated(
    abundance::AbstractVector{<:Real},
    window::AlkaneAbundanceWindow,
    direction::Symbol,
    scancount::Integer,
    maxboundaryapexdistance::Real
)
    isfinite(maxboundaryapexdistance) && maxboundaryapexdistance ≥ 0 ||
        throw(ArgumentError(
            "maxboundaryapexdistance must be finite and nonnegative"))
    1 ≤ window.leftindex ≤ window.rightindex ≤ length(abundance) ||
        return false

    scanindices = window.leftindex:window.rightindex
    peakmodel = max.(Float64.(abundance[scanindices]), 0.0)
    maximum(peakmodel) > 0 || return false
    _, peakoffset = findmax(peakmodel)
    windowmaxscan = first(scanindices) + peakoffset - 1

    if direction ≡ :right
        touchesboundary = window.rightindex ≥ scancount ||
            window.rightstop ≡ :boundary
        boundarydistance = scancount - windowmaxscan
        return touchesboundary &&
            boundarydistance ≤ maxboundaryapexdistance + 1e-9
    elseif direction ≡ :left
        touchesboundary = window.leftindex ≤ 1 ||
            window.leftstop ≡ :boundary
        boundarydistance = windowmaxscan - 1
        return touchesboundary &&
            boundarydistance ≤ maxboundaryapexdistance + 1e-9
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

    apexabundance * exp(-0.5 * abs2(scanerror / positionsigma))
end

function alkane_ladder_position_penalized_cosine_score(
    cosine::Real,
    scanerror::Real,
    positionsigma::Real
)
    isfinite(cosine) && cosine > 0 || return -Inf
    isfinite(scanerror) || return -Inf
    isfinite(positionsigma) && positionsigma > 0 || return -Inf

    z = scanerror / positionsigma

    log(cosine) - 0.5 * abs2(z)
end

function alkane_ladder_addition_spectrum_candidate(
    abundance::AbstractVector{<:Real},
    window::AlkaneAbundanceWindow,
    step::Integer
)
    1 ≤ window.leftindex ≤ window.rightindex ≤ length(abundance) ||
        return nothing

    scanindices = window.leftindex:window.rightindex
    peakmodel = max.(Float64.(abundance[scanindices]), 0.0)
    peakheight, peakoffset = findmax(peakmodel)
    peakheight > 0 && isfinite(peakheight) || return nothing
    peakmodel ./= peakheight

    AlkaneLadderAdditionSpectrumCandidate(
        step,
        window.leftindex,
        first(scanindices) + peakoffset - 1,
        window.rightindex,
        collect(scanindices),
        collect(peakmodel),
        window
    )
end

function alkane_ladder_mass_spectrum_match_context(
    msm::MassScanMatrix,
    variances::AbstractMatrix{<:Real},
    standard::AlkaneStandard,
    enabled::Bool,
    variancefloor::Real
)
    enabled || return nothing
    validate_alkane_series_variances(msm, variances)
    validate_alkane_abundance_variancefloor(variancefloor)

    references = alkane_ladder_mass_spectrum_reference_vectors(msm, standard)
    isempty(references) && return nothing

    AlkaneLadderMassSpectrumMatchContext(
        rawintensities(msm),
        variances,
        Float64(variancefloor),
        references
    )
end

function alkane_ladder_mass_spectrum_reference_vectors(
    msm::MassScanMatrix,
    standard::AlkaneStandard
)
    spectra = alkane_standard_spectra(standard)
    mzbins = alkane_mz_bins(msm)
    indexbymz = Dict(mz => index for (index, mz) in pairs(mzbins))
    references = Dict{Int, Vector{Float64}}()

    for spectrum in spectra
        spectrumattrs = attrs(spectrum)
        :order in keys(spectrumattrs) || continue
        carbon = spectrumattrs.order
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
        references[carbon] = reference
    end

    references
end

function alkane_ladder_candidate_mass_spectrum_match(
    candidate::Union{AlkaneLadderAdditionSpectrumCandidate, AlkaneMolecularIonContrast},
    ::Nothing
)
    AlkaneLadderMassSpectrumMatch(
        NaN,
        NaN,
        0
    )
end

function alkane_ladder_candidate_mass_spectrum_match(
    candidate::Union{AlkaneLadderAdditionSpectrumCandidate, AlkaneMolecularIonContrast},
    context::AlkaneLadderMassSpectrumMatchContext
)
    reference = get(context.references, candidate.ladderstep, nothing)
    isnothing(reference) && return AlkaneLadderMassSpectrumMatch(
        NaN,
        NaN,
        0
    )

    observed = Vector{Float64}(undef, length(reference))
    for mzindex in eachindex(reference)
        fit = alkane_fitted_ion_abundance(
            context.X,
            context.variances,
            candidate.scanindices,
            mzindex,
            candidate.peakmodel,
            context.variancefloor
        )
        observed[mzindex] = max(fit.abundance, 0.0)
    end

    cosine = alkane_ladder_cosine_similarity(observed, reference)
    distance = isfinite(cosine) ? 1.0 - clamp(cosine, 0.0, 1.0) : NaN

    AlkaneLadderMassSpectrumMatch(
        cosine,
        distance,
        Base.count(>(0.0), reference)
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
    expectedscan::Real,
    localstepgap::Real,
    reason::Symbol
)
    AlkaneLadderAdditionDiagnostic(
        :failed,
        reason,
        targetstep,
        source,
        Float64(expectedscan),
        Float64(localstepgap),
        NaN,
        NaN,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        NaN,
        missing,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN
    )
end

function alkane_ladder_addition_diagnostic_reason(
    diagnostic::AlkaneLadderAdditionDiagnostic,
    reason::Symbol
)
    AlkaneLadderAdditionDiagnostic(
        diagnostic.status,
        reason,
        diagnostic.ladderstep,
        diagnostic.source,
        diagnostic.expectedscan,
        diagnostic.localstepgap,
        diagnostic.searchradius,
        diagnostic.positionsigma,
        diagnostic.candidatewindows,
        diagnostic.directionwindows,
        diagnostic.inradiuswindows,
        diagnostic.completewindows,
        diagnostic.validpeakmodelwindows,
        diagnostic.apexconsistentwindows,
        diagnostic.finitecosinewindows,
        diagnostic.passingcosinewindows,
        diagnostic.requiredcosine,
        diagnostic.bestscanindex,
        diagnostic.bestscanerror,
        diagnostic.bestscore,
        diagnostic.bestapexabundance,
        diagnostic.bestapexscanindex,
        diagnostic.bestcosine
    )
end

function alkane_ladder_addition_diagnostic_search(
    diagnostic::AlkaneLadderAdditionDiagnostic,
    searchradius::Real,
    positionsigma::Real
)
    AlkaneLadderAdditionDiagnostic(
        diagnostic.status,
        diagnostic.reason,
        diagnostic.ladderstep,
        diagnostic.source,
        diagnostic.expectedscan,
        diagnostic.localstepgap,
        Float64(searchradius),
        Float64(positionsigma),
        diagnostic.candidatewindows,
        diagnostic.directionwindows,
        diagnostic.inradiuswindows,
        diagnostic.completewindows,
        diagnostic.validpeakmodelwindows,
        diagnostic.apexconsistentwindows,
        diagnostic.finitecosinewindows,
        diagnostic.passingcosinewindows,
        diagnostic.requiredcosine,
        diagnostic.bestscanindex,
        diagnostic.bestscanerror,
        diagnostic.bestscore,
        diagnostic.bestapexabundance,
        diagnostic.bestapexscanindex,
        diagnostic.bestcosine
    )
end

function alkane_ladder_addition_diagnostic_candidate_windows(
    diagnostic::AlkaneLadderAdditionDiagnostic,
    candidatewindows::Integer
)
    AlkaneLadderAdditionDiagnostic(
        diagnostic.status,
        diagnostic.reason,
        diagnostic.ladderstep,
        diagnostic.source,
        diagnostic.expectedscan,
        diagnostic.localstepgap,
        diagnostic.searchradius,
        diagnostic.positionsigma,
        candidatewindows,
        diagnostic.directionwindows,
        diagnostic.inradiuswindows,
        diagnostic.completewindows,
        diagnostic.validpeakmodelwindows,
        diagnostic.apexconsistentwindows,
        diagnostic.finitecosinewindows,
        diagnostic.passingcosinewindows,
        diagnostic.requiredcosine,
        diagnostic.bestscanindex,
        diagnostic.bestscanerror,
        diagnostic.bestscore,
        diagnostic.bestapexabundance,
        diagnostic.bestapexscanindex,
        diagnostic.bestcosine
    )
end

function alkane_ladder_addition_diagnostic_counts(
    diagnostic::AlkaneLadderAdditionDiagnostic,
    directionwindows::Integer,
    inradiuswindows::Integer,
    completewindows::Integer,
    validpeakmodelwindows::Integer,
    apexconsistentwindows::Integer,
    finitecosinewindows::Integer,
    passingcosinewindows::Integer,
    requiredcosine::Real
)
    AlkaneLadderAdditionDiagnostic(
        diagnostic.status,
        diagnostic.reason,
        diagnostic.ladderstep,
        diagnostic.source,
        diagnostic.expectedscan,
        diagnostic.localstepgap,
        diagnostic.searchradius,
        diagnostic.positionsigma,
        diagnostic.candidatewindows,
        directionwindows,
        inradiuswindows,
        completewindows,
        validpeakmodelwindows,
        apexconsistentwindows,
        finitecosinewindows,
        passingcosinewindows,
        Float64(requiredcosine),
        diagnostic.bestscanindex,
        diagnostic.bestscanerror,
        diagnostic.bestscore,
        diagnostic.bestapexabundance,
        diagnostic.bestapexscanindex,
        diagnostic.bestcosine
    )
end

function alkane_ladder_addition_diagnostic_accepted(
    diagnostic::AlkaneLadderAdditionDiagnostic,
    addition::AlkaneLadderAddition
)
    AlkaneLadderAdditionDiagnostic(
        :accepted,
        :passed,
        diagnostic.ladderstep,
        diagnostic.source,
        diagnostic.expectedscan,
        diagnostic.localstepgap,
        diagnostic.searchradius,
        diagnostic.positionsigma,
        diagnostic.candidatewindows,
        diagnostic.directionwindows,
        diagnostic.inradiuswindows,
        diagnostic.completewindows,
        diagnostic.validpeakmodelwindows,
        diagnostic.apexconsistentwindows,
        diagnostic.finitecosinewindows,
        diagnostic.passingcosinewindows,
        diagnostic.requiredcosine,
        addition.scanindex,
        addition.scanerror,
        addition.score,
        addition.apexabundance,
        addition.apexscanindex,
        addition.massspectrumcosine
    )
end
