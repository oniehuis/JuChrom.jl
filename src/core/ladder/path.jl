struct LadderWindowPathState
    pos::Int
    prevpos::Int
    previdx::Int
    curridx::Int
    objective::Float64
    pathlength::Int
    parent::Union{Nothing, Int}
end

struct AlkaneLadderPathCandidate{T<:AlkaneMolecularIonContrast}
    ladderstep::Int
    scanindex::Int
    score::Float64
    z::Float64
    centerz::Float64
    isolationz::Float64
    centervslowerz::Float64
    centervsupperz::Float64
    molecularionscore::Float64
    massspectrumcosine::Float64
    massspectrumdistance::Float64
    massspectrumioncount::Int
    source::Symbol
    misupported::Bool
    gapfilled::Bool
    expectedscan::Float64
    scanerror::Float64
    requiredcosine::Float64
    window::T
end

struct AlkaneLadderPathRunBest{T<:AlkaneLadderPathCandidate}
    objective::Float64
    path::Vector{T}
end

struct AlkaneLadderPathRunResult
    laddersteps::Vector{Int}
    success::Bool
    objective::Union{Missing, Float64}
    pathlength::Int
end

struct AlkaneLadderPathInfo{
    T1<:AlkaneLadderPathCandidate,
    T2<:AbstractDict{Int, <:AbstractVector{<:AlkaneLadderPathCandidate}},
    T3<:AbstractVector{<:AbstractVector{Int}},
    T4<:NamedTuple
}
    success::Bool
    status::Symbol
    reason::Symbol
    failurereason::Union{Nothing, AbstractString}
    message::AbstractString
    laddersteps::Vector{Int}
    scanindices::Vector{Int}
    scores::Vector{Float64}
    zscores::Vector{Float64}
    centerzs::Vector{Float64}
    isolationzs::Vector{Float64}
    centervslowerzs::Vector{Float64}
    centervsupperzs::Vector{Float64}
    massspectrumcosines::Vector{Float64}
    massspectrumdistances::Vector{Float64}
    sources::Vector{Symbol}
    misupported::Vector{Bool}
    gapfilled::Vector{Bool}
    expectedscanindices::Vector{Float64}
    scanerrors::Vector{Float64}
    requiredcosines::Vector{Float64}
    gaps::Vector{Int}
    stepgaps::Vector{Int}
    scanstepgaps::Vector{Float64}
    missingsteps::Int
    objective::Float64
    score::Float64
    windows::Vector{<:AlkaneMolecularIonContrast}
    path::Vector{T1}
    candidates::T2
    candidatesbystep::T2
    candidateruns::T3
    runresults::Vector{AlkaneLadderPathRunResult}
    settings::T4
end

function alkaneladderpath(
    molecularioninfo::AlkaneMolecularIonInfo;
    kwargs...
)
    alkaneladderpath(molecularioninfo.contrasts; kwargs...)
end

function alkaneladderpath(
    contrasts::AbstractDict{Int, <:AbstractVector{<:AlkaneMolecularIonContrast}};
    centerzmin::Real,
    isolationzmin::Real,
    minsteps::Integer,
    stepreward::Real,
    maxcandidatesperstep::Union{Nothing, Integer},
    spacingweight::Real,
    gapincreaseweight::Real,
    maxgapratio::Real,
    maxmissingsteps::Integer,
    missingsteppenalty::Real,
    msm::Union{Nothing, MassScanMatrix},
    variances::Union{Nothing, AbstractMatrix{<:Real}},
    standard::Union{Nothing, AlkaneStandard},
    massspectrummatch::Bool,
    massspectrummatchdistanceweight::Real
)
    validate_alkane_ladder_path_settings(
        centerzmin,
        isolationzmin,
        minsteps,
        stepreward,
        maxcandidatesperstep,
        spacingweight,
        gapincreaseweight,
        maxgapratio,
        maxmissingsteps,
        missingsteppenalty,
        massspectrummatchdistanceweight
    )

    context = alkane_ladder_path_mass_spectrum_match_context(
        msm,
        variances,
        standard,
        massspectrummatch
    )
    settings = (
        centerzmin=centerzmin,
        isolationzmin=isolationzmin,
        minsteps=minsteps,
        stepreward=stepreward,
        maxcandidatesperstep=maxcandidatesperstep,
        spacingweight=spacingweight,
        gapincreaseweight=gapincreaseweight,
        maxgapratio=maxgapratio,
        maxmissingsteps=maxmissingsteps,
        missingsteppenalty=missingsteppenalty,
        massspectrummatch=massspectrummatch,
        massspectrummatchactive=!isnothing(context),
        massspectrummatchdistanceweight=massspectrummatchdistanceweight
    )

    candidatesbystep = alkane_ladder_path_candidates(
        contrasts,
        centerzmin,
        isolationzmin,
        maxcandidatesperstep,
        context,
        massspectrummatchdistanceweight
    )
    candidateruns = alkane_ladder_candidate_runs(
        candidatesbystep,
        maxmissingsteps
    )
    isempty(candidateruns) && return alkane_ladder_path_failure(
        :no_candidates,
        "no molecular-ion candidates pass centerzmin=$(centerzmin) and isolationzmin=$(isolationzmin)",
        candidatesbystep,
        candidateruns,
        AlkaneLadderPathRunResult[],
        settings
    )

    best = nothing
    runresults = AlkaneLadderPathRunResult[]
    for run in candidateruns
        runbest = alkane_best_ladder_path_in_run(
            run,
            candidatesbystep,
            minsteps,
            stepreward,
            spacingweight,
            gapincreaseweight,
            maxgapratio,
            maxmissingsteps,
            missingsteppenalty
        )
        push!(
            runresults,
            AlkaneLadderPathRunResult(
                run,
                !isnothing(runbest),
                isnothing(runbest) ? missing : runbest.objective,
                isnothing(runbest) ? 0 : length(runbest.path)
            )
        )
        if !isnothing(runbest) && (isnothing(best) || runbest.objective > best.objective)
            best = runbest
        end
    end

    isnothing(best) && return alkane_ladder_path_failure(
        :no_valid_path,
        "no valid path with at least $(minsteps) strictly increasing steps",
        candidatesbystep,
        candidateruns,
        runresults,
        settings
    )

    alkane_ladder_path_success(best, candidatesbystep, candidateruns, runresults, settings)
end

function validate_alkane_ladder_path_settings(
    centerzmin::Real,
    isolationzmin::Real,
    minsteps::Integer,
    stepreward::Real,
    maxcandidatesperstep::Union{Nothing, Integer},
    spacingweight::Real,
    gapincreaseweight::Real,
    maxgapratio::Real,
    maxmissingsteps::Integer,
    missingsteppenalty::Real,
    massspectrummatchdistanceweight::Real
)
    isfinite(centerzmin) && centerzmin ≥ 0 || throw(ArgumentError(
        "centerzmin must be finite and nonnegative"))
    isfinite(isolationzmin) && isolationzmin ≥ 0 || throw(ArgumentError(
        "isolationzmin must be finite and nonnegative"))
    minsteps ≥ 1 || throw(ArgumentError("minsteps must be at least 1"))
    isfinite(stepreward) || throw(ArgumentError("stepreward must be finite"))
    isnothing(maxcandidatesperstep) ||
        maxcandidatesperstep ≥ 1 ||
        throw(ArgumentError("maxcandidatesperstep must be positive"))
    isfinite(spacingweight) && spacingweight ≥ 0 || throw(ArgumentError(
        "spacingweight must be finite and nonnegative"))
    isfinite(gapincreaseweight) && gapincreaseweight ≥ 0 || throw(ArgumentError(
        "gapincreaseweight must be finite and nonnegative"))
    isfinite(maxgapratio) && maxgapratio ≥ 1 || throw(ArgumentError(
        "maxgapratio must be finite and at least 1"))
    maxmissingsteps in (0, 1) || throw(ArgumentError(
        "maxmissingsteps must be 0 or 1"))
    isfinite(missingsteppenalty) && missingsteppenalty ≥ 0 || throw(ArgumentError(
        "missingsteppenalty must be finite and nonnegative"))
    isfinite(massspectrummatchdistanceweight) &&
        massspectrummatchdistanceweight ≥ 0 ||
        throw(ArgumentError(
            "mass-spectrum match distance weight must be finite and nonnegative"))

    nothing
end

function alkane_ladder_path_mass_spectrum_match_context(
    msm::Union{Nothing, MassScanMatrix},
    variances::Union{Nothing, AbstractMatrix{<:Real}},
    standard::Union{Nothing, AlkaneStandard},
    enabled::Bool
)
    enabled || return nothing
    (isnothing(msm) || isnothing(variances) || isnothing(standard)) && return nothing

    alkane_ladder_mass_spectrum_match_context(msm, variances, standard, true, 1.0)
end

function alkane_ladder_path_candidates(
    contrasts::AbstractDict{Int, <:AbstractVector{T}},
    centerzmin::Real,
    isolationzmin::Real,
    maxcandidatesperstep::Union{Nothing,Integer},
    massspectrummatchcontext,
    massspectrummatchdistanceweight::Real
) where {T<:AlkaneMolecularIonContrast}

    candidatesbystep = Dict{Int, Vector{AlkaneLadderPathCandidate{T}}}()
    for step in sort(collect(keys(contrasts)))
        candidates = AlkaneLadderPathCandidate{T}[]
        for contrast in contrasts[step]
            alkane_molecular_ion_window_has_candidate_evidence(
                contrast,
                centerzmin,
                isolationzmin
            ) || continue
            match = alkane_ladder_candidate_mass_spectrum_match(
                contrast,
                massspectrummatchcontext
            )
            candidate_score = alkane_ladder_path_candidate_dp_score(
                contrast.molecularionscore,
                match.distance,
                massspectrummatchdistanceweight
            )
            push!(candidates, AlkaneLadderPathCandidate(
                    step,
                    contrast.apexindex,
                    candidate_score,
                    contrast.z,
                    contrast.centerz,
                    contrast.isolationz,
                    contrast.centervslowerz,
                    contrast.centervsupperz,
                    contrast.molecularionscore,
                    match.cosine,
                    match.distance,
                    match.ioncount,
                    :molecular_ion_dp,
                    true,
                    false,
                    NaN,
                    NaN,
                    NaN,
                    contrast
                )
            )
        end
        candidatesbystep[step] = alkane_limit_ladder_path_candidates(
            candidates,
            maxcandidatesperstep
        )
    end

    candidatesbystep
end

function alkane_molecular_ion_window_has_candidate_evidence(
    contrast::AlkaneMolecularIonContrast,
    centerzmin::Real,
    isolationzmin::Real
)
    isfinite(contrast.centerz) && contrast.centerz ≥ centerzmin || return false
    isfinite(contrast.isolationz) && contrast.isolationz ≥ isolationzmin || return false
    isfinite(contrast.molecularionscore) && contrast.molecularionscore > 0 || return false

    true
end

function alkane_ladder_path_candidate_dp_score(
    score::Real,
    distance::Real,
    distanceweight::Real
)
    isfinite(score) && score > 0 || return 0.0
    isfinite(distance) || return score

    score * exp(-distanceweight * clamp(distance, 0.0, 1.0))
end

function alkane_limit_ladder_path_candidates(
    candidates::Vector{<:AlkaneLadderPathCandidate},
    maxcandidates::Union{Nothing,Integer}
)
    sorted = sort(candidates; by=candidate -> candidate.scanindex)
    isnothing(maxcandidates) && return sorted
    length(sorted) ≤ maxcandidates && return sorted

    order = sortperm(sorted; by=candidate -> candidate.score, rev=true)
    limited = sorted[order[1:maxcandidates]]
    sort!(limited; by=candidate -> candidate.scanindex)

    limited
end

function alkane_ladder_candidate_runs(
    candidatesbystep::AbstractDict{Int, <:AbstractVector{<:AlkaneLadderPathCandidate}},
    maxmissingsteps::Integer
)
    maxmissingsteps in (0, 1) || throw(ArgumentError("maxmissingsteps must be 0 or 1"))

    runs = Vector{Vector{Int}}()
    current = Int[]
    for step in sort(collect(keys(candidatesbystep)))
        isempty(candidatesbystep[step]) && continue
        if isempty(current) || step - last(current) ≤ maxmissingsteps + 1
            push!(current, step)
        else
            push!(runs, current)
            current = [step]
        end
    end
    !isempty(current) && push!(runs, current)

    runs
end

function alkane_best_ladder_path_in_run(
    steprun::AbstractVector{<:Integer},
    candidatesbystep::AbstractDict{Int, <:AbstractVector{<:AlkaneLadderPathCandidate}},
    minsteps::Integer,
    stepreward::Real,
    spacingweight::Real,
    gapincreaseweight::Real,
    maxgapratio::Real,
    maxmissingsteps::Integer,
    missingsteppenalty::Real
)
    maxmissingsteps in (0, 1) || throw(ArgumentError("maxmissingsteps must be 0 or 1"))

    candidates = [candidatesbystep[step] for step in steprun]
    length(candidates) < minsteps && return nothing
    candidate_scanindices = [
        [candidate.scanindex for candidate in stepcandidates]
        for stepcandidates in candidates
    ]
    maxgaplogratio = log(maxgapratio)

    bestobjective = -Inf
    beststateid = nothing
    states = LadderWindowPathState[]
    statesbypos = [Dict{Tuple{Int, Int, Int, Int}, Int}() for _ in candidates]

    for pos in eachindex(candidates)
        currentstates = statesbypos[pos]
        current_scanindices = candidate_scanindices[pos]
        for curridx in eachindex(candidates[pos])
            curr = candidates[pos][curridx]
            objective = curr.score + stepreward
            stateid = alkane_push_ladder_path_state!(
                states,
                currentstates,
                pos,
                0,
                0,
                curridx,
                objective,
                1,
                minsteps,
                nothing
            )
            if minsteps ≤ 1 && states[stateid].objective > bestobjective
                bestobjective = states[stateid].objective
                beststateid = stateid
            end
        end

        for prevpos in Iterators.reverse(1:(pos - 1))
            stepgap = steprun[pos] - steprun[prevpos]
            stepgap ≤ maxmissingsteps + 1 || break
            stepgap ≥ 1 || continue
            missingsteps = stepgap - 1

            for stateid in values(statesbypos[prevpos])
                state = states[stateid]
                previous = candidates[prevpos][state.curridx]
                firstcurridx = searchsortedlast(current_scanindices, previous.scanindex) + 1
                for curridx in firstcurridx:length(candidates[pos])
                    curr = candidates[pos][curridx]

                    penalty = missingsteppenalty * missingsteps
                    if state.pathlength >= 2
                        parentid = state.parent
                        isnothing(parentid) && continue
                        parent = states[parentid]
                        prevprev = candidates[parent.pos][parent.curridx]
                        previousstepgap = steprun[state.pos] - steprun[parent.pos]
                        nextstepgap = steprun[pos] - steprun[state.pos]
                        previousgap = (previous.scanindex - prevprev.scanindex) /
                            previousstepgap
                        nextgap = (curr.scanindex - previous.scanindex) / nextstepgap
                        gaplogratio = alkane_ladder_gap_log_ratio(previousgap, nextgap)
                        isfinite(gaplogratio) || continue
                        abs(gaplogratio) ≤ maxgaplogratio || continue
                        penalty += spacingweight * abs2(gaplogratio)
                        penalty += gapincreaseweight * abs2(max(0.0, gaplogratio))
                    end
                    isfinite(penalty) || continue

                    objective = state.objective + curr.score + stepreward - penalty
                    newstateid = alkane_push_ladder_path_state!(
                        states,
                        currentstates,
                        pos,
                        prevpos,
                        state.curridx,
                        curridx,
                        objective,
                        state.pathlength + 1,
                        minsteps,
                        stateid
                    )
                    newstate = states[newstateid]
                    if newstate.pathlength ≥ minsteps && newstate.objective > bestobjective
                        bestobjective = newstate.objective
                        beststateid = newstateid
                    end
                end
            end
        end
    end

    isnothing(beststateid) && return nothing

    AlkaneLadderPathRunBest(
        bestobjective,
        alkane_reconstruct_ladder_path(states, beststateid, candidates)
    )
end

function alkane_push_ladder_path_state!(
    states::Vector{LadderWindowPathState},
    statemap::Dict{Tuple{Int, Int, Int, Int}, Int},
    pos::Int,
    prevpos::Int,
    previdx::Int,
    curridx::Int,
    objective::Float64,
    pathlength::Int,
    lengthcap::Integer,
    parent::Union{Nothing, Int}
)
    key = (prevpos, previdx, curridx, min(pathlength, lengthcap))
    existing = get(statemap, key, nothing)
    if !isnothing(existing) && states[existing].objective ≥ objective
        return existing
    end

    push!(states, LadderWindowPathState(
            pos,
            prevpos,
            previdx,
            curridx,
            objective,
            pathlength,
            parent
        )
    )
    stateid = length(states)
    statemap[key] = stateid

    stateid
end

function alkane_reconstruct_ladder_path(
    states::Vector{LadderWindowPathState},
    stateid::Int,
    candidates::AbstractVector{<:AbstractVector{T}}
) where {T<:AlkaneLadderPathCandidate}

    path = T[]
    alkane_reconstruct_ladder_path!(path, states, stateid, candidates)

    path
end

function alkane_reconstruct_ladder_path!(
    path::Vector{T},
    states::Vector{LadderWindowPathState},
    stateid::Int,
    candidates::AbstractVector{<:AbstractVector{T}}
) where {T<:AlkaneLadderPathCandidate}

    state = states[stateid]
    if !isnothing(state.parent)
        alkane_reconstruct_ladder_path!(path, states, state.parent, candidates)
    end
    push!(path, candidates[state.pos][state.curridx])

    path
end

function alkane_ladder_gap_log_ratio(previousgap::Real, nextgap::Real)
    previousgap > 0 && nextgap > 0 || return Inf
    isfinite(previousgap) && isfinite(nextgap) || return Inf

    log(nextgap / previousgap)
end

alkane_ladder_candidate_window(candidate::AlkaneLadderPathCandidate) = candidate.window

function alkane_ladder_path_success(
    best::AlkaneLadderPathRunBest,
    candidatesbystep::AbstractDict{Int, <:AbstractVector{<:AlkaneLadderPathCandidate}},
    candidateruns::AbstractVector{<:AbstractVector{Int}},
    runresults::AbstractVector{AlkaneLadderPathRunResult},
    settings::NamedTuple
)
    path = best.path
    laddersteps = [candidate.ladderstep for candidate in path]
    scanindices = [candidate.scanindex for candidate in path]
    stepgaps = diff(laddersteps)
    scangaps = diff(scanindices)

    AlkaneLadderPathInfo(
        true,
        :success,
        :success,
        nothing,
        "",
        laddersteps,
        scanindices,
        [candidate.score for candidate in path],
        [candidate.z for candidate in path],
        [candidate.centerz for candidate in path],
        [candidate.isolationz for candidate in path],
        [candidate.centervslowerz for candidate in path],
        [candidate.centervsupperz for candidate in path],
        [candidate.massspectrumcosine for candidate in path],
        [candidate.massspectrumdistance for candidate in path],
        [candidate.source for candidate in path],
        [candidate.misupported for candidate in path],
        [candidate.gapfilled for candidate in path],
        [candidate.expectedscan for candidate in path],
        [candidate.scanerror for candidate in path],
        [candidate.requiredcosine for candidate in path],
        scangaps,
        stepgaps,
        scangaps ./ stepgaps,
        sum(max.(stepgaps .- 1, 0)),
        best.objective,
        best.objective,
        [alkane_ladder_candidate_window(candidate) for candidate in path],
        path,
        candidatesbystep,
        candidatesbystep,
        candidateruns,
        runresults,
        settings
    )
end

function alkane_ladder_path_failure(
    reason::Symbol,
    message::AbstractString,
    candidatesbystep::AbstractDict{Int, <:AbstractVector{<:AlkaneLadderPathCandidate}},
    candidateruns::AbstractVector{<:AbstractVector{Int}},
    runresults::AbstractVector{AlkaneLadderPathRunResult},
    settings::NamedTuple
)
    AlkaneLadderPathInfo(
        false,
        :failed,
        reason,
        message,
        message,
        Int[],
        Int[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Symbol[],
        Bool[],
        Bool[],
        Float64[],
        Float64[],
        Float64[],
        Int[],
        Int[],
        Float64[],
        0,
        -Inf,
        -Inf,
        AlkaneMolecularIonContrast[],
        AlkaneLadderPathCandidate[],
        candidatesbystep,
        candidatesbystep,
        candidateruns,
        runresults,
        settings
    )
end
