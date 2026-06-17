struct LadderWindowPathState
    pos::Int
    prevpos::Int
    previdx::Int
    curridx::Int
    objective::Float64
    pathlength::Int
    parent::Union{Nothing,Int}
end

function alkaneladderpath(
    molecularioninfo;
    kwargs...
)
    hasproperty(molecularioninfo, :contrasts) || throw(ArgumentError(
        "molecularioninfo must contain contrasts"))

    alkaneladderpath(molecularioninfo.contrasts; kwargs...)
end

function alkaneladderpath(
    contrasts::AbstractDict;
    centerzmin::Real=1.645,
    isolationzmin::Real=1.645,
    minsteps::Integer=5,
    stepreward::Real=0.05,
    maxcandidatesperstep::Union{Nothing,Integer}=100,
    spacingweight::Real=25.0,
    gapincreaseweight::Real=5.0,
    maxgapratio::Real=2.5,
    allowsinglemissingstep::Bool=true,
    maxmissingsteps::Integer=allowsinglemissingstep ? 1 : 0,
    missingsteppenalty::Real=2.0,
    msm=nothing,
    variances=nothing,
    standard=defaultalkanestandard(),
    massspectrummatch::Bool=true,
    massspectrummatchdistanceweight::Real=5.0
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
        standard;
        enabled=massspectrummatch
    )
    settings = (
        centerzmin=Float64(centerzmin),
        isolationzmin=Float64(isolationzmin),
        minsteps=Int(minsteps),
        stepreward=Float64(stepreward),
        maxcandidatesperstep=maxcandidatesperstep,
        spacingweight=Float64(spacingweight),
        gapincreaseweight=Float64(gapincreaseweight),
        maxgapratio=Float64(maxgapratio),
        allowsinglemissingstep=Bool(allowsinglemissingstep),
        maxmissingsteps=Int(maxmissingsteps),
        missingsteppenalty=Float64(missingsteppenalty),
        massspectrummatch=Bool(massspectrummatch),
        massspectrummatchactive=!isnothing(context),
        massspectrummatchdistanceweight=Float64(massspectrummatchdistanceweight)
    )

    candidatesbystep = alkane_ladder_path_candidates(
        contrasts;
        centerzmin=settings.centerzmin,
        isolationzmin=settings.isolationzmin,
        maxcandidatesperstep=maxcandidatesperstep,
        massspectrummatchcontext=context,
        massspectrummatchdistanceweight=settings.massspectrummatchdistanceweight
    )
    candidateruns = alkane_ladder_candidate_runs(
        candidatesbystep;
        maxmissingsteps=settings.maxmissingsteps
    )
    isempty(candidateruns) && return alkane_ladder_path_failure(
        :no_candidates,
        "no molecular-ion candidates pass centerzmin=$(settings.centerzmin) and isolationzmin=$(settings.isolationzmin)",
        candidatesbystep,
        candidateruns,
        NamedTuple[],
        settings
    )

    best = nothing
    runresults = NamedTuple[]
    for run in candidateruns
        runbest = alkane_best_ladder_path_in_run(
            run,
            candidatesbystep;
            minsteps=settings.minsteps,
            stepreward=settings.stepreward,
            spacingweight=settings.spacingweight,
            gapincreaseweight=settings.gapincreaseweight,
            maxgapratio=settings.maxgapratio,
            maxmissingsteps=settings.maxmissingsteps,
            missingsteppenalty=settings.missingsteppenalty
        )
        push!(
            runresults,
            (
                laddersteps=collect(Int.(run)),
                success=!isnothing(runbest),
                objective=isnothing(runbest) ? missing : runbest.objective,
                pathlength=isnothing(runbest) ? 0 : length(runbest.path)
            )
        )
        if !isnothing(runbest) && (isnothing(best) || runbest.objective > best.objective)
            best = runbest
        end
    end

    isnothing(best) && return alkane_ladder_path_failure(
        :no_valid_path,
        "no valid path with at least $(settings.minsteps) strictly increasing steps",
        candidatesbystep,
        candidateruns,
        runresults,
        settings
    )

    alkane_ladder_path_success(best, candidatesbystep, candidateruns, runresults, settings)
end

function validate_alkane_ladder_path_settings(
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
    massspectrummatchdistanceweight=5.0
)
    isfinite(centerzmin) && centerzmin >= 0 || throw(ArgumentError(
        "centerzmin must be finite and nonnegative"))
    isfinite(isolationzmin) && isolationzmin >= 0 || throw(ArgumentError(
        "isolationzmin must be finite and nonnegative"))
    minsteps isa Integer || throw(ArgumentError("minsteps must be an integer"))
    minsteps >= 1 || throw(ArgumentError("minsteps must be at least 1"))
    isfinite(stepreward) || throw(ArgumentError("stepreward must be finite"))
    isnothing(maxcandidatesperstep) ||
        maxcandidatesperstep isa Integer ||
        throw(ArgumentError("maxcandidatesperstep must be nothing or an integer"))
    isnothing(maxcandidatesperstep) ||
        maxcandidatesperstep >= 1 ||
        throw(ArgumentError("maxcandidatesperstep must be positive"))
    isfinite(spacingweight) && spacingweight >= 0 || throw(ArgumentError(
        "spacingweight must be finite and nonnegative"))
    isfinite(gapincreaseweight) && gapincreaseweight >= 0 || throw(ArgumentError(
        "gapincreaseweight must be finite and nonnegative"))
    isfinite(maxgapratio) && maxgapratio >= 1 || throw(ArgumentError(
        "maxgapratio must be finite and at least 1"))
    maxmissingsteps isa Integer || throw(ArgumentError(
        "maxmissingsteps must be an integer"))
    maxmissingsteps in (0, 1) || throw(ArgumentError(
        "maxmissingsteps must be 0 or 1"))
    isfinite(missingsteppenalty) && missingsteppenalty >= 0 || throw(ArgumentError(
        "missingsteppenalty must be finite and nonnegative"))
    isfinite(massspectrummatchdistanceweight) &&
        massspectrummatchdistanceweight >= 0 ||
        throw(ArgumentError(
            "mass-spectrum match distance weight must be finite and nonnegative"))

    nothing
end

function alkane_ladder_path_mass_spectrum_match_context(
    msm,
    variances,
    standard;
    enabled::Bool
)
    enabled || return nothing
    (isnothing(msm) || isnothing(variances) || isnothing(standard)) && return nothing

    alkane_ladder_mass_spectrum_match_context(
        msm,
        variances,
        standard;
        enabled=true,
        variancefloor=1.0
    )
end

function alkane_ladder_path_candidates(
    contrasts::AbstractDict;
    centerzmin::Real,
    isolationzmin::Real,
    maxcandidatesperstep::Union{Nothing,Integer},
    massspectrummatchcontext,
    massspectrummatchdistanceweight::Real
)
    candidatesbystep = Dict{Int, Vector{NamedTuple}}()
    for step in sort(Int.(collect(keys(contrasts))))
        candidates = NamedTuple[]
        for contrast in contrasts[step]
            alkane_molecular_ion_window_has_candidate_evidence(
                contrast;
                centerzmin=centerzmin,
                isolationzmin=isolationzmin
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
            push!(
                candidates,
                (
                    ladderstep=step,
                    scanindex=Int(contrast.apexindex),
                    score=candidate_score,
                    z=Float64(contrast.z),
                    centerz=Float64(contrast.centerz),
                    isolationz=Float64(contrast.isolationz),
                    centervslowerz=Float64(contrast.centervslowerz),
                    centervsupperz=Float64(contrast.centervsupperz),
                    molecularionscore=Float64(contrast.molecularionscore),
                    massspectrumcosine=match.cosine,
                    massspectrumdistance=match.distance,
                    massspectrumioncount=match.ioncount,
                    window=contrast
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
    contrast;
    centerzmin::Real,
    isolationzmin::Real
)
    hasproperty(contrast, :centerz) &&
        isfinite(contrast.centerz) &&
        contrast.centerz >= centerzmin ||
        return false
    hasproperty(contrast, :isolationz) &&
        isfinite(contrast.isolationz) &&
        contrast.isolationz >= isolationzmin ||
        return false
    hasproperty(contrast, :molecularionscore) &&
        isfinite(contrast.molecularionscore) &&
        contrast.molecularionscore > 0 ||
        return false

    true
end

function alkane_ladder_path_candidate_dp_score(
    score::Real,
    distance::Real,
    distanceweight::Real
)
    isfinite(score) && score > 0 || return 0.0
    isfinite(distance) || return score

    score * exp(-Float64(distanceweight) * clamp(Float64(distance), 0.0, 1.0))
end

function alkane_limit_ladder_path_candidates(
    candidates::AbstractVector,
    maxcandidates::Union{Nothing,Integer}
)
    sorted = sort(collect(candidates); by=candidate -> candidate.scanindex)
    isnothing(maxcandidates) && return sorted
    length(sorted) <= maxcandidates && return sorted

    order = sortperm(sorted; by=candidate -> candidate.score, rev=true)
    limited = collect(sorted[order[1:maxcandidates]])
    sort!(limited; by=candidate -> candidate.scanindex)

    limited
end

function alkane_ladder_candidate_runs(
    candidatesbystep::AbstractDict{Int};
    maxmissingsteps::Int
)
    maxmissingsteps in (0, 1) || throw(ArgumentError(
        "maxmissingsteps must be 0 or 1"))

    runs = Vector{Vector{Int}}()
    current = Int[]
    for step in sort(collect(keys(candidatesbystep)))
        isempty(candidatesbystep[step]) && continue
        if isempty(current) || step - last(current) <= maxmissingsteps + 1
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
    candidatesbystep::AbstractDict{Int};
    minsteps::Integer,
    stepreward::Float64,
    spacingweight::Float64,
    gapincreaseweight::Float64,
    maxgapratio::Float64,
    maxmissingsteps::Int,
    missingsteppenalty::Float64
)
    maxmissingsteps in (0, 1) || throw(ArgumentError(
        "maxmissingsteps must be 0 or 1"))

    candidates = [candidatesbystep[Int(step)] for step in steprun]
    length(candidates) < minsteps && return nothing

    bestobjective = -Inf
    beststateid = nothing
    states = LadderWindowPathState[]
    statesbypos = [Dict{Tuple{Int,Int,Int,Int}, Int}() for _ in candidates]

    for pos in eachindex(candidates)
        currentstates = statesbypos[pos]
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
            if minsteps <= 1 && states[stateid].objective > bestobjective
                bestobjective = states[stateid].objective
                beststateid = stateid
            end
        end

        for prevpos in 1:(pos - 1)
            stepgap = Int(steprun[pos]) - Int(steprun[prevpos])
            1 <= stepgap <= maxmissingsteps + 1 || continue
            missingsteps = stepgap - 1

            for stateid in values(statesbypos[prevpos])
                state = states[stateid]
                previous = candidates[prevpos][state.curridx]
                for curridx in eachindex(candidates[pos])
                    curr = candidates[pos][curridx]
                    curr.scanindex > previous.scanindex || continue

                    penalty = missingsteppenalty * missingsteps
                    if state.pathlength >= 2
                        parentid = state.parent
                        isnothing(parentid) && continue
                        parent = states[parentid]
                        prevprev = candidates[parent.pos][parent.curridx]
                        previousstepgap =
                            Int(steprun[state.pos]) - Int(steprun[parent.pos])
                        nextstepgap =
                            Int(steprun[pos]) - Int(steprun[state.pos])
                        previousgap =
                            (previous.scanindex - prevprev.scanindex) /
                            previousstepgap
                        nextgap =
                            (curr.scanindex - previous.scanindex) / nextstepgap
                        gaplogratio = alkane_ladder_gap_log_ratio(previousgap, nextgap)
                        isfinite(gaplogratio) || continue
                        abs(gaplogratio) <= log(maxgapratio) || continue
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
                    if newstate.pathlength >= minsteps &&
                            newstate.objective > bestobjective
                        bestobjective = newstate.objective
                        beststateid = newstateid
                    end
                end
            end
        end
    end

    isnothing(beststateid) && return nothing

    (
        objective=bestobjective,
        path=alkane_reconstruct_ladder_path(states, beststateid, candidates)
    )
end

function alkane_push_ladder_path_state!(
    states::Vector{LadderWindowPathState},
    statemap::Dict{Tuple{Int,Int,Int,Int}, Int},
    pos::Int,
    prevpos::Int,
    previdx::Int,
    curridx::Int,
    objective::Float64,
    pathlength::Int,
    lengthcap::Int,
    parent::Union{Nothing,Int}
)
    key = (prevpos, previdx, curridx, min(pathlength, lengthcap))
    existing = get(statemap, key, nothing)
    if !isnothing(existing) && states[existing].objective >= objective
        return existing
    end

    push!(
        states,
        LadderWindowPathState(
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
    candidates::AbstractVector
)
    state = states[stateid]
    if isnothing(state.parent)
        return NamedTuple[candidates[state.pos][state.curridx]]
    end
    path = alkane_reconstruct_ladder_path(states, state.parent, candidates)
    push!(path, candidates[state.pos][state.curridx])

    path
end

function alkane_ladder_gap_log_ratio(previousgap::Real, nextgap::Real)
    previousgap > 0 && nextgap > 0 || return Inf
    isfinite(previousgap) && isfinite(nextgap) || return Inf

    log(nextgap / previousgap)
end

function alkane_ladder_candidate_float_field(candidate, field::Symbol)
    hasproperty(candidate, field) && return Float64(getproperty(candidate, field))
    if hasproperty(candidate, :window) &&
            !isnothing(candidate.window) &&
            hasproperty(candidate.window, field)
        return Float64(getproperty(candidate.window, field))
    end

    NaN
end

function alkane_ladder_path_success(
    best,
    candidatesbystep,
    candidateruns,
    runresults,
    settings
)
    path = best.path
    laddersteps = [candidate.ladderstep for candidate in path]
    scanindices = [candidate.scanindex for candidate in path]
    stepgaps = diff(laddersteps)
    scangaps = diff(scanindices)

    (
        success=true,
        status=:success,
        reason=:success,
        failurereason=nothing,
        message="",
        laddersteps=laddersteps,
        scanindices=scanindices,
        scores=[candidate.score for candidate in path],
        zscores=[candidate.z for candidate in path],
        centerzs=[
            alkane_ladder_candidate_float_field(candidate, :centerz)
            for candidate in path
        ],
        isolationzs=[
            alkane_ladder_candidate_float_field(candidate, :isolationz)
            for candidate in path
        ],
        centervslowerzs=[
            alkane_ladder_candidate_float_field(candidate, :centervslowerz)
            for candidate in path
        ],
        centervsupperzs=[
            alkane_ladder_candidate_float_field(candidate, :centervsupperz)
            for candidate in path
        ],
        massspectrumcosines=[
            alkane_ladder_candidate_float_field(candidate, :massspectrumcosine)
            for candidate in path
        ],
        massspectrumdistances=[
            alkane_ladder_candidate_float_field(candidate, :massspectrumdistance)
            for candidate in path
        ],
        sources=[
            hasproperty(candidate, :source) ? candidate.source : :molecular_ion_dp
            for candidate in path
        ],
        misupported=[
            hasproperty(candidate, :misupported) ? Bool(candidate.misupported) : true
            for candidate in path
        ],
        gapfilled=[
            hasproperty(candidate, :gapfilled) ? Bool(candidate.gapfilled) : false
            for candidate in path
        ],
        expectedscanindices=[
            alkane_ladder_candidate_float_field(candidate, :expectedscan)
            for candidate in path
        ],
        scanerrors=[
            alkane_ladder_candidate_float_field(candidate, :scanerror)
            for candidate in path
        ],
        requiredcosines=[
            alkane_ladder_candidate_float_field(candidate, :requiredcosine)
            for candidate in path
        ],
        gaps=scangaps,
        stepgaps=stepgaps,
        scanstepgaps=scangaps ./ stepgaps,
        missingsteps=sum(max.(stepgaps .- 1, 0)),
        objective=best.objective,
        score=best.objective,
        windows=[candidate.window for candidate in path],
        path=path,
        candidates=candidatesbystep,
        candidatesbystep=candidatesbystep,
        candidateruns=candidateruns,
        runresults=runresults,
        settings=settings
    )
end

function alkane_ladder_path_failure(
    reason::Symbol,
    message,
    candidatesbystep,
    candidateruns,
    runresults,
    settings
)
    (
        success=false,
        status=:failed,
        reason=reason,
        failurereason=message,
        message=message,
        laddersteps=Int[],
        scanindices=Int[],
        scores=Float64[],
        zscores=Float64[],
        centerzs=Float64[],
        isolationzs=Float64[],
        centervslowerzs=Float64[],
        centervsupperzs=Float64[],
        massspectrumcosines=Float64[],
        massspectrumdistances=Float64[],
        sources=Symbol[],
        misupported=Bool[],
        gapfilled=Bool[],
        expectedscanindices=Float64[],
        scanerrors=Float64[],
        requiredcosines=Float64[],
        gaps=Int[],
        stepgaps=Int[],
        scanstepgaps=Float64[],
        missingsteps=0,
        objective=-Inf,
        score=-Inf,
        windows=NamedTuple[],
        path=NamedTuple[],
        candidates=candidatesbystep,
        candidatesbystep=candidatesbystep,
        candidateruns=candidateruns,
        runresults=runresults,
        settings=settings
    )
end
