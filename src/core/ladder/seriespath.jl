"""
    alkaneseriespath(evidencetraces; kwargs...)

Find the best contiguous n-alkane series path through per-carbon evidence traces.

Candidate steps are local maxima in each evidence trace that pass absolute and relative
height thresholds. The dynamic-programming path may start and end at any available carbon
number, but all selected steps inside the path are consecutive carbon numbers and strictly
increase in scan index and retention.
"""
function alkaneseriespath(
    evidencetraces::AbstractDict;
    minsteps::Integer=5,
    stepreward::Real=0.05,
    minrelativepeakheight::Real=1e-4,
    minpeakheight::Real=1e-4,
    maxcandidatespertrace::Union{Nothing, Integer}=100,
    spacingweight::Real=1.0,
    gapincreaseweight::Real=0.25,
)
    validate_series_path_settings(
        minsteps,
        stepreward,
        minrelativepeakheight,
        minpeakheight,
        maxcandidatespertrace,
        spacingweight,
        gapincreaseweight,
    )

    settings = (
        minsteps=Int(minsteps),
        stepreward=Float64(stepreward),
        minrelativepeakheight=Float64(minrelativepeakheight),
        minpeakheight=Float64(minpeakheight),
        maxcandidatespertrace=maxcandidatespertrace,
        spacingweight=Float64(spacingweight),
        gapincreaseweight=Float64(gapincreaseweight),
    )
    candidatesbycarbon = series_path_candidates(
        evidencetraces;
        minrelativepeakheight=settings.minrelativepeakheight,
        minpeakheight=settings.minpeakheight,
        maxcandidatespertrace=maxcandidatespertrace,
    )
    candidateruns = consecutive_candidate_runs(candidatesbycarbon)
    isempty(candidateruns) && return failed_series_path_result(
        "no local maxima found in evidence traces",
        candidatesbycarbon,
        candidateruns,
        NamedTuple[],
        settings,
    )

    best = nothing
    runresults = NamedTuple[]
    for run in candidateruns
        run_best = best_series_path_in_run(
            run,
            candidatesbycarbon;
            minsteps=Int(minsteps),
            stepreward=settings.stepreward,
            spacingweight=settings.spacingweight,
            gapincreaseweight=settings.gapincreaseweight,
        )
        push!(
            runresults,
            (
                carbonrange=collect(Int.(run)),
                success=!isnothing(run_best),
                objective=isnothing(run_best) ? missing : run_best.objective,
                pathlength=isnothing(run_best) ? 0 : length(run_best.path),
            ),
        )
        if !isnothing(run_best) && (isnothing(best) || run_best.objective > best.objective)
            best = run_best
        end
    end

    isnothing(best) && return failed_series_path_result(
        "no valid path with at least $(Int(minsteps)) strictly increasing steps",
        candidatesbycarbon,
        candidateruns,
        runresults,
        settings,
    )

    successful_series_path_result(
        best,
        candidatesbycarbon,
        candidateruns,
        runresults,
        settings,
    )
end

struct SeriesPathState
    pos::Int
    previdx::Int
    curridx::Int
    objective::Float64
    pathlength::Int
    parent::Union{Nothing, Int}
end

function validate_series_path_settings(
    minsteps,
    stepreward,
    minrelativepeakheight,
    minpeakheight,
    maxcandidatespertrace,
    spacingweight,
    gapincreaseweight,
)
    minsteps >= 1 || throw(ArgumentError("minsteps must be at least 1"))
    isfinite(stepreward) || throw(ArgumentError("stepreward must be finite"))
    isfinite(minrelativepeakheight) && 0 <= minrelativepeakheight <= 1 || throw(
        ArgumentError("minrelativepeakheight must be finite and in [0, 1]"))
    isfinite(minpeakheight) && minpeakheight >= 0 || throw(
        ArgumentError("minpeakheight must be finite and nonnegative"))
    isnothing(maxcandidatespertrace) || maxcandidatespertrace >= 1 || throw(
        ArgumentError("maxcandidatespertrace must be nothing or a positive integer"))
    isfinite(spacingweight) && spacingweight >= 0 || throw(
        ArgumentError("spacingweight must be finite and nonnegative"))
    isfinite(gapincreaseweight) && gapincreaseweight >= 0 || throw(
        ArgumentError("gapincreaseweight must be finite and nonnegative"))

    nothing
end

function series_path_candidates(
    evidencetraces::AbstractDict;
    minrelativepeakheight::Float64,
    minpeakheight::Float64,
    maxcandidatespertrace::Union{Nothing, Integer},
)
    candidatesbycarbon = Dict{Int, Vector{NamedTuple}}()
    for carbon in sort(Int.(collect(keys(evidencetraces))))
        trace = evidencetraces[carbon]
        values = Float64.(rawintensities(trace))
        retentions = Float64.(rawretentions(trace))
        length(values) == length(retentions) || throw(DimensionMismatch(
            "evidence trace for C$(carbon) has different intensity and retention lengths"))
        all(isfinite, values) || throw(ArgumentError(
            "evidence trace for C$(carbon) contains nonfinite intensities"))
        all(isfinite, retentions) || throw(ArgumentError(
            "evidence trace for C$(carbon) contains nonfinite retentions"))

        tracemax = isempty(values) ? 0.0 : maximum(values)
        minheight = tracemax > 0 ?
            max(minrelativepeakheight * tracemax, minpeakheight) :
            minpeakheight
        candidates = [
            (
                carbon=carbon,
                scanindex=scanindex,
                retention=retentions[scanindex],
                score=values[scanindex],
            )
            for scanindex in localmaxima(values)
            if values[scanindex] >= minheight
        ]
        candidatesbycarbon[carbon] = limit_series_path_candidates(
            candidates,
            maxcandidatespertrace,
        )
    end

    candidatesbycarbon
end

function limit_series_path_candidates(
    candidates::AbstractVector{<:NamedTuple},
    maxcandidates::Union{Nothing, Integer},
)
    isnothing(maxcandidates) && return collect(candidates)
    length(candidates) <= maxcandidates && return collect(candidates)

    order = sortperm(candidates; by=candidate -> candidate.score, rev=true)
    limited = collect(candidates[order[1:maxcandidates]])
    sort!(limited; by=candidate -> (candidate.scanindex, -candidate.score))
    limited
end

function consecutive_candidate_runs(candidatesbycarbon::AbstractDict{Int})
    runs = Vector{Vector{Int}}()
    current = Int[]
    for carbon in sort(collect(keys(candidatesbycarbon)))
        if isempty(candidatesbycarbon[carbon])
            !isempty(current) && push!(runs, current)
            current = Int[]
            continue
        end

        if isempty(current) || carbon == last(current) + 1
            push!(current, carbon)
        else
            push!(runs, current)
            current = [carbon]
        end
    end
    !isempty(current) && push!(runs, current)

    runs
end

function best_series_path_in_run(
    carbonrun::AbstractVector{<:Integer},
    candidatesbycarbon::AbstractDict{Int};
    minsteps::Integer,
    stepreward::Float64,
    spacingweight::Float64,
    gapincreaseweight::Float64,
)
    candidates = [candidatesbycarbon[Int(carbon)] for carbon in carbonrun]
    length(candidates) < minsteps && return nothing

    bestobjective = -Inf
    beststateid = nothing
    states = SeriesPathState[]
    previousstates = Dict{Tuple{Int, Int, Int}, Int}()

    for pos in eachindex(candidates)
        currentstates = Dict{Tuple{Int, Int, Int}, Int}()

        for curridx in eachindex(candidates[pos])
            curr = candidates[pos][curridx]
            objective = curr.score + stepreward
            stateid = push_series_path_state!(
                states,
                currentstates,
                pos,
                0,
                curridx,
                objective,
                1,
                minsteps,
                nothing,
            )
            if minsteps <= 1 && states[stateid].objective > bestobjective
                bestobjective = states[stateid].objective
                beststateid = stateid
            end
        end

        for stateid in values(previousstates)
            state = states[stateid]
            previous = candidates[pos - 1][state.curridx]

            for curridx in eachindex(candidates[pos])
                curr = candidates[pos][curridx]
                series_step_is_forward(previous, curr) || continue

                penalty = 0.0
                if state.pathlength >= 2
                    prevprev = candidates[pos - 2][state.previdx]
                    previousgap = previous.retention - prevprev.retention
                    nextgap = curr.retention - previous.retention
                    gaplogratio = series_gap_log_ratio(previousgap, nextgap)
                    isfinite(gaplogratio) || continue
                    penalty += spacingweight * abs2(gaplogratio)
                    penalty += gapincreaseweight * abs2(max(0.0, gaplogratio))
                    isfinite(penalty) || continue
                end

                objective = state.objective + curr.score + stepreward - penalty
                newstateid = push_series_path_state!(
                    states,
                    currentstates,
                    pos,
                    state.curridx,
                    curridx,
                    objective,
                    state.pathlength + 1,
                    minsteps,
                    stateid,
                )
                newstate = states[newstateid]
                if newstate.pathlength >= minsteps && newstate.objective > bestobjective
                    bestobjective = newstate.objective
                    beststateid = newstateid
                end
            end
        end

        previousstates = currentstates
    end

    isnothing(beststateid) && return nothing
    (
        objective=bestobjective,
        path=reconstruct_series_path(states, beststateid, candidates),
    )
end

function series_step_is_forward(previouscandidate, nextcandidate)
    nextcandidate.scanindex > previouscandidate.scanindex &&
        nextcandidate.retention > previouscandidate.retention
end

function series_gap_log_ratio(previousgap::Real, nextgap::Real)
    previousgap > 0 && nextgap > 0 || return Inf
    isfinite(previousgap) && isfinite(nextgap) || return Inf
    log(nextgap / previousgap)
end

function push_series_path_state!(
    states::Vector{SeriesPathState},
    statemap::Dict{Tuple{Int, Int, Int}, Int},
    pos::Int,
    previdx::Int,
    curridx::Int,
    objective::Float64,
    pathlength::Int,
    lengthcap::Int,
    parent::Union{Nothing, Int},
)
    key = (previdx, curridx, min(pathlength, lengthcap))
    existing = get(statemap, key, nothing)
    if !isnothing(existing) && states[existing].objective >= objective
        return existing
    end

    push!(states, SeriesPathState(pos, previdx, curridx, objective, pathlength, parent))
    stateid = length(states)
    statemap[key] = stateid
    stateid
end

function reconstruct_series_path(
    states::Vector{SeriesPathState},
    stateid::Int,
    candidates::AbstractVector,
)
    state = states[stateid]
    if isnothing(state.parent)
        return [candidates[state.pos][state.curridx]]
    end

    path = reconstruct_series_path(states, state.parent, candidates)
    push!(path, candidates[state.pos][state.curridx])
    path
end

function successful_series_path_result(
    best,
    candidatesbycarbon,
    candidateruns,
    runresults,
    settings,
)
    path = best.path
    carbonnumbers = [candidate.carbon for candidate in path]
    scanindices = [candidate.scanindex for candidate in path]
    retentions = [candidate.retention for candidate in path]
    scores = [candidate.score for candidate in path]

    (
        success=true,
        failurereason=nothing,
        carbonnumbers=carbonnumbers,
        scanindices=scanindices,
        retentions=retentions,
        scores=scores,
        gaps=diff(retentions),
        objective=best.objective,
        path=path,
        candidatesbycarbon=candidatesbycarbon,
        candidateruns=candidateruns,
        runresults=runresults,
        settings=settings,
    )
end

function failed_series_path_result(
    reason,
    candidatesbycarbon,
    candidateruns,
    runresults,
    settings,
)
    (
        success=false,
        failurereason=reason,
        carbonnumbers=Int[],
        scanindices=Int[],
        retentions=Float64[],
        scores=Float64[],
        gaps=Float64[],
        objective=-Inf,
        path=NamedTuple[],
        candidatesbycarbon=candidatesbycarbon,
        candidateruns=candidateruns,
        runresults=runresults,
        settings=settings,
    )
end
