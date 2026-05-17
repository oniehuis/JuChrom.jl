using CairoMakie
using JuChrom
using LinearAlgebra
using Random

include("parafac2_synthetic_data.jl")

CairoMakie.activate!()

function cosinesimilarity(x::AbstractVector, y::AbstractVector)
    denominator = norm(x) * norm(y)
    denominator == 0 && return 0.0
    dot(x, y) / denominator
end

function componentpermutations(n::Integer)
    n == 1 && return [[1]]
    permutations = Vector{Vector{Int}}()
    for firstcomponent in 1:n
        remaining = collect(1:n)
        deleteat!(remaining, firstcomponent)
        for tail in componentpermutations(n - 1)
            push!(permutations, vcat(firstcomponent, remaining[tail]))
        end
    end
    permutations
end

function scaletoreference(estimate::AbstractVector, reference::AbstractVector)
    denominator = dot(estimate, estimate)
    denominator == 0 && return 1.0
    dot(estimate, reference) / denominator
end

function normalizedforplot(values::AbstractVector)
    scale = maximum(abs, values)
    scale == 0 && return copy(values)
    values ./ scale
end

function spectrumsegments(mzvalues::AbstractVector, intensities::AbstractVector; offset=0.0)
    segments = Point2f[]
    sizehint!(segments, 2 * length(mzvalues))
    for index in eachindex(mzvalues, intensities)
        mz = Float32(mzvalues[index] + offset)
        intensity = Float32(intensities[index])
        push!(segments, Point2f(mz, 0))
        push!(segments, Point2f(mz, intensity))
    end
    segments
end

function matchparafac2components(data, fit)
    truespectra = data.spectra
    estimatedspectra = parafac2spectra(fit)
    ncomponents = size(truespectra, 2)
    scores = parafac2scores(fit)
    weights = parafac2abundances(fit)
    bestpermutation = collect(1:ncomponents)
    bestscore = -Inf

    for permutation in componentpermutations(ncomponents)
        score = sum(
            abs(cosinesimilarity(truespectra[:, component],
                estimatedspectra[:, permutation[component]]))
            for component in 1:ncomponents
        )
        if score > bestscore
            bestscore = score
            bestpermutation = permutation
        end
    end

    spectralsigns = [
        sign(cosinesimilarity(truespectra[:, component],
            estimatedspectra[:, bestpermutation[component]]))
        for component in 1:ncomponents
    ]
    replace!(spectralsigns, 0.0 => 1.0)

    scaledspectra = similar(truespectra)
    spectralcosines = Vector{Float64}(undef, ncomponents)
    scaledweights = similar(data.abundances)
    weightedprofiles = [
        Matrix{Float64}(undef, size(data.profiles[sampleindex]))
        for sampleindex in eachindex(data.profiles)
    ]

    for component in 1:ncomponents
        estimatedcomponent = bestpermutation[component]
        signfactor = spectralsigns[component]

        spectrumestimate = signfactor .* estimatedspectra[:, estimatedcomponent]
        spectralscale = scaletoreference(spectrumestimate, truespectra[:, component])
        scaledspectra[:, component] .= spectralscale .* spectrumestimate
        spectralcosines[component] = cosinesimilarity(
            truespectra[:, component],
            scaledspectra[:, component]
        )

        weightestimate = signfactor .* weights[:, estimatedcomponent]
        weightscale = scaletoreference(weightestimate, data.abundances[:, component])
        scaledweights[:, component] .= weightscale .* weightestimate

        estimatedprofilevectors = Vector{Float64}()
        trueprofilevectors = Vector{Float64}()
        for sampleindex in eachindex(scores)
            estimatedprofile = signfactor .* scores[sampleindex][:, estimatedcomponent] .*
                weights[sampleindex, estimatedcomponent]
            append!(estimatedprofilevectors, estimatedprofile)
            append!(trueprofilevectors, data.profiles[sampleindex][:, component])
        end

        profilescale = scaletoreference(estimatedprofilevectors, trueprofilevectors)
        for sampleindex in eachindex(scores)
            weightedprofiles[sampleindex][:, component] .= profilescale .* signfactor .*
                scores[sampleindex][:, estimatedcomponent] .* weights[sampleindex, estimatedcomponent]
        end
    end

    (
        permutation=bestpermutation,
        spectralcosines=spectralcosines,
        spectra=scaledspectra,
        abundances=scaledweights,
        profiles=weightedprofiles,
    )
end

function plotparafac2overlap(data, fit, matched; samples=(1, length(data.X)), figuresize=(1450, 930))
    colors = Makie.wong_colors()
    fig = Figure(; size=figuresize)
    Label(fig[0, 1:4],
        "PARAFAC2 overlapping synthetic GC/MS benchmark",
        fontsize=22,
        font=:bold,
    )

    for (plotindex, sampleindex) in enumerate(samples)
        retention = data.retentions[sampleindex]
        ax = Axis(fig[1, plotindex],
            title="sample $(sampleindex): synthetic data",
            xlabel="retention",
            ylabel="m/z",
        )
        heatmap!(ax, retention, data.mzvalues, data.X[sampleindex])
    end

    for (plotindex, sampleindex) in enumerate(samples)
        retention = data.retentions[sampleindex]
        ax = Axis(fig[2, plotindex],
            title="sample $(sampleindex): true vs recovered profiles",
            xlabel="retention",
            ylabel="relative intensity",
        )
        for component in axes(data.spectra, 2)
            trueprofile = normalizedforplot(data.profiles[sampleindex][:, component])
            recoveredprofile = normalizedforplot(matched.profiles[sampleindex][:, component])
            lines!(ax, retention, trueprofile;
                color=colors[component],
                linewidth=3,
                label="true C$(component)",
            )
            lines!(ax, retention, recoveredprofile;
                color=colors[component],
                linewidth=2,
                linestyle=:dash,
                label="PARAFAC2 C$(component)",
            )
        end
        axislegend(ax; position=:rt, framevisible=false)
    end

    for component in axes(data.spectra, 2)
        ax = Axis(fig[3, component],
            title="component $(component): mass spectrum",
            xlabel="m/z",
            ylabel="relative intensity",
        )
        mzstep = minimum(diff(data.mzvalues))
        offset = 0.12 * mzstep
        truespectrum = normalizedforplot(data.spectra[:, component])
        recoveredspectrum = normalizedforplot(matched.spectra[:, component])
        linesegments!(ax, spectrumsegments(data.mzvalues, truespectrum; offset=-offset);
            color=(:gray55, 0.7),
            linewidth=3,
            label="true",
        )
        linesegments!(ax, spectrumsegments(data.mzvalues, recoveredspectrum; offset=offset);
            color=(colors[component], 0.95),
            linewidth=1.4,
            label="PARAFAC2",
        )
        scatter!(ax, data.mzvalues .+ offset, recoveredspectrum;
            color=(colors[component], 0.95),
            markersize=3,
        )
        ylims!(ax, -0.03, 1.05)
        axislegend(ax; position=:rt, framevisible=false)
        text!(ax, minimum(data.mzvalues), 0.9;
            text="cos = $(round(matched.spectralcosines[component]; digits=3))",
            align=(:left, :center),
            fontsize=13,
        )
    end

    axabundance = Axis(fig[3, 3],
        title="abundances across samples",
        xlabel="sample",
        ylabel="scaled abundance",
    )
    sampleindices = collect(1:size(data.abundances, 1))
    for component in axes(data.spectra, 2)
        lines!(axabundance, sampleindices, data.abundances[:, component];
            color=colors[component],
            linewidth=3,
            label="true C$(component)",
        )
        scatter!(axabundance, sampleindices, matched.abundances[:, component];
            color=colors[component],
            marker=:xcross,
            markersize=14,
            label="PARAFAC2 C$(component)",
        )
    end
    axislegend(axabundance; position=:rt, framevisible=false)

    residuals = parafac2residuals(fit, data.X)
    sampleindex = first(samples)
    axresidual = Axis(fig[3, 4],
        title="sample $(sampleindex): residual",
        xlabel="retention",
        ylabel="m/z",
    )
    heatmap!(axresidual, data.retentions[sampleindex], data.mzvalues, residuals[sampleindex])

    axloss = Axis(fig[1, 3:4],
        title="ALS loss",
        xlabel="iteration",
        ylabel="RSS",
    )
    lines!(axloss, 0:(length(fit.loss) - 1), fit.loss; color=:black, linewidth=2)
    scatter!(axloss, 0:(length(fit.loss) - 1), fit.loss; color=:black, markersize=6)

    axfit = Axis(fig[2, 3:4],
        title="clean total ion chromatograms",
        xlabel="retention",
        ylabel="intensity",
    )
    for sampleindex in samples
        tic = vec(sum(data.Xclean[sampleindex]; dims=2))
        lines!(axfit, data.retentions[sampleindex], normalizedforplot(tic);
            linewidth=2,
            label="sample $(sampleindex)",
        )
    end
    axislegend(axfit; position=:rt, framevisible=false)

    fig
end

function parafac2overlapdemo(;
    nsamples::Integer=24,
    nretentions::Integer=301,
    mzvalues=defaultparafac2mzvalues(),
    spectra=nothing,
    abundances=nothing,
    retentionrange=(970.0, 1035.0),
    noiselevel::Real=0.01,
    apexes=(1000.0, 1006.0),
    widths=(7.2, 7.8),
    sampleshiftwidth::Real=2.2,
    componentshiftwidth::Real=0.45,
    retentionjitter::Real=0.0,
    nstarts::Integer=20,
    maxiters::Integer=200,
    tol::Real=1e-9,
    compression::Symbol=:cholesky,
    datarng::AbstractRNG=Random.MersenneTwister(2026),
    fitrng::AbstractRNG=Random.MersenneTwister(2027),
    samples=(1, nsamples),
    output=nothing,
)
    data = simulateparafac2overlap(;
        nsamples=nsamples,
        nretentions=nretentions,
        mzvalues=mzvalues,
        spectra=isnothing(spectra) ? defaultparafac2spectra(mzvalues) : spectra,
        abundances=isnothing(abundances) ? defaultparafac2abundances(nsamples) : abundances,
        retentionrange=retentionrange,
        noiselevel=noiselevel,
        apexes=apexes,
        widths=widths,
        sampleshiftwidth=sampleshiftwidth,
        componentshiftwidth=componentshiftwidth,
        retentionjitter=retentionjitter,
        rng=datarng,
    )
    fit = parafac2(
        data.X,
        2;
        retentions=data.retentions,
        mzvalues=data.mzvalues,
        samplelabels=data.samplelabels,
        nstarts=nstarts,
        maxiters=maxiters,
        tol=tol,
        compression=compression,
        rng=fitrng,
    )
    matched = matchparafac2components(data, fit)
    fig = plotparafac2overlap(data, fit, matched; samples=samples)

    if !isnothing(output)
        save(output, fig)
    end

    (
        data=data,
        fit=fit,
        matched=matched,
        figure=fig,
        fitpercent=parafac2fitpercent(fit, data.X),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = parafac2overlapdemo(; output=joinpath(@__DIR__, "parafac2_overlap_demo.png"))
    println("fit percent: ", round(100 * result.fitpercent; digits=2), "%")
    println("spectral cosines: ", round.(result.matched.spectralcosines; digits=3))
    println("saved: ", joinpath(@__DIR__, "parafac2_overlap_demo.png"))
end
