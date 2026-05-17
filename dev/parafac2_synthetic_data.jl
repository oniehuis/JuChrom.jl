using LinearAlgebra
using Random

function gaussianpeak(retentions, apex::Real, width::Real)
    width > 0 || throw(ArgumentError("width must be positive."))
    exp.(-0.5 .* ((retentions .- apex) ./ width) .^ 2)
end

function normalizecolumns(matrix::AbstractMatrix)
    out = float.(copy(matrix))
    for col in axes(out, 2)
        scale = maximum(abs, view(out, :, col))
        if scale > 0
            out[:, col] ./= scale
        end
    end
    out
end

function defaultparafac2abundances(nsamples::Integer)
    nsamples >= 2 || throw(ArgumentError("nsamples must be at least 2."))
    phase = range(0, 2π; length=nsamples)
    abundances = Matrix{Float64}(undef, nsamples, 2)
    abundances[:, 1] .= 1.0 .+ 0.45 .* sin.(phase .+ 0.25)
    abundances[:, 2] .= 1.0 .+ 0.55 .* cos.(phase .* 0.85 .+ 0.9)
    abundances
end

defaultparafac2mzvalues() = collect(50.0:249.0)

function spectralpeak(mzvalues, center::Real, width::Real)
    exp.(-0.5 .* ((mzvalues .- center) ./ width) .^ 2)
end

function defaultparafac2spectra(mzvalues=defaultparafac2mzvalues())
    mzvalues = Float64.(mzvalues)
    centers = [53.0, 57.0, 65.0, 73.0, 81.0, 91.0, 105.0, 119.0,
               133.0, 147.0, 161.0, 175.0, 189.0, 205.0, 221.0, 239.0]
    widths = [0.75, 0.90, 1.15, 0.85, 1.30, 0.95, 1.40, 1.10,
              1.55, 1.20, 1.65, 1.25, 1.45, 1.35, 1.50, 1.10]
    intensities1 = [1.00, 0.48, 0.22, 0.78, 0.16, 0.92, 0.35, 0.58,
                    0.24, 0.44, 0.18, 0.37, 0.29, 0.21, 0.15, 0.12]
    intensities2 = [0.42, 0.98, 0.55, 0.30, 0.76, 0.36, 0.64, 0.23,
                    0.50, 0.18, 0.40, 0.14, 0.34, 0.27, 0.20, 0.16]
    spectra = fill(0.025, length(mzvalues), 2)

    for peakindex in eachindex(centers)
        peak = spectralpeak(mzvalues, centers[peakindex], widths[peakindex])
        spectra[:, 1] .+= intensities1[peakindex] .* peak
        spectra[:, 2] .+= intensities2[peakindex] .* peak
    end

    normalizecolumns(spectra)
end

"""
    simulateparafac2overlap(; kwargs...)

Create a small synthetic GC/MS-like PARAFAC2 benchmark with two strongly overlapping
baseline-free peaks.

The returned named tuple contains vector-of-matrices data `X`, per-sample `retentions`,
shared `mzvalues`, true component `spectra`, true weighted retention `profiles`, sample
`abundances`, clean component matrices, and sample labels. The data matrix convention is
`retentions × m/z`.
"""
function simulateparafac2overlap(;
    nsamples::Integer=24,
    nretentions::Integer=301,
    mzvalues=defaultparafac2mzvalues(),
    spectra=nothing,
    abundances=defaultparafac2abundances(nsamples),
    retentionrange=(970.0, 1035.0),
    apexes=(1000.0, 1006.0),
    widths=(7.2, 7.8),
    sampleshiftwidth::Real=2.2,
    componentshiftwidth::Real=0.45,
    retentionjitter::Real=0.0,
    noiselevel::Real=0.01,
    rng::AbstractRNG=Random.MersenneTwister(2026)
)
    spectra = isnothing(spectra) ? defaultparafac2spectra(mzvalues) : spectra
    size(spectra, 2) == 2 || throw(ArgumentError("spectra must have two columns."))
    size(spectra, 1) == length(mzvalues) || throw(DimensionMismatch(
        "spectra must have one row per m/z value."
    ))
    size(abundances) == (nsamples, 2) || throw(DimensionMismatch(
        "abundances must have size (nsamples, 2)."
    ))
    nretentions >= 3 || throw(ArgumentError("nretentions must be at least 3."))
    first(retentionrange) < last(retentionrange) || throw(ArgumentError(
        "retentionrange must be increasing."
    ))
    noiselevel >= 0 || throw(ArgumentError("noiselevel must be nonnegative."))
    sampleshiftwidth >= 0 || throw(ArgumentError("sampleshiftwidth must be nonnegative."))
    componentshiftwidth >= 0 || throw(ArgumentError(
        "componentshiftwidth must be nonnegative."
    ))
    retentionjitter >= 0 || throw(ArgumentError("retentionjitter must be nonnegative."))

    spectra = normalizecolumns(spectra)
    basegrid = collect(range(first(retentionrange), last(retentionrange); length=nretentions))
    sampleshifts = sampleshiftwidth .* randn(rng, nsamples)
    componentshifts = componentshiftwidth .* randn(rng, nsamples, 2)
    retentions = Vector{Vector{Float64}}(undef, nsamples)
    profiles = Vector{Matrix{Float64}}(undef, nsamples)
    componentmatrices = Vector{Vector{Matrix{Float64}}}(undef, nsamples)
    Xclean = Vector{Matrix{Float64}}(undef, nsamples)

    for sampleindex in 1:nsamples
        jitter = retentionjitter == 0 ? zeros(nretentions) : retentionjitter .* randn(rng, nretentions)
        retention = sort(basegrid .+ sampleshifts[sampleindex] .+ jitter)
        retentions[sampleindex] = retention

        sampleprofiles = Matrix{Float64}(undef, nretentions, 2)
        samplecomponents = Vector{Matrix{Float64}}(undef, 2)
        for component in 1:2
            apex = apexes[component] + sampleshifts[sampleindex] +
                componentshifts[sampleindex, component]
            profile = gaussianpeak(retention, apex, widths[component])
            peakmax = maximum(profile)
            if peakmax > 0
                profile ./= peakmax
            end
            sampleprofiles[:, component] .= abundances[sampleindex, component] .* profile
            samplecomponents[component] =
                sampleprofiles[:, component] * transpose(spectra[:, component])
        end

        profiles[sampleindex] = sampleprofiles
        componentmatrices[sampleindex] = samplecomponents
        Xclean[sampleindex] = samplecomponents[1] .+ samplecomponents[2]
    end

    maxsignal = maximum(maximum, Xclean)
    noisefactor = noiselevel * maxsignal
    X = [Xclean[sampleindex] .+ noisefactor .* randn(rng, size(Xclean[sampleindex]))
         for sampleindex in 1:nsamples]
    samplelabels = ["sample $(sampleindex)" for sampleindex in 1:nsamples]

    (
        X=X,
        Xclean=Xclean,
        retentions=retentions,
        mzvalues=Float64.(mzvalues),
        spectra=spectra,
        profiles=profiles,
        abundances=Matrix{Float64}(abundances),
        componentmatrices=componentmatrices,
        sampleshifts=sampleshifts,
        componentshifts=componentshifts,
        samplelabels=samplelabels,
        parameters=(
            apexes=apexes,
            widths=widths,
            sampleshiftwidth=sampleshiftwidth,
            componentshiftwidth=componentshiftwidth,
            retentionjitter=retentionjitter,
            noiselevel=noiselevel,
        )
    )
end
