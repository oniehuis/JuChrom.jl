const Parafac2StartDiagnostic{T<:Real} = NamedTuple{
    (:start, :iterations, :converged, :stopreason, :loss),
    Tuple{Int, Int, Bool, Symbol, T}
}

"""
    Parafac2Fit

Container for PARAFAC2 decomposition results.

The PARAFAC2 model represented by this type is

    X[k] ≈ bases[k] * core * Diagonal(weights[k, :]) * loadings'

where each `X[k]` is a sample matrix with retentions on the first axis, m/z channels on
the second axis, and `ncomponents` latent components. The retention axis can contain
Kovats indices or any other retention coordinate used upstream.

# Fields
- `ncomponents::Int`: number of PARAFAC2 components.
- `retentioncounts::Vector{Int}`: number of retention positions in each sample matrix.
- `mzcount::Int`: common number of m/z channels in all input sample matrices.
- `retentions`: retention coordinates per sample, stored without units as one vector per
  sample, or `nothing`.
- `retentionunit`: retention coordinate unit, or `nothing` when unitless or absent.
- `mzvalues`: m/z values, stored without units, or `nothing`.
- `mzunit`: m/z unit, or `nothing` when unitless or absent.
- `samplelabels`: sample labels, or `nothing`.
- `loadings`: shared m/z-mode loadings, or `nothing` before fitting.
- `core`: shared PARAFAC2 score core, or `nothing` before fitting.
- `weights`: sample/component weights, or `nothing` before fitting.
- `bases`: sample-specific orthonormal bases, or `nothing` before fitting.
- `loss::Vector`: objective values recorded during fitting.
- `converged::Bool`: whether fitting met its convergence criterion.
- `stopreason::Symbol`: why fitting stopped (`:tol`, `:maxiters`,
  `:objective_increase`, or `:nonfinite`).
- `iterations::Int`: number of fitting iterations performed by the selected start.
- `nstarts::Int`: number of starts attempted during fitting.
- `beststart::Int`: one-based index of the selected start.
- `startdiagnostics`: per-start iteration, convergence, stop reason, and loss summary.
- `compression::Symbol`: requested cross-product compression method.
- `compressed::Vector{Bool}`: whether each sample matrix was actually compressed.
- `nonnegative::Tuple`: nonnegative constraints used during fitting.

At the current implementation stage, [`parafac2`](@ref) fits the direct PARAFAC2
alternating least-squares algorithm with optional nonnegative spectra and intensity
constraints.

See also [`parafac2`](@ref).
"""
struct Parafac2Fit{T<:Real}
    ncomponents::Int
    retentioncounts::Vector{Int}
    mzcount::Int
    retentions::Union{Nothing, Vector{Vector{Float64}}}
    retentionunit::Union{Nothing, Unitful.Units}
    mzvalues::Union{Nothing, Vector{Float64}}
    mzunit::Union{Nothing, Unitful.Units}
    samplelabels::Union{Nothing, Vector{String}}
    loadings::Union{Nothing, Matrix{T}}
    core::Union{Nothing, Matrix{T}}
    weights::Union{Nothing, Matrix{T}}
    bases::Union{Nothing, Vector{Matrix{T}}}
    loss::Vector{T}
    converged::Bool
    stopreason::Symbol
    iterations::Int
    nstarts::Int
    beststart::Int
    startdiagnostics::Vector{Parafac2StartDiagnostic{T}}
    compression::Symbol
    compressed::Vector{Bool}
    nonnegative::Tuple{Vararg{Symbol}}
end

Base.Broadcast.broadcastable(fit::Parafac2Fit) = Base.RefValue(fit)

parafac2totaliterations(fit::Parafac2Fit) =
    isempty(fit.startdiagnostics) ? 0 :
    sum(diagnostic.iterations for diagnostic in fit.startdiagnostics)

parafac2ncompressed(fit::Parafac2Fit) = count(fit.compressed)

function parafac2compressionlabel(fit::Parafac2Fit)
    if fit.compression ≡ :cholesky
        "$(fit.compression) ($(parafac2ncompressed(fit))/$(length(fit.compressed)) slices compressed)"
    else
        string(fit.compression)
    end
end

function parafac2startdiagnosticslabel(fit::Parafac2Fit)
    isempty(fit.startdiagnostics) && return "not stored"
    join(
        (
            "$(diagnostic.start):$(diagnostic.iterations)/$(diagnostic.stopreason)"
            for diagnostic in fit.startdiagnostics
        ),
        ", "
    )
end

function Base.summary(fit::Parafac2Fit)
    "Parafac2Fit (samples=$(length(fit.retentioncounts)), mz=$(fit.mzcount), " *
    "components=$(fit.ncomponents), iterations=$(fit.iterations), " *
    "totaliterations=$(parafac2totaliterations(fit)), " *
    "starts=$(fit.nstarts), beststart=$(fit.beststart), " *
    "compression=$(parafac2compressionlabel(fit)), nonnegative=$(fit.nonnegative), " *
    "converged=$(fit.converged), stopreason=$(fit.stopreason))"
end

Base.show(io::IO, fit::Parafac2Fit) = print(io, summary(fit))

function Base.show(io::IO, ::MIME"text/plain", fit::Parafac2Fit)
    println(io, summary(fit))
    println(io, "  retention counts : ", fit.retentioncounts)
    println(io, "  retentions       : ", isnothing(fit.retentions) ? "not stored" : "stored")
    println(io, "  m/z values       : ", isnothing(fit.mzvalues) ? "not stored" : "stored")
    println(io, "  sample labels    : ", isnothing(fit.samplelabels) ? "not stored" : "stored")
    println(io, "  starts           : ", fit.nstarts, " (best: ", fit.beststart,
        ", total iterations: ", parafac2totaliterations(fit), ")")
    println(io, "  start diagnostics: ", parafac2startdiagnosticslabel(fit))
    println(io, "  compression      : ", parafac2compressionlabel(fit))
    println(io, "  nonnegative      : ", fit.nonnegative)
    println(io, "  stop reason      : ", fit.stopreason)
    println(io, "  factor state     : ", isnothing(fit.loadings) ? "unfitted" : "fitted")
end

function parafac2inputmatrices(
    X::AbstractVector{<:AbstractMatrix{<:Real}},
    ::Type{T}
) where {T<:AbstractFloat}
    [Matrix{T}(Xk) for Xk in X]
end

function parafac2inputscale(
    X::AbstractVector{<:AbstractMatrix{T}}
) where {T<:AbstractFloat}
    scale = zero(T)
    for Xk in X
        scale = hypot(scale, norm(Xk))
    end

    isfinite(scale) && scale > zero(T) ? scale : one(T)
end

function parafac2squarednorms(
    X::AbstractVector{<:AbstractMatrix{T}}
) where {T<:AbstractFloat}
    norms = Vector{T}(undef, length(X))
    for sampleindex in eachindex(X)
        norms[sampleindex] = sum(abs2, X[sampleindex])
    end

    norms
end

const PARAFAC2_THREAD_MINITEMS = 8

parafac2threaded(nitems::Integer) =
    Base.Threads.nthreads() > 1 && nitems >= PARAFAC2_THREAD_MINITEMS

function parafac2compression(compression::Symbol)
    compression in (:none, :cholesky) || throw(ArgumentError(
        "compression must be :none or :cholesky."
    ))
    compression
end

function parafac2nonnegative(nonnegative)
    requested = if nonnegative ≡ true
        (:spectra, :intensities)
    elseif nonnegative ≡ false || isnothing(nonnegative)
        ()
    elseif nonnegative isa Symbol
        (nonnegative,)
    elseif nonnegative isa Tuple || nonnegative isa AbstractVector
        Tuple(nonnegative)
    else
        throw(ArgumentError(
            "nonnegative must be a Bool, Symbol, tuple, or vector of symbols."
        ))
    end

    canonical = Symbol[]
    for item in requested
        item isa Symbol || throw(ArgumentError("nonnegative entries must be symbols."))
        constraint = item
        constraint in (:spectra, :intensities) || throw(ArgumentError(
            "nonnegative entries must be :spectra or :intensities."
        ))
        constraint in canonical || push!(canonical, constraint)
    end

    ordered = Symbol[]
    :spectra in canonical && push!(ordered, :spectra)
    :intensities in canonical && push!(ordered, :intensities)
    Tuple(ordered)
end

parafac2constrainspectra(nonnegative::Tuple{Vararg{Symbol}}) = :spectra in nonnegative
parafac2constrainintensities(nonnegative::Tuple{Vararg{Symbol}}) =
    :intensities in nonnegative

function parafac2compressionmatrices(
    X::AbstractVector{<:AbstractMatrix{T}},
    compression::Symbol
) where {T<:AbstractFloat}
    method = parafac2compression(compression)
    compressed = falses(length(X))
    method ≡ :none && return X, compressed

    matrices = Vector{Matrix{T}}(undef, length(X))
    for sampleindex in eachindex(X)
        Xk = X[sampleindex]
        if size(Xk, 1) > size(Xk, 2)
            crossproduct = transpose(Xk) * Xk
            factorization = cholesky(Symmetric(crossproduct); check=false)
            if issuccess(factorization)
                Hk = Matrix{T}(factorization.U)
                all(isfinite, Hk) || throw(ArgumentError(
                    "Cholesky compression produced non-finite values for X[$(sampleindex)]."
                ))
                matrices[sampleindex] = Hk
                compressed[sampleindex] = true
            else
                matrices[sampleindex] = Xk
            end
        else
            matrices[sampleindex] = Xk
        end
    end

    matrices, compressed
end

function parafac2initialloadings(
    X::AbstractVector{<:AbstractMatrix{T}},
    ncomponents::Integer
) where {T<:AbstractFloat}
    mzcount = size(first(X), 2)
    crossproduct = zeros(T, mzcount, mzcount)
    for Xk in X
        crossproduct .+= transpose(Xk) * Xk
    end

    decomposition = eigen(Symmetric(crossproduct))
    order = sortperm(decomposition.values; rev=true)[1:ncomponents]
    Matrix{T}(decomposition.vectors[:, order])
end

function parafac2basis(
    X::AbstractMatrix{T},
    coefficient::AbstractMatrix{T}
) where {T<:AbstractFloat}
    ncomponents = size(coefficient, 1)
    decomposition = svd(coefficient * transpose(X); full=false)
    Matrix{T}(
        decomposition.V[:, 1:ncomponents] *
        transpose(decomposition.U[:, 1:ncomponents])
    )
end

function parafac2basis!(
    basis::AbstractMatrix{T},
    X::AbstractMatrix{T},
    coefficient::AbstractMatrix{T}
) where {T<:AbstractFloat}
    ncomponents = size(coefficient, 1)
    decomposition = svd(coefficient * transpose(X); full=false)
    @views mul!(
        basis,
        decomposition.V[:, 1:ncomponents],
        transpose(decomposition.U[:, 1:ncomponents])
    )
    basis
end

function parafac2basisworkspace(
    X::AbstractVector{<:AbstractMatrix{T}},
    ncomponents::Integer
) where {T<:AbstractFloat}
    [Matrix{T}(undef, size(Xk, 1), ncomponents) for Xk in X]
end

function parafac2updatebases!(
    bases::AbstractVector{<:AbstractMatrix{T}},
    X::AbstractVector{<:AbstractMatrix{T}},
    loadings::AbstractMatrix{T},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T}
) where {T<:AbstractFloat}
    length(bases) == length(X) || throw(DimensionMismatch(
        "bases and X must have the same length."
    ))

    if parafac2threaded(length(X))
        @threads :static for sampleindex in eachindex(X)
            coefficient = core * Diagonal(view(weights, sampleindex, :)) * transpose(loadings)
            parafac2basis!(bases[sampleindex], X[sampleindex], coefficient)
        end
    else
        for sampleindex in eachindex(X)
            coefficient = core * Diagonal(view(weights, sampleindex, :)) * transpose(loadings)
            parafac2basis!(bases[sampleindex], X[sampleindex], coefficient)
        end
    end

    bases
end

function parafac2updatebases(
    X::AbstractVector{<:AbstractMatrix{T}},
    loadings::AbstractMatrix{T},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T}
) where {T<:AbstractFloat}
    parafac2updatebases!(
        parafac2basisworkspace(X, size(core, 1)),
        X,
        loadings,
        core,
        weights
    )
end

function parafac2transformedworkspace(
    X::AbstractVector{<:AbstractMatrix{T}},
    bases::AbstractVector{<:AbstractMatrix{T}}
) where {T<:AbstractFloat}
    length(bases) == length(X) || throw(DimensionMismatch(
        "bases and X must have the same length."
    ))
    [Matrix{T}(undef, size(bases[sampleindex], 2), size(X[sampleindex], 2))
        for sampleindex in eachindex(X)]
end

function parafac2transformeddata!(
    transformed::AbstractVector{<:AbstractMatrix{T}},
    X::AbstractVector{<:AbstractMatrix{T}},
    bases::AbstractVector{<:AbstractMatrix{T}}
) where {T<:AbstractFloat}
    length(transformed) == length(X) == length(bases) || throw(DimensionMismatch(
        "transformed data, X, and bases must have the same length."
    ))

    if parafac2threaded(length(X))
        @threads :static for sampleindex in eachindex(X)
            mul!(transformed[sampleindex], transpose(bases[sampleindex]), X[sampleindex])
        end
    else
        for sampleindex in eachindex(X)
            mul!(transformed[sampleindex], transpose(bases[sampleindex]), X[sampleindex])
        end
    end

    transformed
end

function parafac2transformeddata(
    X::AbstractVector{<:AbstractMatrix{T}},
    bases::AbstractVector{<:AbstractMatrix{T}}
) where {T<:AbstractFloat}
    parafac2transformeddata!(
        parafac2transformedworkspace(X, bases),
        X,
        bases
    )
end

function khatrirao(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(A, 2) == size(B, 2) || throw(DimensionMismatch(
        "A and B must have the same number of columns."
    ))

    rows = size(A, 1) * size(B, 1)
    cols = size(A, 2)
    out = Matrix{T}(undef, rows, cols)

    for (outcol, acol, bcol) in zip(axes(out, 2), axes(A, 2), axes(B, 2))
        row = firstindex(out, 1)
        for a in axes(A, 1), b in axes(B, 1)
            out[row, outcol] = A[a, acol] * B[b, bcol]
            row += 1
        end
    end

    out
end

function nonnegativeleastsquares(
    design::AbstractMatrix{T},
    response::AbstractVector{T}
) where {T<:AbstractFloat}
    size(design, 1) == length(response) || throw(DimensionMismatch(
        "design rows must match response length."
    ))

    ncoef = size(design, 2)
    x = zeros(T, ncoef)
    z = zeros(T, ncoef)
    passive = falses(ncoef)
    tolerance = sqrt(eps(T)) * max(norm(design), one(T)) * max(norm(response), one(T))
    maxiters = max(30, 30 * ncoef * ncoef)
    iteration = 0
    gradient = transpose(design) * (response - design * x)

    while iteration < maxiters
        candidate = 0
        candidategradient = T(tolerance)
        for index in axes(passive, 1)
            if !passive[index] && gradient[index] > candidategradient
                candidate = index
                candidategradient = gradient[index]
            end
        end
        candidate == 0 && break

        passive[candidate] = true
        while true
            active = findall(passive)
            z .= zero(T)
            if !isempty(active)
                z[active] .= design[:, active] \ response
            end

            all(index -> z[index] > T(tolerance), active) && break

            alpha = one(T)
            for index in active
                if z[index] ≤ T(tolerance)
                    denominator = x[index] - z[index]
                    if denominator > eps(T)
                        alpha = min(alpha, x[index] / denominator)
                    else
                        alpha = zero(T)
                    end
                end
            end

            x .+= alpha .* (z .- x)
            for index in active
                if x[index] ≤ T(tolerance)
                    passive[index] = false
                    x[index] = zero(T)
                end
            end

            iteration += 1
            iteration ≥ maxiters && break
        end

        x .= max.(z, zero(T))
        gradient = transpose(design) * (response - design * x)
        iteration += 1
    end

    max.(x, zero(T))
end

function nonnegativeleastsquaresgram!(
    x::AbstractVector{T},
    gram::AbstractMatrix{T},
    cross::AbstractVector{T},
    tolerance::T
) where {T<:AbstractFloat}
    ncoef = length(x)
    z = zeros(T, ncoef)
    passive = falses(ncoef)
    gradient = Vector{T}(cross)
    maxiters = max(30, 30 * ncoef * ncoef)
    iteration = 0
    fill!(x, zero(T))

    while iteration < maxiters
        candidate = 0
        candidategradient = tolerance
        for index in axes(passive, 1)
            if !passive[index] && gradient[index] > candidategradient
                candidate = index
                candidategradient = gradient[index]
            end
        end
        candidate == 0 && return true

        passive[candidate] = true
        while true
            active = findall(passive)
            z .= zero(T)
            if !isempty(active)
                active_solution = try
                    gram[active, active] \ cross[active]
                catch err
                    err isa InterruptException && rethrow()
                    return false
                end
                all(isfinite, active_solution) || return false
                z[active] .= active_solution
            end

            all(index -> z[index] > tolerance, active) && break

            alpha = one(T)
            for index in active
                if z[index] ≤ tolerance
                    denominator = x[index] - z[index]
                    if denominator > eps(T)
                        alpha = min(alpha, x[index] / denominator)
                    else
                        alpha = zero(T)
                    end
                end
            end

            x .+= alpha .* (z .- x)
            for index in active
                if x[index] ≤ tolerance
                    passive[index] = false
                    x[index] = zero(T)
                end
            end

            iteration += 1
            iteration ≥ maxiters && break
        end

        x .= max.(z, zero(T))
        gradient .= cross
        mul!(gradient, gram, x, -one(T), one(T))
        iteration += 1
    end

    true
end

function nonnegativeleastsquares_batch_rhs!(
    coefficients::AbstractMatrix{T},
    design::AbstractMatrix{T},
    responses::AbstractMatrix{T},
    gram::AbstractMatrix{T},
    cross::AbstractMatrix{T},
    designnorm::T,
    rhsindex::Integer
) where {T<:AbstractFloat}
    response = view(responses, :, rhsindex)
    tolerance = sqrt(eps(T)) * max(designnorm, one(T)) * max(norm(response), one(T))
    coefficient = view(coefficients, :, rhsindex)
    success = nonnegativeleastsquaresgram!(
        coefficient,
        gram,
        view(cross, :, rhsindex),
        T(tolerance)
    )
    if !success
        coefficient .= nonnegativeleastsquares(design, response)
    end

    coefficients
end

function nonnegativeleastsquares_batch(
    design::AbstractMatrix{T},
    responses::AbstractMatrix{T}
) where {T<:AbstractFloat}
    size(design, 1) == size(responses, 1) || throw(DimensionMismatch(
        "design rows must match response rows."
    ))

    ncoef = size(design, 2)
    nrhs = size(responses, 2)
    coefficients = Matrix{T}(undef, ncoef, nrhs)
    gram = transpose(design) * design
    cross = transpose(design) * responses
    designnorm = norm(design)

    if parafac2threaded(nrhs)
        @threads :static for rhsindex in axes(responses, 2)
            nonnegativeleastsquares_batch_rhs!(
                coefficients,
                design,
                responses,
                gram,
                cross,
                designnorm,
                rhsindex
            )
        end
    else
        for rhsindex in axes(responses, 2)
            nonnegativeleastsquares_batch_rhs!(
                coefficients,
                design,
                responses,
                gram,
                cross,
                designnorm,
                rhsindex
            )
        end
    end

    coefficients
end

function leastsquaressolution(
    design::AbstractMatrix{T},
    response::AbstractVector{T},
    nonnegative::Bool
) where {T<:AbstractFloat}
    nonnegative ? nonnegativeleastsquares(design, response) : design \ response
end

function parafac2updatecore(
    transformed::AbstractVector{<:AbstractMatrix{T}},
    loadings::AbstractMatrix{T},
    weights::AbstractMatrix{T}
) where {T<:AbstractFloat}
    ncomponents = size(weights, 2)
    mzcount = size(loadings, 1)
    nsamples = length(transformed)
    Ywide = Matrix{T}(undef, ncomponents, mzcount * nsamples)
    Bwide = Matrix{T}(undef, ncomponents, mzcount * nsamples)

    if parafac2threaded(nsamples)
        @threads :static for sampleindex in 1:nsamples
            cols = ((sampleindex - 1) * mzcount + 1):(sampleindex * mzcount)
            Ywide[:, cols] .= transformed[sampleindex]
            Bwide[:, cols] .= Diagonal(view(weights, sampleindex, :)) * transpose(loadings)
        end
    else
        for sampleindex in 1:nsamples
            cols = ((sampleindex - 1) * mzcount + 1):(sampleindex * mzcount)
            Ywide[:, cols] .= transformed[sampleindex]
            Bwide[:, cols] .= Diagonal(view(weights, sampleindex, :)) * transpose(loadings)
        end
    end

    Matrix{T}(transpose(transpose(Bwide) \ transpose(Ywide)))
end

function parafac2updateloadings(
    transformed::AbstractVector{<:AbstractMatrix{T}},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T},
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    ncomponents = size(weights, 2)
    mzcount = size(first(transformed), 2)
    nsamples = length(transformed)
    Ywide = Matrix{T}(undef, mzcount, ncomponents * nsamples)
    Bwide = Matrix{T}(undef, ncomponents, ncomponents * nsamples)

    if parafac2threaded(nsamples)
        @threads :static for sampleindex in 1:nsamples
            cols = ((sampleindex - 1) * ncomponents + 1):(sampleindex * ncomponents)
            Ywide[:, cols] .= transpose(transformed[sampleindex])
            Bwide[:, cols] .= Diagonal(view(weights, sampleindex, :)) * transpose(core)
        end
    else
        for sampleindex in 1:nsamples
            cols = ((sampleindex - 1) * ncomponents + 1):(sampleindex * ncomponents)
            Ywide[:, cols] .= transpose(transformed[sampleindex])
            Bwide[:, cols] .= Diagonal(view(weights, sampleindex, :)) * transpose(core)
        end
    end

    design = transpose(Bwide)
    response = transpose(Ywide)

    if parafac2constrainspectra(nonnegative)
        coefficients = nonnegativeleastsquares_batch(design, response)
        return Matrix{T}(transpose(coefficients))
    end

    Matrix{T}(transpose(design \ response))
end

function parafac2updateweights(
    transformed::AbstractVector{<:AbstractMatrix{T}},
    core::AbstractMatrix{T},
    loadings::AbstractMatrix{T},
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    nsamples = length(transformed)
    ncomponents = size(core, 2)
    weights = Matrix{T}(undef, nsamples, ncomponents)
    design = khatrirao(loadings, core)

    if parafac2constrainintensities(nonnegative)
        responses = Matrix{T}(undef, size(design, 1), nsamples)
        if parafac2threaded(nsamples)
            @threads :static for sampleindex in 1:nsamples
                responses[:, sampleindex] .= vec(transformed[sampleindex])
            end
        else
            for sampleindex in 1:nsamples
                responses[:, sampleindex] .= vec(transformed[sampleindex])
            end
        end
        coefficients = nonnegativeleastsquares_batch(design, responses)
        weights .= transpose(coefficients)
    else
        if parafac2threaded(nsamples)
            @threads :static for sampleindex in 1:nsamples
                weights[sampleindex, :] .= design \ vec(transformed[sampleindex])
            end
        else
            for sampleindex in 1:nsamples
                weights[sampleindex, :] .= design \ vec(transformed[sampleindex])
            end
        end
    end

    weights
end

function parafac2normalize!(
    core::AbstractMatrix{T},
    loadings::AbstractMatrix{T},
    weights::AbstractMatrix{T}
) where {T<:AbstractFloat}
    for component in axes(loadings, 2)
        corenorm = norm(view(core, :, component))
        if isfinite(corenorm) && corenorm > eps(T)
            core[:, component] ./= corenorm
            weights[:, component] .*= corenorm
        end

        loadingnorm = norm(view(loadings, :, component))
        if isfinite(loadingnorm) && loadingnorm > eps(T)
            loadings[:, component] ./= loadingnorm
            weights[:, component] .*= loadingnorm
        end
    end

    core, loadings, weights
end

function parafac2cpcycle(
    transformed::AbstractVector{<:AbstractMatrix{T}},
    core::AbstractMatrix{T},
    loadings::AbstractMatrix{T},
    weights::AbstractMatrix{T},
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    core_new = parafac2updatecore(transformed, loadings, weights)
    loadings_new = parafac2updateloadings(transformed, core_new, weights, nonnegative)
    weights_new = parafac2updateweights(transformed, core_new, loadings_new, nonnegative)
    parafac2normalize!(core_new, loadings_new, weights_new)

    core_new, loadings_new, weights_new
end

function parafac2loss(
    X::AbstractVector{<:AbstractMatrix{T}},
    bases::AbstractVector{<:AbstractMatrix{T}},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T},
    loadings::AbstractMatrix{T}
) where {T<:AbstractFloat}
    if parafac2threaded(length(X))
        sampleindices = collect(eachindex(X))
        nchunks = min(Base.Threads.nthreads(), length(sampleindices))
        partials = zeros(T, nchunks)
        @threads :static for chunkindex in 1:nchunks
            chunkloss = zero(T)
            firstposition = div((chunkindex - 1) * length(sampleindices), nchunks) + 1
            lastposition = div(chunkindex * length(sampleindices), nchunks)
            for position in firstposition:lastposition
                sampleindex = sampleindices[position]
                fitted = bases[sampleindex] *
                    core *
                    Diagonal(view(weights, sampleindex, :)) *
                    transpose(loadings)
                chunkloss += sum(abs2, X[sampleindex] .- fitted)
            end
            partials[chunkindex] = chunkloss
        end

        return sum(partials)
    end

    loss = zero(T)
    for sampleindex in eachindex(X)
        fitted = bases[sampleindex] *
            core *
            Diagonal(view(weights, sampleindex, :)) *
            transpose(loadings)
        loss += sum(abs2, X[sampleindex] .- fitted)
    end

    loss
end

parafac2usestransformedloss(::AbstractVector{<:Matrix}) = true
parafac2usestransformedloss(::AbstractVector) = false

function parafac2smallresidualloss(
    transformed::AbstractMatrix{T},
    core::AbstractMatrix{T},
    weights::AbstractVector{T},
    loadings::AbstractMatrix{T}
) where {T<:AbstractFloat}
    loss = zero(T)
    for mzindex in axes(loadings, 1)
        for component in axes(core, 1)
            fitted = zero(T)
            for latent in axes(core, 2)
                fitted += core[component, latent] *
                    weights[latent] *
                    loadings[mzindex, latent]
            end
            residual = transformed[component, mzindex] - fitted
            loss += residual * residual
        end
    end

    loss
end

function parafac2transformedloss_sample(
    transformed::AbstractMatrix{T},
    datanorm::T,
    core::AbstractMatrix{T},
    weights::AbstractVector{T},
    loadings::AbstractMatrix{T}
) where {T<:AbstractFloat}
    datanorm - sum(abs2, transformed) +
        parafac2smallresidualloss(transformed, core, weights, loadings)
end

function parafac2loss(
    transformed::AbstractVector{<:AbstractMatrix{T}},
    datanorms::AbstractVector{T},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T},
    loadings::AbstractMatrix{T}
) where {T<:AbstractFloat}
    length(transformed) == length(datanorms) || throw(DimensionMismatch(
        "transformed data and data norms must have the same length."
    ))

    if parafac2threaded(length(transformed))
        nsamples = length(transformed)
        nchunks = min(Base.Threads.nthreads(), nsamples)
        partials = zeros(T, nchunks)
        @threads :static for chunkindex in 1:nchunks
            chunkloss = zero(T)
            firstposition = div((chunkindex - 1) * nsamples, nchunks) + 1
            lastposition = div(chunkindex * nsamples, nchunks)
            for sampleindex in firstposition:lastposition
                chunkloss += parafac2transformedloss_sample(
                    transformed[sampleindex],
                    datanorms[sampleindex],
                    core,
                    view(weights, sampleindex, :),
                    loadings
                )
            end
            partials[chunkindex] = chunkloss
        end

        return sum(partials)
    end

    loss = zero(T)
    for sampleindex in eachindex(transformed)
        loss += parafac2transformedloss_sample(
            transformed[sampleindex],
            datanorms[sampleindex],
            core,
            view(weights, sampleindex, :),
            loadings
        )
    end

    loss
end

function parafac2transformedlossunstable(
    loss::T,
    datanorms::AbstractVector{T}
) where {T<:AbstractFloat}
    !isfinite(loss) && return true
    loss ≤ sqrt(eps(T)) * max(sum(datanorms), one(T))
end

function parafac2checkedtransformedloss(
    X::AbstractVector{<:AbstractMatrix{T}},
    bases::AbstractVector{<:AbstractMatrix{T}},
    transformed::AbstractVector{<:AbstractMatrix{T}},
    datanorms::AbstractVector{T},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T},
    loadings::AbstractMatrix{T}
) where {T<:AbstractFloat}
    loss = parafac2loss(transformed, datanorms, core, weights, loadings)
    if parafac2transformedlossunstable(loss, datanorms)
        return parafac2loss(X, bases, core, weights, loadings)
    end

    loss
end

function parafac2requirefitted(fit::Parafac2Fit)
    if isnothing(fit.loadings) ||
            isnothing(fit.core) ||
            isnothing(fit.weights) ||
            isnothing(fit.bases)
        throw(ArgumentError("fit does not contain fitted PARAFAC2 factors."))
    end
    nothing
end

function parafac2checksampleindex(fit::Parafac2Fit, sampleindex::Integer)
    1 ≤ sampleindex ≤ length(fit.retentioncounts) || throw(BoundsError(
        fit.retentioncounts,
        sampleindex
    ))
    Int(sampleindex)
end

"""
    parafac2scores(fit::Parafac2Fit)
    parafac2scores(fit::Parafac2Fit, sampleindex::Integer)

Return sample-specific PARAFAC2 score matrices.

For sample `k`, scores are computed as

    fit.bases[k] * fit.core

The one-argument method returns one score matrix per sample. The two-argument method
returns only the requested sample's score matrix.

See also [`parafac2reconstruct`](@ref).
"""
function parafac2scores(fit::Parafac2Fit)
    parafac2requirefitted(fit)
    [fit.bases[sampleindex] * fit.core for sampleindex in eachindex(fit.bases)]
end

function parafac2scores(fit::Parafac2Fit, sampleindex::Integer)
    parafac2requirefitted(fit)
    checkedindex = parafac2checksampleindex(fit, sampleindex)
    fit.bases[checkedindex] * fit.core
end

"""
    parafac2profileminima(fit::Parafac2Fit)

Return the minimum fitted profile value for each sample/component as a
`samples × components` matrix.

Negative entries indicate that the corresponding fitted chromatographic profile contains
negative values.

See also [`parafac2scores`](@ref), [`parafac2profilediagnostics`](@ref).
"""
function parafac2profileminima(fit::Parafac2Fit)
    profiles = parafac2scores(fit)
    minima = Matrix{eltype(fit.core)}(undef, length(profiles), fit.ncomponents)

    for sampleindex in eachindex(profiles)
        for component in axes(profiles[sampleindex], 2)
            minima[sampleindex, component] = minimum(view(profiles[sampleindex], :, component))
        end
    end

    minima
end

"""
    parafac2profilediagnostics(fit::Parafac2Fit)

Return profile sign diagnostics as a named tuple.

The returned fields are `minima`, `minimum`, `hasnegative`, and `negativecounts`.
`minima` contains the minimum profile value for each sample/component, and
`negativecounts` contains the number of negative profile entries for each
sample/component.

See also [`parafac2profileminima`](@ref), [`parafac2scores`](@ref).
"""
function parafac2profilediagnostics(fit::Parafac2Fit)
    profiles = parafac2scores(fit)
    minima = Matrix{eltype(fit.core)}(undef, length(profiles), fit.ncomponents)
    negativecounts = Matrix{Int}(undef, length(profiles), fit.ncomponents)

    for sampleindex in eachindex(profiles)
        for component in axes(profiles[sampleindex], 2)
            profile = view(profiles[sampleindex], :, component)
            minima[sampleindex, component] = minimum(profile)
            negativecounts[sampleindex, component] =
                Base.count(<(zero(eltype(fit.core))), profile)
        end
    end

    (
        minima=minima,
        minimum=minimum(minima),
        hasnegative=any(<(zero(eltype(fit.core))), minima),
        negativecounts=negativecounts,
    )
end

"""
    parafac2spectra(fit::Parafac2Fit; metadata=false, unit=nothing)

Return component spectral loadings as an `m/z × components` matrix.

By default this returns a copy of `fit.loadings`. If `metadata=true`, return a named
tuple `(values=spectra, mzvalues=mzvalues)`, where `mzvalues` is `nothing` when no m/z
axis metadata was stored. Pass `unit` with `metadata=true` to convert stored unitful m/z
metadata.

Fits use nonnegative spectral loadings by default. If the fit was created with
`nonnegative=()` or without `:spectra`, these loadings can be signed.

See also [`parafac2intensities`](@ref), [`parafac2apexes`](@ref).
"""
function parafac2spectra(
    fit::Parafac2Fit;
    metadata::Bool=false,
    unit::Union{Nothing, Unitful.Units}=nothing
)
    parafac2requirefitted(fit)
    metadata || isnothing(unit) || throw(ArgumentError(
        "unit can only be used when metadata=true."
    ))
    spectra = copy(fit.loadings)
    metadata || return spectra

    (values=spectra, mzvalues=mzvalues(fit; unit=unit))
end

"""
    parafac2intensities(fit::Parafac2Fit)

Return sample/component weights as a `samples × components` matrix.

Fits use nonnegative intensity weights by default. If the fit was created with
`nonnegative=()` or without `:intensities`, these values can be signed and retain the
scale convention of the selected solution.

See also [`parafac2spectra`](@ref), [`parafac2apexes`](@ref).
"""
function parafac2intensities(fit::Parafac2Fit)
    parafac2requirefitted(fit)
    copy(fit.weights)
end

function parafac2positivearea(
    values::AbstractVector{<:Real},
    axis::Union{Nothing, AbstractVector{<:Real}}
)
    if !isnothing(axis)
        length(axis) == length(values) || throw(DimensionMismatch(
            "integration axis length must match signal length."
        ))
    end

    area = 0.0
    for index in 1:(length(values) - 1)
        left = Float64(values[index])
        right = Float64(values[index + 1])
        dx = isnothing(axis) ? 1.0 : Float64(axis[index + 1]) - Float64(axis[index])

        if left > 0 && right > 0
            area += 0.5 * (left + right) * dx
        elseif left > 0
            area += 0.5 * left * dx * left / (left - right)
        elseif right > 0
            area += 0.5 * right * dx * right / (right - left)
        end
    end

    area
end

"""
    parafac2areas(fit::Parafac2Fit)

Return positive component peak areas as a `samples × components` matrix.

For sample `k` and component `r`, the integrated signal is

    parafac2intensities(fit)[k, r] .* parafac2scores(fit, k)[:, r]

Areas are trapezoidal integrals over the stored retention/Kovats axis when retention
metadata are available, and over unit-spaced scan indices otherwise. Only the positive part
of each linearly interpolated signal is integrated, so sign changes are handled by linear
zero-crossing interpolation rather than by clipping sampled values.

See also [`parafac2scores`](@ref), [`parafac2intensities`](@ref).
"""
function parafac2areas(fit::Parafac2Fit)
    parafac2requirefitted(fit)

    weights = parafac2intensities(fit)
    retentionvectors = rawretentions(fit)
    nsamples = length(fit.retentioncounts)
    areas = Matrix{Float64}(undef, nsamples, fit.ncomponents)

    for sampleindex in 1:nsamples
        profiles = parafac2scores(fit, sampleindex)
        axis = isnothing(retentionvectors) ? nothing : retentionvectors[sampleindex]
        for component in axes(profiles, 2)
            signal = weights[sampleindex, component] .* view(profiles, :, component)
            areas[sampleindex, component] = parafac2positivearea(signal, axis)
        end
    end

    areas
end

"""
    parafac2apexes(fit::Parafac2Fit; unit=nothing)

Return apex information for each sample/component retention profile.

The return value is a named tuple `(indices=indices, retentions=retentions,
values=values)`, where each matrix has dimensions `samples × components`. `indices`
contains the retention index of the maximum value in `parafac2scores(fit, sample)[:, component]`.
`values` contains the corresponding profile values. `retentions` is `nothing` when the
fit has no stored retention metadata; otherwise it contains the apex retention coordinate,
converted to `unit` when requested.

The apex is the maximum fitted profile value, not the maximum absolute value. If the fit
was created without nonnegative intensity weights, profile signs are model-dependent.

See also [`parafac2scores`](@ref), [`parafac2spectra`](@ref).
"""
function parafac2apexes(
    fit::Parafac2Fit;
    unit::Union{Nothing, Unitful.Units}=nothing
)
    parafac2requirefitted(fit)
    profiles = parafac2scores(fit)
    nsamples = length(profiles)
    ncomponents = fit.ncomponents
    indices = Matrix{Int}(undef, nsamples, ncomponents)
    values = Matrix{eltype(fit.core)}(undef, nsamples, ncomponents)
    retentionvectors = retentions(fit; unit=unit)
    apexretentions = if isnothing(retentionvectors)
        nothing
    else
        Matrix{typeof(retentionvectors[1][1])}(undef, nsamples, ncomponents)
    end

    for sampleindex in eachindex(profiles)
        profilematrix = profiles[sampleindex]
        for component in axes(profilematrix, 2)
            profile = view(profilematrix, :, component)
            apexindex = argmax(profile)
            indices[sampleindex, component] = apexindex
            values[sampleindex, component] = profile[apexindex]
            if !isnothing(apexretentions)
                apexretentions[sampleindex, component] = retentionvectors[sampleindex][apexindex]
            end
        end
    end

    (indices=indices, retentions=apexretentions, values=values)
end

"""
    parafac2reconstruct(fit::Parafac2Fit)
    parafac2reconstruct(fit::Parafac2Fit, sampleindex::Integer)

Reconstruct fitted sample matrix or matrices from a PARAFAC2 fit.

For sample `k`, the reconstruction is

    parafac2scores(fit, k) * Diagonal(fit.weights[k, :]) * fit.loadings'

The one-argument method returns one reconstructed `retentions × m/z` matrix per sample.
The two-argument method returns only the requested sample's reconstruction.
"""
function parafac2reconstruct(fit::Parafac2Fit, sampleindex::Integer)
    parafac2requirefitted(fit)
    checkedindex = parafac2checksampleindex(fit, sampleindex)
    parafac2scores(fit, checkedindex) *
        Diagonal(view(fit.weights, checkedindex, :)) *
        transpose(fit.loadings)
end

function parafac2reconstruct(fit::Parafac2Fit)
    parafac2requirefitted(fit)
    [parafac2reconstruct(fit, sampleindex) for sampleindex in eachindex(fit.bases)]
end

function parafac2diagnosticmatrices(
    fit::Parafac2Fit,
    X::AbstractVector{<:AbstractMatrix{<:Real}}
)
    parafac2requirefitted(fit)
    length(X) == length(fit.retentioncounts) || throw(DimensionMismatch(
        "X must contain $(length(fit.retentioncounts)) sample matrices, got $(length(X))."
    ))

    T = float(promote_type(eltype(fit.loadings), map(eltype, X)...))
    matrices = Vector{Matrix{T}}(undef, length(X))
    for sampleindex in eachindex(X)
        expectedsize = (fit.retentioncounts[sampleindex], fit.mzcount)
        size(X[sampleindex]) == expectedsize || throw(DimensionMismatch(
            "X[$(sampleindex)] must have size $(expectedsize), got $(size(X[sampleindex]))."
        ))
        all(isfinite, X[sampleindex]) || throw(ArgumentError(
            "X[$(sampleindex)] must contain only finite values."
        ))
        matrices[sampleindex] = Matrix{T}(X[sampleindex])
    end

    matrices
end

function parafac2diagnosticmatrices(
    fit::Parafac2Fit,
    X::AbstractArray{<:Real, 3}
)
    samples = [view(X, k, :, :) for k in axes(X, 1)]
    parafac2diagnosticmatrices(fit, samples)
end

"""
    parafac2residuals(fit::Parafac2Fit, X)

Return residual matrices `X[k] - parafac2reconstruct(fit, k)` for each sample.

`X` can be either the original vector of `retentions × m/z` matrices or a
`samples × retentions × m/z` tensor. The dimensions must match the fitted data layout.

See also [`parafac2loss`](@ref), [`parafac2fitpercent`](@ref).
"""
function parafac2residuals(
    fit::Parafac2Fit,
    X::Union{AbstractVector{<:AbstractMatrix{<:Real}}, AbstractArray{<:Real, 3}}
)
    matrices = parafac2diagnosticmatrices(fit, X)
    reconstructions = parafac2reconstruct(fit)
    [matrices[sampleindex] .- reconstructions[sampleindex] for sampleindex in eachindex(matrices)]
end

"""
    parafac2loss(fit::Parafac2Fit, X)

Return the residual sum of squares for `fit` against data `X`.

`X` can be either a vector of sample matrices or a `samples × retentions × m/z` tensor.

See also [`parafac2residuals`](@ref), [`parafac2fitpercent`](@ref).
"""
function parafac2loss(
    fit::Parafac2Fit,
    X::Union{AbstractVector{<:AbstractMatrix{<:Real}}, AbstractArray{<:Real, 3}}
)
    residuals = parafac2residuals(fit, X)
    sum(sum(abs2, residual) for residual in residuals)
end

"""
    parafac2fitpercent(fit::Parafac2Fit, X)

Return the fraction of total sum of squares explained by the PARAFAC2 fit.

The returned value is

    1 - parafac2loss(fit, X) / sum(abs2, X)

If the total sum of squares is zero, the return value is `NaN`.

See also [`parafac2loss`](@ref), [`parafac2reconstruct`](@ref).
"""
function parafac2fitpercent(
    fit::Parafac2Fit,
    X::Union{AbstractVector{<:AbstractMatrix{<:Real}}, AbstractArray{<:Real, 3}}
)
    matrices = parafac2diagnosticmatrices(fit, X)
    totalsumsquares = sum(sum(abs2, matrix) for matrix in matrices)
    totalsumsquares == 0 && return NaN
    1 - parafac2loss(fit, matrices) / totalsumsquares
end

function parafac2rationalstart(
    X::AbstractVector{<:AbstractMatrix{T}},
    ncomponents::Integer,
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    loadings = parafac2initialloadings(X, ncomponents)
    if parafac2constrainspectra(nonnegative)
        loadings .= abs.(loadings)
    end
    core = Matrix{T}(I, ncomponents, ncomponents)
    weights = ones(T, length(X), ncomponents)
    parafac2normalize!(core, loadings, weights)

    loadings, core, weights
end

function parafac2randomstart(
    X::AbstractVector{<:AbstractMatrix{T}},
    ncomponents::Integer,
    rng::Random.AbstractRNG,
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    mzcount = size(first(X), 2)
    decomposition = svd(Random.randn(rng, T, mzcount, ncomponents); full=false)
    loadings = Matrix{T}(decomposition.U[:, 1:ncomponents])
    if parafac2constrainspectra(nonnegative)
        loadings .= abs.(loadings)
    end
    core = Random.randn(rng, T, ncomponents, ncomponents)
    weights = Random.randn(rng, T, length(X), ncomponents)
    if parafac2constrainintensities(nonnegative)
        weights .= abs.(weights)
    end
    parafac2normalize!(core, loadings, weights)

    loadings, core, weights
end

function parafac2fitstart(
    X::AbstractVector{<:AbstractMatrix{T}},
    loadings::AbstractMatrix{T},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T},
    maxiters::Integer,
    tol::Real,
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    parafac2fitstart(
        X,
        parafac2squarednorms(X),
        loadings,
        core,
        weights,
        maxiters,
        tol,
        nonnegative
    )
end

function parafac2fitstart(
    X::AbstractVector{<:AbstractMatrix{T}},
    datanorms::AbstractVector{T},
    loadings::AbstractMatrix{T},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T},
    maxiters::Integer,
    tol::Real,
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    if parafac2usestransformedloss(X)
        return parafac2fitstart_transformedloss(
            X,
            datanorms,
            loadings,
            core,
            weights,
            maxiters,
            tol,
            nonnegative
        )
    end

    parafac2fitstart_fullloss(
        X,
        loadings,
        core,
        weights,
        maxiters,
        tol,
        nonnegative
    )
end

function parafac2fitstart_fullloss(
    X::AbstractVector{<:AbstractMatrix{T}},
    loadings::AbstractMatrix{T},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T},
    maxiters::Integer,
    tol::Real,
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    bases = parafac2updatebases(X, loadings, core, weights)
    losses = T[parafac2loss(X, bases, core, weights, loadings)]
    converged = false
    stopreason = :maxiters
    iterations = 0

    for iteration in 1:maxiters
        previousloss = last(losses)

        previousloadings = loadings
        previouscore = core
        previousweights = weights
        previousbases = bases

        transformed = parafac2transformeddata(X, bases)
        core, loadings, weights = parafac2cpcycle(
            transformed,
            core,
            loadings,
            weights,
            nonnegative
        )
        bases = parafac2updatebases(X, loadings, core, weights)
        currentloss = parafac2loss(X, bases, core, weights, loadings)

        increasetolerance = sqrt(eps(T)) * max(one(T), previousloss)
        if !isfinite(currentloss) || currentloss > previousloss + increasetolerance
            loadings = previousloadings
            core = previouscore
            weights = previousweights
            bases = previousbases
            stopreason = isfinite(currentloss) ? :objective_increase : :nonfinite
            break
        elseif currentloss > previousloss
            loadings = previousloadings
            core = previouscore
            weights = previousweights
            bases = previousbases
            converged = true
            stopreason = :tol
            break
        end

        push!(losses, currentloss)
        iterations = iteration

        improvement = previousloss - currentloss
        threshold = T(tol) * max(previousloss, eps(T))
        if improvement ≤ threshold
            converged = true
            stopreason = :tol
            break
        end
    end

    (
        loadings=loadings,
        core=core,
        weights=weights,
        bases=bases,
        losses=losses,
        converged=converged,
        stopreason=stopreason,
        iterations=iterations
    )
end

function parafac2fitstart_transformedloss(
    X::AbstractVector{<:AbstractMatrix{T}},
    datanorms::AbstractVector{T},
    loadings::AbstractMatrix{T},
    core::AbstractMatrix{T},
    weights::AbstractMatrix{T},
    maxiters::Integer,
    tol::Real,
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    bases = parafac2updatebases(X, loadings, core, weights)
    candidatebases = parafac2basisworkspace(X, size(core, 1))
    transformed = parafac2transformeddata(X, bases)
    candidatetransformed = parafac2transformedworkspace(X, bases)
    losses = T[
        parafac2checkedtransformedloss(
            X,
            bases,
            transformed,
            datanorms,
            core,
            weights,
            loadings
        )
    ]
    converged = false
    stopreason = :maxiters
    iterations = 0

    for iteration in 1:maxiters
        previousloss = last(losses)

        previousloadings = loadings
        previouscore = core
        previousweights = weights

        core, loadings, weights = parafac2cpcycle(
            transformed,
            core,
            loadings,
            weights,
            nonnegative
        )
        parafac2updatebases!(candidatebases, X, loadings, core, weights)
        parafac2transformeddata!(candidatetransformed, X, candidatebases)
        currentloss = parafac2checkedtransformedloss(
            X,
            candidatebases,
            candidatetransformed,
            datanorms,
            core,
            weights,
            loadings
        )

        increasetolerance = sqrt(eps(T)) * max(one(T), previousloss)
        if !isfinite(currentloss) || currentloss > previousloss + increasetolerance
            loadings = previousloadings
            core = previouscore
            weights = previousweights
            stopreason = isfinite(currentloss) ? :objective_increase : :nonfinite
            break
        elseif currentloss > previousloss
            loadings = previousloadings
            core = previouscore
            weights = previousweights
            converged = true
            stopreason = :tol
            break
        end

        bases, candidatebases = candidatebases, bases
        transformed, candidatetransformed = candidatetransformed, transformed
        push!(losses, currentloss)
        iterations = iteration

        improvement = previousloss - currentloss
        threshold = T(tol) * max(previousloss, eps(T))
        if improvement ≤ threshold
            converged = true
            stopreason = :tol
            break
        end
    end

    (
        loadings=loadings,
        core=core,
        weights=weights,
        bases=bases,
        losses=losses,
        converged=converged,
        stopreason=stopreason,
        iterations=iterations
    )
end

function parafac2startdiagnostic(
    startindex::Integer,
    fit
)
    (
        start=Int(startindex),
        iterations=Int(fit.iterations),
        converged=Bool(fit.converged),
        stopreason=fit.stopreason,
        loss=last(fit.losses)
    )
end

function parafac2fit(
    X::AbstractVector{<:AbstractMatrix{T}},
    ncomponents::Integer,
    maxiters::Integer,
    tol::Real,
    nstarts::Integer,
    rng::Random.AbstractRNG,
    nonnegative::Tuple{Vararg{Symbol}}
) where {T<:AbstractFloat}
    nstarts ≥ 1 || throw(ArgumentError("nstarts must be at least 1."))
    datanorms = parafac2squarednorms(X)

    loadings, core, weights = parafac2rationalstart(X, ncomponents, nonnegative)
    bestfit = parafac2fitstart(
        X,
        datanorms,
        loadings,
        core,
        weights,
        maxiters,
        tol,
        nonnegative
    )
    bestloss = last(bestfit.losses)
    beststart = 1
    startdiagnostics = Vector{Parafac2StartDiagnostic{T}}(undef, Int(nstarts))
    startdiagnostics[1] = parafac2startdiagnostic(1, bestfit)

    for startindex in 2:nstarts
        loadings, core, weights = parafac2randomstart(X, ncomponents, rng, nonnegative)
        currentfit = parafac2fitstart(
            X,
            datanorms,
            loadings,
            core,
            weights,
            maxiters,
            tol,
            nonnegative
        )
        currentloss = last(currentfit.losses)
        startdiagnostics[startindex] = parafac2startdiagnostic(startindex, currentfit)
        if currentloss < bestloss
            bestfit = currentfit
            bestloss = currentloss
            beststart = startindex
        end
    end

    (
        bestfit.loadings,
        bestfit.core,
        bestfit.weights,
        bestfit.bases,
        bestfit.losses,
        bestfit.converged,
        bestfit.stopreason,
        bestfit.iterations,
        beststart,
        startdiagnostics
    )
end

function parafac2coordinatevector(
    values,
    expectedlength::Integer,
    name::AbstractString;
    targetunit::Union{Nothing, Unitful.Units}=nothing,
    requirepositive::Bool=false
)
    length(values) == expectedlength || throw(DimensionMismatch(
        "$name must have length $(expectedlength), got $(length(values))."
    ))

    hasunits = any(isunitful, values)
    rawvalues = Vector{Float64}(undef, length(values))
    valueunit = targetunit

    if hasunits
        all(isunitful, values) || throw(ArgumentError(
            "$name must not mix unitless and unitful values."
        ))
        if isnothing(valueunit)
            valueunit = unit(first(values))
        end
        for (i, value) in enumerate(values)
            try
                rawvalues[i] = Float64(ustrip(valueunit, value))
            catch err
                err isa Unitful.DimensionError || rethrow()
                throw(ArgumentError("$name contains values incompatible with $(valueunit)."))
            end
        end
    else
        isnothing(targetunit) || throw(ArgumentError(
            "$name must have units compatible with $(targetunit)."
        ))
        all(value -> value isa Real, values) || throw(ArgumentError(
            "$name must contain real values or Unitful quantities."
        ))
        rawvalues .= Float64.(values)
    end

    all(isfinite, rawvalues) || throw(ArgumentError("$name must contain only finite values."))
    all(diff(rawvalues) .> 0) || throw(ArgumentError("$name must be strictly increasing."))
    if requirepositive
        all(rawvalues .> 0) || throw(ArgumentError("$name must be positive."))
    end

    rawvalues, valueunit
end

function parafac2retentionmetadata(retentions, retentioncounts::AbstractVector{<:Integer})
    nsamples = length(retentioncounts)
    isnothing(retentions) && return nothing, nothing

    retentionvectors = Vector{Vector{Float64}}(undef, nsamples)
    retentionunit = nothing

    if retentions isa AbstractMatrix
        size(retentions, 1) == nsamples || throw(DimensionMismatch(
            "retentions must have one row per sample ($(nsamples)), got $(size(retentions, 1))."
        ))
        for sampleindex in 1:nsamples
            values, sampleunit = parafac2coordinatevector(
                view(retentions, sampleindex, :),
                retentioncounts[sampleindex],
                "retentions for sample $(sampleindex)";
                targetunit=retentionunit
            )
            if sampleindex > 1 && isnothing(retentionunit) != isnothing(sampleunit)
                throw(ArgumentError(
                    "retentions must be consistently unitless or consistently unitful."
                ))
            end
            retentionunit = sampleindex == 1 ? sampleunit : retentionunit
            retentionvectors[sampleindex] = values
        end
        return retentionvectors, retentionunit
    end

    if retentions isa AbstractVector && all(value -> value isa AbstractVector, retentions)
        length(retentions) == nsamples || throw(DimensionMismatch(
            "retentions must contain one vector per sample ($(nsamples)), got $(length(retentions))."
        ))
        for sampleindex in 1:nsamples
            values, sampleunit = parafac2coordinatevector(
                retentions[sampleindex],
                retentioncounts[sampleindex],
                "retentions for sample $(sampleindex)";
                targetunit=retentionunit
            )
            if sampleindex > 1 && isnothing(retentionunit) != isnothing(sampleunit)
                throw(ArgumentError(
                    "retentions must be consistently unitless or consistently unitful."
                ))
            end
            retentionunit = sampleindex == 1 ? sampleunit : retentionunit
            retentionvectors[sampleindex] = values
        end
        return retentionvectors, retentionunit
    end

    sharedvalues, retentionunit = parafac2coordinatevector(
        retentions,
        first(retentioncounts),
        "retentions"
    )
    for sampleindex in 1:nsamples
        length(sharedvalues) == retentioncounts[sampleindex] || throw(DimensionMismatch(
            "shared retentions have length $(length(sharedvalues)), but sample " *
            "$(sampleindex) has $(retentioncounts[sampleindex]) retention positions."
        ))
        retentionvectors[sampleindex] = copy(sharedvalues)
    end

    retentionvectors, retentionunit
end

function parafac2mzmetadata(mzvalues, mzcount::Integer)
    isnothing(mzvalues) && return nothing, nothing
    parafac2coordinatevector(mzvalues, mzcount, "m/z values"; requirepositive=true)
end

function parafac2samplelabelmetadata(samplelabels, nsamples::Integer)
    isnothing(samplelabels) && return nothing
    length(samplelabels) == nsamples || throw(DimensionMismatch(
        "samplelabels must have length $(nsamples), got $(length(samplelabels))."
    ))

    labels = string.(collect(samplelabels))
    all(!isempty, labels) || throw(ArgumentError("samplelabels must not contain empty labels."))
    length(unique(labels)) == length(labels) || throw(ArgumentError(
        "samplelabels must be unique."
    ))

    labels
end

function convertparafac2rawvalues(
    ::Nothing,
    ::Union{Nothing, Unitful.Units},
    ::AbstractString,
    ::Union{Nothing, Unitful.Units}
)
    nothing
end

function convertparafac2rawvalues(
    values::Union{Nothing, Vector{Float64}},
    valueunit::Union{Nothing, Unitful.Units},
    name::AbstractString,
    targetunit::Union{Nothing, Unitful.Units}
)
    isnothing(values) && return nothing
    if isnothing(valueunit)
        isnothing(targetunit) || throw(ArgumentError("Cannot convert unitless $name to a unit."))
        return copy(values)
    end

    outputunit = isnothing(targetunit) ? valueunit : targetunit
    Float64.(ustrip.(Base.RefValue(outputunit), values .* valueunit))
end

function convertparafac2rawvalues(
    values::Union{Nothing, Vector{Vector{Float64}}},
    valueunit::Union{Nothing, Unitful.Units},
    name::AbstractString,
    targetunit::Union{Nothing, Unitful.Units}
)
    isnothing(values) && return nothing
    if isnothing(valueunit)
        isnothing(targetunit) || throw(ArgumentError("Cannot convert unitless $name to a unit."))
        return copy.(values)
    end

    outputunit = isnothing(targetunit) ? valueunit : targetunit
    [Float64.(ustrip.(Base.RefValue(outputunit), value .* valueunit)) for value in values]
end

function convertparafac2values(
    ::Nothing,
    ::Union{Nothing, Unitful.Units},
    ::AbstractString,
    ::Union{Nothing, Unitful.Units}
)
    nothing
end

function convertparafac2values(
    values::Union{Nothing, Vector{Float64}},
    valueunit::Union{Nothing, Unitful.Units},
    name::AbstractString,
    targetunit::Union{Nothing, Unitful.Units}
)
    isnothing(values) && return nothing
    if isnothing(valueunit)
        isnothing(targetunit) || throw(ArgumentError("Cannot convert unitless $name to a unit."))
        return copy(values)
    end

    outputunit = isnothing(targetunit) ? valueunit : targetunit
    uconvert.(Base.RefValue(outputunit), values .* valueunit)
end

function convertparafac2values(
    values::Union{Nothing, Vector{Vector{Float64}}},
    valueunit::Union{Nothing, Unitful.Units},
    name::AbstractString,
    targetunit::Union{Nothing, Unitful.Units}
)
    isnothing(values) && return nothing
    if isnothing(valueunit)
        isnothing(targetunit) || throw(ArgumentError("Cannot convert unitless $name to a unit."))
        return copy.(values)
    end

    outputunit = isnothing(targetunit) ? valueunit : targetunit
    [uconvert.(Base.RefValue(outputunit), value .* valueunit) for value in values]
end

rawretentions(fit::Parafac2Fit; unit::Union{Nothing, Unitful.Units}=nothing) =
    convertparafac2rawvalues(fit.retentions, fit.retentionunit, "retentions", unit)

retentions(fit::Parafac2Fit; unit::Union{Nothing, Unitful.Units}=nothing) =
    convertparafac2values(fit.retentions, fit.retentionunit, "retentions", unit)

retentionunit(fit::Parafac2Fit) = fit.retentionunit

rawmzvalues(fit::Parafac2Fit; unit::Union{Nothing, Unitful.Units}=nothing) =
    convertparafac2rawvalues(fit.mzvalues, fit.mzunit, "m/z values", unit)

mzvalues(fit::Parafac2Fit; unit::Union{Nothing, Unitful.Units}=nothing) =
    convertparafac2values(fit.mzvalues, fit.mzunit, "m/z values", unit)

mzunit(fit::Parafac2Fit) = fit.mzunit

function parafac2dimensions(
    X::AbstractVector{<:AbstractMatrix{<:Real}},
    ncomponents::Integer
)
    isempty(X) && throw(ArgumentError("X must contain at least one sample matrix."))
    ncomponents ≥ 1 || throw(ArgumentError("ncomponents must be at least 1."))

    mzcount = size(first(X), 2)
    retentioncounts = Vector{Int}(undef, length(X))
    smallestmode = typemax(Int)

    for (sampleindex, Xk) in enumerate(X)
        nretentions, nmz = size(Xk)
        retentioncounts[sampleindex] = nretentions

        if nmz != mzcount
            throw(DimensionMismatch(
                "all sample matrices must have the same number of m/z channels " *
                "(sample 1 has $(mzcount), sample $(sampleindex) has $(nmz))."
            ))
        end

        all(isfinite, Xk) || throw(ArgumentError(
            "X[$(sampleindex)] must contain only finite values."
        ))

        smallestmode = min(smallestmode, nretentions, nmz)
    end

    if ncomponents > smallestmode
        throw(ArgumentError(
            "ncomponents must be between 1 and min(size(X[k])...)=$(smallestmode)."
        ))
    end

    retentioncounts, mzcount
end

"""
    parafac2(X::AbstractVector{<:AbstractMatrix{<:Real}}, ncomponents::Integer;
        retentions=nothing, mzvalues=nothing, samplelabels=nothing, maxiters=100,
        tol=1e-8, nstarts=1, rng=Random.default_rng(), compression=:cholesky,
        nonnegative=(:spectra, :intensities))
    parafac2(X::AbstractArray{<:Real, 3}, ncomponents::Integer;
        retentions=nothing, mzvalues=nothing, samplelabels=nothing, maxiters=100,
        tol=1e-8, nstarts=1, rng=Random.default_rng(), compression=:cholesky,
        nonnegative=(:spectra, :intensities))

Fit a native PARAFAC2 decomposition and return a [`Parafac2Fit`](@ref). The implemented
algorithm is described by Kiers et al. (1999).

For the vector method, each `X[k]` is one `retentions × m/z` sample matrix. All matrices
must have the same number of m/z channels, may have different numbers of retention
positions, and must contain only finite real values. For the tensor method, `X` is
interpreted as `samples × retentions × m/z` and each sample is converted internally to a
`retentions × m/z` matrix view.

Optional `retentions` can be a shared retention coordinate vector, a vector of retention
coordinate vectors (one per sample), or a `samples × retentions` matrix. A shared vector
is only valid when all samples have the same number of retention positions; it is expanded
internally to one copied retention vector per sample. A vector of vectors is the preferred
form when samples have their own retention axes, for example when data have not been binned
onto a common grid. A retention matrix provides one row per sample and is therefore only
rectangular.

The returned [`Parafac2Fit`](@ref) stores retention metadata per sample regardless of input
form. Consequently, `rawretentions(fit)` and `retentions(fit)` return `nothing` or a vector
whose `k`th entry is the retention axis for sample `k`. For a `samples × retentions × m/z`
tensor, a single `retentions` vector should be understood as the common retention grid for
the tensor. If retention axes differ between samples, prefer the vector-of-matrices method
and pass `retentions` as a vector of vectors.

Optional `mzvalues` is a shared m/z coordinate vector. Retention and m/z coordinates may be
unitless real values or `Unitful.AbstractQuantity` values; units are stripped and stored
separately in the returned fit. Optional `samplelabels` must contain one unique label per
sample.

The direct fitting algorithm follows Kiers, Ten Berge, and Bro's alternating scheme:
initialize `loadings` from PCA of `sum(X[k]' * X[k])`, initialize `core = I` and
`weights .= 1`, update sample-specific `bases` by the SVD/Procrustes step, then apply one
PARAFAC1/CP ALS cycle to the small tensor with frontal slices `bases[k]' * X[k]`.
Start 1 is this deterministic rational initialization. Additional starts use random
factors, and the returned fit stores the start with the smallest final objective value.

If `compression=:cholesky`, each sample matrix with more retention rows than m/z columns
and a positive-definite cross-product
is replaced during fitting by a smaller matrix `H[k]` satisfying `H[k]' * H[k] == X[k]' * X[k]`,
following the paper's cross-product sufficiency step. Final sample-specific bases are then
computed from the original matrices, so scores, reconstructions, residuals, and apexes keep
the original retention dimensions. The default `compression=:cholesky` therefore applies
this acceleration only when a sample matrix is tall enough to benefit from it. Use
`compression=:none` to keep the original matrices throughout fitting.

For numerical conditioning, fitting is performed internally after dividing all sample
matrices by their shared global Frobenius norm. The returned fit is transformed back to the
original intensity scale, so weights, reconstructions, losses, areas, and diagnostics use
the same units as the input data.

By default, spectral loadings and sample/component intensity weights are constrained
nonnegative. Pass `nonnegative=()` or `nonnegative=false` to run the paper's original
unconstrained direct ALS update. `nonnegative` may also be `:spectra`, `:intensities`,
or a tuple/vector containing those symbols.

`ncomponents` must satisfy `1 ≤ ncomponents ≤ min(size(X[k])...)` across all sample
matrices. `maxiters` controls the number of ALS cycles and may be zero to return the
rational initialization for `nstarts=1`. `tol` is the relative objective-decrease
threshold used for convergence. `nstarts` must be at least one; pass `rng` to make random
starts reproducible. `compression` must be `:none` or `:cholesky`.

Reference
Kiers HA, ten Berge JMF, Bro R (1999): PARAFAC2—Part I. A direct fitting algorithm for the 
PARAFAC2 model. — Journal of Chemometrics 13_: 275–294.
"""
function parafac2(
    X::AbstractVector{<:AbstractMatrix{<:Real}},
    ncomponents::Integer;
    retentions=nothing,
    mzvalues=nothing,
    samplelabels=nothing,
    maxiters::Integer=100,
    tol::Real=1e-8,
    nstarts::Integer=1,
    rng::Random.AbstractRNG=Random.default_rng(),
    compression::Symbol=:cholesky,
    nonnegative=(:spectra, :intensities)
)
    retentioncounts, mzcount = parafac2dimensions(X, ncomponents)
    maxiters ≥ 0 || throw(ArgumentError("maxiters must be nonnegative."))
    isfinite(tol) && tol ≥ 0 || throw(ArgumentError("tol must be finite and nonnegative."))
    nstarts ≥ 1 || throw(ArgumentError("nstarts must be at least 1."))
    compression = parafac2compression(compression)
    nonnegative = parafac2nonnegative(nonnegative)
    retentionmetadata, retentionunit = parafac2retentionmetadata(retentions, retentioncounts)
    mzmetadata, mzunit = parafac2mzmetadata(mzvalues, mzcount)
    labels = parafac2samplelabelmetadata(samplelabels, length(X))
    T = float(promote_type(map(eltype, X)...))
    Xoriginal = parafac2inputmatrices(X, T)
    inputscale = parafac2inputscale(Xoriginal)
    Xscaled = inputscale == one(T) ? Xoriginal : [Xk ./ inputscale for Xk in Xoriginal]
    Xfit, compressed = parafac2compressionmatrices(Xscaled, compression)
    loadings, core, weights, bases, losses, converged, stopreason, iterations, beststart,
    startdiagnostics = parafac2fit(
            Xfit,
            ncomponents,
            maxiters,
            tol,
            nstarts,
            rng,
            nonnegative
        )
    if any(compressed)
        bases = parafac2updatebases(
            Xscaled,
            loadings,
            core,
            weights
        )
        losses = copy(losses)
        losses[end] = parafac2loss(Xscaled, bases, core, weights, loadings)
        startdiagnostics = copy(startdiagnostics)
        bestdiagnostic = startdiagnostics[beststart]
        startdiagnostics[beststart] = (
            start=bestdiagnostic.start,
            iterations=bestdiagnostic.iterations,
            converged=bestdiagnostic.converged,
            stopreason=bestdiagnostic.stopreason,
            loss=last(losses)
        )
    end

    if inputscale != one(T)
        weights = copy(weights)
        weights .*= inputscale
        losses = copy(losses)
        lossscale = inputscale^2
        losses .*= lossscale
        startdiagnostics = Parafac2StartDiagnostic{T}[
            (
                start=diagnostic.start,
                iterations=diagnostic.iterations,
                converged=diagnostic.converged,
                stopreason=diagnostic.stopreason,
                loss=diagnostic.loss * lossscale
            )
            for diagnostic in startdiagnostics
        ]
    end

    Parafac2Fit{T}(
        Int(ncomponents),
        retentioncounts,
        mzcount,
        retentionmetadata,
        retentionunit,
        mzmetadata,
        mzunit,
        labels,
        loadings,
        core,
        weights,
        bases,
        losses,
        converged,
        stopreason,
        iterations,
        Int(nstarts),
        beststart,
        startdiagnostics,
        compression,
        collect(compressed),
        nonnegative
    )
end

function parafac2(
    X::AbstractArray{<:Real, 3},
    ncomponents::Integer;
    retentions=nothing,
    mzvalues=nothing,
    samplelabels=nothing,
    maxiters::Integer=100,
    tol::Real=1e-8,
    nstarts::Integer=1,
    rng::Random.AbstractRNG=Random.default_rng(),
    compression::Symbol=:cholesky,
    nonnegative=(:spectra, :intensities)
)
    samples = [view(X, k, :, :) for k in axes(X, 1)]
    parafac2(samples, ncomponents;
        retentions=retentions,
        mzvalues=mzvalues,
        samplelabels=samplelabels,
        maxiters=maxiters,
        tol=tol,
        nstarts=nstarts,
        rng=rng,
        compression=compression,
        nonnegative=nonnegative
    )
end
