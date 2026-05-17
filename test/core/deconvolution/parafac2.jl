module TestParafac2

using Test
using JuChrom
using LinearAlgebra
using Random
using Unitful

struct Parafac2LossProbeMatrix <: AbstractMatrix{Float64}
    data::Matrix{Float64}
end

Base.size(matrix::Parafac2LossProbeMatrix) = size(matrix.data)
Base.getindex(matrix::Parafac2LossProbeMatrix, row::Int, col::Int) =
    matrix.data[row, col]
Base.IndexStyle(::Type{<:Parafac2LossProbeMatrix}) = IndexCartesian()

const PARAFAC2_LOSS_PROBE_CALLS = Ref(0)

function JuChrom.parafac2loss(
    X::Vector{Parafac2LossProbeMatrix},
    bases::AbstractVector{<:AbstractMatrix{Float64}},
    core::AbstractMatrix{Float64},
    weights::AbstractMatrix{Float64},
    loadings::AbstractMatrix{Float64}
)
    PARAFAC2_LOSS_PROBE_CALLS[] += 1
    PARAFAC2_LOSS_PROBE_CALLS[] == 1 ? 1.0 : Inf
end

@testset "parafac2 vector input validation" begin
    X = [
        reshape(collect(1.0:15.0), 5, 3),
        reshape(collect(1.0:21.0), 7, 3)
    ]

    fit = parafac2(X, 2; maxiters=0)

    @test fit isa Parafac2Fit{Float64}
    @test fit.ncomponents == 2
    @test fit.retentioncounts == [5, 7]
    @test fit.mzcount == 3
    @test fit.retentions === nothing
    @test fit.retentionunit === nothing
    @test fit.mzvalues === nothing
    @test fit.mzunit === nothing
    @test fit.samplelabels === nothing
    @test size(fit.loadings) == (3, 2)
    @test size(fit.core) == (2, 2)
    @test size(fit.weights) == (2, 2)
    @test length(fit.bases) == 2
    @test size(fit.bases[1]) == (5, 2)
    @test size(fit.bases[2]) == (7, 2)
    @test length(fit.loss) == 1
    @test isfinite(only(fit.loss))
    @test fit.converged == false
    @test fit.stopreason == :maxiters
    @test fit.iterations == 0
    @test fit.nstarts == 1
    @test fit.beststart == 1
    @test fit.compression == :none
    @test fit.nonnegative == (:spectra, :abundances)
    @test all(fit.loadings .>= 0)
    @test all(fit.weights .>= 0)
    @test fit.bases[1]' * fit.bases[1] ≈ I(2) atol=1e-6
    @test fit.bases[2]' * fit.bases[2] ≈ I(2) atol=1e-6
end

@testset "parafac2 tensor input validation" begin
    X = reshape(collect(Float32, 1:60), 3, 5, 4)
    fit = parafac2(X, 3; maxiters=0)

    @test fit isa Parafac2Fit{Float32}
    @test fit.ncomponents == 3
    @test fit.retentioncounts == [5, 5, 5]
    @test fit.mzcount == 4
    @test fit.nstarts == 1
    @test fit.beststart == 1
    @test fit.compression == :none
    @test fit.nonnegative == (:spectra, :abundances)
end

@testset "parafac2 unitless metadata" begin
    X = [
        reshape(collect(1.0:12.0), 4, 3),
        reshape(collect(1.0:12.0), 4, 3)
    ]
    fit = parafac2(
        X,
        2;
        retentions=[100.0, 110.0, 120.0, 130.0],
        mzvalues=[50.0, 75.0, 100.0],
        samplelabels=[:blank, :sample],
        maxiters=0
    )

    @test rawretentions(fit) == [[100.0, 110.0, 120.0, 130.0],
                                 [100.0, 110.0, 120.0, 130.0]]
    @test retentions(fit) == rawretentions(fit)
    @test retentionunit(fit) === nothing
    @test rawmzvalues(fit) == [50.0, 75.0, 100.0]
    @test mzvalues(fit) == rawmzvalues(fit)
    @test mzunit(fit) === nothing
    @test fit.samplelabels == ["blank", "sample"]
end

@testset "parafac2 unitful metadata" begin
    X = [
        reshape(collect(1.0:12.0), 4, 3),
        reshape(collect(1.0:15.0), 5, 3)
    ]
    retention_axes = [
        [1.0, 2.0, 3.0, 4.0]u"s",
        [0.1, 0.2, 0.3, 0.4, 0.5]u"minute"
    ]
    fit = parafac2(
        X,
        2;
        retentions=retention_axes,
        mzvalues=[100.0, 125.0, 150.0]u"Th",
        samplelabels=["a", "b"],
        maxiters=0
    )

    @test retentionunit(fit) == u"s"
    @test rawretentions(fit)[1] ≈ [1.0, 2.0, 3.0, 4.0]
    @test rawretentions(fit)[2] ≈ [6.0, 12.0, 18.0, 24.0, 30.0]
    @test rawretentions(fit; unit=u"minute")[1] ≈ [1.0, 2.0, 3.0, 4.0] ./ 60.0
    @test retentions(fit)[1] == [1.0, 2.0, 3.0, 4.0]u"s"
    @test retentions(fit; unit=u"minute")[2] == [0.1, 0.2, 0.3, 0.4, 0.5]u"minute"

    @test mzunit(fit) == u"Th"
    @test rawmzvalues(fit) ≈ [100.0, 125.0, 150.0]
    @test rawmzvalues(fit; unit=u"kTh") ≈ [0.1, 0.125, 0.15]
    @test mzvalues(fit) == [100.0, 125.0, 150.0]u"Th"
end

@testset "parafac2 tensor metadata" begin
    X = reshape(collect(Float64, 1:24), 2, 4, 3)
    retention_matrix = [1.0 2.0 3.0 4.0;
                        1.5 2.5 3.5 4.5]

    fit = parafac2(
        X,
        2;
        retentions=retention_matrix,
        mzvalues=[40.0, 41.0, 42.0],
        samplelabels=["s1", "s2"],
        maxiters=0
    )

    @test fit.retentioncounts == [4, 4]
    @test rawretentions(fit) == [[1.0, 2.0, 3.0, 4.0],
                                 [1.5, 2.5, 3.5, 4.5]]
    @test rawmzvalues(fit) == [40.0, 41.0, 42.0]
    @test fit.samplelabels == ["s1", "s2"]
end

@testset "parafac2 input validation errors" begin
    @test_throws ArgumentError parafac2(Matrix{Float64}[], 1)
    @test_throws ArgumentError parafac2([ones(3, 2)], 0)
    @test_throws ArgumentError parafac2([ones(3, 2)], 3)
    @test_throws ArgumentError parafac2([ones(0, 2)], 1)
    @test_throws ArgumentError parafac2(zeros(0, 3, 2), 1)

    @test_throws DimensionMismatch parafac2([ones(3, 2), ones(4, 3)], 1)

    Xnan = [1.0 NaN; 2.0 3.0]
    Xinf = [1.0 Inf; 2.0 3.0]
    @test_throws ArgumentError parafac2([Xnan], 1)
    @test_throws ArgumentError parafac2([Xinf], 1)
    @test_throws ArgumentError parafac2(reshape([1.0, NaN, 2.0, 3.0], 2, 2, 1), 1)
end

@testset "parafac2 unconstrained ALS exact data" begin
    ncomponents = 2
    mzcount = 4
    loadings_true = Matrix(qr([
        1.0 0.2
        0.1 1.0
        0.5 0.3
        0.2 0.4
    ]).Q[:, 1:ncomponents])
    core_true = [1.0 0.2; -0.1 0.8]
    weights_true = [1.0 2.0;
                    2.0 1.0;
                    1.5 0.7]

    bases_true = [
        Matrix(qr([1.0 0.0; 0.0 1.0; 0.2 0.3; 0.4 -0.1; -0.2 0.5]).Q[:, 1:ncomponents]),
        Matrix(qr([0.3 0.1; 1.0 0.2; -0.1 0.8; 0.4 0.0; 0.2 -0.4; -0.3 0.6]).Q[:, 1:ncomponents]),
        Matrix(qr([0.4 0.3; -0.2 1.0; 0.1 -0.2; 1.0 0.0; 0.3 0.2; -0.1 0.5; 0.2 -0.3]).Q[:, 1:ncomponents])
    ]

    X = [
        bases_true[k] *
        core_true *
        Diagonal(weights_true[k, :]) *
        loadings_true'
        for k in eachindex(bases_true)
    ]

    fit = parafac2(X, ncomponents; maxiters=100, tol=1e-10, nonnegative=())
    reconstructed = parafac2reconstruct(fit)
    residual = sum(sum(abs2, X[k] .- reconstructed[k]) for k in eachindex(X))
    total = sum(sum(abs2, Xk) for Xk in X)

    @test fit.iterations <= 100
    @test fit.nonnegative == ()
    @test length(fit.loss) == fit.iterations + 1
    @test all(isfinite, fit.loss)
    @test fit.loss[end] <= fit.loss[1] + 1e-8
    @test all(diff(fit.loss) .<= 1e-6)
    @test residual / total < 1e-6
    @test fit.bases[1]' * fit.bases[1] ≈ I(ncomponents) atol=1e-7
    @test fit.bases[2]' * fit.bases[2] ≈ I(ncomponents) atol=1e-7
    @test fit.bases[3]' * fit.bases[3] ≈ I(ncomponents) atol=1e-7
    @test parafac2loss(fit, X) ≈ fit.loss[end] rtol=1e-10 atol=1e-10
    @test parafac2fitpercent(fit, X) > 1 - 1e-6
end

@testset "parafac2 rank-one exact convergence" begin
    ncomponents = 1
    loadings_true = reshape([1.0, 0.5, 0.2, 0.1] ./ norm([1.0, 0.5, 0.2, 0.1]), :, 1)
    core_true = reshape([1.3], 1, 1)
    weights_true = reshape([1.0, 2.0, 0.7], 3, 1)
    bases_true = [
        reshape([1.0, 0.2, 0.3, 0.4] ./ norm([1.0, 0.2, 0.3, 0.4]), :, 1),
        reshape([0.1, 1.0, 0.2, 0.3, 0.5] ./ norm([0.1, 1.0, 0.2, 0.3, 0.5]), :, 1),
        reshape([0.4, 0.1, 0.6, 0.2, 0.3, 0.7] ./ norm([0.4, 0.1, 0.6, 0.2, 0.3, 0.7]), :, 1)
    ]
    X = [
        bases_true[k] *
        core_true *
        Diagonal(weights_true[k, :]) *
        loadings_true'
        for k in eachindex(bases_true)
    ]

    fit = parafac2(X, ncomponents; maxiters=20, tol=1e-12)
    reconstructed = parafac2reconstruct(fit)
    residual = sum(sum(abs2, X[k] .- reconstructed[k]) for k in eachindex(X))
    total = sum(sum(abs2, Xk) for Xk in X)

    @test fit.converged
    @test fit.stopreason == :tol
    @test residual / total < 1e-20
end

@testset "parafac2 nonnegative spectra and abundances" begin
    X = [
        [1.0 0.2 0.4 1.3; 0.7 1.1 0.3 0.5; 0.2 0.5 1.3 0.7; 1.4 0.8 0.6 0.3],
        [0.8 1.0 0.2 0.4; 1.3 0.4 0.5 1.1; 0.6 0.9 1.4 0.2; 1.5 0.7 0.8 0.6],
        [0.9 0.3 0.7 1.2; 0.4 1.5 0.9 0.5; 0.2 0.8 1.7 0.9; 1.2 0.3 1.0 0.4]
    ]

    fit = parafac2(X, 2; maxiters=5, tol=0.0)
    @test fit.nonnegative == (:spectra, :abundances)
    @test all(parafac2spectra(fit) .>= 0)
    @test all(parafac2abundances(fit) .>= 0)

    fitbool = parafac2(X, 2; maxiters=0, nonnegative=true)
    @test fitbool.nonnegative == (:spectra, :abundances)

    fitspectra = parafac2(X, 2; maxiters=2, tol=0.0, nonnegative=:spectra)
    @test fitspectra.nonnegative == (:spectra,)
    @test all(parafac2spectra(fitspectra) .>= 0)

    fitweights = parafac2(X, 2; maxiters=2, tol=0.0, nonnegative=:weights)
    @test fitweights.nonnegative == (:abundances,)
    @test all(parafac2abundances(fitweights) .>= 0)

    fitunconstrained = parafac2(X, 2; maxiters=0, nonnegative=false)
    @test fitunconstrained.nonnegative == ()
end

@testset "parafac2 cross-product compression" begin
    X = [
        [1.0 0.2 0.4; 0.7 1.1 0.3; 0.2 0.5 1.3; 1.4 0.8 0.6;
         0.9 1.6 0.5; 0.3 0.7 1.8; 1.1 0.4 0.9; 0.5 1.2 1.0],
        [0.8 1.0 0.2; 1.3 0.4 0.5; 0.6 0.9 1.4; 1.5 0.7 0.8;
         0.4 1.5 0.9; 0.2 0.8 1.7; 1.2 0.3 1.0]
    ]

    fitnone = parafac2(X, 2; maxiters=0, compression=:none)
    fitchol = parafac2(X, 2; maxiters=0, compression=:cholesky)

    @test fitnone.compression == :none
    @test fitchol.compression == :cholesky
    @test fitchol.retentioncounts == [8, 7]
    @test size(fitchol.bases[1]) == (8, 2)
    @test size(fitchol.bases[2]) == (7, 2)
    @test size(parafac2scores(fitchol, 1)) == (8, 2)
    @test size(parafac2scores(fitchol, 2)) == (7, 2)
    @test abs.(fitchol.loadings) ≈ abs.(fitnone.loadings)
    @test fitchol.core ≈ fitnone.core
    @test fitchol.weights ≈ fitnone.weights
    @test parafac2loss(fitchol, X) ≈ last(fitchol.loss) rtol=1e-10 atol=1e-10

    fititer = parafac2(X, 2; maxiters=5, tol=0.0, compression=:cholesky)
    @test fititer.compression == :cholesky
    @test fititer.bases[1]' * fititer.bases[1] ≈ I(2) atol=1e-6
    @test fititer.bases[2]' * fititer.bases[2] ≈ I(2) atol=1e-6
    @test parafac2loss(fititer, X) ≈ last(fititer.loss) rtol=1e-8 atol=1e-8

    tensor = Array{Float64}(undef, 2, 8, 3)
    tensor[1, :, :] .= X[1]
    tensor[2, 1:7, :] .= X[2]
    tensor[2, 8, :] .= [1.0, 1.1, 1.2]
    tensorfit = parafac2(tensor, 2; maxiters=1, tol=0.0, compression=:cholesky)
    @test tensorfit.compression == :cholesky
    @test tensorfit.retentioncounts == [8, 8]
    @test size(tensorfit.bases[1]) == (8, 2)
    @test size(tensorfit.bases[2]) == (8, 2)

    Xwide = [
        [1.0 0.2 0.4; 0.7 1.1 0.3],
        [0.8 1.0 0.2; 1.3 0.4 0.5]
    ]
    fitwide = parafac2(Xwide, 2; maxiters=0, compression=:cholesky)
    fitwidenone = parafac2(Xwide, 2; maxiters=0, compression=:none)
    @test fitwide.compression == :cholesky
    @test fitwide.retentioncounts == [2, 2]
    @test size(fitwide.bases[1]) == (2, 2)
    @test fitwide.loss ≈ fitwidenone.loss
end

@testset "parafac2 least-squares helpers" begin
    design_prune = [-3.0 -3.0; -3.0 -2.0]
    response_prune = [-3.0, -2.0]
    @test JuChrom.nonnegativeleastsquares(design_prune, response_prune) ≈ [0.0, 1.0]

    design_zero_alpha = [
        -1.8024720054589953 -1.0100699865483473
         0.9500573649054608 -0.9140623330009525
        -0.4653449227441103 -1.0488926824878395
        -2.089590332600559   0.057505351867115514
        -0.948734830619434   0.058121993966906614
    ]
    response_zero_alpha = [
        -1.256565821949869,
         0.5431532750982856,
        -1.6748079068851502,
         1.1427848002128,
         0.06600532091593743
    ]
    solution = JuChrom.nonnegativeleastsquares(design_zero_alpha, response_zero_alpha)
    @test solution ≈ [9.880553004028577e-8, 0.8772628742396302] rtol=1e-10
    @test all(value -> value >= 0, solution)

    design = Matrix{Float64}(I, 2, 2)
    response = [-1.0, 2.0]
    @test JuChrom.leastsquaressolution(design, response, false) == response
    @test JuChrom.leastsquaressolution(design, response, true) == [0.0, 2.0]
end

@testset "parafac2 multiple starts" begin
    X = [
        [1.0 0.2 1.3 0.7; 0.8 0.4 1.0 0.6; 0.5 0.7 0.8 0.4; 0.3 0.9 0.5 0.2;
         0.2 1.1 0.3 0.1],
        [0.9 1.1 0.2 0.6; 1.0 0.8 0.4 0.7; 1.2 0.5 0.6 0.8; 1.4 0.3 0.7 0.9;
         1.1 0.2 0.8 1.0; 0.8 0.1 0.9 1.1],
        [0.4 0.8 1.1 0.3; 0.5 0.9 1.0 0.4; 0.7 1.0 0.8 0.5; 0.9 1.1 0.6 0.6;
         1.0 1.0 0.4 0.7]
    ]

    single = parafac2(X, 2; maxiters=6, tol=0.0, nstarts=1)
    multi = parafac2(X, 2; maxiters=6, tol=0.0, nstarts=4,
        rng=Random.MersenneTwister(123))

    @test single.nstarts == 1
    @test single.beststart == 1
    @test multi.nstarts == 4
    @test 1 <= multi.beststart <= 4
    @test length(multi.loss) == multi.iterations + 1
    @test last(multi.loss) <= last(single.loss) + 1e-8
    @test parafac2loss(multi, X) ≈ last(multi.loss) rtol=1e-10 atol=1e-10

    fit1 = parafac2(X, 2; maxiters=5, tol=0.0, nstarts=3,
        rng=Random.MersenneTwister(456))
    fit2 = parafac2(X, 2; maxiters=5, tol=0.0, nstarts=3,
        rng=Random.MersenneTwister(456))
    @test fit1.beststart == fit2.beststart
    @test fit1.loss ≈ fit2.loss
    @test fit1.loadings ≈ fit2.loadings
    @test fit1.core ≈ fit2.core
    @test fit1.weights ≈ fit2.weights

    tensor = reshape(collect(1.0:60.0), 3, 5, 4)
    tensorfit = parafac2(tensor, 2; maxiters=2, tol=0.0, nstarts=2,
        rng=Random.MersenneTwister(789))
    @test tensorfit.nstarts == 2
    @test 1 <= tensorfit.beststart <= 2
end

@testset "parafac2 public reconstruction helpers" begin
    ncomponents = 2
    X = [
        [1.0 2.0 0.5; 1.5 2.5 0.7; 2.0 2.1 0.9; 1.8 1.7 0.6],
        [0.8 1.2 1.9; 1.0 1.5 1.7; 1.3 1.7 1.4; 1.6 1.8 1.2; 1.9 2.0 1.0]
    ]
    fit = parafac2(X, ncomponents; maxiters=5, tol=0.0)

    scores = parafac2scores(fit)
    reconstructed = parafac2reconstruct(fit)

    @test length(scores) == length(X)
    @test size(scores[1]) == (4, ncomponents)
    @test size(scores[2]) == (5, ncomponents)
    @test parafac2scores(fit, 1) == scores[1]
    @test length(reconstructed) == length(X)
    @test size(reconstructed[1]) == size(X[1])
    @test size(reconstructed[2]) == size(X[2])
    @test parafac2reconstruct(fit, 2) == reconstructed[2]
    @test_throws BoundsError parafac2scores(fit, 0)
    @test_throws BoundsError parafac2reconstruct(fit, 3)

    unfitted = Parafac2Fit{Float64}(
        1,
        [2],
        2,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        Float64[],
        false,
        :maxiters,
        0,
        1,
        1,
        :none,
        ()
    )
    @test_throws ArgumentError parafac2scores(unfitted)
    @test_throws ArgumentError parafac2reconstruct(unfitted)
    @test_throws ArgumentError parafac2spectra(unfitted)
end

@testset "parafac2 component summaries" begin
    X = [
        [1.0 2.0 0.5; 1.5 2.5 0.7; 2.0 2.1 0.9; 1.8 1.7 0.6],
        [0.8 1.2 1.9; 1.0 1.5 1.7; 1.3 1.7 1.4; 1.6 1.8 1.2; 1.9 2.0 1.0]
    ]
    fit = parafac2(
        X,
        2;
        retentions=[[10.0, 20.0, 30.0, 40.0],
                    [11.0, 22.0, 33.0, 44.0, 55.0]],
        mzvalues=[50.0, 75.0, 100.0],
        maxiters=5,
        tol=0.0
    )

    spectra = parafac2spectra(fit)
    spectrametadata = parafac2spectra(fit; metadata=true)
    abundances = parafac2abundances(fit)

    @test spectra == fit.loadings
    @test spectra !== fit.loadings
    @test spectrametadata.values == fit.loadings
    @test spectrametadata.values !== fit.loadings
    @test spectrametadata.mzvalues == [50.0, 75.0, 100.0]
    @test abundances == fit.weights
    @test abundances !== fit.weights

    apexes = parafac2apexes(fit)
    scores = parafac2scores(fit)
    storedretentions = retentions(fit)

    @test size(apexes.indices) == (2, 2)
    @test size(apexes.retentions) == (2, 2)
    @test size(apexes.values) == (2, 2)

    for sampleindex in eachindex(scores), component in axes(scores[sampleindex], 2)
        profile = view(scores[sampleindex], :, component)
        apexindex = argmax(profile)
        @test apexes.indices[sampleindex, component] == apexindex
        @test apexes.values[sampleindex, component] == profile[apexindex]
        @test apexes.retentions[sampleindex, component] == storedretentions[sampleindex][apexindex]
    end

    fitwithoutmetadata = parafac2(X, 2; maxiters=0)
    @test rawretentions(fitwithoutmetadata) === nothing
    @test rawmzvalues(fitwithoutmetadata) === nothing
    @test parafac2spectra(fitwithoutmetadata; metadata=true).mzvalues === nothing
    @test parafac2apexes(fitwithoutmetadata).retentions === nothing

    Xunit = [
        reshape(collect(1.0:6.0), 3, 2),
        reshape(collect(1.0:8.0), 4, 2)
    ]
    fitunit = parafac2(
        Xunit,
        1;
        retentions=[[1.0, 2.0, 3.0]u"s",
                    [0.1, 0.2, 0.3, 0.4]u"minute"],
        mzvalues=[100.0, 200.0]u"Th",
        maxiters=0
    )
    unitmetadata = parafac2spectra(fitunit; metadata=true, unit=u"kTh")
    unitapexes = parafac2apexes(fitunit; unit=u"minute")
    unitretentions = retentions(fitunit; unit=u"minute")

    @test unitmetadata.mzvalues ≈ [0.1, 0.2]u"kTh"
    for sampleindex in 1:2
        @test unitapexes.retentions[sampleindex, 1] ==
            unitretentions[sampleindex][unitapexes.indices[sampleindex, 1]]
    end

    @test_throws ArgumentError parafac2spectra(fit; metadata=true, unit=u"Th")
    @test_throws ArgumentError parafac2spectra(fit; unit=u"Th")
    @test_throws ArgumentError parafac2apexes(fit; unit=u"s")
end

@testset "parafac2 public diagnostics" begin
    X = [
        [1.0 0.2 0.3; 0.8 0.4 0.5; 0.6 0.7 0.9; 0.4 1.0 1.1],
        [0.3 1.2 0.4; 0.5 1.0 0.7; 0.8 0.8 1.0; 1.0 0.5 1.2]
    ]
    fit = parafac2(X, 2; maxiters=10, tol=0.0)
    residuals = parafac2residuals(fit, X)
    reconstructed = parafac2reconstruct(fit)
    tensor = Array{Float64}(undef, 2, 4, 3)
    tensor[1, :, :] .= X[1]
    tensor[2, :, :] .= X[2]

    @test length(residuals) == 2
    @test residuals[1] ≈ X[1] .- reconstructed[1]
    @test residuals[2] ≈ X[2] .- reconstructed[2]
    @test parafac2loss(fit, X) ≈ sum(sum(abs2, residual) for residual in residuals)
    @test parafac2loss(fit, X) ≈ fit.loss[end] rtol=1e-10 atol=1e-10
    @test parafac2loss(fit, tensor) ≈ parafac2loss(fit, X)
    @test parafac2residuals(fit, tensor)[1] ≈ residuals[1]
    @test 0 <= parafac2fitpercent(fit, X) <= 1
    @test parafac2fitpercent(fit, tensor) ≈ parafac2fitpercent(fit, X)
    profileminima = parafac2profileminima(fit)
    profilediagnostics = parafac2profilediagnostics(fit)
    scores = parafac2scores(fit)
    @test size(profileminima) == (2, 2)
    @test profilediagnostics.minima == profileminima
    @test profilediagnostics.minimum == minimum(profileminima)
    @test profilediagnostics.hasnegative == any(<(0), profileminima)
    @test size(profilediagnostics.negativecounts) == (2, 2)
    for sampleindex in eachindex(scores), component in axes(scores[sampleindex], 2)
        profile = view(scores[sampleindex], :, component)
        @test profileminima[sampleindex, component] == minimum(profile)
        @test profilediagnostics.negativecounts[sampleindex, component] == count(<(0), profile)
    end

    @test_throws DimensionMismatch parafac2residuals(fit, [X[1]])
    @test_throws DimensionMismatch parafac2loss(fit, [ones(3, 3), X[2]])
    @test_throws DimensionMismatch parafac2fitpercent(fit, zeros(3, 4, 3))
    Xbad = deepcopy(X)
    Xbad[1][1, 1] = NaN
    @test_throws ArgumentError parafac2residuals(fit, Xbad)

    zero_fit = parafac2([zeros(3, 2)], 1; maxiters=0)
    @test isnan(parafac2fitpercent(zero_fit, [zeros(3, 2)]))
end

@testset "parafac2 ALS keyword validation" begin
    X = [ones(4, 3), ones(5, 3)]

    @test_throws ArgumentError parafac2(X, 2; maxiters=-1)
    @test_throws ArgumentError parafac2(X, 2; tol=-1.0)
    @test_throws ArgumentError parafac2(X, 2; tol=NaN)
    @test_throws ArgumentError parafac2(X, 2; nstarts=0)
    @test_throws ArgumentError parafac2(X, 2; compression=:svd)
    @test_throws ArgumentError parafac2(X, 2; nonnegative=:profiles)
    @test_throws ArgumentError parafac2(X, 2; nonnegative=:scores)
    @test_throws ArgumentError parafac2(X, 2; nonnegative=(:spectra, 1))
end

@testset "parafac2 metadata validation errors" begin
    X = [ones(4, 3), ones(5, 3)]

    @test_throws DimensionMismatch parafac2(X, 2; retentions=[1.0, 2.0, 3.0])
    @test_throws DimensionMismatch parafac2(X, 2; retentions=[[1.0, 2.0, 3.0, 4.0]])
    @test_throws DimensionMismatch parafac2(X, 2; retentions=[1.0 2.0 3.0 4.0])
    @test_throws ArgumentError parafac2(X, 2; retentions=[[1.0, 3.0, 2.0, 4.0],
                                                          [1.0, 2.0, 3.0, 4.0, 5.0]])
    @test_throws ArgumentError parafac2(X, 2; retentions=[[1.0, 2.0, 3.0, 4.0],
                                                          [1.0u"s", 2.0u"s", 3.0u"s", 4.0u"s", 5.0u"s"]])

    Xsame = [ones(4, 3), ones(4, 3)]
    mixed_retention_matrix = Any[
        1.0       2.0       3.0       4.0
        1.0u"s"   2.0u"s"   3.0u"s"   4.0u"s"
    ]
    @test_throws ArgumentError parafac2(Xsame, 2; retentions=mixed_retention_matrix)

    @test_throws DimensionMismatch parafac2(X, 2; mzvalues=[1.0, 2.0])
    @test_throws ArgumentError parafac2(X, 2; mzvalues=[1.0, 1.0, 2.0])
    @test_throws ArgumentError parafac2(X, 2; mzvalues=[0.0, 1.0, 2.0])
    @test_throws ArgumentError parafac2(X, 2; mzvalues=[1.0, NaN, 2.0])
    @test_throws ArgumentError parafac2(X, 2; mzvalues=Any[1.0u"Th", 2.0u"s", 3.0u"Th"])

    @test_throws DimensionMismatch parafac2(X, 2; samplelabels=["a"])
    @test_throws ArgumentError parafac2(X, 2; samplelabels=["a", "a"])
    @test_throws ArgumentError parafac2(X, 2; samplelabels=["a", ""])

    fit = parafac2(X, 2; retentions=[[1.0, 2.0, 3.0, 4.0],
                                     [1.0, 2.0, 3.0, 4.0, 5.0]],
                         mzvalues=[1.0, 2.0, 3.0],
                         maxiters=0)
    @test_throws ArgumentError rawretentions(fit; unit=u"s")
    @test_throws ArgumentError rawmzvalues(fit; unit=u"Th")
end

@testset "parafac2 ALS loss safeguard" begin
    X = [
        Parafac2LossProbeMatrix([1.0 0.2; 0.3 0.4]),
        Parafac2LossProbeMatrix([0.5 0.6; 0.7 0.8; 0.9 1.0])
    ]
    loadings = reshape([1.0, 0.0], :, 1)
    core = ones(1, 1)
    weights = ones(2, 1)

    PARAFAC2_LOSS_PROBE_CALLS[] = 0
    fitstart = JuChrom.parafac2fitstart(X, loadings, core, weights, 1, 0.0, ())

    @test PARAFAC2_LOSS_PROBE_CALLS[] == 2
    @test fitstart.stopreason == :nonfinite
    @test fitstart.iterations == 0
    @test fitstart.losses == [1.0]
    @test fitstart.loadings === loadings
    @test fitstart.core === core
    @test fitstart.weights === weights
end

@testset "parafac2 display and broadcasting" begin
    fit = parafac2([ones(3, 2), ones(4, 2)], 1; samplelabels=["a", "b"], maxiters=0)

    @test sprint(show, fit) == summary(fit)
    pretty = sprint(io -> show(io, MIME"text/plain"(), fit))
    @test occursin("Parafac2Fit", pretty)
    @test occursin("fitted", pretty)
    @test occursin("sample labels", pretty)
    @test occursin("stop reason", pretty)
    @test identity.(fit) === fit
    @test (fit .=== fit) === true
end

end
