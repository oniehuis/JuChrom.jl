module TestDeconvolution

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./deconvolution./unimodalfit.jl
# ─────────────────────────────────────────────────────────────────────────────

using Test
using JuChrom
using JuChrom: nnlspenalized

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests
# ─────────────────────────────────────────────────────────────────────────────

# ── nnlspenalized ────────────────────────────────────────────────────────────

@testset "nnlspenalized branches" begin
    F = [1.0 0.0; 0.0 1.0]
    y = [1.0, 2.0]

    a_default = nnlspenalized(F, y)
    @test isapprox(a_default, y; atol=1e-6)

    a_weighted = nnlspenalized(F, y; w=[1.0, 4.0])
    @test isapprox(a_weighted, y; atol=1e-6)

    @test_throws DimensionMismatch nnlspenalized(F, y; w=[1.0])
    @test_throws ArgumentError nnlspenalized(F, y; w=[1.0, -1.0])

    μ = [0.5, 0.5]
    λ = [1.0, 2.0]
    a_prior = nnlspenalized(F, y; μ=μ, λ=λ)
    @test all(a_prior .>= 0.0)

    @test_throws DimensionMismatch nnlspenalized(F, y; μ=[0.1], λ=[1.0, 2.0])
    @test_throws DimensionMismatch nnlspenalized(F, y; μ=[0.1, 0.2], λ=[1.0])

    a_mu_only = nnlspenalized(F, y; μ=μ)
    @test isapprox(a_mu_only, y; atol=1e-6)

    a_lambda_only = nnlspenalized(F, y; λ=λ)
    @test isapprox(a_lambda_only, y; atol=1e-6)
end


# ── unimodalfit ──────────────────────────────────────────────────────────────

@testset "unimodalfit basic outputs" begin
    n_mz = 2
    n_scans = 5
    t_actual = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    Y = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    t0_times = [1.5, 3.0]

    A_hat, basis, C, (tgrid, Fhat_grid), Yfit = unimodalfit(
        Y,
        t_actual,
        t0_times;
        K=6,
        tgrid_n=20,
        iters=1,
        ridge=1e-6,
        nnls_ridge=1e-8
    )

    @test size(A_hat) == (n_mz, length(t0_times))
    @test size(C, 2) == length(t0_times)
    @test length(tgrid) == 20
    @test size(Fhat_grid) == (20, length(t0_times))
    @test size(Yfit) == size(Y)
    @test all(isfinite.(A_hat))
    @test all(isfinite.(Yfit))
end

@testset "unimodalfit sigma and masks" begin
    n_mz = 2
    n_scans = 5
    t_actual = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    Y = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    t0_times = [1.5, 3.0]
    σ = fill(1.0, n_mz, n_scans)
    A_init = fill(0.1, n_mz, length(t0_times))
    lock_Azeros = falses(n_mz, length(t0_times))
    lock_Azeros[1, 2] = true

    A_hat, basis, C, (tgrid, Fhat_grid), Yfit = unimodalfit(
        Y,
        t_actual,
        t0_times;
        σ=σ,
        A_init=A_init,
        lock_Azeros=lock_Azeros,
        K=6,
        tgrid_n=20,
        iters=1
    )

    @test A_hat[1, 2] == 0.0
    @test all(A_hat .>= 0.0)
    @test all(isfinite.(Yfit))
end

@testset "unimodalfit priors" begin
    n_mz = 2
    n_scans = 5
    t_actual = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    Y = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    t0_times = [1.5, 3.0]
    A_prior = fill(0.2, n_mz, length(t0_times))
    lambda_peaks = [0.1, 0.2]

    A_hat, basis, C, (tgrid, Fhat_grid), Yfit = unimodalfit(
        Y,
        t_actual,
        t0_times;
        A_prior=A_prior,
        lambda_peaks=lambda_peaks,
        K=6,
        tgrid_n=20,
        iters=1
    )

    @test all(A_hat .>= 0.0)
    @test all(isfinite.(A_hat))
    @test all(isfinite.(Yfit))
end

@testset "unimodalfit input validation" begin
    n_mz = 2
    n_scans = 5
    t_actual = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    Y = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    t0_times = [1.5, 3.0]

    t_bad = zeros(n_mz, n_scans + 1)
    @test_throws DimensionMismatch unimodalfit(Y, t_bad, t0_times)

    σ_bad = [1.0 1.0 1.0 1.0 0.0;
             1.0 1.0 1.0 1.0 1.0]
    @test_throws ArgumentError unimodalfit(Y, t_actual, t0_times; σ=σ_bad)

    @test_throws DimensionMismatch unimodalfit(Y, t_actual, t0_times; A_init=ones(n_mz + 1, length(t0_times)))
    @test_throws DimensionMismatch unimodalfit(Y, t_actual, t0_times; lock_Azeros=falses(n_mz, length(t0_times) + 1))
    @test_throws DimensionMismatch unimodalfit(Y, t_actual, t0_times; A_prior=ones(n_mz + 1, length(t0_times)))
    @test_throws ArgumentError unimodalfit(Y, t_actual, t0_times; lambda_peaks=[0.1])
    @test_throws ArgumentError unimodalfit(Y, t_actual, t0_times; lambda_peaks=[0.1, -0.2])
    @test_throws ArgumentError unimodalfit(Y, t_actual, t0_times; A_prior=ones(n_mz, length(t0_times)))
end


# ── unimodalfit_t0 ───────────────────────────────────────────────────────────

@testset "unimodalfit_t0 basic outputs" begin
    n_mz = 2
    n_scans = 5
    t_actual = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    Y = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    t0_guess = [1.4, 3.1]

    t0_best, sse, A_hat, basis, C, (tgrid, Fhat_grid), Yfit = unimodalfit_t0(
        Y,
        t_actual;
        t0_guess=t0_guess,
        strategy=:single,
        half_width=0.5,
        ngrid=5,
        coord_sweeps=1,
        max_iter=5,
        K=6,
        tgrid_n=20,
        iters=1,
        ridge=1e-6,
        nnls_ridge=1e-8,
        verbose=false
    )

    @test length(t0_best) == length(t0_guess)
    @test isfinite(sse)
    @test size(A_hat) == (n_mz, length(t0_guess))
    @test size(C, 2) == length(t0_guess)
    @test length(tgrid) == 20
    @test size(Fhat_grid) == (20, length(t0_guess))
    @test size(Yfit) == size(Y)
end

@testset "unimodalfit_t0 iterative stages" begin
    n_mz = 2
    n_scans = 5
    t_actual = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    Y = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    t0_guess = [1.4, 3.1]
    stages = [(half_width=0.6, ngrid=5), (half_width=0.3, ngrid=7)]

    t0_best, sse, A_hat, basis, C, (tgrid, Fhat_grid), Yfit = unimodalfit_t0(
        Y,
        t_actual;
        t0_guess=t0_guess,
        strategy=:iterative,
        stages=stages,
        coord_sweeps=1,
        max_iter=5,
        K=6,
        tgrid_n=20,
        iters=1,
        verbose=false
    )

    @test length(t0_best) == length(t0_guess)
    @test isfinite(sse)
    @test size(Yfit) == size(Y)
end

@testset "unimodalfit_t0 input validation" begin
    n_mz = 2
    n_scans = 5
    t_actual = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    Y = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]

    @test_throws ArgumentError unimodalfit_t0(Y, t_actual; t0_guess=Float64[])
    @test_throws ArgumentError unimodalfit_t0(Y, t_actual; t0_guess=[1.0], half_width=0.0)
    @test_throws ArgumentError unimodalfit_t0(Y, t_actual; t0_guess=[1.0], ngrid=2)
    @test_throws ArgumentError unimodalfit_t0(Y, t_actual; t0_guess=[1.0], coord_sweeps=0)
    @test_throws ArgumentError unimodalfit_t0(Y, t_actual; t0_guess=[1.0], max_iter=0)
    @test_throws ArgumentError unimodalfit_t0(Y, t_actual; t0_guess=[1.0], tol_sse=0.0)
    @test_throws ArgumentError unimodalfit_t0(Y, t_actual; t0_guess=[1.0, 1.05], min_sep=0.2)
end

end  # module TestDeconvolution
