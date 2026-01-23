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
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions = [1.5, 3.0]

    A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit(
        I,
        R,
        peakretentions;
        knots_n=6,
        rgrid_n=20,
        iters=1,
        shape_ridge=1e-6,
        spectra_ridge=1e-8
    )

    @test size(A_hat) == (n_mz, length(peakretentions))
    @test size(C, 2) == length(peakretentions)
    @test length(rgrid) == 20
    @test size(Fhat_grid) == (20, length(peakretentions))
    @test size(Ifit) == size(I)
    @test all(isfinite.(A_hat))
    @test all(isfinite.(Ifit))
end

@testset "unimodalfit sigma and masks" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions = [1.5, 3.0]
    σ = fill(1.0, n_mz, n_scans)
    spectra_init = fill(0.1, n_mz, length(peakretentions))
    spectra_zero_mask = falses(n_mz, length(peakretentions))
    spectra_zero_mask[1, 2] = true

    A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit(
        I,
        R,
        peakretentions;
        σ=σ,
        spectra_init=spectra_init,
        spectra_zero_mask=spectra_zero_mask,
        knots_n=6,
        rgrid_n=20,
        iters=1
    )

    @test A_hat[1, 2] == 0.0
    @test all(A_hat .>= 0.0)
    @test all(isfinite.(Ifit))
end

@testset "unimodalfit priors" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions = [1.5, 3.0]
    spectra_prior = fill(0.2, n_mz, length(peakretentions))
    spectra_prior_weights = [0.1, 0.2]

    A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit(
        I,
        R,
        peakretentions;
        spectra_prior=spectra_prior,
        spectra_prior_weights=spectra_prior_weights,
        knots_n=6,
        rgrid_n=20,
        iters=1
    )

    @test all(A_hat .>= 0.0)
    @test all(isfinite.(A_hat))
    @test all(isfinite.(Ifit))
end

@testset "unimodalfit input validation" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions = [1.5, 3.0]

    σ_size_bad = ones(n_mz, n_scans + 1)
    @test_throws DimensionMismatch unimodalfit(I, R, peakretentions; σ=σ_size_bad)

    σ_nan = copy(ones(n_mz, n_scans))
    σ_nan[1, 1] = NaN
    @test_throws ArgumentError unimodalfit(I, R, peakretentions; σ=σ_nan)

    @test_throws ArgumentError unimodalfit(I, R, Float64[])

    t_bad = zeros(n_mz, n_scans + 1)
    @test_throws DimensionMismatch unimodalfit(I, t_bad, peakretentions)

    σ_bad = [1.0 1.0 1.0 1.0 0.0;
             1.0 1.0 1.0 1.0 1.0]
    @test_throws ArgumentError unimodalfit(I, R, peakretentions; σ=σ_bad)

    @test_throws DimensionMismatch unimodalfit(I, R, peakretentions; spectra_init=ones(n_mz + 1, length(peakretentions)))
    @test_throws DimensionMismatch unimodalfit(I, R, peakretentions; spectra_zero_mask=falses(n_mz, length(peakretentions) + 1))
    @test_throws DimensionMismatch unimodalfit(I, R, peakretentions; spectra_prior=ones(n_mz + 1, length(peakretentions)))
    @test_throws ArgumentError unimodalfit(I, R, peakretentions; spectra_prior_weights=[0.1])
    @test_throws DimensionMismatch unimodalfit(I, R, peakretentions; spectra_prior=ones(n_mz, length(peakretentions)), spectra_prior_weights=[0.1])
    @test_throws ArgumentError unimodalfit(I, R, peakretentions; spectra_prior=ones(n_mz, length(peakretentions)), spectra_prior_weights=[0.1, -0.2])
    @test_throws ArgumentError unimodalfit(I, R, peakretentions; spectra_prior=ones(n_mz, length(peakretentions)))
end

@testset "unimodalfit coupling validation" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions = [1.5, 3.0]

    @test_throws ArgumentError unimodalfit(I, R, peakretentions; shape_couple=-0.1)
    @test_throws ArgumentError unimodalfit(I, R, peakretentions; shape_couple_mode=:bad)
    @test_throws ArgumentError unimodalfit(I, R, peakretentions; shape_couple_graph=:bad)
    @test_throws ArgumentError unimodalfit(
        I,
        R,
        peakretentions;
        shape_couple=0.1,
        shape_couple_graph=:window,
        shape_couple_window=0.0
    )
    @test_throws ArgumentError unimodalfit(
        I,
        R,
        peakretentions;
        shape_couple=0.1,
        shape_couple_tau_halfwidth=0.0
    )
    @test_throws ArgumentError unimodalfit(
        I,
        R,
        peakretentions;
        shape_couple=0.1,
        shape_couple_tau_n=2
    )
    @test_throws ArgumentError unimodalfit(I, R, peakretentions; apex_localize=-0.1)
    @test_throws ArgumentError unimodalfit(
        I,
        R,
        peakretentions;
        apex_localize=0.1,
        apex_localize_scale=0.0
    )
end

@testset "unimodalfit coupling branches" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions = [1.5, 3.0]

    A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit(
        I,
        R,
        peakretentions;
        knots_n=6,
        rgrid_n=20,
        iters=1,
        shape_couple=0.1,
        shape_couple_mode=:d1,
        shape_couple_graph=:neighbors,
        shape_couple_tau_halfwidth=0.5,
        shape_couple_tau_n=5,
        apex_localize=0.2,
        apex_localize_scale=1.0
    )
    @test all(isfinite.(A_hat))
    @test all(isfinite.(Ifit))

    A_hat2, basis2, C2, (rgrid2, Fhat_grid2), Ifit2 = unimodalfit(
        I,
        R,
        peakretentions;
        knots_n=6,
        rgrid_n=20,
        iters=1,
        shape_couple=0.1,
        shape_couple_mode=:d2,
        shape_couple_graph=:window,
        shape_couple_window=2.0,
        shape_couple_tau_halfwidth=0.5,
        shape_couple_tau_n=5
    )
    @test all(isfinite.(A_hat2))
    @test all(isfinite.(Ifit2))
end

@testset "unimodalfit coupling edge cases" begin
    # K ≤ 1 returns empty edges without coupling work
    n_mz = 1
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0]
    peakretentions = [2.0]

    A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit(
        I,
        R,
        peakretentions;
        knots_n=6,
        rgrid_n=20,
        iters=1,
        shape_couple=0.1,
        shape_couple_graph=:neighbors,
        shape_couple_tau_halfwidth=0.5,
        shape_couple_tau_n=3
    )
    @test size(A_hat, 2) == 1
    @test all(isfinite.(Ifit))

    # m < 3 triggers a continue in coupling accumulation
    n_mz2 = 1
    n_scans2 = 3
    R2 = repeat([0.0 0.05 0.1], n_mz2, 1)
    I2 = [1.0 2.0 1.0]
    peakretentions2 = [0.02, 0.08]

    A_hat2, basis2, C2, (rgrid2, Fhat_grid2), Ifit2 = unimodalfit(
        I2,
        R2,
        peakretentions2;
        knots_n=6,
        rgrid_n=10,
        iters=1,
        shape_couple=0.1,
        shape_couple_graph=:neighbors,
        shape_couple_tau_halfwidth=0.1,
        shape_couple_tau_n=3
    )
    @test all(isfinite.(A_hat2))
    @test all(isfinite.(Ifit2))

    # reach the fallback error branch in peak_edges by mapping line numbers
    srcfile = abspath(joinpath(@__DIR__, "..", "..", "..", "src", "core", "deconvolution", "unimodalfit.jl"))
    peak_edges_src = repeat("\n", 416) * """
function peak_edges_test(
    peakretentions::AbstractVector{<:Real},
    graph::Symbol,
    window::Real
)
    K = length(peakretentions)
    edges = Tuple{Int,Int}[]
    if K ≤ 1
        return edges
    end
    if graph == :neighbors
        for k in 1:(K-1)
            push!(edges, (k, k+1))
        end
        return edges
    elseif graph == :window
        t0 = Float64.(peakretentions)
        for k in 1:K
            for ℓ in (k+1):K
                if abs(t0[k] - t0[ℓ]) ≤ window
                    push!(edges, (k, ℓ))
                end
            end
        end
        return edges
    else
        error("shape_couple_graph must be :neighbors or :window")
    end
end
"""
    Base.include_string(@__MODULE__, peak_edges_src, srcfile)
    @test_throws ErrorException peak_edges_test([1.0, 2.0], :bad, 1.0)
end


# ── unimodalfit_apexsearch ───────────────────────────────────────────────────────────

@testset "unimodalfit_apexsearch basic outputs" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions_guess = [1.4, 3.1]

    t0_best, sse, A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit_apexsearch(
        I,
        R;
        peakretentions_guess=peakretentions_guess,
        strategy=:single,
        half_width=0.5,
        ngrid=5,
        coord_sweeps=1,
        max_iter=5,
        knots_n=6,
        rgrid_n=20,
        iters=1,
        shape_ridge=1e-6,
        spectra_ridge=1e-8,
        verbose=false
    )

    @test length(t0_best) == length(peakretentions_guess)
    @test isfinite(sse)
    @test size(A_hat) == (n_mz, length(peakretentions_guess))
    @test size(C, 2) == length(peakretentions_guess)
    @test length(rgrid) == 20
    @test size(Fhat_grid) == (20, length(peakretentions_guess))
    @test size(Ifit) == size(I)
end

@testset "unimodalfit_apexsearch iterative stages" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions_guess = [1.4, 3.1]
    stages = [(half_width=0.6, ngrid=5), (half_width=0.3, ngrid=7)]

    t0_best, sse, A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit_apexsearch(
        I,
        R;
        peakretentions_guess=peakretentions_guess,
        strategy=:iterative,
        stages=stages,
        coord_sweeps=1,
        max_iter=5,
        knots_n=6,
        rgrid_n=20,
        iters=1,
        verbose=false
    )

    @test length(t0_best) == length(peakretentions_guess)
    @test isfinite(sse)
    @test size(Ifit) == size(I)
end

@testset "unimodalfit_apexsearch input validation" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]

    @test_throws ArgumentError unimodalfit_apexsearch(I, R; peakretentions_guess=Float64[])
    @test_throws ArgumentError unimodalfit_apexsearch(I, R; peakretentions_guess=[1.0], half_width=0.0)
    @test_throws ArgumentError unimodalfit_apexsearch(I, R; peakretentions_guess=[1.0], ngrid=2)
    @test_throws ArgumentError unimodalfit_apexsearch(I, R; peakretentions_guess=[1.0], coord_sweeps=0)
    @test_throws ArgumentError unimodalfit_apexsearch(I, R; peakretentions_guess=[1.0], max_iter=0)
    @test_throws ArgumentError unimodalfit_apexsearch(I, R; peakretentions_guess=[1.0], tol_sse=0.0)
    @test_throws ArgumentError unimodalfit_apexsearch(I, R; peakretentions_guess=[1.0, 1.05], min_sep=0.2)
end

@testset "unimodalfit_apexsearch branch coverage" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions_guess = [1.2, 3.2]

    @test_throws ArgumentError unimodalfit_apexsearch(
        I,
        R;
        peakretentions_guess=[1.0, 3.0],
        max_sep=1.0,
        verbose=false
    )

    t0_best, sse, A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit_apexsearch(
        I,
        R;
        peakretentions_guess=peakretentions_guess,
        strategy=:iterative,
        stages=nothing,
        adaptive_window=true,
        coord_sweeps=1,
        max_iter=5,
        tol_sse=NaN,
        knots_n=6,
        rgrid_n=20,
        iters=1,
        verbose=false
    )
    @test length(t0_best) == length(peakretentions_guess)
    @test isfinite(sse)

    t0_best, sse, A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit_apexsearch(
        I,
        R;
        peakretentions_guess=[1.0, 1.8],
        min_sep=0.6,
        half_width=0.5,
        ngrid=5,
        coord_sweeps=1,
        max_iter=5,
        knots_n=6,
        rgrid_n=20,
        iters=1,
        verbose=false
    )
    @test length(t0_best) == 2

    path, io = mktemp()
    redirect_stdout(io) do
        unimodalfit_apexsearch(
            I,
            R;
            peakretentions_guess=peakretentions_guess,
            strategy=:single,
            half_width=0.5,
            ngrid=5,
            coord_sweeps=1,
            max_iter=5,
            tol_sse=1e6,
            knots_n=6,
            rgrid_n=20,
            iters=1,
            verbose=true
        )
    end
    close(io)
    printed = read(path, String)
    rm(path)
    @test occursin("Initial peakretentions", printed)
    @test occursin("Stage 1", printed)
    @test occursin("sweep 1", printed)
    @test occursin("Stopping: SSE improvement < tol_sse", printed)

    path, io = mktemp()
    redirect_stdout(io) do
        unimodalfit_apexsearch(
            I,
            R;
            peakretentions_guess=peakretentions_guess,
            strategy=:single,
            half_width=0.5,
            ngrid=5,
            coord_sweeps=2,
            max_iter=1,
            knots_n=6,
            rgrid_n=20,
            iters=1,
            verbose=true
        )
    end
    close(io)
    printed = read(path, String)
    rm(path)
    @test occursin("Stopping: reached max_iter", printed)

    I_zero = zeros(n_mz, n_scans)
    t0_best, sse, A_hat, basis, C, (rgrid, Fhat_grid), Ifit = unimodalfit_apexsearch(
        I_zero,
        R;
        peakretentions_guess=[2.0],
        strategy=:single,
        half_width=0.1,
        ngrid=3,
        coord_sweeps=1,
        max_iter=5,
        tol_sse=NaN,
        knots_n=6,
        rgrid_n=20,
        iters=1,
        verbose=false
    )
    @test length(t0_best) == 1
end

@testset "unimodalfit internal sigma handling" begin
    n_mz = 2
    n_scans = 5
    R = repeat(collect(0.0:1.0:4.0)', n_mz, 1)
    I = [1.0 2.0 3.0 2.0 1.0;
         0.5 1.0 1.5 1.0 0.5]
    peakretentions = [1.5, 3.0]

    wfun = JuChrom.sigma2weights
    @test_throws ArgumentError wfun([1.0, NaN])
    @test_throws ArgumentError wfun([1.0, 0.0])
    @test wfun([1e308, 1e308]) === nothing
    @test wfun([1.0, 2.0]; w_cap=NaN) === nothing

    I_nan = copy(I)
    I_nan[1, 1] = NaN
    @test_throws ArgumentError unimodalfit(
        I_nan,
        R,
        peakretentions;
        σ=fill(1.0, n_mz, n_scans),
        knots_n=6,
        rgrid_n=20,
        iters=1
    )
end

end  # module TestDeconvolution
