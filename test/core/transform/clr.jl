module TestClr

using Statistics: mean
using Test
using JuChrom

# ─────────────────────────────────────────────────────────────────────────────
# clr
# ─────────────────────────────────────────────────────────────────────────────

@testset "clr(x::AbstractArray{<:Real})" begin
    @testset "basic correctness & zero-sum" begin
        x = [1.0, 2.0, 3.0]
        y = JuChrom.clr(x)
        # manual expectation: log.(x) .- mean(log.(x))
        exp = log.(x) .- mean(log.(x))
        @test y ≈ exp
        @test isapprox(sum(y), 0.0; atol=1e-12)
        @test size(y) == size(x)
    end

    @testset "scale invariance" begin
        x = [0.5, 1.5, 3.0, 6.0]
        y1 = JuChrom.clr(x)
        y2 = JuChrom.clr(5.0 .* x)  # scaling should not change clr
        @test y1 ≈ y2
    end

    @testset "matrix input preserves shape and zero-sum over all elements" begin
        X = [1.0  2.0
             3.0  6.0]
        Y = JuChrom.clr(X)
        @test size(Y) == size(X)
        @test isapprox(sum(Y), 0.0; atol=1e-12)
        @test Y ≈ (log.(X) .- mean(log.(X)))
    end

    @testset "errors on non-positive components" begin
        @test_throws DomainError JuChrom.clr([0.0, 1.0, 2.0])
        @test_throws DomainError JuChrom.clr([-1.0, 1.0, 2.0])
    end
end

@testset "clr(x, variances)" begin
    @testset "dimension checks" begin
        x = [1.0, 2.0, 3.0]
        v_bad = [0.1, 0.2]  # wrong length
        @test_throws DimensionMismatch JuChrom.clr(x, v_bad)
    end

    @testset "domain checks" begin
        x_bad = [0.0, 1.0, 2.0]
        v = [0.01, 0.01, 0.04]
        @test_throws DomainError JuChrom.clr(x_bad, v)
    end

    @testset "correctness of transform and variance propagation (vector)" begin
        x = [1.0, 2.0, 4.0]
        v = [0.01, 0.04, 0.09]  # variances (std = [0.1, 0.2, 0.3])
        y, v_clr = JuChrom.clr(x, v)

        # Expected CLR
        y_exp = log.(x) .- mean(log.(x))
        @test y ≈ y_exp
        @test isapprox(sum(y), 0.0; atol=1e-12)

        # Expected variance propagation (exact finite-sample correction)
        N = length(x)
        σ_log = sqrt.(v) ./ x
        σ2_log = σ_log .^ 2
        Σσ2 = sum(σ2_log)
        v_exp = @. σ2_log * (1 - 2/N) + Σσ2 / N^2

        @test v_clr ≈ v_exp
        @test all(>=(0), v_clr)  # non-negative
        @test size(v_clr) == size(x)
    end

    @testset "matrix input (shape preserved)" begin
        X = [1.0 2.0; 3.0 4.0]
        V = [0.01 0.04; 0.09 0.16]
        Y, Vc = JuChrom.clr(X, V)

        # transform
        @test size(Y) == size(X)
        @test Y ≈ (log.(X) .- mean(log.(X)))
        @test isapprox(sum(Y), 0.0; atol=1e-12)

        # variances
        N = length(X)
        σ_log = sqrt.(V) ./ X
        σ2_log = σ_log .^ 2
        Σσ2 = sum(σ2_log)
        Vexp = @. σ2_log * (1 - 2/N) + Σσ2 / N^2

        @test size(Vc) == size(X)
        @test Vc ≈ Vexp
        @test all(>=(0), Vc)
    end
end

end  # module TestClr
