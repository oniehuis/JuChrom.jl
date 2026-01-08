module TestWhitening

using Test
using JuChrom

# ─────────────────────────────────────────────────────────────────────────────
# whiten
# ─────────────────────────────────────────────────────────────────────────────

@testset "whiten(x::AbstractArray{<:Real}, σ²::AbstractArray{<:Real}, σ²_floor::Real)" begin
    # --- basic correctness (vector) ---
    @testset "vector correctness & flooring" begin
        x   = [1.0, 2.0, 3.0, 4.0]
        σ²  = [1.0, 0.25, 1e-12, 9.0]  # includes very small variance to trigger floor
        f   = 1e-6                      # positive floor
        y   = JuChrom.whiten(x, σ², f)

        σ²safe = max.(σ², f)
        @test y ≈ x ./ sqrt.(σ²safe)

        # inputs are not mutated
        @test x == [1.0, 2.0, 3.0, 4.0]
        @test σ² == [1.0, 0.25, 1e-12, 9.0]
        # element with near-zero variance uses floor
        @test y[3] ≈ x[3] / sqrt(f)
    end

    # --- matrix correctness ---
    @testset "matrix correctness" begin
        X  = [1 2 3;
              4 5 6] .+ 0.0           # make it Float64 explicitly
        V  = [1 4 9;
              1 1 1e-10]
        f  = 1e-6
        Y  = JuChrom.whiten(X, V, f)
        @test Y ≈ X ./ sqrt.(max.(V, f))
        @test size(Y) == size(X)
    end

    # --- dimension mismatch ---
    @testset "dimension checks" begin
        x  = [1.0, 2.0, 3.0]
        v1 = [1.0, 1.0]             # wrong length
        @test_throws DimensionMismatch JuChrom.whiten(x, v1, 1e-6)

        X  = ones(2, 3)
        V  = ones(3, 2)             # wrong shape
        @test_throws DimensionMismatch JuChrom.whiten(X, V, 1e-6)
    end

    # --- floor validation ---
    @testset "σ²_floor must be positive" begin
        x = [1.0, 2.0]
        v = [1.0, 1.0]
        @test_throws ArgumentError JuChrom.whiten(x, v, 0.0)
        @test_throws ArgumentError JuChrom.whiten(x, v, -1.0)
    end

    # --- integer inputs promote to float ---
    @testset "integer inputs promote" begin
        x  = [1, 2, 3]                      # Int
        v  = [1, 4, 9]                      # Int
        f  = 1                              # Int floor (positive)
        y  = JuChrom.whiten(x, v, f)
        @test eltype(y) <: AbstractFloat
        @test y ≈ x ./ sqrt.(max.(v, f))
    end
end


end  # module TestWhitening
