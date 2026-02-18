module TestMath

using Test
using Unitful
using JuChrom

# ── Constants ────────────────────────────────────────────────────────────────
const ATOL = 1e-8  # Tolerance for floating point comparisons

# ── cosdis ───────────────────────────────────────────────────────────────────

@testset "cosdis(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, positive_only::Bool=true)" begin
    # Orthogonal vectors: similarity 0 → distance 1
    x = [1.0, 0.0]; y = [0.0, 1.0]
    @test isapprox(cosdis(x, y, false), 1.0; atol=ATOL)
    @test isapprox(cosdis(x, y, true),  1.0; atol=ATOL)

    # Identical vectors: similarity 1 → distance 0
    x = [1.0, 2.0, 3.0]
    @test isapprox(cosdis(x, x, false), 0.0; atol=ATOL)
    @test isapprox(cosdis(x, x, true),  0.0; atol=ATOL)

    # Opposite vectors:
    #   positive_only=false → sim = -1 → dist = 2
    #   positive_only=true  → sim clamped to 0 → dist = 1
    x = [1.0, 0.0]; y = [-1.0, 0.0]
    @test isapprox(cosdis(x, y, false), 2.0; atol=ATOL)
    @test isapprox(cosdis(x, y, true),  1.0; atol=ATOL)

    # Zero vector → NaN
    x = [0.0, 0.0]; y = [1.0, 2.0]
    @test isnan(cosdis(x, y, false))
    @test isnan(cosdis(x, y, true))

    # Range checks
    x = [1.0, -1.0]; y = [-1.0, 1.0]
    d  = cosdis(x, y, false)
    dp = cosdis(x, y, true)
    @test 0.0 ≤ d  ≤ 2.0
    @test 0.0 ≤ dp ≤ 1.0

    # Works with integers too
    @test isapprox(cosdis([1,2,3], [1,2,3], false), 0.0; atol=ATOL)

    # Weighted variant dispatches on weight vector
    x = [1.0, 2.0]; y = [3.0, 4.0]; w = [1.0, 0.0]
    @test isapprox(cosdis(x, y, w, false), 0.0; atol=ATOL)
    @test isapprox(cosdis(x, y, ones(2), false), cosdis(x, y, false); atol=ATOL)

    # Length checks
    @test_throws DimensionMismatch cosdis([1.0, 2.0], [1.0], false)
    @test_throws DimensionMismatch cosdis([1.0, 2.0], [1.0, 2.0], [1.0], false)
end

# ── cossim ───────────────────────────────────────────────────────────────────

@testset "cossim(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, positive_only::Bool=true)" begin
    # Orthogonal vectors
    x = [1.0, 0.0]; y = [0.0, 1.0]
    @test isapprox(cossim(x, y, false), 0.0; atol=ATOL)
    @test cossim(x, y, true) == 0.0  # clamped

    # Identical vectors
    x = [1.0, 2.0, 3.0]
    @test isapprox(cossim(x, x, false), 1.0; atol=ATOL)
    @test isapprox(cossim(x, x, true),  1.0; atol=ATOL)

    # Opposite vectors
    x = [1.0, 0.0]; y = [-1.0, 0.0]
    @test isapprox(cossim(x, y, false), -1.0; atol=ATOL)
    @test cossim(x, y, true) == 0.0   # clamped

    # Zero vector → NaN
    x = [0.0, 0.0]; y = [1.0, 2.0]
    @test isnan(cossim(x, y, false))
    @test isnan(cossim(x, y, true))

    # Range checks
    x = [1.0, -1.0]; y = [-1.0, 1.0]
    s  = cossim(x, y, false)
    sp = cossim(x, y, true)
    @test -1.0 ≤ s  ≤ 1.0
    @test  0.0 ≤ sp ≤ 1.0

    # Works with integers too
    @test isapprox(cossim([1,2,3], [1,2,3], false), 1.0; atol=ATOL)

    # Weighted variant dispatches on weight vector
    x = [1.0, 2.0]; y = [3.0, 4.0]; w = [1.0, 0.0]
    @test isapprox(cossim(x, y, w, false), 1.0; atol=ATOL)
    @test isapprox(cossim(x, y, ones(2), false), cossim(x, y, false); atol=ATOL)

    # Length checks
    @test_throws DimensionMismatch cossim([1.0, 2.0], [1.0], false)
    @test_throws DimensionMismatch cossim([1.0, 2.0], [1.0, 2.0], [1.0], false)
end

# ─────────────────────────────────────────────────────────────────────────────
# localmaxima (renamed from local_maxima_indices)
# ─────────────────────────────────────────────────────────────────────────────
@testset "localmaxima(values::AbstractVector{<:Real})" begin
    @test JuChrom.localmaxima([1, 3, 2]) == [2]
    @test JuChrom.localmaxima([1, 2, 3, 2, 1]) == [3]
    @test JuChrom.localmaxima([5, 4, 3, 2, 1]) == Int[]
    @test JuChrom.localmaxima([1, 2, 3, 4, 5]) == Int[]

    # Multiple peaks
    @test JuChrom.localmaxima([0, 2, 1, 3, 1, 4, 2]) == [2, 4, 6]

    # Flat (plateau) and equal regions
    @test JuChrom.localmaxima([1, 2, 2, 1]) == Int[]
    @test JuChrom.localmaxima([1, 1, 1]) == Int[]
    @test JuChrom.localmaxima([2, 2, 2, 2, 2]) == Int[]

    # Short input
    @test JuChrom.localmaxima(Real[]) == Int[]
    @test JuChrom.localmaxima(Float64[]) == Int[]
    @test JuChrom.localmaxima(Int[]) == Int[]
    @test JuChrom.localmaxima([1]) == Int[]
    @test JuChrom.localmaxima([1, 2]) == Int[]

    # Noise and dips
    @test JuChrom.localmaxima([1, 3, 2, 4, 3, 5, 4]) == [2, 4, 6]

    # Negative values
    @test JuChrom.localmaxima([-2, -1, -3]) == [2]
    @test JuChrom.localmaxima([-3, -1, -2, 0, -1]) == [2, 4]

    # Positive and negative Float64 values
    @test JuChrom.localmaxima([0.0, 2.1, 1.2, 3.1, 1.8, 4.0, 2.5]) == [2, 4, 6]
    @test JuChrom.localmaxima([-3.0, 1.1, -2.1, 0.0, -1.9]) == [2, 4]
end

# ─────────────────────────────────────────────────────────────────────────────
# minmax_scale / inverse_minmax_scale
# ─────────────────────────────────────────────────────────────────────────────
@testset "minmax_scale & inverse_minmax_scale" begin
    @test JuChrom.minmax_scale(15, 10, 20) == 0.5
    @test JuChrom.minmax_scale(5, 5, 15) == 0.0
    @test JuChrom.minmax_scale(15, 5, 15) == 1.0

    @test JuChrom.inverse_minmax_scale(0.5, 10, 20) == 15.0
    @test JuChrom.inverse_minmax_scale(0.0, 5, 15) == 5.0
    @test JuChrom.inverse_minmax_scale(1.0, 5, 15) == 15.0

    # round-trip
    for v in (10.0, 12.5, 20.0)
        n = JuChrom.minmax_scale(v, 10.0, 20.0)
        @test JuChrom.inverse_minmax_scale(n, 10.0, 20.0) ≈ v
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# maxdecimals
# ─────────────────────────────────────────────────────────────────────────────
@testset "maxdecimals" begin
    @test JuChrom.maxdecimals([1.23, 4.567, 8.9]) == 3
    @test JuChrom.maxdecimals([1.0, 2.0, 3.0]) == 0
    @test JuChrom.maxdecimals(Int[]) == 0
    @test JuChrom.maxdecimals([10, 2.5, 3.14]) == 2
end

end # module TestMath
