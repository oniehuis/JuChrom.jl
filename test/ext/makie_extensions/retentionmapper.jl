module TestRetentionsMapper

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./ext/makie_extensions/retentionsmapper.jl
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

using Test
using CairoMakie
using Unitful

# ── Internal Project Imports ─────────────────────────────────────────────────

using JuChrom
using JuChrom: RetentionMapper

CairoMakie.activate!()  # headless backend for tests

@testset "RetentionMapper Makie plot (unitless)" begin
    rm = fitmap([1.0, 2.0, 4.0], [10.0, 20.0, 40.0])
    fig = Makie.plot(rm; size=(400, 300), reverse=true)

    @test fig isa Makie.Figure
    axes = [c for c in fig.content if c isa Makie.Axis]
    @test length(axes) == 4  # forward + reverse axes when reverse=true
    titles = [ax.title[] for ax in axes]
    @test "Retention A → Retention B" in titles
    @test "Retention B → Retention A" in titles
end

@testset "RetentionMapper plot reverse toggle" begin
    rm = fitmap([0.0, 1.0, 2.0], [0.0, 2.0, 4.0])

    fig_forward = Makie.plot(rm; size=(300, 200), reverse=false)
    axes_fwd = [c for c in fig_forward.content if c isa Makie.Axis]
    @test length(axes_fwd) == 2

    fig_both = Makie.plot(rm; size=(300, 200), reverse=true)
    axes_both = [c for c in fig_both.content if c isa Makie.Axis]
    @test length(axes_both) == 4
end

@testset "RetentionMapper plot data plumbing" begin
    rm = fitmap([0.0, 1.0, 2.0], [0.0, 2.0, 4.0])
    fig = Makie.plot(rm; size=(300, 200), reverse=true)
    axes = Dict(ax.title[] => ax for ax in fig.content if ax isa Makie.Axis)

    # Forward mapping axis: scatter positions should match raw retentions
    fwd_ax = axes["Retention A → Retention B"]
    elements = fwd_ax.scene.plots
    scatters = [p for p in elements if p isa Makie.Scatter]
    @test !isempty(scatters)
    positions = scatters[1].attributes[:positions][]
    xs = Float64[p[1] for p in positions]
    ys = Float64[p[2] for p in positions]
    @test xs ≈ rawretentions_A(rm)
    @test ys ≈ rawretentions_B(rm)

    # Reverse mapping axis: line should follow rawinvmap
    rev_ax = axes["Retention B → Retention A"]
    lines = [p for p in rev_ax.scene.plots if p isa Makie.Lines]
    @test length(lines) == 1
end

@testset "RetentionMapper Makie plot (unitful)" begin
    rm = fitmap([1.0, 2.0, 4.0]u"s", [10.0, 20.0, 40.0]u"s")
    fig = Makie.plot(rm; size=(300, 200), reverse=false)
    @test fig isa Makie.Figure
    axes = [c for c in fig.content if c isa Makie.Axis]
    @test length(axes) >= 2
end

end  # module TestRetentionsMapper
