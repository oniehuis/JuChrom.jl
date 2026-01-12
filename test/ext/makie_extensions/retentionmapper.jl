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

struct DummyRetentionMapper <: JuChrom.AbstractRetentionMapper{Nothing, Nothing}
    rA::Vector{Float64}
    rA_unit::Nothing
end

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

@testset "RetentionMapper residual annotation colors" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    if ext === nothing
        Base.require_extension(JuChrom, :MakieExtension)
        ext = Base.get_extension(JuChrom, :MakieExtension)
    end
    ext === nothing && error("JuChrom MakieExtension could not be loaded; ensure Makie is available.")

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])

    ext.add_single_residual_annotation!(ax, 1.0, 0.0, 1.0, 0.5, 10.0)
    ext.add_single_residual_annotation!(ax, 2.0, 1.0, -1.0, 0.5, 10.0)

    texts = [p for p in ax.scene.plots if p isa Makie.Text]
    @test length(texts) ≥ 2
    text_values = collect(Iterators.flatten([t.attributes[:text][] for t in texts]))
    @test any(t -> occursin("1.0", t), text_values)
    @test any(t -> occursin("-1.0", t), text_values)
    text_colors = [t.attributes[:color][] for t in texts]
    @test any(c -> c == Makie.to_color(:black), text_colors)
    @test any(c -> c == Makie.to_color(:red), text_colors)
end

@testset "RetentionMapper text positioning helper" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    if ext === nothing
        Base.require_extension(JuChrom, :MakieExtension)
        ext = Base.get_extension(JuChrom, :MakieExtension)
    end
    ext === nothing && error("JuChrom MakieExtension could not be loaded; ensure Makie is available.")

    align_low, offset_low = ext.determine_text_position(0.0, 1.0, 10.0)
    @test align_low == (:left, :center)
    @test offset_low == (0, 10.0)

    align_high, offset_high = ext.determine_text_position(2.0, 1.0, 10.0)
    @test align_high == (:right, :center)
    @test offset_high == (0, -10.0)
end

@testset "RetentionMapper inverse derivative unit string" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    if ext === nothing
        Base.require_extension(JuChrom, :MakieExtension)
        ext = Base.get_extension(JuChrom, :MakieExtension)
    end
    ext === nothing && error("JuChrom MakieExtension could not be loaded; ensure Makie is available.")

    @test ext.format_inverse_derivative_unit_string(u"s", nothing) == " [s]"
    norm_unit = s -> replace(s, "⁻¹" => "^-1")
    @test norm_unit(ext.format_inverse_derivative_unit_string(nothing, u"s")) == " [s^-1]"
    @test norm_unit(ext.format_inverse_derivative_unit_string(u"minute", u"s")) == " [minute × s^-1]"
end

@testset "RetentionMapper derivative unit string" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    if ext === nothing
        Base.require_extension(JuChrom, :MakieExtension)
        ext = Base.get_extension(JuChrom, :MakieExtension)
    end
    ext === nothing && error("JuChrom MakieExtension could not be loaded; ensure Makie is available.")

    @test ext.format_derivative_unit_string(nothing, u"s") == " [s]"
    norm_unit = s -> replace(s, "⁻¹" => "^-1")
    @test norm_unit(ext.format_derivative_unit_string(u"s", nothing)) == " [s^-1]"
end

@testset "RetentionMapper validation" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    if ext === nothing
        Base.require_extension(JuChrom, :MakieExtension)
        ext = Base.get_extension(JuChrom, :MakieExtension)
    end
    ext === nothing && error("JuChrom MakieExtension could not be loaded; ensure Makie is available.")

    rm_empty = DummyRetentionMapper(Float64[], nothing)
    @test_throws ArgumentError ext.validate_retention_mapper(rm_empty)

    rm_single = DummyRetentionMapper([1.0], nothing)
    @test_logs (:warn, r"Only one data point available") ext.validate_retention_mapper(rm_single)
end

@testset "RetentionMapper residual annotation loop" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    if ext === nothing
        Base.require_extension(JuChrom, :MakieExtension)
        ext = Base.get_extension(JuChrom, :MakieExtension)
    end
    ext === nothing && error("JuChrom MakieExtension could not be loaded; ensure Makie is available.")

    rm = fitmap([1.0, 2.0, 3.0], [1.0, 2.0, 4.0])
    ret_pred = retentions_B(rm) .+ [0.0, 0.1, -0.2]

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])
    ext.add_residual_annotations!(ax, rm, ret_pred, nothing, nothing, 10.0, 1)

    texts = [p for p in ax.scene.plots if p isa Makie.Text]
    text_values = collect(Iterators.flatten([p.attributes[:text][] for p in texts]))
    @test any(t -> occursin("0.1", t), text_values)
    @test any(t -> occursin("-0.2", t), text_values)
end

end  # module TestRetentionsMapper
