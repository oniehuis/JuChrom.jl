module TestMassSpectrumMakie

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./ext/makie_extensions/massspectrum.jl
# ─────────────────────────────────────────────────────────────────────────────

using Test
using CairoMakie
using Makie
using JuChrom
using Unitful

# Force-load the Makie extension so the recipe is available
let ext = Base.get_extension(JuChrom, :MakieExtension)
    if ext ≡ nothing
        Base.require_extension(JuChrom, :MakieExtension)
    end
end

CairoMakie.activate!()  # headless backend for tests

@testset "massspectrum argument names" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    @test ext ≢ nothing
    @test ext.argument_names(ext.MassSpectrum) == (:mzvalues, :intensities)
end

@testset "massspectrum Makie plot (unitless)" begin
    mz = [100.0, 200.0, 300.0]
    intensities = [10.0, 5.0, 8.0]

    plt = massspectrum(mz, intensities)
    @test plt isa Makie.FigureAxisPlot

    lineseg = plt.plot.linesegments[]
    @test length(lineseg) == 2 * length(mz)  # two points per stick
    labels = plt.plot.mz_labels[]
    @test length(labels) == length(mz)
    @test all(!isempty, string.(labels))
    ax = only([c for c in plt.figure.content if c isa Makie.Axis])
    @test ax.xlabel[] == "m/z"
    @test ax.ylabel[] == "Intensity [unitless]"
end

@testset "massspectrum Makie plot from MassSpectrum container" begin
    mz = [100.0, 200.0, 300.0]
    intensities = [10.0, 5.0, 8.0]
    ms = MassSpectrum(mz, intensities; attrs=(source=:test,))

    plt = massspectrum(ms)
    @test plt isa Makie.FigureAxisPlot
    @test plt.plot.mzvalues[] == mz
    @test plt.plot.intensities[] == intensities

    lineseg = plt.plot.linesegments[]
    @test length(lineseg) == 2 * length(mz)
end

@testset "massspectrum Makie plot from unitful MassSpectrum container" begin
    mz = [100.0, 200.0, 300.0]
    intensities = [10.0, 5.0, 8.0]
    ms = MassSpectrum(mz, nothing, intensities, u"ms^-1"; attrs=(source=:test,))

    plt = massspectrum(ms)
    @test plt isa Makie.FigureAxisPlot
    @test plt.plot.intensities[] == intensities .* u"ms^-1"

    lineseg = plt.plot.linesegments[]
    @test lineseg[2][2] ≈ 10.0
    @test !(lineseg[2][2] isa Unitful.AbstractQuantity)
    ax = only([c for c in plt.figure.content if c isa Makie.Axis])
    @test ax.ylabel[] == "Intensity [" * string(u"ms^-1") * "]"
end

@testset "massspectrum Makie plot (unitful m/z)" begin
    mz = [100.0, 150.0]u"Th"
    intensities = [2.0, 3.0]

    plt = massspectrum(mz, intensities)
    @test plt isa Makie.FigureAxisPlot

    labels = plt.plot.mz_labels[]
    stripped = ustrip.(mz)
    @test all(val -> any(lbl -> occursin(string(val), string(lbl)), labels), stripped)
    ax = only([c for c in plt.figure.content if c isa Makie.Axis])
    @test ax.xlabel[] == "m/z [" * string(u"Th") * "]"
end

@testset "massspectrum Makie plot (unitful intensities)" begin
    mz = [100.0, 150.0]
    intensities = [2.0, 3.0]u"ms^-1"

    plt = massspectrum(mz, intensities)
    @test plt isa Makie.FigureAxisPlot

    lineseg = plt.plot.linesegments[]
    @test lineseg[2][2] ≈ 2.0
    @test lineseg[4][2] ≈ 3.0
    @test !(lineseg[2][2] isa Unitful.AbstractQuantity)
    ax = only([c for c in plt.figure.content if c isa Makie.Axis])
    limits = ax.finallimits[]
    ymax = limits.origin[2] + limits.widths[2]
    @test ymax > maximum(ustrip.(intensities))
    @test ax.ylabel[] == "Intensity [" * string(u"ms^-1") * "]"

    fig = Figure()
    target_ax = Axis(fig[1, 1])
    massspectrum!(target_ax, mz, intensities)
    @test target_ax.ylabel[] == "Intensity [" * string(u"ms^-1") * "]"
end

@testset "massspectrum intensity threshold suppresses labels" begin
    mz = [100.0, 200.0]
    intensities = [0.5, 10.0]
    plt = massspectrum(mz, intensities; intensity_threshold=1.0)
    labels = plt.plot.mz_labels[]
    @test length(labels) == 1
    val = parse(Float64, string(labels[1]))
    @test val ≈ 200.0
end

@testset "massspectrum unitful intensity threshold suppresses labels" begin
    mz = [100.0, 200.0]
    intensities = [0.5, 10.0]u"ms^-1"
    plt = massspectrum(mz, intensities; intensity_threshold=1.0u"ms^-1")
    labels = plt.plot.mz_labels[]
    @test length(labels) == 1
    val = parse(Float64, string(labels[1]))
    @test val ≈ 200.0
end

@testset "massspectrum y-limits pad above max intensity" begin
    mz = [50.0, 75.0]
    intensities = [1.0, 4.0]
    plt = massspectrum(
        mz,
        intensities;
        mz_label_fontsize=24,
        mz_label_offset=(0, 6),
        mz_label_padding=3
    )
    ax = only([c for c in plt.figure.content if c isa Makie.Axis])
    limits = ax.finallimits[]
    viewport = ax.scene.viewport[]
    ymin = limits.origin[2]
    ymax = limits.origin[2] + limits.widths[2]
    reservepixels = (1 - maximum(intensities) / ymax) * viewport.widths[2]
    ext = Base.get_extension(JuChrom, :MakieExtension)
    expectedreserve = ext.massspectrum_label_headroom_pixels(
        plt.plot,
        mz,
        intensities,
        nothing
    )
    box = Makie.text_bb(
        string(mz[2]),
        ext.massspectrum_measurement_font(plt.plot.mz_label_font[]),
        plt.plot.mz_label_fontsize[]
    )
    offset_y = plt.plot.mz_label_offset[][2]
    visibletopgap = reservepixels - offset_y - Float64(box.widths[2])

    @test ymax > maximum(intensities)
    @test ymin ≤ 0
    @test reservepixels ≈ expectedreserve atol = 3
    @test visibletopgap ≈ 1.5 * offset_y atol = 3
end

@testset "massspectrum label headroom stays fixed in pixels" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    mz = [50.0, 75.0]
    intensities = [1.0, 4.0]

    function reserve_pixels(height)
        fig = Figure(size=(400, height))
        ax = Axis(fig[1, 1])
        plot = massspectrum!(
            ax,
            mz,
            intensities;
            mz_label_fontsize=24,
            mz_label_offset=(0, 6),
            mz_label_padding=3
        )
        viewport = ax.scene.viewport[]
        limits = ax.finallimits[]
        ylimit = limits.origin[2] + limits.widths[2]
        reserve = (1 - maximum(intensities) / ylimit) * viewport.widths[2]
        expected = ext.massspectrum_label_headroom_pixels(plot, mz, intensities, nothing)

        reserve, expected
    end

    reserve_small, expected_small = reserve_pixels(260)
    reserve_tall, expected_tall = reserve_pixels(620)

    @test reserve_small ≈ expected_small atol = 3
    @test reserve_tall ≈ expected_tall atol = 3
    @test reserve_small ≈ reserve_tall atol = 3
end

@testset "massspectrum label headroom updates after resize" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    mz = [50.0, 75.0]
    intensities = [1.0, 4.0]
    fig = Figure(size=(400, 260))
    ax = Axis(fig[1, 1])
    plot = massspectrum!(
        ax,
        mz,
        intensities;
        mz_label_fontsize=24,
        mz_label_offset=(0, 6),
        mz_label_padding=3
    )

    function reserve_pixels()
        viewport = ax.scene.viewport[]
        limits = ax.finallimits[]
        ylimit = limits.origin[2] + limits.widths[2]
        (1 - maximum(intensities) / ylimit) * viewport.widths[2]
    end

    expected = ext.massspectrum_label_headroom_pixels(plot, mz, intensities, nothing)
    @test reserve_pixels() ≈ expected atol = 3
    resize!(fig, 400, 620)
    @test reserve_pixels() ≈ expected atol = 3
end

@testset "massspectrum y-limits fallback for nonpositive max intensity" begin
    mz = [10.0, 20.0]
    intensities = [0.0, 0.0]
    plt = massspectrum(mz, intensities)
    ax = only([c for c in plt.figure.content if c isa Makie.Axis])
    limits = ax.finallimits[]
    ymin = limits.origin[2]
    ymax = limits.origin[2] + limits.widths[2]
    @test isapprox(ymin, 0.0; atol=1e-6)
    @test isapprox(ymax, 1.0; atol=1e-6)
end

end  # module
