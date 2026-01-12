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
    if ext === nothing
        Base.require_extension(JuChrom, :MakieExtension)
    end
end

CairoMakie.activate!()  # headless backend for tests

@testset "massspectrum argument names" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    @test ext !== nothing
    @test Makie.argument_names(ext.MassSpectrum, 2) == (:mzvalues, :intensities)
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
end

@testset "massspectrum Makie plot (unitful m/z)" begin
    mz = [100.0, 150.0]u"Th"
    intensities = [2.0, 3.0]

    plt = massspectrum(mz, intensities)
    @test plt isa Makie.FigureAxisPlot

    labels = plt.plot.mz_labels[]
    stripped = ustrip.(mz)
    @test all(val -> any(lbl -> occursin(string(val), string(lbl)), labels), stripped)
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

@testset "massspectrum y-limits pad above max intensity" begin
    mz = [50.0, 75.0]
    intensities = [1.0, 4.0]
    plt = massspectrum(mz, intensities)
    ax = only([c for c in plt.figure.content if c isa Makie.Axis])
    limits = ax.finallimits[]
    ymin = limits.origin[2]
    ymax = limits.origin[2] + limits.widths[2]
    @test ymax > maximum(intensities)
    @test ymin ≤ 0
end

end  # module
