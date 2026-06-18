module TestTicTraceMakie

using Test
using CairoMakie
using Makie
using JuChrom

let ext = Base.get_extension(JuChrom, :MakieExtension)
    if ext ≡ nothing
        Base.require_extension(JuChrom, :MakieExtension)
    end
end

CairoMakie.activate!()

function synthetic_ladder_result()
    peakmodel = [0.0, 0.5, 1.0, 0.5, 0.0]
    X = zeros(5, 5)
    X[:, 2] .= 100.0 .* peakmodel
    X[:, 3] .= 20.0 .* peakmodel
    X[:, 4] .= 5.0 .* peakmodel
    msm = MassScanMatrix(collect(1.0:5.0), [100.0, 114.0, 115.0, 128.0, 129.0], X)
    standard = AlkaneStandard(
        "test C8",
        [JuChrom.alkane_reference_spectrum(8, "octane", 800.0, [114], [100.0])],
        NamedTuple()
    )
    varianceestimate = CountVarianceEstimate(
        ones(size(X)),
        nothing,
        1.0,
        1.0,
        1.0,
        1.0,
        1,
        1,
        5,
        1,
        1.0,
        0.0,
        1.0
    )
    result = findalkaneseries(
        msm,
        varianceestimate;
        standard=standard,
        carbonrange=8:8,
        thresholdfraction=0.05,
        pathminsteps=1,
        apexminioncount=1,
        apexmzscanorder=:ascending,
        apexmzretentionkwargs=(
            retentionref=:middle,
            scaninterval=1e-9,
            mzcount=mzcount(msm),
            dwellref=:middle,
            dwell=:homogeneous
        )
    )

    msm, result
end

function only_axis(fig)
    only([content for content in fig.content if content isa Makie.Axis])
end

@testset "tictrace overlays ladder steps from AlkaneSeriesResult" begin
    msm, result = synthetic_ladder_result()
    ext = Base.get_extension(JuChrom, :MakieExtension)

    fig = tictrace(msm, result)
    ax = only_axis(fig)

    @test ext.argument_names(ext.TicTrace) == (
        :retentions,
        :intensities,
        :baseline_intensities,
        :step_retentions,
        :step_numbers
    )
    @test length(alkaneladdersteps(result)) == 1
    @test length(ax.scene.plots) == 1
    @test ax.xlabel[] == "Retention"
    @test ax.ylabel[] == "TIC [unitless]"
    @test only(ax.scene.plots).labels[] == true
end

@testset "tictrace supports retention unit conversion" begin
    msm, result = synthetic_ladder_result()
    ext = Base.get_extension(JuChrom, :MakieExtension)
    unitful_msm = MassScanMatrix(
        rawretentions(msm),
        u"ms",
        rawmzvalues(msm),
        nothing,
        rawintensities(msm),
        nothing
    )
    unitful_result = AlkaneSeriesResult(
        result.standard,
        result.variances,
        result.varianceinfo,
        result.baselineinfo,
        result.channelinfo,
        result.abundanceinfo,
        result.molecularioninfo,
        result.pathinfo,
        result.apexinfo,
        result.additioninfo,
        JuChrom.alkane_series_datainfo(unitful_msm, unitful_msm, result.variances),
        u"ms"
    )

    fig = tictrace(unitful_msm, unitful_result; unit=u"minute")
    ax = only_axis(fig)
    limits = ax.finallimits[]

    @test ax.xlabel[] == "Retention [minute]"
    @test limits.origin[1] ≈ minimum(rawretentions(unitful_msm; unit=u"minute"))
    @test limits.origin[1] + limits.widths[1] ≈
        maximum(rawretentions(unitful_msm; unit=u"minute"))

    step = only(alkaneladdersteps(result))
    @test_throws ArgumentError ext.tictrace_step_retention(step, nothing, u"minute")
    @test ext.tictrace_step_retention(step, u"ms", nothing) == step.apexretention
    font = ext.tictrace_measurement_font(:regular)
    @test ext.tictrace_measurement_font(font) === Makie.to_font(font)
end

@testset "tictrace! overlays ladder steps on an existing axis" begin
    msm, result = synthetic_ladder_result()
    fig = Figure()
    ax = Axis(fig[1, 1])

    out = tictrace!(ax, msm, result)

    @test out === ax
    @test length(ax.scene.plots) == 1
end

@testset "tictrace validates AlkaneSeriesResult checksum against plotted matrix" begin
    msm, result = synthetic_ladder_result()
    mismatched_msm = MassScanMatrix(
        rawretentions(msm),
        rawmzvalues(msm),
        rawintensities(msm) .+ 1.0
    )

    @test tictrace(msm, result) isa Makie.Figure
    @test tictrace(VarianceMassScanMatrix(msm, result.variances), result) isa Makie.Figure
    @test_throws ArgumentError tictrace(mismatched_msm, result)

    fig = Figure()
    ax = Axis(fig[1, 1])
    @test_throws ArgumentError tictrace!(ax, mismatched_msm, result)
end

@testset "tictrace ladder labels are thinned when dense" begin
    ext = Base.get_extension(JuChrom, :MakieExtension)
    fig = Figure(size=(300, 300))
    ax = Axis(fig[1, 1])
    xlims!(ax, 1, 20)
    ylims!(ax, 0, 1)

    positions, labels = ext.tictrace_visible_label_data(
        ax,
        collect(1.0:20.0),
        collect(1:20),
        true,
        :bold,
        11,
        1,
        4.0,
        0.99,
        ax.scene.viewport[],
        ax.finallimits[]
    )

    @test !isempty(positions)
    @test length(positions) == length(labels)
    @test length(labels) < 20
end

@testset "tictrace accepts integer TIC with floating baseline" begin
    msm, result = synthetic_ladder_result()
    int_msm = MassScanMatrix(
        rawretentions(msm),
        rawmzvalues(msm),
        round.(Int, rawintensities(msm))
    )
    float_baselines = MassScanMatrix(
        rawretentions(msm),
        rawmzvalues(msm),
        fill(0.1, size(rawintensities(msm)))
    )
    patched_result = AlkaneSeriesResult(
        result.standard,
        result.variances,
        result.varianceinfo,
        JuChrom.AlkaneBaselineInfo(
            float_baselines,
            :test,
            1.0,
            true,
            10.0,
            0.5,
            1.0,
            0.2
        ),
        result.channelinfo,
        result.abundanceinfo,
        result.molecularioninfo,
        result.pathinfo,
        result.apexinfo,
        result.additioninfo,
        JuChrom.alkane_series_datainfo(int_msm, int_msm, result.variances),
        result.retentionunit
    )

    fig = tictrace(int_msm, patched_result)

    @test fig isa Makie.Figure
end

end
