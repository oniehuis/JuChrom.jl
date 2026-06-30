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

    sized = tictrace(
        msm,
        result;
        figure=(; size=(320, 240)),
        axis=(; title="Annotated ladder", xgridvisible=true)
    )
    sized_ax = only_axis(sized)

    @test sized.scene.viewport[].widths == [320, 240]
    @test sized_ax.title[] == "Annotated ladder"
    @test sized_ax.xgridvisible[] == true
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
        result.success,
        result.status,
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

    @test_throws MethodError tictrace(unitful_msm, unitful_result; unit=u"minute")

    fig = tictrace(unitful_msm, unitful_result; retentionunit=u"minute")
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

    unitful_signal_msm = MassScanMatrix(
        rawretentions(msm),
        u"ms",
        rawmzvalues(msm),
        nothing,
        rawintensities(msm),
        u"pA"
    )
    baseline_msm = MassScanMatrix(
        rawretentions(msm),
        u"ms",
        rawmzvalues(msm),
        nothing,
        fill(1000.0, size(rawintensities(msm))),
        u"pA"
    )
    unitful_signal_result = AlkaneSeriesResult(
        result.success,
        result.status,
        result.standard,
        result.variances,
        result.varianceinfo,
        JuChrom.AlkaneBaselineInfo(
            baseline_msm,
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
        JuChrom.alkane_series_datainfo(
            unitful_signal_msm,
            unitful_signal_msm,
            result.variances
        ),
        u"ms"
    )

    converted = tictrace(
        unitful_signal_msm,
        unitful_signal_result;
        retentionunit=u"minute",
        intensityunit=u"nA"
    )
    converted_ax = only_axis(converted)
    converted_plot = only(converted_ax.scene.plots)

    @test converted_ax.xlabel[] == "Retention [minute]"
    @test converted_ax.ylabel[] == "TIC [nA]"
    @test converted_plot.retentions[] ≈ rawretentions(unitful_signal_msm; unit=u"minute")
    @test converted_plot.intensities[] ≈
        vec(sum(rawintensities(unitful_signal_msm; unit=u"nA"); dims=2))
    @test converted_plot.baseline_intensities[] ≈ fill(5.0, scancount(unitful_signal_msm))
end

@testset "tictrace plots plain TIC from JuChrom containers" begin
    X = [1000.0 2000.0; 3000.0 4000.0]
    msm = MassScanMatrix(
        [60.0, 120.0],
        u"s",
        [100.0, 101.0],
        nothing,
        X,
        u"pA"
    )

    fig = tictrace(
        msm;
        figure=(; size=(320, 240)),
        retentionunit=u"minute",
        intensityunit=u"nA",
        color=:red,
        linewidth=2
    )
    ax = only_axis(fig)
    plt = only(ax.scene.plots)

    @test fig isa Makie.Figure
    @test fig.scene.viewport[].widths == [320, 240]
    @test ax.xlabel[] == "Retention [minute]"
    @test ax.ylabel[] == "Intensity [nA]"
    @test ax.title[] == ""
    @test plt isa Makie.Lines
    @test plt.arg1[] ≈ [1.0, 2.0]
    @test plt.arg2[] ≈ [3.0, 7.0]
    @test plt.linewidth[] ≈ 2

    titled = tictrace(
        msm;
        title="Signal TIC",
        retentionunit=u"minute",
        intensityunit=u"nA"
    )
    @test only_axis(titled).title[] == "Signal TIC"

    fig_existing = Figure()
    ax_existing = Axis(fig_existing[1, 1])
    plt_existing = tictrace!(
        ax_existing,
        msm;
        title="Existing axis TIC",
        retentionunit=u"minute",
        intensityunit=u"nA",
        linewidth=3
    )
    @test plt_existing isa Makie.Lines
    @test ax_existing.xlabel[] == "Retention [minute]"
    @test ax_existing.ylabel[] == "Intensity [nA]"
    @test ax_existing.title[] == "Existing axis TIC"
    @test plt_existing.arg1[] ≈ [1.0, 2.0]
    @test plt_existing.arg2[] ≈ [3.0, 7.0]
    @test plt_existing.linewidth[] ≈ 3

    vmsm = VarianceMassScanMatrix(msm, ones(size(X)))
    vfig = tictrace(vmsm; retentionunit=u"minute", intensityunit=u"nA")
    vplt = only(only_axis(vfig).scene.plots)
    @test vplt.arg1[] ≈ [1.0, 2.0]
    @test vplt.arg2[] ≈ [3.0, 7.0]

    mss = MassScanSeries([
        MassScan(60.0u"s", [100.0, 101.0], [1000.0, 2000.0]u"pA"),
        MassScan(120.0u"s", [100.0, 101.0], [3000.0, 4000.0]u"pA")
    ])
    mfig = tictrace(mss; retentionunit=u"minute", intensityunit=u"nA")
    mplt = only(only_axis(mfig).scene.plots)
    @test mplt.arg1[] ≈ [1.0, 2.0]
    @test mplt.arg2[] ≈ [3.0, 7.0]

    unitless = MassScanMatrix([1.0, 2.0], [100.0], reshape([5.0, 7.0], 2, 1))
    unitless_ax = only_axis(tictrace(unitless))
    @test unitless_ax.xlabel[] == "Retention [unitless]"
    @test unitless_ax.ylabel[] == "Intensity [unitless]"
    ext = Base.get_extension(JuChrom, :MakieExtension)
    @test ext.tictrace_axis_unit_label(nothing) == "unitless"

    layout_fig = Figure()
    axplot = tictrace(
        layout_fig[1, 1],
        msm;
        retentionunit=u"minute",
        intensityunit=u"nA",
        linewidth=4
    )
    @test axplot isa Makie.AxisPlot
    @test axplot.axis.xlabel[] == "Retention [minute]"
    @test axplot.axis.ylabel[] == "Intensity [nA]"
    @test axplot.plot.arg1[] ≈ [1.0, 2.0]
    @test axplot.plot.arg2[] ≈ [3.0, 7.0]
    @test axplot.plot.linewidth[] ≈ 4
    @test_throws ArgumentError tictrace(
        Figure()[1, 1],
        msm;
        figure=(; size=(320, 240))
    )
end

@testset "lines! plots TIC from JuChrom containers with unit conversion" begin
    X = [1000.0 2000.0; 3000.0 4000.0]
    msm = MassScanMatrix(
        [60.0, 120.0],
        u"s",
        [100.0, 101.0],
        nothing,
        X,
        u"pA"
    )

    fig = Figure()
    ax = Axis(fig[1, 1])
    plt = lines!(
        ax,
        msm;
        retentionunit=u"minute",
        intensityunit=u"nA",
        color=:red,
        linewidth=2
    )

    @test plt isa Makie.Lines
    @test plt.arg1[] ≈ [1.0, 2.0]
    @test plt.arg2[] ≈ [3.0, 7.0]
    @test plt.linewidth[] ≈ 2

    vmsm = VarianceMassScanMatrix(msm, ones(size(X)))
    plt_v = lines!(
        ax,
        vmsm;
        retentionunit=u"minute",
        intensityunit=u"nA"
    )
    @test plt_v isa Makie.Lines
    @test plt_v.arg1[] ≈ [1.0, 2.0]
    @test plt_v.arg2[] ≈ [3.0, 7.0]

    mss = MassScanSeries([
        MassScan(60.0u"s", [100.0, 101.0], [1000.0, 2000.0]u"pA"),
        MassScan(120.0u"s", [100.0, 101.0], [3000.0, 4000.0]u"pA")
    ])
    plt_mss = lines!(
        ax,
        mss;
        retentionunit=u"minute",
        intensityunit=u"nA"
    )
    @test plt_mss isa Makie.Lines
    @test plt_mss.arg1[] ≈ [1.0, 2.0]
    @test plt_mss.arg2[] ≈ [3.0, 7.0]

    css = ChromScanSeries([
        ChromScan(60.0u"s", 1000.0u"pA"),
        ChromScan(120.0u"s", 2000.0u"pA")
    ])
    plt_css = lines!(
        ax,
        css;
        retentionunit=u"minute",
        intensityunit=u"nA"
    )
    @test plt_css isa Makie.Lines
    @test plt_css.arg1[] ≈ [1.0, 2.0]
    @test plt_css.arg2[] ≈ [1.0, 2.0]

    fap = lines(
        msm;
        figure=(; size=(320, 240)),
        axis=(; xlabel="Retention [minute]", ylabel="TIC [nA]"),
        retentionunit=u"minute",
        intensityunit=u"nA",
        color=:blue,
        linewidth=3
    )
    @test fap isa Makie.FigureAxisPlot
    @test fap.plot isa Makie.Lines
    @test fap.figure.scene.viewport[].widths == [320, 240]
    @test fap.axis.xlabel[] == "Retention [minute]"
    @test fap.axis.ylabel[] == "TIC [nA]"
    @test fap.plot.arg1[] ≈ [1.0, 2.0]
    @test fap.plot.arg2[] ≈ [3.0, 7.0]
    @test fap.plot.linewidth[] ≈ 3

    layout_fig = Figure()
    axplot = lines(
        layout_fig[1, 1],
        msm;
        axis=(; xlabel="Retention [minute]"),
        retentionunit=u"minute",
        intensityunit=u"nA",
        linewidth=4
    )
    @test axplot isa Makie.AxisPlot
    @test axplot.axis.xlabel[] == "Retention [minute]"
    @test axplot.plot.arg1[] ≈ [1.0, 2.0]
    @test axplot.plot.arg2[] ≈ [3.0, 7.0]
    @test axplot.plot.linewidth[] ≈ 4
    @test_throws ArgumentError lines(
        Figure()[1, 1],
        msm;
        figure=(; size=(320, 240))
    )
end

@testset "tictrace! overlays ladder steps on an existing axis" begin
    msm, result = synthetic_ladder_result()
    fig = Figure()
    ax = Axis(fig[1, 1])

    out = tictrace!(
        ax,
        msm,
        result;
        title="Existing annotated ladder",
        axis=(; xlabel="RT")
    )

    @test out === ax
    @test length(ax.scene.plots) == 1
    @test ax.title[] == "Existing annotated ladder"
    @test ax.xlabel[] == "RT"
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
        result.success,
        result.status,
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
