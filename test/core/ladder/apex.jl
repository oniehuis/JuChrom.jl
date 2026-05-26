using Test
using JuChrom

const APEX_MZRETENTION_KWARGS = (
    retentionref=:start,
    scaninterval=0.2,
    mzcount=4,
    order=:ascending,
    dwellref=:middle,
    dwell=:homogeneous,
)

function apex_reference_standard()
    spectrum = JuChrom.alkane_reference_spectrum(
        8,
        "octane synthetic",
        800.0,
        [43, 57, 71, 85],
        [0.18, 1.0, 0.42, 0.09],
    )
    (spectra=[spectrum],)
end

function apex_test_msm(; apexretention=5.35, mode=:quadratic)
    retentions = collect(1.0:10.0)
    mzs = [43.0, 57.0, 71.0, 85.0]
    amplitudes = [90.0, 340.0, 170.0, 55.0]
    intensities = Matrix{Float64}(undef, length(retentions), length(mzs))
    for scanindex in eachindex(retentions), mzindex in eachindex(mzs)
        obsretention = mzretention(
            retentions[scanindex];
            APEX_MZRETENTION_KWARGS...,
            mzindex=mzindex,
        )
        intensities[scanindex, mzindex] = if mode == :quadratic
            amplitudes[mzindex] * exp(-0.85 * abs2(obsretention - apexretention))
        elseif mode == :monotone
            amplitudes[mzindex] * exp(0.15 * obsretention)
        else
            error("unsupported test mode")
        end
    end

    MassScanMatrix(retentions, mzs, intensities)
end

function apex_channelinfo(msm)
    JuChrom.alkane_mz_channels(
        msm;
        standard=apex_reference_standard(),
        carbonrange=8:8,
        minrelativeintensity=0.05,
    )
end

@testset "ladder peak apex uses all selected reference ions" begin
    msm = apex_test_msm(; apexretention=5.35)
    variances = fill(1.0, size(rawintensities(msm)))
    channelinfo = apex_channelinfo(msm)
    peak = (carbon=8, scanindex=5, retention=5.0, score=1.0)

    apex = JuChrom.ladder_peak_apex_twopass(
        msm,
        peak,
        variances;
        channelinfo=channelinfo,
        mzretentionkwargs=APEX_MZRETENTION_KWARGS,
        scanwindow=2,
        logfloorfraction=1e-9,
        variancefloor=1e-6,
    )

    @test apex.success
    @test apex.continuous_refinement
    @test !apex.fallback_to_snapshot
    @test apex.apex_retention ≈ 5.35 atol=0.02
    @test apex.apex_scan_index ≈ 5.35 atol=0.02
    @test apex.mz_indices == [1, 2, 3, 4]
    @test apex.mz_values == mzvalues(msm)
    @test apex.reference_intensities ≈ [0.18, 1.0, 0.42, 0.09]
    @test apex.variance_weighted
    @test JuChrom.ladder_peak_apex(
        msm,
        5,
        variances;
        carbon=8,
        channelinfo=channelinfo,
        mzretentionkwargs=APEX_MZRETENTION_KWARGS,
        logfloorfraction=1e-9,
        variancefloor=1e-6,
    ).apex_retention ≈ 5.35 atol=0.02
    @test_throws ArgumentError JuChrom.ladder_peak_apex(
        msm,
        5,
        variances;
        channelinfo=channelinfo,
        mzretentionkwargs=APEX_MZRETENTION_KWARGS,
    )
end

@testset "two-pass apex refinement shifts an off-center scan" begin
    msm = apex_test_msm(; apexretention=6.35)
    variances = fill(1.0, size(rawintensities(msm)))
    channelinfo = apex_channelinfo(msm)
    peak = (carbon=8, scanindex=5, retention=5.0, score=1.0)

    apex = JuChrom.ladder_peak_apex_twopass(
        msm,
        peak,
        variances;
        channelinfo=channelinfo,
        mzretentionkwargs=APEX_MZRETENTION_KWARGS,
        scanwindow=2,
        maxapexoffsetscans=1.0,
        logfloorfraction=1e-9,
        variancefloor=1e-6,
    )

    @test apex.success
    @test apex.shift_attempted
    @test apex.shift_direction == 1
    @test apex.shift_count >= 1
    @test apex.input_scan_index == 6
    @test apex.ladder_scan_index == 5
    @test apex.apex_retention ≈ 6.35 atol=0.02
    @test abs(apex.apex_scan_index - apex.input_scan_index) <
        abs(apex.initial_apex.apex_scan_index - apex.initial_apex.input_scan_index)
end

@testset "invalid continuous apex falls back to ladder snapshot" begin
    msm = apex_test_msm(; mode=:monotone)
    variances = fill(1.0, size(rawintensities(msm)))
    channelinfo = apex_channelinfo(msm)
    peak = (carbon=8, scanindex=5, retention=5.0, score=1.0)

    apex = JuChrom.ladder_peak_apex_twopass(
        msm,
        peak,
        variances;
        channelinfo=channelinfo,
        mzretentionkwargs=APEX_MZRETENTION_KWARGS,
        scanwindow=2,
        logfloorfraction=1e-9,
        variancefloor=1e-6,
    )

    @test apex.success
    @test !apex.continuous_refinement
    @test apex.fallback_to_snapshot
    @test apex.apex_retention == 5.0
    @test apex.apex_scan_index == 5.0
end

@testset "series apex refinement reports failed paths without fitting" begin
    msm = apex_test_msm()
    variances = fill(1.0, size(rawintensities(msm)))
    channelinfo = apex_channelinfo(msm)
    failedpath = (
        success=false,
        failurereason="no local maxima found in evidence traces",
        path=NamedTuple[],
    )

    apexes = JuChrom.refine_alkane_series_apexes(
        msm,
        variances,
        failedpath;
        channelinfo=channelinfo,
    )

    @test !apexes.success
    @test occursin("series path failed", apexes.failurereason)
    @test isempty(apexes.results)
    @test isempty(apexes.bycarbon)
end
