using Test
using JuChrom

function abundance_test_inputs()
    X = [
        20.0 0.0 10.0 0.0
        12.0 0.0 3.0 0.0
        -4.0 0.0 -2.0 0.0
    ]
    variance = [
        1.0 10.0 1.0 10.0
        1.0 10.0 4.0 10.0
        2.0 10.0 2.0 10.0
    ]
    mzindices = [1, 3]
    referenceintensities = [2.0, 1.0]

    X, variance, mzindices, referenceintensities
end

function abundance_test_channelinfo()
    references = [
        AlkaneReferenceChannels(8, NamedTuple(), [1, 3], [29.0, 57.0], [2.0, 1.0]),
        AlkaneReferenceChannels(9, NamedTuple(), [2], [43.0], [4.0])
    ]
    AlkaneChannelInfo(
        [1, 2, 3],
        [29.0, 43.0, 57.0],
        references,
        [8, 9],
        0.0,
        nothing
    )
end

@testset "alkane_mz_channels returns concrete channel info" begin
    X = zeros(2, 3)
    msm = MassScanMatrix([1.0, 2.0], [29.0, 43.0, 57.0], X)
    standard = AlkaneStandard(
        "test",
        [
            JuChrom.alkane_reference_spectrum(
                8,
                "octane",
                800.0,
                [29, 57],
                [10.0, 5.0]
            )
        ],
        NamedTuple()
    )

    channelinfo = JuChrom.alkane_mz_channels(
        msm;
        standard=standard,
        carbonrange=8:8
    )

    @test channelinfo isa AlkaneChannelInfo
    @test channelinfo.mzindices == [1, 3]
    @test channelinfo.mzvalues == [29.0, 57.0]
    @test channelinfo.carbonrange == [8]
    @test channelinfo.minrelativeintensity == 0.0
    @test only(channelinfo.references) isa AlkaneReferenceChannels
    @test only(channelinfo.references).carbon == 8
    @test only(channelinfo.references).mzindices == [1, 3]
    @test only(channelinfo.references).referenceintensities == [10.0, 5.0]

    intmsm = MassScanMatrix([1.0, 2.0], [29, 43, 57], X)
    intchannelinfo = JuChrom.alkane_mz_channels(
        intmsm;
        standard=standard,
        carbonrange=8:8
    )

    @test intchannelinfo isa AlkaneChannelInfo
    @test intchannelinfo.mzindices == [1, 3]
    @test intchannelinfo.mzvalues == [29, 57]
    @test only(intchannelinfo.references) isa AlkaneReferenceChannels
    @test only(intchannelinfo.references).mzvalues == [29, 57]
end

@testset "alkane_abundance weighted least-squares estimate" begin
    X, variance, mzindices, referenceintensities = abundance_test_inputs()

    abundance = JuChrom.alkane_abundance(
        X,
        variance,
        mzindices,
        referenceintensities,
        1.0,
        false,
        8
    )

    @test abundance ≈ [10.0, 99 / 17, -2.0]

    floored = JuChrom.alkane_abundance(
        X,
        variance,
        mzindices,
        referenceintensities,
        4.0,
        false,
        8
    )
    @test floored ≈ [10.0, 27 / 5, -2.0]

    clipped = JuChrom.alkane_abundance(
        X,
        variance,
        mzindices,
        referenceintensities,
        1.0,
        true,
        8
    )
    @test clipped ≈ [10.0, 99 / 17, 0.0]

    float32_reference = JuChrom.alkane_abundance(
        X,
        variance,
        mzindices,
        Float32[2.0, 1.0],
        1.0,
        false,
        8
    )
    @test float32_reference ≈ abundance
end

@testset "alkane_abundance_variance formula" begin
    _, variance, mzindices, referenceintensities = abundance_test_inputs()

    abundancevariance = JuChrom.alkane_abundance_variance(
        variance,
        mzindices,
        referenceintensities,
        1.0,
        8
    )

    @test abundancevariance ≈ [1 / 5, 4 / 17, 2 / 5]

    floored = JuChrom.alkane_abundance_variance(
        variance,
        mzindices,
        referenceintensities,
        4.0,
        8
    )
    @test floored ≈ fill(4 / 5, 3)
end

@testset "alkane_abundance_fit returns abundance and variance in one pass" begin
    X, variance, mzindices, referenceintensities = abundance_test_inputs()

    fit = JuChrom.alkane_abundance_fit(
        X,
        variance,
        mzindices,
        referenceintensities,
        1.0,
        false,
        8
    )

    @test fit isa JuChrom.AlkaneAbundanceFit
    @test fit.abundance ≈ [10.0, 99 / 17, -2.0]
    @test fit.abundancevariance ≈ [1 / 5, 4 / 17, 2 / 5]
end

@testset "alkane_abundances and alkane_abundance_variances wrappers" begin
    X, variance, mzindices, referenceintensities = abundance_test_inputs()
    msm = MassScanMatrix([1.0, 2.0, 3.0], [29.0, 43.0, 57.0, 71.0], X)
    channelinfo = abundance_test_channelinfo()

    abundances = JuChrom.alkane_abundances(
        msm,
        variance,
        channelinfo,
        1.0,
        false
    )
    abundancevariances = JuChrom.alkane_abundance_variances(
        variance,
        channelinfo,
        1.0
    )

    @test sort(collect(keys(abundances))) == [8, 9]
    @test abundances[8] ≈ JuChrom.alkane_abundance(
        X,
        variance,
        mzindices,
        referenceintensities,
        1.0,
        false,
        8
    )
    @test abundances[9] == zeros(3)
    @test abundancevariances[8] ≈ [1 / 5, 4 / 17, 2 / 5]
    @test abundancevariances[9] ≈ fill(10 / 16, 3)

    @test_throws DimensionMismatch JuChrom.alkane_abundances(
        msm,
        ones(2, 4),
        channelinfo,
        1.0,
        false
    )
end

@testset "alkane_abundance_info combines abundance pipeline" begin
    X, variance, _, _ = abundance_test_inputs()
    msm = MassScanMatrix([1.0, 2.0, 3.0], [29.0, 43.0, 57.0, 71.0], X)
    channelinfo = abundance_test_channelinfo()

    abundanceinfo = JuChrom.alkane_abundance_info(
        msm,
        variance,
        channelinfo,
        4.0,
        true,
        0.25,
        5.0
    )

    @test abundanceinfo isa AlkaneAbundanceInfo
    @test sort(collect(keys(abundanceinfo.abundances))) == [8, 9]
    @test sort(collect(keys(abundanceinfo.abundancevariances))) == [8, 9]
    @test sort(collect(keys(abundanceinfo.windows))) == [8, 9]
    @test abundanceinfo.abundances[8] ≈ [10.0, 27 / 5, 0.0]
    @test abundanceinfo.abundancevariances[8] ≈ fill(4 / 5, 3)
    @test abundanceinfo.settings == (
        variancefloor=4.0,
        nonnegative=true,
        thresholdfraction=0.25,
        minrisez=5.0
    )
end

@testset "alkane reference abundance intensity validation" begin
    reference = AlkaneReferenceChannels(
        8,
        NamedTuple(),
        [1, 2],
        [29.0, 43.0],
        [1.0, 0.5]
    )
    intensities = JuChrom.alkane_reference_abundance_intensities(reference)

    @test intensities == [1.0, 0.5]
    @test eltype(intensities) == Float64
    @test_throws ArgumentError JuChrom.alkane_reference_abundance_intensities(
        AlkaneReferenceChannels(
            8,
            NamedTuple(),
            [1, 2],
            [29.0, 43.0],
            [1.0, Inf]
        )
    )
    @test_throws ArgumentError JuChrom.alkane_reference_abundance_intensities(
        AlkaneReferenceChannels(
            8,
            NamedTuple(),
            [1, 2],
            [29.0, 43.0],
            [1.0, -0.5]
        )
    )
end

@testset "abundance input validation helpers" begin
    @test JuChrom.validate_alkane_abundance_variancefloor(1.0) ≡ nothing
    @test_throws ArgumentError JuChrom.validate_alkane_abundance_variancefloor(0.0)
    @test_throws ArgumentError JuChrom.validate_alkane_abundance_variancefloor(Inf)

    @test JuChrom.validate_alkane_abundance_variances_matrix([1.0 2.0]) ≡ nothing
    @test_throws ArgumentError JuChrom.validate_alkane_abundance_variances_matrix([1.0 Inf])
    @test_throws ArgumentError JuChrom.validate_alkane_abundance_variances_matrix([1.0 -1.0])
    @test_throws MethodError JuChrom.validate_alkane_abundance_variances_matrix([1.0, 2.0])

    @test JuChrom.validate_alkane_reference_abundance_channels(
        [1, 2],
        [1.0, 0.5],
        2,
        8
    ) ≡ nothing
    @test_throws DimensionMismatch JuChrom.validate_alkane_reference_abundance_channels(
        [1],
        [1.0, 0.5],
        2,
        8
    )
    @test_throws ArgumentError JuChrom.validate_alkane_reference_abundance_channels(
        Int[],
        Float64[],
        2,
        8
    )
    @test_throws MethodError JuChrom.validate_alkane_reference_abundance_channels(
        [1.5],
        [1.0],
        2,
        8
    )
    @test_throws ArgumentError JuChrom.validate_alkane_reference_abundance_channels(
        [3],
        [1.0],
        2,
        8
    )
    @test_throws ArgumentError JuChrom.validate_alkane_reference_abundance_channels(
        [1, 2],
        [0.0, 0.0],
        2,
        8
    )
    @test_throws ArgumentError JuChrom.validate_alkane_reference_abundance_channels(
        [1, 2],
        [1.0, Inf],
        2,
        8
    )
    @test_throws ArgumentError JuChrom.validate_alkane_reference_abundance_channels(
        [1, 2],
        [1.0, -1.0],
        2,
        8
    )

    X, variance, _, referenceintensities = abundance_test_inputs()
    @test_throws DimensionMismatch JuChrom.alkane_abundance(
        X,
        variance[:, 1:3],
        [1, 3],
        referenceintensities,
        1.0,
        false,
        8
    )
end

@testset "alkane abundance window validation helpers" begin
    @test JuChrom.validate_alkane_abundance_window_settings(0.5, 10.0) ≡ nothing
    @test_throws ArgumentError JuChrom.validate_alkane_abundance_window_settings(-0.1, 10.0)
    @test_throws ArgumentError JuChrom.validate_alkane_abundance_window_settings(1.0, 10.0)
    @test_throws ArgumentError JuChrom.validate_alkane_abundance_window_settings(0.5, -1.0)

    @test JuChrom.alkane_abundance_values([1, 2, 3], 8) == [1.0, 2.0, 3.0]
    @test_throws ArgumentError JuChrom.alkane_abundance_values([1.0, Inf], 8)
    @test_throws MethodError JuChrom.alkane_abundance_values(["bad"], 8)

    variance = Dict(8 => [0.1, 0.2, 0.3])
    @test JuChrom.alkane_abundance_window_variances(variance, 8, [1.0, 2.0, 3.0]) ==
        [0.1, 0.2, 0.3]
    @test_throws ArgumentError JuChrom.alkane_abundance_window_variances(
        Dict{Int, Vector{Float64}}(),
        8,
        [1.0]
    )
    @test_throws DimensionMismatch JuChrom.alkane_abundance_window_variances(
        Dict(8 => [0.1, 0.2]),
        8,
        [1.0]
    )
    @test_throws ArgumentError JuChrom.alkane_abundance_window_variances(
        Dict(8 => [0.1, Inf]),
        8,
        [1.0, 2.0]
    )
    @test_throws ArgumentError JuChrom.alkane_abundance_window_variances(
        Dict(8 => [0.1, -0.2]),
        8,
        [1.0, 2.0]
    )

    windows = [
        AlkaneAbundanceWindow(8, 2, 3, 4, 0.0, 1.0, 0.0, 0.0, :test, :test),
        AlkaneAbundanceWindow(8, 7, 8, 9, 0.0, 1.0, 0.0, 0.0, :test, :test)
    ]
    @test JuChrom.alkane_abundance_index_is_in_window(3, windows)
    @test !JuChrom.alkane_abundance_index_is_in_window(5, windows)
end

@testset "alkane_abundance_windows peak boundaries" begin
    abundance = [0.0, 1.0, 0.1, 2.0, 0.2, 0.0]
    windows = JuChrom.alkane_abundance_windows(
        Dict(8 => abundance),
        nothing,
        0.25,
        10.0
    )

    @test length(windows[8]) == 2
    @test windows[8][1] isa AlkaneAbundanceWindow
    @test windows[8][1].ladderstep == 8
    @test windows[8][1].leftindex == 3
    @test windows[8][1].apexindex == 4
    @test windows[8][1].rightindex == 5
    @test windows[8][1].leftabundance == 0.1
    @test windows[8][1].apexabundance == 2.0
    @test windows[8][1].rightabundance == 0.2
    @test windows[8][1].threshold == 0.5
    @test windows[8][1].leftstop ≡ :localminimum
    @test windows[8][1].rightstop ≡ :threshold
    @test windows[8][2].apexindex == 2
    @test windows[8][2].leftindex == 1
    @test windows[8][2].rightindex == 3

    variance_limited = JuChrom.alkane_abundance_windows(
        Dict(8 => abundance),
        Dict(8 => fill(0.01, length(abundance))),
        0.01,
        10.0
    )
    @test length(variance_limited[8]) == 1
    @test variance_limited[8][1].leftindex == 1
    @test variance_limited[8][1].rightindex == 6
    @test_throws MethodError JuChrom.alkane_abundance_windows(
        Dict(8 => abundance),
        [0.1, 0.2],
        0.05,
        10.0
    )
end

@testset "alkane abundance peak-window helpers" begin
    abundance = [0.0, 1.0, 0.1, 2.0, 0.2, 0.0]
    maxima = JuChrom.localmaxima(abundance)
    minima = Set(JuChrom.localmaxima(-abundance))

    peakwindow = JuChrom.alkane_abundance_peak_window(
        8,
        abundance,
        minima,
        maxima,
        4,
        0.25,
        nothing,
        10.0
    )
    @test peakwindow isa AlkaneAbundanceWindow
    @test peakwindow.leftindex == 3
    @test peakwindow.rightindex == 5
    @test peakwindow.leftstop ≡ :localminimum
    @test peakwindow.rightstop ≡ :threshold

    leftindex, leftstop = JuChrom.alkane_abundance_peak_window_side(
        abundance,
        nothing,
        minima,
        maxima,
        4,
        0.5,
        10.0,
        :left
    )
    rightindex, rightstop = JuChrom.alkane_abundance_peak_window_side(
        abundance,
        nothing,
        minima,
        maxima,
        4,
        0.5,
        10.0,
        :right
    )

    @test (leftindex, leftstop) == (3, :localminimum)
    @test (rightindex, rightstop) == (5, :threshold)

    @test JuChrom.alkane_abundance_peak_window_side(
        [1.0, 2.0, 1.0],
        nothing,
        Set{Int}(),
        [2],
        2,
        0.0,
        10.0,
        :left
    ) == (1, :boundary)
    @test JuChrom.alkane_abundance_peak_window_side(
        [1.0, 2.0, 1.0],
        nothing,
        Set{Int}(),
        [2],
        2,
        0.0,
        10.0,
        :right
    ) == (3, :boundary)
    @test_throws ArgumentError JuChrom.alkane_abundance_peak_window_side(
        abundance,
        nothing,
        minima,
        maxima,
        4,
        0.5,
        10.0,
        :bad
    )
end

@testset "alkane abundance local-minimum stop logic" begin
    abundance = [0.0, 1.0, 0.1, 2.0, 0.2, 0.0]
    maxima = [2, 4]

    @test JuChrom.alkane_abundance_local_minimum_stops_window(
        abundance,
        nothing,
        maxima,
        4,
        3,
        :left,
        10.0
    )
    @test !JuChrom.alkane_abundance_local_minimum_stops_window(
        abundance,
        fill(0.01, length(abundance)),
        [4],
        4,
        3,
        :left,
        10.0
    )
    @test !JuChrom.alkane_abundance_local_minimum_stops_window(
        [0.0, 1.0, 1.0, 2.0],
        fill(0.01, 4),
        [2, 4],
        4,
        3,
        :left,
        10.0
    )
    @test JuChrom.alkane_abundance_local_minimum_stops_window(
        abundance,
        zeros(length(abundance)),
        maxima,
        4,
        3,
        :left,
        10.0
    )
    @test JuChrom.alkane_abundance_local_minimum_stops_window(
        abundance,
        fill(0.01, length(abundance)),
        maxima,
        4,
        3,
        :left,
        5.0
    )
    @test !JuChrom.alkane_abundance_local_minimum_stops_window(
        abundance,
        fill(0.01, length(abundance)),
        maxima,
        4,
        3,
        :left,
        10.0
    )

    @test JuChrom.alkane_neighboring_peak_across_minimum(maxima, 4, 3, :left) == 2
    @test JuChrom.alkane_neighboring_peak_across_minimum(maxima, 2, 3, :right) == 4
    @test JuChrom.alkane_neighboring_peak_across_minimum([4], 4, 3, :left) ≡ nothing
    @test_throws ArgumentError JuChrom.alkane_neighboring_peak_across_minimum(
        maxima,
        4,
        3,
        :bad
    )
end
