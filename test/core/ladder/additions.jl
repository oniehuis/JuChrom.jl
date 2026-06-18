using Test
using JuChrom

function test_ladder_addition_abundance(peaks)
    abundance = fill(0.01, 60)
    scans = collect(1.0:60.0)
    for (apex, height) in peaks
        abundance .+= height .* exp.(-0.5 .* abs2.((scans .- apex) ./ 1.5))
    end

    abundance
end

function ladder_addition_test_inputs()
    abundances = Dict(
        8 => test_ladder_addition_abundance([(10, 20.0)]),
        9 => test_ladder_addition_abundance([(20, 30.0), (26, 100.0)]),
        10 => test_ladder_addition_abundance([(30, 40.0)]),
        11 => test_ladder_addition_abundance([(40, 25.0)])
    )
    windows = Dict(
        8 => [
            AlkaneAbundanceWindow(
                8, 8, 10, 12, 0.01, 20.0, 0.01, 0.01, :threshold, :threshold
            )
        ],
        9 => [
            AlkaneAbundanceWindow(
                9, 18, 20, 22, 0.01, 30.0, 0.01, 0.01, :threshold, :threshold
            ),
            AlkaneAbundanceWindow(
                9, 24, 26, 28, 0.01, 100.0, 0.01, 0.01, :threshold, :threshold
            )
        ],
        10 => [
            AlkaneAbundanceWindow(
                10, 28, 30, 32, 0.01, 40.0, 0.01, 0.01, :threshold, :threshold
            )
        ],
        11 => [
            AlkaneAbundanceWindow(
                11, 38, 40, 42, 0.01, 25.0, 0.01, 0.01, :threshold, :threshold
            )
        ]
    )
    abundanceinfo = AlkaneAbundanceInfo(
        abundances,
        Dict(step => ones(length(abundance)) for (step, abundance) in abundances),
        windows,
        NamedTuple()
    )
    msm = ladder_addition_test_msm(abundanceinfo)
    variances = ones(60, 1)
    pathinfo = (path=NamedTuple[],)
    msm, variances, abundanceinfo, pathinfo
end

function ladder_addition_test_msm(abundanceinfo)
    X = zeros(60, 1)
    for abundance in values(abundanceinfo.abundances)
        X[:, 1] .+= abundance
    end

    MassScanMatrix(collect(1.0:60.0), [100.0], X)
end

function ladder_addition_test_mzretentionkwargs(msm)
    (
        retentionref=:middle,
        scaninterval=1e-9,
        mzcount=mzcount(msm),
        order=:ascending,
        dwellref=:middle,
        dwell=:homogeneous
    )
end

function ladder_addition_test_settings(
    msm,
    variances,
    abundanceinfo,
    apexinfo;
    minradius=5.0,
    radiusfraction=0.15,
    positionsigmafraction=0.05,
    maxextensionsteps=100,
    massspectrummatch=false,
    gapmincosinefloor=0.85,
    gapcosinetolerance=0.03,
    edgemaxanchors=6,
    edgeminradius=5.0,
    edgeradiusfraction=0.2,
    edgemincosinefloor=0.9,
    edgecosinetolerance=0.03,
    edgecosineanchorcount=3,
    edgepositionsigmafraction=0.1,
    massspectrumvariancefloor=1.0,
    apexscanwindow=2,
    apexvariancefloor=1.0,
    apexlogfloorfraction=1e-3,
    apexionexcludemzvalues=JuChrom.DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
    apexionmzvalues=[100.0],
    apexionminrelativeintensity=0.1,
    apexminioncount=1,
    apexmzretentionkwargs=ladder_addition_test_mzretentionkwargs(msm),
    apexmaxshiftfromguess=3.0,
    apexcenteredscantolerance=0.25,
    apexfitqualityoutlierz=3.0,
    apexfitqualityminsteps=6,
    carbonrange=nothing
)
    JuChrom.alkane_ladder_addition_settings(
        msm,
        variances,
        abundanceinfo,
        apexinfo,
        minradius,
        radiusfraction,
        positionsigmafraction,
        maxextensionsteps,
        massspectrummatch,
        gapmincosinefloor,
        gapcosinetolerance,
        edgemaxanchors,
        edgeminradius,
        edgeradiusfraction,
        edgemincosinefloor,
        edgecosinetolerance,
        edgecosineanchorcount,
        edgepositionsigmafraction,
        massspectrumvariancefloor,
        apexscanwindow,
        apexvariancefloor,
        apexlogfloorfraction,
        apexionexcludemzvalues,
        apexionmzvalues,
        apexionminrelativeintensity,
        apexminioncount,
        apexmzretentionkwargs,
        apexmaxshiftfromguess,
        apexcenteredscantolerance,
        apexfitqualityoutlierz,
        apexfitqualityminsteps,
        carbonrange
    )
end

function ladder_additions(
    msm,
    variances,
    abundanceinfo,
    pathinfo,
    apexinfo;
    standard=defaultalkanestandard(),
    kwargs...
)
    settings = ladder_addition_test_settings(
        msm,
        variances,
        abundanceinfo,
        apexinfo;
        kwargs...
    )

    alkaneladderadditions(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexinfo,
        settings,
        standard
    )
end

function validate_ladder_addition_settings(
    minradius,
    radiusfraction,
    positionsigmafraction,
    maxextensionsteps=100
)
    JuChrom.validate_alkane_ladder_addition_settings(
        minradius,
        radiusfraction,
        positionsigmafraction,
        maxextensionsteps,
        0.85,
        0.03,
        6,
        5.0,
        0.2,
        0.9,
        0.03,
        3,
        0.1,
        1.0
    )
end

@testset "alkane ladder edge extension cap follows current anchors and carbon range" begin
    anchors = [
        JuChrom.AlkaneLadderAdditionAnchor(12, 10.0, NaN, :molecularion, NaN, NaN),
        JuChrom.AlkaneLadderAdditionAnchor(15, 20.0, NaN, :molecularion, NaN, NaN)
    ]

    @test JuChrom.alkane_ladder_edge_extension_step_cap(anchors, :left, 8:40) == 4
    @test JuChrom.alkane_ladder_edge_extension_step_cap(anchors, :right, 8:40) == 25
    @test JuChrom.alkane_ladder_edge_extension_step_cap(anchors, :left, 12:15) == 0
    @test JuChrom.alkane_ladder_edge_extension_step_cap(anchors, :right, 12:15) == 0
    @test_throws ArgumentError JuChrom.alkane_ladder_edge_extension_step_cap(
        anchors,
        :bad,
        8:40
    )
end

@testset "alkaneladderadditions proposes single-step gap fills" begin
    msm, variances, abundanceinfo, pathinfo = ladder_addition_test_inputs()
    apexinfo = (
        apexes=[
            (success=true, ladderstep=8, apexscanindex=10.0),
            (success=true, ladderstep=10, apexscanindex=30.0)
        ],
    )

    additions = ladder_additions(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexinfo;
        minradius=2.0,
        radiusfraction=0.5,
        positionsigmafraction=0.5,
        massspectrummatch=false,
        apexionmzvalues=[100.0],
        apexminioncount=1,
        apexmzretentionkwargs=ladder_addition_test_mzretentionkwargs(msm)
    )
    fill = only(additions.gapfilled)

    @test additions.status ≡ :success
    @test fill.ladderstep == 9
    @test fill.source ≡ :gapfilled
    @test fill.scanindex == 20
    @test fill.expectedscan == 20.0
    @test fill.scanerror == 0.0
    @test fill.score == 30.0
    @test isempty(additions.leftextended)
    @test only(additions.rightextended).ladderstep == 11
end

@testset "alkaneladderadditions gates additions by mass-spectrum cosine" begin
    msm, variances, abundanceinfo, pathinfo = ladder_addition_test_inputs()
    X = rawintensities(msm)
    X[:, 1] .= abundanceinfo.abundances[9]
    msm = MassScanMatrix(rawretentions(msm), rawmzvalues(msm), X)
    standard = AlkaneStandard(
        "test C9",
        [JuChrom.alkane_reference_spectrum(9, "nonane", 900.0, [100], [100.0])],
        NamedTuple()
    )
    apexinfo = (
        apexes=[
            (success=true, ladderstep=8, apexscanindex=10.0),
            (success=true, ladderstep=10, apexscanindex=30.0)
        ],
    )

    additions = ladder_additions(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexinfo;
        standard=standard,
        massspectrummatch=true,
        gapmincosinefloor=0.99,
        minradius=2.0,
        radiusfraction=0.5,
        positionsigmafraction=0.5,
        apexionmzvalues=[100.0],
        apexminioncount=1,
        apexmzretentionkwargs=ladder_addition_test_mzretentionkwargs(msm)
    )
    fill = only(additions.gapfilled)

    @test fill.ladderstep == 9
    @test fill.massspectrumcosine ≈ 1.0
    @test fill.requiredcosine == 0.99
end

@testset "alkaneladderadditions proposes one-step edge extensions" begin
    msm, variances, abundanceinfo, pathinfo = ladder_addition_test_inputs()
    apexinfo = (
        apexes=[
            (success=true, ladderstep=9, apexscanindex=20.0),
            (success=true, ladderstep=10, apexscanindex=30.0)
        ],
    )

    additions = ladder_additions(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexinfo;
        minradius=2.0,
        radiusfraction=0.5,
        positionsigmafraction=0.5,
        massspectrummatch=false,
        apexionmzvalues=[100.0],
        apexminioncount=1,
        apexmzretentionkwargs=ladder_addition_test_mzretentionkwargs(msm)
    )
    left = only(additions.leftextended)
    right = only(additions.rightextended)

    @test left.ladderstep == 8
    @test left.source ≡ :leftextended
    @test left.scanindex == 10
    @test left.expectedscan == 10.0
    @test right.ladderstep == 11
    @test right.source ≡ :rightextended
    @test right.scanindex == 40
    @test right.expectedscan == 40.0
    @test [a.ladderstep for a in additions.additions] == [8, 11]
end

@testset "alkaneladderadditions iterates edge extensions" begin
    msm, variances, abundanceinfo, pathinfo = ladder_addition_test_inputs()
    abundanceinfo.abundances[12] = test_ladder_addition_abundance([(50, 15.0)])
    abundanceinfo.windows[12] = [
        AlkaneAbundanceWindow(
            12, 48, 50, 52, 0.01, 15.0, 0.01, 0.01, :threshold, :threshold
        )
    ]
    msm = ladder_addition_test_msm(abundanceinfo)
    apexinfo = (
        apexes=[
            (success=true, ladderstep=9, apexscanindex=20.0),
            (success=true, ladderstep=10, apexscanindex=30.0)
        ],
    )

    additions = ladder_additions(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexinfo;
        minradius=2.0,
        radiusfraction=0.5,
        positionsigmafraction=0.5,
        massspectrummatch=false,
        apexionmzvalues=[100.0],
        apexminioncount=1,
        apexmzretentionkwargs=ladder_addition_test_mzretentionkwargs(msm)
    )

    @test [a.ladderstep for a in additions.rightextended] == [11, 12]
    @test [a.source for a in additions.rightextended] ==
        [:rightextended, :rightextended]
    @test all(!hasproperty(a, :edgeanchor) for a in additions.rightextended)
    @test all(!hasproperty(a, :edgeanchors) for a in additions.rightextended)
    @test length(additions.diagnostics.rightextended) == 2
    @test additions.diagnostics.rightextended[1].status ≡ :accepted
    @test additions.diagnostics.rightextended[2].status ≡ :accepted
end

@testset "alkaneladderadditions stops extension when added apex fails" begin
    msm, variances, abundanceinfo, pathinfo = ladder_addition_test_inputs()
    abundanceinfo.abundances[12] = test_ladder_addition_abundance([(50, 15.0)])
    abundanceinfo.windows[12] = [
        AlkaneAbundanceWindow(
            12, 50, 50, 50, 0.01, 15.0, 0.01, 0.01, :threshold, :threshold
        )
    ]
    apexinfo = (
        apexes=[
            (success=true, ladderstep=9, apexscanindex=20.0),
            (success=true, ladderstep=10, apexscanindex=30.0)
        ],
    )

    additions = ladder_additions(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexinfo;
        minradius=2.0,
        radiusfraction=0.5,
        positionsigmafraction=0.5,
        massspectrummatch=false,
        apexionmzvalues=[100.0],
        apexminioncount=1,
        apexmzretentionkwargs=ladder_addition_test_mzretentionkwargs(msm)
    )

    @test [a.ladderstep for a in additions.rightextended] == [11]
    @test length(additions.diagnostics.rightextended) == 2
    @test additions.diagnostics.rightextended[1].status ≡ :accepted
    @test additions.diagnostics.rightextended[2].status ≡ :failed
end

@testset "alkaneladderadditions returns empty result without two refined anchors" begin
    msm, variances, abundanceinfo, pathinfo = ladder_addition_test_inputs()
    apexinfo = (
        apexes=[
            (success=true, ladderstep=8, apexscanindex=10.0),
            (success=false, ladderstep=9, apexscanindex=20.0)
        ],
    )

    additions = ladder_additions(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        apexinfo;
        massspectrummatch=false,
        apexionmzvalues=[100.0],
        apexminioncount=1,
        apexmzretentionkwargs=ladder_addition_test_mzretentionkwargs(msm)
    )

    @test additions.status ≡ :empty
    @test isempty(additions.additions)
    @test isempty(additions.gapfilled)
    @test isempty(additions.leftextended)
    @test isempty(additions.rightextended)
end

@testset "alkane ladder addition helper coverage" begin
    window = AlkaneAbundanceWindow(
        9, 5, 10, 15, 0.0, 10.0, 0.0, 0.0, :threshold, :threshold
    )
    candidate = JuChrom.AlkaneLadderAdditionCandidate(
        9,
        5,
        10,
        15,
        collect(5:15),
        ones(11),
        :gapfilled,
        10,
        10.0,
        0.0,
        5.0,
        1.0,
        10.0,
        10.0,
        10.0,
        NaN,
        NaN,
        0,
        NaN,
        window
    )

    @test JuChrom.alkane_ladder_fallback_scanindices(
        collect(1.0:30.0),
        10,
        8,
        candidate
    ) == collect(5:15)

    inferredkwargs = (order=:ascending,)
    apexinfo = (settings=(mzretentionkwargs=inferredkwargs,),)
    @test JuChrom.alkane_ladder_addition_apex_mzretentionkwargs(
        apexinfo,
        nothing
    ) == inferredkwargs
    @test JuChrom.alkane_ladder_addition_apex_mzretentionkwargs(
        apexinfo,
        (order=:descending,)
    ) == (order=:descending,)

    anchors = [
        JuChrom.AlkaneLadderAdditionAnchor(8, 10.0, 10.0, :molecularion, 0.95, 0.9),
        JuChrom.AlkaneLadderAdditionAnchor(9, 22.0, 22.0, :molecularion, 0.90, 0.9),
        JuChrom.AlkaneLadderAdditionAnchor(10, 40.0, 40.0, :molecularion, 0.85, 0.9)
    ]
    @test JuChrom.alkane_ladder_edge_extension_anchor_results(
        anchors,
        :left,
        2
    ) == anchors[1:2]
    @test JuChrom.alkane_ladder_edge_extension_anchor_results(
        anchors,
        :right,
        2
    ) == anchors[2:3]
    @test_throws ArgumentError JuChrom.alkane_ladder_edge_extension_anchor_results(
        anchors,
        :bad,
        2
    )

    @test isnan(JuChrom.alkane_ladder_anchor_mass_spectrum_cosine(nothing))
    @test JuChrom.alkane_single_gap_required_cosine(
        anchors[1],
        anchors[2],
        0.8,
        0.03
    ) ≈ 0.87
    @test JuChrom.alkane_single_gap_required_cosine(
        nothing,
        nothing,
        0.8,
        0.03
    ) == 0.8
    @test JuChrom.alkane_edge_extension_required_cosine(
        anchors,
        :right,
        0.8,
        0.03,
        2
    ) ≈ 0.82

    rightprediction = JuChrom.alkane_ladder_edge_extension_scan_prediction(
        anchors,
        11
    )
    @test rightprediction.expectedscan ≈ 67.0
    @test rightprediction.localstepgap ≈ 27.0

    leftprediction = JuChrom.alkane_ladder_edge_extension_scan_prediction(
        anchors,
        7
    )
    @test leftprediction.expectedscan ≈ 2.0
    @test leftprediction.localstepgap ≈ 8.0

    @test JuChrom.alkane_ladder_edge_extension_scan_prediction(anchors, 9) ≡
        nothing

    boundarywindow = AlkaneAbundanceWindow(
        8, 1, 1, 3, 0.0, 10.0, 0.0, 0.0, :boundary, :threshold
    )
    @test JuChrom.alkane_ladder_edge_window_apex_is_boundary_truncated(
        [10.0, 5.0, 1.0],
        boundarywindow,
        :left,
        3,
        3.0
    )
    @test_throws ArgumentError JuChrom.alkane_ladder_edge_window_apex_is_boundary_truncated(
        [10.0, 5.0, 1.0],
        boundarywindow,
        :bad,
        3,
        3.0
    )
end

@testset "alkane ladder edge addition diagnostics cover rejected windows" begin
    msm, variances, abundanceinfo, pathinfo = ladder_addition_test_inputs()
    apexinfo = (
        apexes=[
            (success=true, ladderstep=8, apexscanindex=10.0),
            (success=true, ladderstep=10, apexscanindex=30.0)
        ],
    )
    settings = ladder_addition_test_settings(
        msm,
        variances,
        abundanceinfo,
        apexinfo;
        apexionmzvalues=[100.0],
        apexminioncount=1,
        apexmzretentionkwargs=ladder_addition_test_mzretentionkwargs(msm)
    )
    settings_with_one_edge_anchor = JuChrom.AlkaneLadderAdditionSettings(
        settings.minradius,
        settings.radiusfraction,
        settings.positionsigmafraction,
        settings.maxextensionsteps,
        settings.massspectrummatch,
        settings.gapmincosinefloor,
        settings.gapcosinetolerance,
        1,
        settings.edgeminradius,
        settings.edgeradiusfraction,
        settings.edgemincosinefloor,
        settings.edgecosinetolerance,
        settings.edgecosineanchorcount,
        settings.edgepositionsigmafraction,
        settings.massspectrumvariancefloor,
        settings.apexscanwindow,
        settings.apexvariancefloor,
        settings.apexlogfloorfraction,
        settings.apexionexcludemzvalues,
        settings.apexionmzvalues,
        settings.apexionminrelativeintensity,
        settings.apexminioncount,
        settings.apexmzretentionkwargs,
        settings.apexmaxshiftfromguess,
        settings.apexcenteredscantolerance,
        settings.apexfitqualityoutlierz,
        settings.apexfitqualityminsteps,
        settings.carbonrange
    )
    anchors = [
        JuChrom.AlkaneLadderAdditionAnchor(8, 10.0, 10.0, :molecularion, NaN, NaN),
        JuChrom.AlkaneLadderAdditionAnchor(10, 30.0, 30.0, :molecularion, NaN, NaN)
    ]
    leftadditions, leftdiagnostic = JuChrom.alkane_edge_addition(
        msm,
        variances,
        abundanceinfo,
        anchors,
        :left,
        scancount(msm),
        settings_with_one_edge_anchor,
        nothing,
        defaultalkanestandard()
    )
    rightadditions, rightdiagnostic = JuChrom.alkane_edge_addition(
        msm,
        variances,
        abundanceinfo,
        anchors,
        :right,
        scancount(msm),
        settings_with_one_edge_anchor,
        nothing,
        defaultalkanestandard()
    )

    @test isempty(leftadditions)
    @test isempty(rightadditions)
    @test leftdiagnostic.reason ≡ :insufficient_refined_anchors
    @test leftdiagnostic.ladderstep == 7
    @test rightdiagnostic.reason ≡ :insufficient_refined_anchors
    @test rightdiagnostic.ladderstep == 11
    @test_throws ArgumentError JuChrom.alkane_edge_addition(
        msm,
        variances,
        abundanceinfo,
        anchors,
        :bad,
        scancount(msm),
        settings,
        nothing,
        defaultalkanestandard()
    )

    truncatedabundances = Dict(7 => [10.0, 5.0, 1.0, zeros(57)...])
    truncatedwindow = AlkaneAbundanceWindow(
        7, 1, 1, 3, 0.0, 10.0, 0.0, 0.0, :boundary, :threshold
    )
    truncatedinfo = AlkaneAbundanceInfo(
        truncatedabundances,
        Dict(7 => ones(60)),
        Dict(7 => [truncatedwindow]),
        NamedTuple()
    )
    truncatedevaluation = JuChrom.alkane_best_ladder_addition_window(
        msm,
        variances,
        truncatedinfo,
        7,
        :leftextended,
        1.0,
        10.0,
        60,
        settings,
        nothing,
        defaultalkanestandard(),
        :left,
        nothing,
        nothing,
        anchors[1],
        anchors
    )
    @test isnothing(truncatedevaluation.addition)
    @test truncatedevaluation.diagnostic.reason ≡
        :abundance_window_truncated_by_run_boundary

    mismatchabundance = zeros(60)
    mismatchabundance[6] = 1.0
    mismatchabundance[7] = 10.0
    mismatchabundance[8] = 1.0
    mismatchwindow = AlkaneAbundanceWindow(
        11, 6, 6, 8, 0.0, 1.0, 0.0, 0.0, :threshold, :threshold
    )
    mismatchanchor = JuChrom.AlkaneLadderAdditionAnchor(
        10,
        5.0,
        5.0,
        :molecularion,
        NaN,
        NaN
    )
    mismatchinfo = AlkaneAbundanceInfo(
        Dict(11 => mismatchabundance),
        Dict(11 => ones(60)),
        Dict(11 => [mismatchwindow]),
        NamedTuple()
    )
    mismatchevaluation = JuChrom.alkane_best_ladder_addition_window(
        msm,
        variances,
        mismatchinfo,
        11,
        :rightextended,
        6.0,
        10.0,
        60,
        settings,
        nothing,
        defaultalkanestandard(),
        :right,
        nothing,
        nothing,
        mismatchanchor,
        [mismatchanchor, anchors[2]]
    )
    @test isnothing(mismatchevaluation.addition)
    @test mismatchevaluation.diagnostic.reason ≡ :abundance_window_apex_mismatch
end

@testset "alkane ladder addition validation" begin
    @test validate_ladder_addition_settings(2.0, 0.5, 0.5) ≡ nothing
    @test_throws ArgumentError validate_ladder_addition_settings(
        -1.0,
        0.5,
        0.5
    )
    @test_throws ArgumentError validate_ladder_addition_settings(
        2.0,
        -0.5,
        0.5
    )
    @test_throws ArgumentError validate_ladder_addition_settings(
        2.0,
        0.5,
        0.0
    )
    @test_throws ArgumentError validate_ladder_addition_settings(
        2.0,
        0.5,
        0.5,
        -1
    )
end
