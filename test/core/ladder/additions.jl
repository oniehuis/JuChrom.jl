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

@testset "alkane ladder edge extension cap follows current anchors and carbon range" begin
    anchors = [
        (ladderstep=12,),
        (ladderstep=15,)
    ]

    @test JuChrom.alkane_ladder_edge_extension_step_cap(anchors, :left, 8:40) == 4
    @test JuChrom.alkane_ladder_edge_extension_step_cap(anchors, :right, 8:40) == 25
    @test JuChrom.alkane_ladder_edge_extension_step_cap(anchors, :left, 12:15) == 0
    @test JuChrom.alkane_ladder_edge_extension_step_cap(anchors, :right, 12:15) == 0
end

@testset "alkaneladderadditions proposes single-step gap fills" begin
    msm, variances, abundanceinfo, pathinfo = ladder_addition_test_inputs()
    apexinfo = (
        apexes=[
            (success=true, ladderstep=8, apexscanindex=10.0),
            (success=true, ladderstep=10, apexscanindex=30.0)
        ],
    )

    additions = alkaneladderadditions(
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

    @test additions.status == :success
    @test fill.ladderstep == 9
    @test fill.source == :gapfilled
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

    additions = alkaneladderadditions(
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

    additions = alkaneladderadditions(
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
    @test left.source == :leftextended
    @test left.scanindex == 10
    @test left.expectedscan == 10.0
    @test right.ladderstep == 11
    @test right.source == :rightextended
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

    additions = alkaneladderadditions(
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
    @test additions.diagnostics.rightextended[1].status == :accepted
    @test additions.diagnostics.rightextended[2].status == :accepted
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

    additions = alkaneladderadditions(
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
    @test additions.diagnostics.rightextended[1].status == :accepted
    @test additions.diagnostics.rightextended[2].status == :failed
end

@testset "alkaneladderadditions returns empty result without two refined anchors" begin
    msm, variances, abundanceinfo, pathinfo = ladder_addition_test_inputs()
    apexinfo = (
        apexes=[
            (success=true, ladderstep=8, apexscanindex=10.0),
            (success=false, ladderstep=9, apexscanindex=20.0)
        ],
    )

    additions = alkaneladderadditions(
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

    @test additions.status == :empty
    @test isempty(additions.additions)
    @test isempty(additions.gapfilled)
    @test isempty(additions.leftextended)
    @test isempty(additions.rightextended)
end

@testset "alkane ladder addition validation" begin
    @test JuChrom.validate_alkane_ladder_addition_settings(2.0, 0.5, 0.5) === nothing
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_addition_settings(
        -1.0,
        0.5,
        0.5
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_addition_settings(
        2.0,
        -0.5,
        0.5
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_addition_settings(
        2.0,
        0.5,
        0.0
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_addition_settings(
        2.0,
        0.5,
        0.5,
        -1
    )
end
