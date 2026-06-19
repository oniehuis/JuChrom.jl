using Test
using JuChrom
using Unitful

function test_findalkanes_msm()
    n = 25
    retentions = collect(1.0:n)
    peak = exp.(-0.5 .* abs2.((retentions .- 13) ./ 2))
    noise = mod.(collect(1:n), 4)
    X = zeros(n, 5)
    X[:, 1] .= 5 .+ noise
    X[:, 2] .= 20 .+ noise .+ 100 .* peak
    X[:, 3] .= 10 .+ noise .+ 20 .* peak
    X[:, 4] .= 5 .+ noise .+ 5 .* peak
    X[:, 5] .= 3 .+ noise

    MassScanMatrix(
        retentions,
        [100.0, 114.0, 115.0, 128.0, 129.0],
        X
    )
end

function test_findalkanes_standard()
    AlkaneStandard(
        "test C8",
        [JuChrom.alkane_reference_spectrum(8, "octane", 800.0, [114], [100.0])],
        NamedTuple()
    )
end

function test_findalkanes_mzretentionkwargs(msm)
    (
        retentionref=:middle,
        scaninterval=1e-9,
        mzcount=mzcount(msm),
        dwellref=:middle,
        dwell=:homogeneous
    )
end

function test_ladder_apex_settings()
    JuChrom.AlkaneLadderApexSettings(
        defaultalkanestandard(),
        2,
        1.0,
        1e-3,
        (
            retentionref=:middle,
            scaninterval=1e-9,
            mzcount=1,
            order=:ascending,
            dwellref=:middle,
            dwell=:homogeneous
        ),
        :ascending,
        JuChrom.DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
        nothing,
        0.1,
        1,
        3.0,
        0.25,
        3.0,
        6,
        3,
        5,
        2.0,
        3,
        14,
        5,
        8
    )
end

function test_ladder_apex(c, scan; source=:molecularion, cosine=NaN, required=NaN,
        good=true, success=true)
    retention = Float64(scan) * 10.0
    JuChrom.AlkaneLadderApex(
        success,
        success ? :success : :failed,
        success ? nothing : "failed",
        success,
        !success,
        c,
        scan,
        retention,
        Float64(scan),
        retention,
        Float64(scan),
        retention,
        Float64(scan),
        retention,
        0.0,
        0.0,
        1.0,
        [scan],
        [1.0],
        1.0,
        0.0,
        Float64[],
        1e-3,
        source,
        source ≡ :gapfilled,
        source in (:leftextended, :rightextended),
        Float64(cosine),
        Float64(required),
        1.0,
        0.0,
        0.0,
        false,
        3.0,
        NaN,
        NaN,
        0,
        false,
        false,
        nothing,
        good,
        (ladderstep=c, scanindex=scan),
        nothing
    )
end

function test_ladder_apex_info(apexes)
    settings = test_ladder_apex_settings()
    JuChrom.AlkaneLadderApexInfo(
        :success,
        :success,
        apexes,
        Dict(apex.ladderstep => apex for apex in apexes),
        settings,
        JuChrom.alkane_ladder_no_scan_order_info(:test),
        [apex.calibration_excluded for apex in apexes],
        [apex.good_for_calibration for apex in apexes],
        [apex.apex_fit_quality_score for apex in apexes],
        [apex.apex_fit_quality_zscore for apex in apexes]
    )
end

function test_ladder_addition_settings()
    JuChrom.AlkaneLadderAdditionSettings(
        5.0,
        0.15,
        0.05,
        100,
        false,
        0.85,
        0.03,
        6,
        5.0,
        0.2,
        0.9,
        0.03,
        3,
        0.1,
        1.0,
        2,
        1.0,
        1e-3,
        JuChrom.DEFAULT_ALKANE_APEX_EXCLUDED_MZVALUES,
        nothing,
        0.1,
        1,
        test_ladder_apex_settings().mzretentionkwargs,
        3.0,
        0.25,
        3.0,
        6,
        collect(8:12)
    )
end

function test_ladder_addition(apex)
    scan = Int(round(apex.apexscanindex))
    window = AlkaneAbundanceWindow(
        apex.ladderstep,
        scan,
        scan,
        scan,
        1.0,
        1.0,
        1.0,
        0.05,
        :threshold,
        :threshold
    )
    JuChrom.AlkaneLadderAddition(
        apex.ladderstep,
        scan,
        scan,
        scan,
        [scan],
        [1.0],
        apex.source,
        scan,
        Float64(scan),
        0.0,
        5.0,
        1.0,
        10.0,
        1.0,
        1.0,
        apex.mass_spectrum_cosine,
        1.0 - apex.mass_spectrum_cosine,
        1,
        apex.required_cosine,
        window,
        apex,
        apex.success,
        apex.apexscanindex,
        apex.apexretention,
        apex.reason
    )
end

function test_ladder_addition_info(gapfilled, leftextended, rightextended)
    additions = vcat(gapfilled, rightextended)
    T = eltype(additions)
    typed_gapfilled = T[gapfilled...]
    typed_leftextended = T[leftextended...]
    typed_rightextended = T[rightextended...]
    JuChrom.AlkaneLadderAdditionInfo(
        :success,
        additions,
        typed_gapfilled,
        typed_leftextended,
        typed_rightextended,
        JuChrom.AlkaneLadderAdditionDiagnostics(
            JuChrom.AlkaneLadderAdditionDiagnostic[],
            JuChrom.AlkaneLadderAdditionDiagnostic[],
            JuChrom.AlkaneLadderAdditionDiagnostic[]
        ),
        test_ladder_addition_settings()
    )
end

function test_ladder_step_result(; standard=defaultalkanestandard(), retentionunit=nothing)
    molecular_apexes = [
        test_ladder_apex(10, 30; source=:molecularion, good=true),
        test_ladder_apex(8, 10; source=:molecularion, good=false),
        test_ladder_apex(99, 99; source=:molecularion, success=false)
    ]
    gapfilled = [
        test_ladder_addition(test_ladder_apex(
            9,
            20;
            source=:gapfilled,
            cosine=0.95,
            required=0.85
        ))
    ]
    rightextended = [
        test_ladder_addition(test_ladder_apex(
            11,
            40;
            source=:rightextended,
            cosine=0.93,
            required=0.9
        ))
    ]

    AlkaneSeriesResult(
        standard,
        ones(1, 1),
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        test_ladder_apex_info(molecular_apexes),
        test_ladder_addition_info(gapfilled, [], rightextended),
        nothing,
        retentionunit
    )
end

@testset "alkaneladdersteps returns merged refined ladder step view" begin
    result = test_ladder_step_result()
    compact = AlkaneSeriesResult(
        result.standard,
        result.variances,
        result.varianceinfo,
        result.baselineinfo,
        result.channelinfo,
        result.abundanceinfo,
        result.molecularioninfo,
        result.pathinfo,
        result.apexinfo,
        result.additioninfo
    )

    steps = alkaneladdersteps(result)

    @test compact.datainfo === nothing
    @test retentionunit(compact) === nothing
    @test steps isa Vector{AlkaneLadderStep}
    @test [step.ladderstep for step in steps] == [8, 9, 10, 11]
    @test [step.source for step in steps] ==
        [:molecularion, :gapfilled, :molecularion, :rightextended]
    @test steps[2].massspectrumcosine == 0.95
    @test steps[2].requiredcosine == 0.85
    @test steps[1].goodforcalibration == false
    @test steps[2].apex.ladderstep == 9
end

@testset "alkaneladdersteps filters by evidence source" begin
    result = test_ladder_step_result()

    mi = alkaneladdersteps(result; gapfilled=false, edgeextended=false)
    gap = alkaneladdersteps(result; molecularion=false, edgeextended=false)
    edge = alkaneladdersteps(result; molecularion=false, gapfilled=false)

    @test [step.ladderstep for step in mi] == [8, 10]
    @test [step.ladderstep for step in gap] == [9]
    @test [step.ladderstep for step in edge] == [11]
end

@testset "alkaneladdercalibrationpoints returns curated RT RI anchors" begin
    result = test_ladder_step_result(; retentionunit=u"minute")

    points = alkaneladdercalibrationpoints(result)

    @test points isa Vector{AlkaneLadderCalibrationPoint}
    @test [point.ladderstep for point in points] == [9, 10, 11]
    @test [point.retention for point in points] == [200.0, 300.0, 400.0]
    @test [point.retentionindex for point in points] == [900.0, 1000.0, 1100.0]
    @test all(point.retentionunit == u"minute" for point in points)
    @test [point.source for point in points] ==
        [:gapfilled, :molecularion, :rightextended]

    with_bad_step = alkaneladdercalibrationpoints(result; include=[8])
    @test [point.ladderstep for point in with_bad_step] == [8, 9, 10, 11]

    without_gap = alkaneladdercalibrationpoints(result; gapfilled=false)
    @test [point.ladderstep for point in without_gap] == [10, 11]

    excluded = alkaneladdercalibrationpoints(result; exclude=[10])
    @test [point.ladderstep for point in excluded] == [9, 11]

    replacement = AlkaneLadderCalibrationPoint(10, 305.0, u"minute", 1000.0, :manual, true)
    replaced = alkaneladdercalibrationpoints(result; exclude=[10], extra=[replacement])
    @test [point.ladderstep for point in replaced] == [9, 10, 11]
    @test replaced[2].retention == 305.0
    @test replaced[2].source ≡ :manual

    @test_throws ArgumentError alkaneladdercalibrationpoints(result; extra=[replacement])
    @test_throws ArgumentError alkaneladdercalibrationpoints(
        test_ladder_step_result(; standard=nothing)
    )
end

@testset "fitmap accepts alkane ladder calibration points and results" begin
    result = test_ladder_step_result(; retentionunit=u"minute")

    points = alkaneladdercalibrationpoints(result)
    mapper_from_points = fitmap(points)
    mapper_from_result = fitmap(result)
    mapper_with_curation = fitmap(result; goodforcalibration=false, exclude=[8])

    @test mapper_from_points isa RetentionMapper
    @test mapper_from_result isa RetentionMapper
    @test mapper_with_curation isa RetentionMapper
    @test mapper_from_points.lambda == 3e-9
    @test mapper_from_result.lambda == 3e-9
    @test fitmap(points; λ=1e-8).lambda == 1e-8
    @test retentionunit_A(mapper_from_points) == u"minute"
    @test rawretentions_A(mapper_from_points) ≈ [200.0, 300.0, 400.0]
    @test rawretentions_B(mapper_from_points) ≈ [900.0, 1000.0, 1100.0]
    @test rawretentions_A(mapper_from_result) ≈ rawretentions_A(mapper_from_points)
end

@testset "alkane m/z channel matching reports missing standard entries" begin
    msm = test_findalkanes_msm()
    standard = AlkaneStandard(
        "missing C8",
        [JuChrom.alkane_reference_spectrum(9, "nonane", 900.0, [128], [100.0])],
        NamedTuple()
    )

    @test_throws ArgumentError JuChrom.alkane_mz_channels(msm, standard, 8:8, 0.0)
end

@testset "findalkanes preprocessing and variance wrappers" begin
    msm = test_findalkanes_msm()
    standard = test_findalkanes_standard()
    mzretentionkwargs = test_findalkanes_mzretentionkwargs(msm)

    @test_throws ArgumentError findalkanes(
        msm;
        standard=standard,
        subtractbaseline=false
    )

    result = findalkanes(
        msm;
        standard=standard,
        carbonrange=8:8,
        pathminsteps=1,
        thresholdfraction=0.05,
        variancewindowsize=5,
        variancemintransitioncount=1,
        baselineλ=1e3,
        apexminioncount=1,
        apexmzscanorder=:ascending,
        apexmzretentionkwargs=mzretentionkwargs
    )

    @test result.varianceinfo isa CountVarianceEstimate
    @test result.baselineinfo isa JuChrom.AlkaneBaselineInfo
    @test result.baselineinfo.estimator ≡ :arpls
    @test result.pathinfo.status ≡ :success

    vmsm = VarianceMassScanMatrix(msm, ones(size(rawintensities(msm))))
    direct = findalkaneseries(
        vmsm;
        standard=standard,
        carbonrange=8:8,
        pathminsteps=1,
        thresholdfraction=0.05,
        apexminioncount=1,
        apexmzscanorder=:ascending,
        apexmzretentionkwargs=mzretentionkwargs
    )

    @test direct.varianceinfo === vmsm
    @test direct.pathinfo.status ≡ :success
end
