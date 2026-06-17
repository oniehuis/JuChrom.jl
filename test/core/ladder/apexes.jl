using Test
using JuChrom

function test_alkane_ion_apex_inputs(; apexretention=3.25, centerindex=3)
    retentions = collect(1.0:7.0)
    mzs = [100.0, 101.0, 102.0]
    ionheights = [100.0, 60.0, 30.0]
    profile = exp.(-0.5 .* abs2.((retentions .- apexretention) ./ 1.0))
    X = profile * ionheights'
    msm = MassScanMatrix(retentions, mzs, X)
    variances = ones(size(X))
    abundance = vec(sum(X; dims=2))
    abundanceinfo = (
        abundances=Dict(8 => abundance),
        abundancevariances=Dict(8 => ones(length(retentions))),
    )
    candidate = (
        ladderstep=8,
        scanindex=centerindex,
        window=(leftindex=1, rightindex=length(retentions)),
    )
    mzkwargs = (
        retentionref=:middle,
        scaninterval=1e-9,
        mzcount=length(mzs),
        order=:ascending,
        dwellref=:middle,
        dwell=:homogeneous,
    )

    msm, variances, abundanceinfo, candidate, mzkwargs
end

function test_mzretention_timing_kwargs(mzkwargs)
    (
        retentionref=mzkwargs.retentionref,
        scaninterval=mzkwargs.scaninterval,
        mzcount=mzkwargs.mzcount,
        dwellref=mzkwargs.dwellref,
        dwell=mzkwargs.dwell,
    )
end

function test_scan_order_inputs(; trueorder=:ascending, carbons=8:11)
    retentions = collect(1.0:(10.0 + 5.0 * length(carbons)))
    mzs = collect(43.0:14.0:211.0)
    ionheights = [100.0, 92.0, 84.0, 76.0, 68.0, 60.0, 52.0,
        44.0, 36.0, 28.0, 20.0, 14.0, 10.0]
    spectra = [
        MassSpectrum(
            mzs,
            ionheights;
            attrs=(order=carbon, label="C$(carbon)", ri=100.0 * carbon),
        )
        for carbon in carbons
    ]
    standard = AlkaneStandard("synthetic scan-order standard", spectra, NamedTuple())
    mzkwargs = (
        retentionref=:start,
        scaninterval=1.0,
        mzcount=length(mzs),
        order=trueorder,
        dwellref=:middle,
        dwell=:homogeneous,
    )

    X = zeros(length(retentions), length(mzs))
    path = NamedTuple[]
    for (offset, carbon) in enumerate(carbons)
        apexretention = 6.0 + 5.0 * (offset - 1)
        for mzindex in eachindex(mzs), scanindex in eachindex(retentions)
            observation = JuChrom.alkane_ladder_observation_retention(
                MassScanMatrix(retentions, mzs, zeros(length(retentions), length(mzs))),
                retentions,
                scanindex,
                mzindex,
                mzkwargs,
            )
            X[scanindex, mzindex] += ionheights[mzindex] *
                exp(-0.5 * abs2((observation - apexretention) / 0.65))
        end
        inputscan = round(Int, apexretention)
        push!(
            path,
            (
                ladderstep=carbon,
                scanindex=inputscan,
                apexabundance=maximum(X[:, offset]),
                window=(
                    leftindex=max(1, inputscan - 4),
                    rightindex=min(length(retentions), inputscan + 4),
                    apexabundance=maximum(X[:, offset]),
                ),
            ),
        )
    end

    msm = MassScanMatrix(retentions, mzs, X)
    variances = ones(size(X))
    timingkwargs = test_mzretention_timing_kwargs(mzkwargs)

    msm, variances, (path=path,), standard, timingkwargs
end

@testset "alkane scan-order contrast ion selection uses paired m/z extremes" begin
    mzs = vcat(collect(43.0:14.0:225.0), 226.0)
    msm = MassScanMatrix(collect(1.0:3.0), mzs, zeros(3, length(mzs)))
    spectrum = MassSpectrum(
        mzs,
        collect(range(100.0, 20.0; length=length(mzs)));
        attrs=(order=16, label="C16", ri=1600.0),
    )
    standard = AlkaneStandard("test C16", [spectrum], NamedTuple())

    selection = JuChrom.alkane_ladder_scan_order_extreme_apex_mzvalues(
        msm,
        standard,
        16;
        extremeioncount=5,
    )

    @test selection.selection == :reference_alkane_series_extreme_grid_ions
    @test selection.mzvalues == [43.0, 57.0, 71.0, 85.0, 99.0,
        169.0, 183.0, 197.0, 225.0, 226.0]
    @test 211.0 in selection.excluded_mzvalues
    @test iseven(length(selection.mzvalues))
end

@testset "alkaneladderapexes refines selected path apex from selected ions" begin
    msm, variances, abundanceinfo, candidate, mzkwargs = test_alkane_ion_apex_inputs()
    pathinfo = (path=[candidate],)

    apexinfo = alkaneladderapexes(
        msm,
        variances,
        abundanceinfo,
        pathinfo;
        scanwindow=2,
        apexionmzvalues=[100.0, 101.0, 102.0],
        minioncount=3,
        mzretentionkwargs=test_mzretention_timing_kwargs(mzkwargs),
        mzscanorder=:ascending,
        variancefloor=1.0,
        logfloorfraction=1e-9,
    )
    apex = only(apexinfo.apexes)

    @test apexinfo.status == :success
    @test apex.success
    @test apex.ladderstep == 8
    @test apex.scanindex == 3
    @test apex.apexretention ≈ 3.25 atol = 1e-5
    @test apex.apexscanindex ≈ 3.25 atol = 1e-5
    @test apex.apexoffsetscans ≈ 0.25 atol = 1e-5
    @test apex.apex_ion_selection == :explicit_mzvalues
    @test apex.mz_indices == [1, 2, 3]
    @test apex.recenter_attempts == 0
    @test apexinfo.bycarbon[8].apexretention ≈ apex.apexretention
end

@testset "alkaneladderapex recenters an off-center initial ion fit" begin
    msm, variances, abundanceinfo, candidate, mzkwargs =
        test_alkane_ion_apex_inputs(; apexretention=4.2, centerindex=3)

    apex = alkaneladderapex(
        msm,
        variances,
        abundanceinfo,
        candidate;
        scanwindow=1,
        apexionmzvalues=[100.0, 101.0, 102.0],
        minioncount=3,
        mzretentionkwargs=mzkwargs,
        variancefloor=1.0,
        logfloorfraction=1e-9,
        maxapexshiftfromguess=3.0,
    )

    @test apex.success
    @test apex.apexretention ≈ 4.2 atol = 1e-5
    @test apex.recenter_attempts > 0
    @test apex.recenter_used
    @test any(a -> a.success && abs(a.apex_offset_from_fit_center_scans) <= 0.25,
        apex.all_apex_attempts)
    @test apex.candidate == candidate
end

@testset "alkaneladderapexes returns empty result for failed path" begin
    msm, variances, abundanceinfo, _, _ = test_alkane_ion_apex_inputs()
    pathinfo = (path=NamedTuple[],)

    apexinfo = alkaneladderapexes(msm, variances, abundanceinfo, pathinfo)

    @test apexinfo.status == :failed
    @test apexinfo.reason == :no_path
    @test isempty(apexinfo.apexes)
    @test apexinfo.scanorderinfo.status == :no_path
    @test apexinfo.scanorderinfo.selected_order == :unknown
    @test isnothing(apexinfo.scanorderinfo.evidence_score)
end

@testset "alkaneladderapex reports structured failure when too few ions are available" begin
    msm, variances, abundanceinfo, candidate, mzkwargs = test_alkane_ion_apex_inputs()

    apex = alkaneladderapex(
        msm,
        variances,
        abundanceinfo,
        candidate;
        scanwindow=2,
        apexionmzvalues=[100.0],
        minioncount=3,
        mzretentionkwargs=mzkwargs,
    )

    @test !apex.success
    @test apex.reason == :apex_fit_failed
    @test occursin("fewer than 3 apex ion", apex.failurereason)
    @test apex.apexscanindex == 3.0
    @test apex.apexretention == 3.0
end

@testset "alkane ladder apex validation" begin
    @test JuChrom.validate_alkane_ladder_apex_settings(2, 1.0, 1e-3) === nothing
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_apex_settings(
        0,
        1.0,
        1e-3,
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_apex_settings(
        2,
        0.0,
        1e-3,
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_apex_settings(
        2,
        1.0,
        0.0,
    )
end

@testset "alkaneladderscanorder keeps explicit user order" begin
    msm, variances, pathinfo, standard, timingkwargs = test_scan_order_inputs()
    provided = merge(timingkwargs, (order=:ascending,))

    scanorder = alkaneladderscanorder(
        msm,
        variances,
        pathinfo;
        standard=standard,
        mzretentionkwargs=provided,
    )

    @test scanorder.mzretentionkwargs.order == :ascending
    @test scanorder.info.status == :provided
    @test scanorder.info.selected_order == :ascending
    @test isnothing(scanorder.info.evidence_score)

    simultaneous = merge(timingkwargs, (order=:simultaneous,))
    scanorder = alkaneladderscanorder(
        msm,
        variances,
        pathinfo;
        standard=standard,
        mzretentionkwargs=simultaneous,
    )

    @test scanorder.mzretentionkwargs.order == :simultaneous
    @test scanorder.info.status == :provided
    @test scanorder.info.selected_order == :simultaneous

    @test_throws ArgumentError alkaneladderscanorder(
        msm,
        variances,
        pathinfo;
        standard=standard,
        mzretentionkwargs=provided,
        mzscanorder=:descending,
    )
    @test_throws ArgumentError alkaneladderscanorder(
        msm,
        variances,
        pathinfo;
        standard=standard,
        mzscanorder=:unknown,
    )
end

@testset "alkaneladderscanorder orders limited peaks across the path" begin
    pathinfo = (
        path=[
            (
                ladderstep=step,
                scanindex=step,
                window=(leftindex=step, rightindex=step),
            )
            for step in 8:40
        ],
    )

    selected = JuChrom.alkane_ladder_scan_order_candidates(pathinfo; maxpeaks=5)
    @test [candidate.ladderstep for candidate in selected[1:5]] == [8, 16, 24, 32, 40]
    @test sort([candidate.ladderstep for candidate in selected]) == collect(8:40)

    selected = JuChrom.alkane_ladder_scan_order_candidates(pathinfo; maxpeaks=6)
    @test [candidate.ladderstep for candidate in selected[1:6]] == [8, 14, 21, 27, 34, 40]

    selected = JuChrom.alkane_ladder_scan_order_candidates(pathinfo; maxpeaks=nothing)
    @test [candidate.ladderstep for candidate in selected] == collect(8:40)

    shortpath = (path=pathinfo.path[1:5],)
    selected = JuChrom.alkane_ladder_scan_order_candidates(shortpath; maxpeaks=4)
    @test [candidate.ladderstep for candidate in selected[1:4]] == [8, 9, 11, 12]

    @test_throws ArgumentError JuChrom.validate_alkane_ladder_scan_order_settings(
        0,
        4,
        1.25,
        3,
        14,
        5,
        8,
    )
end

@testset "alkaneladderscanorder expands ambiguous subset to all peaks" begin
    msm, variances, pathinfo, standard, timingkwargs =
        test_scan_order_inputs(; trueorder=:ascending, carbons=8:15)

    scanorder = alkaneladderscanorder(
        msm,
        variances,
        pathinfo;
        standard=standard,
        mzretentionkwargs=timingkwargs,
        scanwindow=2,
        logfloorfraction=1e-9,
        minpeaks=3,
        minapexvarianceratio=1e12,
        minioncount=8,
    )

    @test scanorder.info.status == :ambiguous
    @test scanorder.info.expanded_for_ambiguity
    @test scanorder.info.initial_n_peaks_tried == 3
    @test scanorder.info.trials.ascending.n_tried == length(pathinfo.path)
    @test scanorder.info.trials.descending.n_tried == length(pathinfo.path)
end

@testset "alkaneladderscanorder infers sequential scan direction symmetrically" begin
    for trueorder in (:ascending, :descending)
        msm, variances, pathinfo, standard, timingkwargs =
            test_scan_order_inputs(; trueorder=trueorder)

        scanorder = alkaneladderscanorder(
            msm,
            variances,
            pathinfo;
            standard=standard,
            mzretentionkwargs=timingkwargs,
            scanwindow=2,
            logfloorfraction=1e-9,
            maxpeaks=1,
            minpeaks=3,
            minioncount=8,
        )

        @test scanorder.mzretentionkwargs.order == trueorder
        @test scanorder.info.selected_order == trueorder
        @test scanorder.info.requested_order == :inferdirection
        @test scanorder.info.status == :accepted
        @test scanorder.info.n_peaks_used == 3
        @test scanorder.info.trials.ascending.n_tried >= 3
        @test scanorder.info.trials.descending.n_tried >= 3
    end
end
