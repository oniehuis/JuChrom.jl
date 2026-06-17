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
    abundanceinfo = AlkaneAbundanceInfo(
        Dict(8 => abundance),
        Dict(8 => ones(length(retentions))),
        Dict{Int, Vector{AlkaneAbundanceWindow}}(),
        NamedTuple()
    )
    candidate = (
        ladderstep=8,
        scanindex=centerindex,
        window=(leftindex=1, rightindex=length(retentions))
    )
    mzkwargs = (
        retentionref=:middle,
        scaninterval=1e-9,
        mzcount=length(mzs),
        order=:ascending,
        dwellref=:middle,
        dwell=:homogeneous
    )

    msm, variances, abundanceinfo, candidate, mzkwargs
end

function test_apex_path_contrast(candidate)
    scanindex = candidate.scanindex
    leftindex = candidate.window.leftindex
    rightindex = candidate.window.rightindex
    AlkaneMolecularIonContrast(
        candidate.ladderstep,
        leftindex,
        scanindex,
        rightindex,
        collect(leftindex:rightindex),
        ones(rightindex - leftindex + 1),
        1.0,
        114,
        Int[],
        Int[],
        Int[],
        Float64[],
        Float64[],
        Float64[],
        false,
        false,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        2.0,
        0.0,
        0.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0
    )
end

function test_apex_path_candidate(candidate)
    contrast = test_apex_path_contrast(candidate)
    JuChrom.AlkaneLadderPathCandidate(
        candidate.ladderstep,
        candidate.scanindex,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        NaN,
        NaN,
        0,
        :molecular_ion_dp,
        true,
        false,
        Float64(candidate.scanindex),
        0.0,
        NaN,
        contrast
    )
end

function test_apex_pathinfo(candidate)
    test_apex_pathinfo_from_candidates([test_apex_path_candidate(candidate)])
end

function test_apex_pathinfo_from_candidates(pathcandidates)
    candidatesbystep = Dict(
        candidate.ladderstep => [candidate] for candidate in pathcandidates
    )
    JuChrom.AlkaneLadderPathInfo(
        true,
        :success,
        :success,
        nothing,
        "",
        [candidate.ladderstep for candidate in pathcandidates],
        [candidate.scanindex for candidate in pathcandidates],
        [candidate.score for candidate in pathcandidates],
        [candidate.z for candidate in pathcandidates],
        [candidate.centerz for candidate in pathcandidates],
        [candidate.isolationz for candidate in pathcandidates],
        [candidate.centervslowerz for candidate in pathcandidates],
        [candidate.centervsupperz for candidate in pathcandidates],
        [candidate.massspectrumcosine for candidate in pathcandidates],
        [candidate.massspectrumdistance for candidate in pathcandidates],
        [candidate.source for candidate in pathcandidates],
        [candidate.misupported for candidate in pathcandidates],
        [candidate.gapfilled for candidate in pathcandidates],
        [candidate.expectedscan for candidate in pathcandidates],
        [candidate.scanerror for candidate in pathcandidates],
        [candidate.requiredcosine for candidate in pathcandidates],
        Int[],
        Int[],
        Float64[],
        0,
        sum(candidate.score for candidate in pathcandidates),
        sum(candidate.score for candidate in pathcandidates),
        [candidate.window for candidate in pathcandidates],
        pathcandidates,
        candidatesbystep,
        candidatesbystep,
        Vector{Int}[],
        JuChrom.AlkaneLadderPathRunResult[],
        NamedTuple()
    )
end

function test_empty_apex_pathinfo()
    candidatesbystep = Dict{Int, Vector{JuChrom.AlkaneLadderPathCandidate{
        AlkaneMolecularIonContrast{Float64}
    }}}()
    JuChrom.AlkaneLadderPathInfo(
        false,
        :failed,
        :no_candidates,
        "no candidates",
        "no candidates",
        Int[],
        Int[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Symbol[],
        Bool[],
        Bool[],
        Float64[],
        Float64[],
        Float64[],
        Int[],
        Int[],
        Float64[],
        0,
        -Inf,
        -Inf,
        AlkaneMolecularIonContrast[],
        JuChrom.AlkaneLadderPathCandidate[],
        candidatesbystep,
        candidatesbystep,
        Vector{Int}[],
        JuChrom.AlkaneLadderPathRunResult[],
        NamedTuple()
    )
end

function test_mzretention_timing_kwargs(mzkwargs)
    (
        retentionref=mzkwargs.retentionref,
        scaninterval=mzkwargs.scaninterval,
        mzcount=mzkwargs.mzcount,
        dwellref=mzkwargs.dwellref,
        dwell=mzkwargs.dwell
    )
end

function test_apex_settings(;
    standard=defaultalkanestandard(),
    scanwindow=2,
    variancefloor=1.0,
    logfloorfraction=1e-3,
    mzretentionkwargs=nothing,
    mzscanorder=:inferdirection,
    apexionexcludemzvalues=(),
    apexionmzvalues=nothing,
    apexionminrelativeintensity=0.1,
    minioncount=3,
    maxapexshiftfromguess=3.0,
    apexfitqualityoutlierz=3.0,
    mzscanordermaxpeaks=3,
    mzscanorderminpeaks=5,
    mzscanorderminapexvarianceratio=2.0,
    mzscanordershapeioncount=3,
    mzscanordershapemzspacing=14,
    mzscanorderextremeioncount=5,
    mzscanorderminioncount=8
)
    JuChrom.AlkaneLadderApexSettings(
        standard,
        scanwindow,
        variancefloor,
        logfloorfraction,
        mzretentionkwargs,
        mzscanorder,
        apexionexcludemzvalues,
        apexionmzvalues,
        apexionminrelativeintensity,
        minioncount,
        maxapexshiftfromguess,
        apexfitqualityoutlierz,
        mzscanordermaxpeaks,
        mzscanorderminpeaks,
        mzscanorderminapexvarianceratio,
        mzscanordershapeioncount,
        mzscanordershapemzspacing,
        mzscanorderextremeioncount,
        mzscanorderminioncount
    )
end

function test_scan_order_inputs(; trueorder=:ascending, carbons=8:11)
    retentions = collect(1.0:(10.0 + 5.0 * length(carbons)))
    mzs = collect(43:14:211)
    ionheights = [100.0, 92.0, 84.0, 76.0, 68.0, 60.0, 52.0,
        44.0, 36.0, 28.0, 20.0, 14.0, 10.0]
    spectra = [
        MassSpectrum(
            mzs,
            ionheights;
            attrs=(order=carbon, label="C$(carbon)", ri=100.0 * carbon)
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
        dwell=:homogeneous
    )

    X = zeros(length(retentions), length(mzs))
    path = NamedTuple[]
    for (offset, carbon) in enumerate(carbons)
        apexretention = 6.0 + 5.0 * (offset - 1)
        for mzindex in eachindex(mzs), scanindex in eachindex(retentions)
            observation = JuChrom.alkane_ladder_observation_retention(
                retentions,
                scanindex,
                mzindex,
                mzkwargs
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
                    apexabundance=maximum(X[:, offset])
                )
            )
        )
    end

    msm = MassScanMatrix(retentions, mzs, X)
    variances = ones(size(X))
    timingkwargs = test_mzretention_timing_kwargs(mzkwargs)
    pathinfo = test_apex_pathinfo_from_candidates([
        test_apex_path_candidate(candidate) for candidate in path
    ])

    msm, variances, pathinfo, standard, timingkwargs
end

@testset "alkane scan-order contrast ion selection uses paired m/z extremes" begin
    mzs = vcat(collect(43.0:14.0:225.0), 226.0)
    msm = MassScanMatrix(collect(1.0:3.0), mzs, zeros(3, length(mzs)))
    spectrum = MassSpectrum(
        mzs,
        collect(range(100.0, 20.0; length=length(mzs)));
        attrs=(order=16, label="C16", ri=1600.0)
    )
    standard = AlkaneStandard("test C16", [spectrum], NamedTuple())

    selection = JuChrom.alkane_ladder_scan_order_extreme_apex_mzvalues(
        msm,
        standard,
        16,
        5
    )

    @test selection.selection == :reference_alkane_series_extreme_grid_ions
    @test selection.mzvalues == [43.0, 57.0, 71.0, 85.0, 99.0,
        169.0, 183.0, 197.0, 225.0, 226.0]
    @test 211.0 in selection.excluded_mzvalues
    @test iseven(length(selection.mzvalues))
end

@testset "alkaneladderapexes refines selected path apex from selected ions" begin
    msm, variances, abundanceinfo, candidate, mzkwargs = test_alkane_ion_apex_inputs()
    pathinfo = test_apex_pathinfo(candidate)

    settings = test_apex_settings(
        scanwindow=2,
        apexionmzvalues=[100.0, 101.0, 102.0],
        minioncount=3,
        mzretentionkwargs=test_mzretention_timing_kwargs(mzkwargs),
        mzscanorder=:ascending,
        variancefloor=1.0,
        logfloorfraction=1e-9
    )
    apexinfo = alkaneladderapexes(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        settings
    )
    apex = only(apexinfo.apexes)

    @test apexinfo.status == :success
    @test apex.success
    @test apex.ladderstep == 8
    @test apex.scanindex == 3
    @test apex.apexretention ≈ 3.25 atol = 1e-5
    @test apex.apexscanindex ≈ 3.25 atol = 1e-5
    @test apex.apexoffsetscans ≈ 0.25 atol = 1e-5
    @test apex.fit.fit.apex_ion_selection == :explicit_mzvalues
    @test apex.fit.fit.mz_indices == [1, 2, 3]
    @test apex.fit.recenter_attempts == 0
    @test apexinfo.bycarbon[8].apexretention ≈ apex.apexretention
end

@testset "alkaneladderapex recenters an off-center initial ion fit" begin
    msm, variances, abundanceinfo, candidate, mzkwargs =
        test_alkane_ion_apex_inputs(; apexretention=4.2, centerindex=3)

    settings = test_apex_settings(
        scanwindow=1,
        apexionmzvalues=[100.0, 101.0, 102.0],
        minioncount=3,
        mzretentionkwargs=mzkwargs,
        variancefloor=1.0,
        logfloorfraction=1e-9,
        maxapexshiftfromguess=3.0
    )
    apex = alkaneladderapex(
        msm,
        variances,
        abundanceinfo,
        candidate,
        settings
    )

    @test apex.success
    @test apex.apexretention ≈ 4.2 atol = 1e-5
    @test apex.fit.recenter_attempts > 0
    @test apex.fit.recenter_used
    @test any(a -> a.success && abs(a.apex_offset_from_fit_center_scans) <= 0.25,
        apex.fit.all_apex_attempts)
    @test apex.candidate == candidate
end

@testset "alkaneladderapexes returns empty result for failed path" begin
    msm, variances, abundanceinfo, _, _ = test_alkane_ion_apex_inputs()
    pathinfo = test_empty_apex_pathinfo()

    apexinfo = alkaneladderapexes(
        msm,
        variances,
        abundanceinfo,
        pathinfo,
        test_apex_settings()
    )

    @test apexinfo.status == :failed
    @test apexinfo.reason == :no_path
    @test isempty(apexinfo.apexes)
    @test apexinfo.scanorderinfo.status == :no_path
    @test apexinfo.scanorderinfo.selected_order == :unknown
    @test isnothing(apexinfo.scanorderinfo.evidence_score)
end

@testset "alkaneladderapex reports structured failure when too few ions are available" begin
    msm, variances, abundanceinfo, candidate, mzkwargs = test_alkane_ion_apex_inputs()

    settings = test_apex_settings(
        scanwindow=2,
        apexionmzvalues=[100.0],
        minioncount=3,
        mzretentionkwargs=mzkwargs
    )
    apex = alkaneladderapex(
        msm,
        variances,
        abundanceinfo,
        candidate,
        settings
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
        1e-3
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_apex_settings(
        2,
        0.0,
        1e-3
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_apex_settings(
        2,
        1.0,
        0.0
    )
end

@testset "alkaneladderscanorder keeps explicit user order" begin
    msm, variances, pathinfo, standard, timingkwargs = test_scan_order_inputs()
    provided = merge(timingkwargs, (order=:ascending,))

    scanorder = alkaneladderscanorder(
        msm,
        variances,
        pathinfo,
        test_apex_settings(; standard=standard, mzretentionkwargs=provided)
    )

    @test scanorder.mzretentionkwargs.order == :ascending
    @test scanorder.info.status == :provided
    @test scanorder.info.selected_order == :ascending
    @test isnothing(scanorder.info.evidence_score)

    simultaneous = merge(timingkwargs, (order=:simultaneous,))
    scanorder = alkaneladderscanorder(
        msm,
        variances,
        pathinfo,
        test_apex_settings(; standard=standard, mzretentionkwargs=simultaneous)
    )

    @test scanorder.mzretentionkwargs.order == :simultaneous
    @test scanorder.info.status == :provided
    @test scanorder.info.selected_order == :simultaneous

    @test_throws ArgumentError alkaneladderscanorder(
        msm,
        variances,
        pathinfo,
        test_apex_settings(;
            standard=standard,
            mzretentionkwargs=provided,
            mzscanorder=:descending
        )
    )
    @test_throws ArgumentError alkaneladderscanorder(
        msm,
        variances,
        pathinfo,
        test_apex_settings(; standard=standard, mzscanorder=:unknown)
    )
end

@testset "alkaneladderscanorder orders limited peaks across the path" begin
    pathinfo = test_apex_pathinfo_from_candidates([
        test_apex_path_candidate((
            ladderstep=step,
            scanindex=step,
            window=(leftindex=step, rightindex=step)
        ))
        for step in 8:40
    ])

    selected = JuChrom.alkane_ladder_scan_order_candidates(pathinfo, 5)
    @test [candidate.ladderstep for candidate in selected[1:5]] == [8, 16, 24, 32, 40]
    @test sort([candidate.ladderstep for candidate in selected]) == collect(8:40)

    selected = JuChrom.alkane_ladder_scan_order_candidates(pathinfo, 6)
    @test [candidate.ladderstep for candidate in selected[1:6]] == [8, 14, 21, 27, 34, 40]

    selected = JuChrom.alkane_ladder_scan_order_candidates(pathinfo, nothing)
    @test [candidate.ladderstep for candidate in selected] == collect(8:40)

    shortpath = test_apex_pathinfo_from_candidates(pathinfo.path[1:5])
    selected = JuChrom.alkane_ladder_scan_order_candidates(shortpath, 4)
    @test [candidate.ladderstep for candidate in selected[1:4]] == [8, 9, 11, 12]

    @test_throws ArgumentError JuChrom.validate_alkane_ladder_scan_order_settings(
        0,
        4,
        1.25,
        3,
        14,
        5,
        8
    )
end

@testset "alkaneladderscanorder expands ambiguous subset to all peaks" begin
    msm, variances, pathinfo, standard, timingkwargs =
        test_scan_order_inputs(; trueorder=:ascending, carbons=8:15)

    scanorder = alkaneladderscanorder(
        msm,
        variances,
        pathinfo,
        test_apex_settings(;
            standard=standard,
            mzretentionkwargs=timingkwargs,
            scanwindow=2,
            logfloorfraction=1e-9,
            mzscanorderminpeaks=3,
            mzscanorderminapexvarianceratio=1e12,
            mzscanorderminioncount=8
        )
    )

    @test scanorder.info.status == :ambiguous
    @test scanorder.info.expanded_for_ambiguity
    @test scanorder.info.expanded_ions
    @test scanorder.info.initial_n_peaks_tried == 3
    @test scanorder.info.trials.ascending.n_tried == length(pathinfo.path)
    @test scanorder.info.trials.descending.n_tried == length(pathinfo.path)
    @test maximum(scanorder.info.trials.ascending.ion_counts) > 10
    @test maximum(scanorder.info.trials.descending.ion_counts) > 10
end

@testset "alkaneladderscanorder infers sequential scan direction symmetrically" begin
    for trueorder in (:ascending, :descending)
        msm, variances, pathinfo, standard, timingkwargs =
            test_scan_order_inputs(; trueorder=trueorder)

        scanorder = alkaneladderscanorder(
            msm,
            variances,
            pathinfo,
            test_apex_settings(;
                standard=standard,
                mzretentionkwargs=timingkwargs,
                scanwindow=2,
                logfloorfraction=1e-9,
                mzscanordermaxpeaks=1,
                mzscanorderminpeaks=3,
                mzscanorderminioncount=8
            )
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
