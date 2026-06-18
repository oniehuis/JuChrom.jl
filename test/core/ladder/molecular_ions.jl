using Test
using JuChrom

function molecular_ion_test_inputs()
    peakmodel = [0.0, 0.5, 1.0, 0.5, 0.0]
    X = zeros(5, 5)
    X[:, 1] .= 3.0 .* peakmodel
    X[:, 2] .= 10.0 .* peakmodel
    X[:, 3] .= 2.0 .* peakmodel
    X[:, 4] .= -5.0 .* peakmodel
    msm = MassScanMatrix(collect(1.0:5.0), [100.0, 114.0, 115.0, 128.0, 129.0], X)
    variances = ones(size(X))
    abundance = [0.0, 1.0, 2.0, 1.0, 0.0]
    window = AlkaneAbundanceWindow(
        8,
        1,
        3,
        5,
        0.0,
        2.0,
        0.0,
        0.1,
        :threshold,
        :threshold
    )
    abundanceinfo = AlkaneAbundanceInfo(
        Dict(8 => abundance),
        Dict(8 => fill(0.1, length(abundance))),
        Dict(8 => [window]),
        NamedTuple()
    )

    msm, variances, abundanceinfo, window, peakmodel
end

function test_molecular_ion_contrast(leftindex, rightindex, z, molecularionscore)
    AlkaneMolecularIonContrast(
        8,
        leftindex,
        leftindex,
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
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        molecularionscore,
        z
    )
end

@testset "alkane_molecular_ion_info combines molecular-ion pipeline" begin
    msm, variances, abundanceinfo, _, _ = molecular_ion_test_inputs()

    molecularioninfo = JuChrom.alkane_molecular_ion_info(
        msm,
        variances,
        abundanceinfo,
        1,
        14,
        1.0,
        1.645,
        1.645
    )

    @test molecularioninfo isa AlkaneMolecularIonInfo
    @test sort(collect(keys(molecularioninfo.contrasts))) == [8]
    @test length(molecularioninfo.contrasts[8]) == 1
    @test only(molecularioninfo.contrasts[8]) isa AlkaneMolecularIonContrast
    @test molecularioninfo.zscorevectors[8] ≈
        fill(only(molecularioninfo.contrasts[8]).z, 5)
    @test molecularioninfo.settings == (
        ionwindow=1,
        stepmass=14,
        variancefloor=1.0,
        centerzmin=1.645,
        isolationzmin=1.645
    )
end

@testset "alkane_molecular_ion_window_contrasts fits center and control ions" begin
    msm, variances, abundanceinfo, _, peakmodel = molecular_ion_test_inputs()

    contrasts = JuChrom.alkane_molecular_ion_window_contrasts(
        msm,
        variances,
        abundanceinfo.abundances,
        abundanceinfo.windows,
        1,
        14,
        1.0,
        1.645,
        1.645
    )
    contrast = only(contrasts[8])

    @test contrast.ladderstep == 8
    @test contrast.molecularion == 114
    @test contrast.centerindices == [2, 3]
    @test contrast.lowerindices == [1]
    @test contrast.upperindices == [4, 5]
    @test contrast.centerions == [114.0, 115.0]
    @test contrast.lowerions == [100.0]
    @test contrast.upperions == [128.0, 129.0]
    @test contrast.scanindices == collect(1:5)
    @test contrast.peakmodel ≈ peakmodel
    @test contrast.apexindex == 3
    @test contrast.centerabundance ≈ 12.0
    @test contrast.lowerabundance ≈ 3.0
    @test contrast.upperabundance ≈ -5.0
    @test contrast.controlpenalty ≈ 3.0
    @test contrast.signedcontrastabundance ≈ 14.0
    @test contrast.contrastabundance ≈ 9.0
    @test contrast.centerstderr ≈ sqrt(4 / 3)
    @test contrast.contrastvariance ≈ 10 / 3
    @test contrast.z ≈ 9 / sqrt(10 / 3)
    @test contrast.centerz ≈ 12 / sqrt(4 / 3)
    @test contrast.centervslowerz ≈ 9 / sqrt(2)
    @test contrast.centervsupperz ≈ 17 / sqrt(8 / 3)
    @test contrast.isolationz ≈ 9 / sqrt(2)
    @test contrast.molecularionscore ≈ contrast.z

    gatedcontrasts = JuChrom.alkane_molecular_ion_window_contrasts(
        msm,
        variances,
        abundanceinfo.abundances,
        abundanceinfo.windows,
        1,
        14,
        1.0,
        100.0,
        1.645
    )
    gatedcontrast = only(gatedcontrasts[8])
    @test gatedcontrast.z ≈ contrast.z
    @test gatedcontrast.molecularionscore == 0.0

    noabundance = copy(abundanceinfo.abundances)
    delete!(noabundance, 8)
    @test_throws ArgumentError JuChrom.alkane_molecular_ion_window_contrasts(
        msm,
        variances,
        noabundance,
        abundanceinfo.windows,
        1,
        14,
        1.0,
        1.645,
        1.645
    )

    badabundances = Dict(8 => [1.0, 2.0])
    @test_throws DimensionMismatch JuChrom.alkane_molecular_ion_window_contrasts(
        msm,
        variances,
        badabundances,
        abundanceinfo.windows,
        1,
        14,
        1.0,
        1.645,
        1.645
    )
end

@testset "alkane_molecular_ion_zscore_vectors maps window z-scores to scans" begin
    abundances = Dict(8 => zeros(6))
    contrasts = Dict(8 => [
        test_molecular_ion_contrast(2, 4, 3.0, 2.0),
        test_molecular_ion_contrast(3, 5, 5.0, 4.0),
        test_molecular_ion_contrast(6, 6, NaN, NaN)
    ])

    zscores = JuChrom.alkane_molecular_ion_zscore_vectors(
        abundances,
        contrasts,
        1.0
    )

    @test zscores[8] == [1.0, 2.0, 4.0, 4.0, 4.0, 1.0]
    @test_throws ArgumentError JuChrom.alkane_molecular_ion_zscore_vectors(
        abundances,
        contrasts,
        Inf
    )
    @test_throws DimensionMismatch JuChrom.alkane_molecular_ion_zscore_vectors(
        abundances,
        Dict(8 => [test_molecular_ion_contrast(1, 7, 2.0, 2.0)]),
        0.0
    )
end

@testset "alkane molecular-ion helpers" begin
    @test JuChrom.validate_alkane_molecular_ion_settings(
        1,
        14,
        1.0,
        1.645,
        1.645
    ) ≡ nothing
    @test_throws ArgumentError JuChrom.validate_alkane_molecular_ion_settings(
        -1,
        14,
        1.0,
        1.645,
        1.645
    )
    @test_throws ArgumentError JuChrom.validate_alkane_molecular_ion_settings(
        1,
        0,
        1.0,
        1.645,
        1.645
    )
    @test_throws ArgumentError JuChrom.validate_alkane_molecular_ion_settings(
        1,
        14,
        0.0,
        1.645,
        1.645
    )
    @test_throws ArgumentError JuChrom.validate_alkane_molecular_ion_settings(
        1,
        14,
        1.0,
        Inf,
        1.645
    )
    @test_throws ArgumentError JuChrom.validate_alkane_molecular_ion_settings(
        1,
        14,
        1.0,
        1.645,
        -1.0
    )

    mzs = [100, 114, 115, 128, 129]
    @test JuChrom.alkane_molecular_ion(8) == 114
    @test JuChrom.alkane_mz_window_indices(mzs, 128, 1) == [4, 5]
    @test JuChrom.alkane_molecular_ion_group_indices(mzs, 114, 1) == [2, 3]
    @test JuChrom.alkane_molecular_ion_group_is_annotatable(mzs, 114, [2, 3])
    @test !JuChrom.alkane_molecular_ion_group_is_annotatable(mzs, 114, [2])
    @test !JuChrom.alkane_molecular_ion_group_is_annotatable(mzs, 116, Int[])

    positivewindow = AlkaneAbundanceWindow(
        8,
        1,
        2,
        3,
        0.0,
        1.0,
        0.0,
        0.0,
        :test,
        :test
    )
    badwindow = AlkaneAbundanceWindow(
        8,
        1,
        2,
        3,
        0.0,
        0.0,
        0.0,
        0.0,
        :test,
        :test
    )

    @test JuChrom.alkane_abundance_window_has_positive_peak_model(
        [0.0, 1.0, 0.0],
        positivewindow
    )
    @test !JuChrom.alkane_abundance_window_has_positive_peak_model(
        [0.0, 0.0, 0.0],
        badwindow
    )
    @test_throws DimensionMismatch JuChrom.alkane_abundance_window_has_positive_peak_model(
        [0.0, 1.0],
        positivewindow
    )
end

@testset "alkane molecular-ion fit kernels" begin
    _, variances, _, _, peakmodel = molecular_ion_test_inputs()
    X = hcat(10.0 .* peakmodel)

    fit = JuChrom.alkane_fitted_ion_abundance(
        X,
        variances[:, 1:1],
        1:5,
        1,
        peakmodel,
        1.0
    )

    @test fit isa JuChrom.AlkaneFittedIonAbundance
    @test fit.abundance ≈ 10.0
    @test fit.variance ≈ 2 / 3
    zerofit = JuChrom.alkane_fitted_ion_abundance(
        X,
        variances[:, 1:1],
        1:5,
        1,
        zeros(5),
        1.0
    )
    @test zerofit isa JuChrom.AlkaneFittedIonAbundance
    @test zerofit.abundance == 0.0
    @test zerofit.variance == 0.0
    @test_throws DimensionMismatch JuChrom.alkane_fitted_ion_abundance(
        X,
        variances[:, 1:1],
        1:5,
        1,
        ones(4),
        1.0
    )

    groupfit = JuChrom.alkane_fitted_ion_group_abundance(
        hcat(10.0 .* peakmodel, 2.0 .* peakmodel),
        ones(5, 2),
        1:5,
        [1, 2],
        peakmodel,
        1.0
    )
    @test groupfit.abundance ≈ 12.0
    @test groupfit.variance ≈ 4 / 3

    @test JuChrom.alkane_molecular_ion_control_veto_z(
        JuChrom.AlkaneFittedIonAbundance(3.0, 1.0),
        JuChrom.AlkaneFittedIonAbundance(5.0, 3.0)
    ) ≈ -1.0
    @test JuChrom.alkane_molecular_ion_center_vs_control_z(
        JuChrom.AlkaneFittedIonAbundance(5.0, 3.0),
        JuChrom.AlkaneFittedIonAbundance(3.0, 1.0)
    ) ≈ 1.0
    @test JuChrom.alkane_molecular_ion_isolation_z(2.0, 3.0) == 2.0
    @test isnan(JuChrom.alkane_molecular_ion_isolation_z(NaN, NaN))
end

@testset "alkane ladder fails fast without molecular-ion channels" begin
    msm = MassScanMatrix(collect(1.0:3.0), [43.0, 57.0, 71.0], zeros(3, 3))

    err = try
        findalkanes(
            msm;
            carbonrange=8:12,
            pathminsteps=3
        )
        nothing
    catch err
        err
    end
    @test err isa ArgumentError
    @test occursin("molecular-ion channels", sprint(showerror, err))

    @test_throws ArgumentError findalkaneseries(
        msm,
        ones(size(rawintensities(msm)));
        carbonrange=8:12,
        pathminsteps=3
    )
end

@testset "findalkaneseries stores molecular-ion info" begin
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

    result = findalkaneseries(
        msm,
        ones(size(X));
        standard=standard,
        carbonrange=8:8,
        thresholdfraction=0.05,
        pathminsteps=1,
        apexmzscanorder=:ascending,
        apexmzretentionkwargs=(
            retentionref=:middle,
            scaninterval=1e-9,
            mzcount=mzcount(msm),
            dwellref=:middle,
            dwell=:homogeneous
        )
    )

    @test hasproperty(result, :molecularioninfo)
    @test sort(collect(keys(result.molecularioninfo.contrasts))) == [8]
    @test length(result.molecularioninfo.contrasts[8]) == 1
    @test result.molecularioninfo.settings.ionwindow == 1
    @test result.molecularioninfo.settings.stepmass == 14
    @test result.molecularioninfo.settings.centerzmin == 1.645
    @test result.molecularioninfo.settings.isolationzmin == 1.645
    @test hasproperty(result, :pathinfo)
    @test result.pathinfo.status ≡ :success
    @test result.pathinfo.laddersteps == [8]
    @test hasproperty(result, :apexinfo)
    @test result.apexinfo.status ≡ :success
    @test only(result.apexinfo.apexes).ladderstep == 8
    @test hasproperty(result, :additioninfo)
    @test result.additioninfo.status ≡ :empty
end
