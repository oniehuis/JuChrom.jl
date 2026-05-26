using Test
using JuChrom

const SPECTRA_MZRETENTION_KWARGS = (
    retentionref=:start,
    scaninterval=0.2,
    mzcount=6,
    order=:ascending,
    dwellref=:middle,
    dwell=:homogeneous,
)

function spectra_reference_standard()
    spectrum = JuChrom.alkane_reference_spectrum(
        8,
        "octane synthetic",
        800.0,
        [43, 57, 71, 85],
        [0.18, 1.0, 0.42, 0.09],
    )
    (spectra=[spectrum],)
end

function spectra_test_msm(; apexretention=5.35)
    retentions = collect(1.0:10.0)
    mzs = [43.0, 57.0, 71.0, 85.0, 99.0, 113.0]
    heights = [90.0, 340.0, 170.0, 55.0, 22.0, -15.0]
    intensities = Matrix{Float64}(undef, length(retentions), length(mzs))
    for scanindex in eachindex(retentions), mzindex in eachindex(mzs)
        obsretention = mzretention(
            retentions[scanindex];
            SPECTRA_MZRETENTION_KWARGS...,
            mzindex=mzindex,
        )
        intensities[scanindex, mzindex] =
            heights[mzindex] * exp(-0.85 * abs2(obsretention - apexretention))
    end

    MassScanMatrix(retentions, mzs, intensities)
end

function spectra_channelinfo(msm)
    JuChrom.alkane_mz_channels(
        msm;
        standard=spectra_reference_standard(),
        carbonrange=8:8,
        minrelativeintensity=0.05,
    )
end

function spectra_test_apex(msm, variances, channelinfo)
    JuChrom.ladder_peak_apex_twopass(
        msm,
        (carbon=8, scanindex=5, retention=5.0, score=1.0),
        variances;
        channelinfo=channelinfo,
        mzretentionkwargs=SPECTRA_MZRETENTION_KWARGS,
        scanwindow=2,
        logfloorfraction=1e-9,
        variancefloor=1e-6,
    )
end

@testset "ladder step mass spectrum reuses apex ions and fits remaining ions" begin
    msm = spectra_test_msm()
    variances = fill(1.0, size(rawintensities(msm)))
    channelinfo = spectra_channelinfo(msm)
    apex = spectra_test_apex(msm, variances, channelinfo)

    spectrum = JuChrom.ladder_step_mass_spectrum(
        msm,
        apex,
        variances;
        scanwindow=2,
        mzretentionkwargs=SPECTRA_MZRETENTION_KWARGS,
        variancefloor=1e-6,
        allownegative=true,
    )

    @test spectrum isa MassSpectrum
    @test mzvalues(spectrum) == mzvalues(msm)
    @test intensities(spectrum) ≈ [90.0, 340.0, 170.0, 55.0, 22.0, -15.0] atol=0.15

    spectrum_attrs = attrs(spectrum)
    @test spectrum_attrs.carbon == 8
    @test spectrum_attrs.model == :fixed_apex_peak_shape
    @test spectrum_attrs.reference_mz_indices == [1, 2, 3, 4]
    @test spectrum_attrs.height_sources[1:4] == fill(:apex_model, 4)
    @test spectrum_attrs.height_sources[5:6] == fill(:fixed_shape_fit, 2)
    @test all(spectrum_attrs.fit_success)
    @test spectrum_attrs.n_observations[1:4] == fill(size(apex.fit_intensities, 1), 4)
    @test all(==(5), spectrum_attrs.n_observations[5:6])

    nonnegative = JuChrom.ladder_step_mass_spectrum(
        msm,
        apex,
        variances;
        scanwindow=2,
        mzretentionkwargs=SPECTRA_MZRETENTION_KWARGS,
        variancefloor=1e-6,
        allownegative=false,
    )
    @test intensities(nonnegative)[1:5] ≈ [90.0, 340.0, 170.0, 55.0, 22.0] atol=0.15
    @test intensities(nonnegative)[6] == 0.0
end

@testset "alkane series step spectra wraps per-step spectra" begin
    msm = spectra_test_msm()
    variances = fill(1.0, size(rawintensities(msm)))
    channelinfo = spectra_channelinfo(msm)
    apex = spectra_test_apex(msm, variances, channelinfo)
    apex9 = merge(apex, (carbon=9,))
    apex10 = merge(apex, (carbon=10,))
    apexes = (
        success=true,
        failurereason=nothing,
        results=[
            (
                carbon=8,
                apex=apex,
                success=true,
            ),
            (
                carbon=9,
                apex=apex9,
                success=true,
            ),
            (
                carbon=10,
                apex=apex10,
                success=true,
            ),
        ],
    )

    stepspectra = JuChrom.alkane_series_step_spectra(
        msm,
        variances,
        apexes;
        scanwindow=2,
        mzretentionkwargs=SPECTRA_MZRETENTION_KWARGS,
        variancefloor=1e-6,
    )

    @test stepspectra.success
    @test stepspectra.carbonnumbers == [8, 9, 10]
    @test sort(collect(keys(stepspectra.spectra))) == [8, 9, 10]
    @test stepspectra.spectra[8] isa MassSpectrum
    @test intensities(stepspectra.spectra[8])[5:6] ≈ [22.0, -15.0] atol=0.15
    @test stepspectra.settings.carbons === :all

    single = JuChrom.alkane_series_step_spectra(
        msm,
        variances,
        apexes;
        carbons=9,
        scanwindow=2,
        mzretentionkwargs=SPECTRA_MZRETENTION_KWARGS,
        variancefloor=1e-6,
    )
    @test single.success
    @test single.carbonnumbers == [9]
    @test only(keys(single.spectra)) == 9
    @test single.settings.carbons == [9]

    selected = JuChrom.alkane_series_step_spectra(
        msm,
        variances,
        apexes;
        carbons=[8, 10],
        scanwindow=2,
        mzretentionkwargs=SPECTRA_MZRETENTION_KWARGS,
        variancefloor=1e-6,
    )
    @test selected.success
    @test selected.carbonnumbers == [8, 10]
    @test sort(collect(keys(selected.spectra))) == [8, 10]
    @test selected.settings.carbons == [8, 10]

    ranged = JuChrom.alkane_series_step_spectra(
        msm,
        variances,
        apexes;
        carbons=9:10,
        scanwindow=2,
        mzretentionkwargs=SPECTRA_MZRETENTION_KWARGS,
        variancefloor=1e-6,
    )
    @test ranged.success
    @test ranged.carbonnumbers == [9, 10]
    @test sort(collect(keys(ranged.spectra))) == [9, 10]
    @test ranged.settings.carbons == [9, 10]

    absent = JuChrom.alkane_series_step_spectra(
        msm,
        variances,
        apexes;
        carbons=11,
        scanwindow=2,
    )
    @test !absent.success
    @test isempty(absent.spectra)
    @test absent.settings.carbons == [11]
    @test occursin("requested carbon numbers", absent.failurereason)

    @test_throws ArgumentError JuChrom.alkane_series_step_spectra(
        msm,
        variances,
        apexes;
        carbons=Int[],
    )
    @test_throws ArgumentError JuChrom.alkane_series_step_spectra(
        msm,
        variances,
        apexes;
        carbons=[8, 9.5],
    )

    failed = JuChrom.alkane_series_step_spectra(
        msm,
        variances,
        (
            success=false,
            failurereason="series path failed",
            results=NamedTuple[],
        );
        scanwindow=2,
    )
    @test !failed.success
    @test isempty(failed.spectra)
end
