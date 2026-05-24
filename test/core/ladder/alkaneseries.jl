using Test
using JuChrom

function alkaneseries_test_matrix()
    validwindow(base, step) = base .+ [0.0, step, -step, step, -step]
    vcat(
        [0.0 10.0 10.0;
         10.0 10.0 10.0;
         10.0 10.0 10.0;
         10.0 10.0 10.0;
         10.0 10.0 10.0],
        hcat(validwindow(30.0, 3.0), validwindow(50.0, 5.0), validwindow(70.0, 7.0)),
        hcat(validwindow(35.0, 3.5), validwindow(55.0, 5.5), validwindow(75.0, 7.5)),
        hcat(validwindow(40.0, 4.0), validwindow(60.0, 6.0), validwindow(80.0, 8.0)),
    )
end

@testset "alkane reference spectrum helper" begin
    mzs = Int16[43, 57, 71]
    relative_abundances = Float32[0.25, 1.0, 0.5]
    label = SubString("pentadecane reference", 1, 11)

    spectrum = JuChrom.alkane_reference_spectrum(
        Int8(15),
        label,
        1500//1,
        mzs,
        relative_abundances,
    )

    @test spectrum isa MassSpectrum
    @test mzvalues(spectrum) == [43, 57, 71]
    @test eltype(mzvalues(spectrum)) == Int
    @test intensities(spectrum) == [0.25, 1.0, 0.5]
    @test eltype(intensities(spectrum)) == Float64

    spectrum_attrs = attrs(spectrum)
    @test spectrum_attrs.ri === 1500.0
    @test spectrum_attrs.order === 15
    @test spectrum_attrs.label == "pentadecane"
    @test spectrum_attrs.label isa String
    @test spectrum_attrs.source == :experimental_ladder_spectra
    @test spectrum_attrs.normalization == :base_peak
    @test spectrum_attrs.relative_abundance_threshold ==
        JuChrom.ALKANE_REFERENCE_RELATIVE_ABUNDANCE_THRESHOLD

    mzs[1] = 41
    relative_abundances[1] = 0.75
    @test mzvalues(spectrum) == [43, 57, 71]
    @test intensities(spectrum) == [0.25, 1.0, 0.5]

    @test_throws DimensionMismatch JuChrom.alkane_reference_spectrum(
        15,
        "invalid",
        1500.0,
        [43, 57],
        [1.0],
    )
end

@testset "AlkaneSeriesResult convenience constructors" begin
    standard = (name="test standard",)
    variances = [1.0 2.0; 3.0 4.0]
    varianceinfo = (estimator=:unit,)
    baselineinfo = (estimator=:baseline,)
    channelinfo = (mzindices=[1, 2],)
    traces = (evidence=Dict{Int, ChromScanSeries}(),)

    with_varianceinfo = AlkaneSeriesResult(standard, variances, varianceinfo)
    @test with_varianceinfo.standard === standard
    @test with_varianceinfo.variances === variances
    @test with_varianceinfo.varianceinfo === varianceinfo
    @test with_varianceinfo.baselineinfo === nothing
    @test with_varianceinfo.channelinfo === nothing
    @test with_varianceinfo.traces === nothing
    @test with_varianceinfo.seriespath === nothing

    with_baselineinfo = AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        baselineinfo,
    )
    @test with_baselineinfo.baselineinfo === baselineinfo
    @test with_baselineinfo.channelinfo === nothing
    @test with_baselineinfo.traces === nothing
    @test with_baselineinfo.seriespath === nothing

    with_channelinfo = AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        baselineinfo,
        channelinfo,
    )
    @test with_channelinfo.channelinfo === channelinfo
    @test with_channelinfo.traces === nothing
    @test with_channelinfo.seriespath === nothing

    with_traces = AlkaneSeriesResult(
        standard,
        variances,
        varianceinfo,
        baselineinfo,
        channelinfo,
        traces,
    )
    @test with_traces.traces === traces
    @test with_traces.seriespath === nothing
end

@testset "alkane series preprocessing helpers" begin
    varianceestimate = [1.0 2.0; 3.0 4.0]
    extracted_variances, extracted_info =
        JuChrom.extract_alkane_series_variances(varianceestimate)

    @test extracted_variances === varianceestimate
    @test extracted_info === nothing
end

@testset "findalkaneseries preprocessing boundary" begin
    rawcounts = alkaneseries_test_matrix()
    msm = MassScanMatrix(collect(1.0:size(rawcounts, 1)), [50.0, 60.0, 70.0], rawcounts)
    standard = defaultalkanestandard()
    variancekwargs = (
        variancewindowsize=5,
        variancemintransitioncount=1,
        variancezerothresholdquantile=1.0,
    )

    @testset "built-in reference spectra" begin
        spectra = alkanereferencespectra()
        standard_default = defaultalkanestandard()

        @test standard_default isa AlkaneStandard
        @test standard_default.name == "n-alkane C8-C40 reference spectra"
        @test standard_default.spectra !== spectra
        @test length(spectra) == 33
        @test [attrs(spectrum).order for spectrum in spectra] == collect(8:40)
        @test all(spectrum -> spectrum isa MassSpectrum, spectra)

        for spectrum in spectra
            spectrum_attrs = attrs(spectrum)
            @test spectrum_attrs.source == :experimental_ladder_spectra
            @test spectrum_attrs.normalization == :base_peak
            @test spectrum_attrs.relative_abundance_threshold == 0.01
            @test maximum(intensities(spectrum)) ≈ 1.0
            @test all(>(0.01), intensities(spectrum))
            @test all(>(0), mzvalues(spectrum))
            @test all(>(0), diff(mzvalues(spectrum)))
        end

        c8 = alkanereferencespectrum(8)
        @test attrs(c8).label == "octane"
        @test first(mzvalues(c8)) == 29

        intensities(c8)[1] = 99.0
        @test first(intensities(alkanereferencespectrum(8))) != 99.0
        @test_throws ArgumentError alkanereferencespectrum(7)
    end

    @testset "reference m/z channels are matched to msm channels" begin
        channels = JuChrom.alkane_mz_channels(
            msm;
            standard=standard,
            carbonrange=8:10,
            minrelativeintensity=0.05,
        )

        @test channels.mzindices == [3]
        @test channels.mzvalues == [70.0]
        @test channels.carbonrange == [8, 9, 10]
        @test channels.minrelativeintensity == 0.05
        @test [reference.carbon for reference in channels.references] == [8, 9, 10]
        @test all(reference -> reference.mzindices == [3], channels.references)
        @test all(reference -> reference.mzvalues == [70.0], channels.references)

        vector_standard_channels = JuChrom.alkane_mz_channels(
            msm;
            standard=alkanereferencespectra(),
            carbonrange=8:8,
            minrelativeintensity=0.05,
        )
        @test only(vector_standard_channels.references).carbon == 8
        @test only(vector_standard_channels.references).referenceattrs.label == "octane"

        c8 = first(channels.references)
        @test c8.referenceattrs.label == "octane"
        @test c8.referenceintensities ≈ [0.21725361]

        unitful_msm = MassScanMatrix(
            collect(1.0:size(rawcounts, 1)),
            [0.05, 0.06, 0.07]u"kTh",
            rawcounts,
        )
        unitful_channels = JuChrom.alkane_mz_channels(
            unitful_msm;
            standard=standard,
            carbonrange=8:8,
            minrelativeintensity=0.05,
        )
        @test unitful_channels.mzindices == [3]
        @test unitful_channels.mzvalues == [0.07]u"kTh"
        @test unitful_channels.mzunit == u"kTh"
        @test only(unitful_channels.references).mzvalues == [0.07]u"kTh"

        noninteger_msm = MassScanMatrix(
            collect(1.0:size(rawcounts, 1)),
            [69.95],
            rawcounts[:, 1:1],
        )
        @test_throws ArgumentError JuChrom.alkane_mz_channels(
            noninteger_msm;
            standard=standard,
            carbonrange=8:8,
        )
        incompatible_mz_unit_msm = MassScanMatrix(
            collect(1.0:size(rawcounts, 1)),
            [50.0, 60.0, 70.0]u"s",
            rawcounts,
        )
        @test_throws ArgumentError JuChrom.alkane_mz_channels(
            incompatible_mz_unit_msm;
            standard=standard,
            carbonrange=8:8,
        )

        @test_throws ArgumentError JuChrom.alkane_mz_channels(
            MassScanMatrix(collect(1.0:size(rawcounts, 1)), [500.0], rawcounts[:, 1:1]);
            standard=standard,
            carbonrange=8:8,
        )
        @test_throws ArgumentError JuChrom.alkane_mz_channels(
            msm;
            standard=(name="missing spectra",),
            carbonrange=8:8,
        )
        @test_throws ArgumentError JuChrom.alkane_mz_channels(
            msm;
            standard=standard,
            carbonrange=8:8,
            minrelativeintensity=1.0,
        )
    end

    @testset "core analysis requires provided variances" begin
        σ² = fill(2.0, size(rawcounts))
        result = findalkaneseries(msm, σ²; standard=standard)
        defaultresult = findalkaneseries(msm, σ²)

        @test result isa AlkaneSeriesResult
        @test result.standard === standard
        @test result.variances === σ²
        @test result.varianceinfo === nothing
        @test result.baselineinfo === nothing
        @test result.channelinfo.mzindices == [3]
        @test length(result.channelinfo.references) == 33
        @test result.traces.carbonrange == collect(8:40)
        @test sort(collect(keys(result.traces.match))) == collect(8:40)
        @test sort(collect(keys(result.traces.molecularion))) == collect(8:40)
        @test sort(collect(keys(result.traces.mzpeakdistance))) == collect(8:40)
        @test sort(collect(keys(result.traces.evidence))) == collect(8:40)
        @test result.traces.match[8] isa ChromScanSeries
        @test result.traces.molecularion[8] isa ChromScanSeries
        @test result.traces.mzpeakdistance[8] isa ChromScanSeries
        @test result.traces.evidence[8] isa ChromScanSeries
        @test result.seriespath !== nothing
        @test !result.seriespath.success
        @test result.seriespath.failurereason ==
            "no local maxima found in evidence traces"
        @test sort(collect(keys(result.seriespath.candidatesbycarbon))) == collect(8:40)
        @test all(isempty, values(result.seriespath.candidatesbycarbon))
        @test rawretentions(result.traces.match[8]) == rawretentions(msm)
        @test rawretentions(result.traces.molecularion[8]) == rawretentions(msm)
        @test rawretentions(result.traces.mzpeakdistance[8]) == rawretentions(msm)
        @test rawretentions(result.traces.evidence[8]) == rawretentions(msm)
        @test extras(result.traces.evidence[8])["distance_evidence_included"]
        @test rawintensities(result.traces.evidence[8]) ≈
            rawintensities(result.traces.molecularion[8]) .*
            rawintensities(result.traces.match[8]) .*
            rawintensities(result.traces.mzpeakdistance[8])
        @test defaultresult.standard isa AlkaneStandard
        @test defaultresult.channelinfo.mzindices == [3]
        @test_throws MethodError findalkaneseries(msm; standard=standard)
        @test_throws MethodError findalkaneseries(
            msm;
            standard=standard,
            variances=σ²,
        )
    end

    @testset "core can skip distance traces" begin
        σ² = fill(2.0, size(rawcounts))
        result = findalkaneseries(
            msm,
            σ²;
            standard=standard,
            carbonrange=8:8,
            includedistance=false,
        )

        @test result.traces.carbonrange == [8]
        @test sort(collect(keys(result.traces.match))) == [8]
        @test sort(collect(keys(result.traces.molecularion))) == [8]
        @test isempty(result.traces.mzpeakdistance)
        @test sort(collect(keys(result.traces.evidence))) == [8]
        @test !result.seriespath.success
        @test !extras(result.traces.evidence[8])["distance_evidence_included"]
        @test rawintensities(result.traces.evidence[8]) ≈
            rawintensities(result.traces.molecularion[8]) .*
            rawintensities(result.traces.match[8])
    end

    @testset "core stores upstream baseline info" begin
        σ² = fill(2.0, size(rawcounts))
        baselines = MassScanMatrix(
            collect(1.0:size(rawcounts, 1)),
            [50.0, 60.0, 70.0],
            fill(1.0, size(rawcounts)),
        )
        signal = msm - baselines
        varianceinfo = (source=:external,)
        baselineinfo = (source=:external, baselines=baselines)
        result = findalkaneseries(
            signal,
            σ²;
            standard=standard,
            varianceinfo=varianceinfo,
            baselineinfo=baselineinfo,
        )

        @test result.variances === σ²
        @test result.varianceinfo === varianceinfo
        @test result.baselineinfo === baselineinfo
    end

    @testset "raw-count wrapper estimates variances and stores baseline info" begin
        result = findalkanes(
            msm;
            standard=standard,
            baselineλ=1e3,
            variancekwargs...,
        )

        @test result isa AlkaneSeriesResult
        @test result.standard === standard
        @test result.varianceinfo !== nothing
        @test result.variances == result.varianceinfo.variances
        @test result.baselineinfo.estimator == :arpls
        @test result.baselineinfo.baselines isa MassScanMatrix
        @test size(rawintensities(result.baselineinfo.baselines)) == size(rawcounts)
        @test all(isfinite, rawintensities(result.baselineinfo.baselines))
    end

    @testset "raw-count wrapper requires variances when baseline subtraction is disabled" begin
        σ² = fill(2.0, size(rawcounts))
        result = findalkanes(
            msm;
            standard=standard,
            variances=σ²,
            subtractbaseline=false,
        )

        @test result.variances === σ²
        @test result.varianceinfo === nothing
        @test result.baselineinfo === nothing
        @test_throws ArgumentError findalkanes(
            msm;
            standard=standard,
            subtractbaseline=false,
        )
    end

    @testset "raw-count wrapper can use provided variances for baseline subtraction" begin
        σ² = fill(3.0, size(rawcounts))
        result = findalkanes(
            msm;
            standard=standard,
            variances=σ²,
            baselineλ=1e3,
        )

        @test result.variances === σ²
        @test result.varianceinfo === nothing
        @test result.baselineinfo.baselines isa MassScanMatrix
    end

    @testset "variance validation" begin
        @test_throws DimensionMismatch findalkanes(
            msm;
            standard=standard,
            variances=ones(size(rawcounts, 1), size(rawcounts, 2) + 1),
        )

        badvariances = fill(1.0, size(rawcounts))
        badvariances[1, 1] = -1.0
        @test_throws ArgumentError findalkanes(
            msm;
            standard=standard,
            variances=badvariances,
        )
    end
end
