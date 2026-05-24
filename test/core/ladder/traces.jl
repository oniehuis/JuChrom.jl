using Test
using JuChrom

@testset "alkane trace inference" begin
    @testset "alkanematchtrace" begin
        reference = [1.0, 0.74324805, 0.21725361]
        rawcounts = [
            reference';
            0.0 0.0 0.0;
            0.0 0.0 1.0;
        ]
        msm = MassScanMatrix(
            [1.0, 2.0, 3.0],
            [43.0, 57.0, 70.0],
            rawcounts;
            sample=(name="trace test",),
            extras=Dict("source" => "unit"),
        )

        trace = JuChrom.alkanematchtrace(
            msm,
            8;
            spectralpower=1.0,
            smoothing=0,
        )
        expected = [
            1.0,
            0.0,
            cossim([0.0, 0.0, 1.0], reference),
        ]

        @test trace isa ChromScanSeries
        @test rawretentions(trace) == [1.0, 2.0, 3.0]
        @test rawintensities(trace) ≈ expected
        @test intensityunit(trace) === nothing
        @test sample(trace) == sample(msm)
        @test extras(trace)["source"] == "unit"
        @test extras(trace)["trace_type"] == "alkane_match"
        @test extras(trace)["carbon"] == 8

        traceattrs = attrs(first(trace))
        @test traceattrs.carbon == 8
        @test traceattrs.referenceattrs.label == "octane"
        @test traceattrs.matchedmz == [43.0, 57.0, 70.0]
        @test traceattrs.mzunit === nothing
        @test traceattrs.minrelativeintensity == 0.05
        @test traceattrs.spectralpower == 1.0
        @test traceattrs.positiveonly
        @test traceattrs.smoothing == 0

        channelinfo = JuChrom.alkane_mz_channels(msm; carbonrange=8:8)
        reused = JuChrom.alkanematchtrace(
            msm,
            8;
            channelinfo=channelinfo,
            spectralpower=1.0,
            smoothing=0,
        )
        @test rawintensities(reused) ≈ rawintensities(trace)
        @test_throws ArgumentError JuChrom.alkanematchtrace(
            msm,
            9;
            channelinfo=channelinfo,
            spectralpower=1.0,
            smoothing=0,
        )

        unitful = MassScanMatrix(
            [1.0, 2.0, 3.0],
            [0.043, 0.057, 0.070]u"kTh",
            rawcounts .* u"pA",
        )
        unitful_trace = JuChrom.alkanematchtrace(
            unitful,
            8;
            spectralpower=1.0,
            smoothing=0,
        )
        @test rawintensities(unitful_trace) ≈ expected
        @test intensityunit(unitful_trace) === nothing
        @test attrs(first(unitful_trace)).matchedmz == [0.043, 0.057, 0.070]u"kTh"
        @test attrs(first(unitful_trace)).mzunit == u"kTh"
        @test extras(unitful_trace)["mz_unit"] == u"kTh"

        noninteger = MassScanMatrix([1.0], [43.2], reshape([1.0], 1, 1))
        @test_throws ArgumentError JuChrom.alkanematchtrace(noninteger, 8)
        incompatible_mz_unit = MassScanMatrix(
            [1.0],
            [43.0]u"s",
            reshape([1.0], 1, 1),
        )
        @test_throws ArgumentError JuChrom.alkanematchtrace(incompatible_mz_unit, 8)
        @test_throws ArgumentError JuChrom.alkanematchtrace(msm, 8; spectralpower=0.0)
        @test_throws ArgumentError JuChrom.alkanematchtrace(msm, 8; smoothing=-1)
    end

    @testset "alkane_molecular_ion_trace" begin
        rawcounts = [
            1.0 10.0 2.0;
            4.0 5.0 1.0;
            0.0 0.0 0.0;
        ]
        msm = MassScanMatrix(
            [1.0, 2.0, 3.0],
            [100.0, 114.0, 128.0],
            rawcounts;
            extras=Dict("source" => "unit"),
        )

        trace = JuChrom.alkane_molecular_ion_trace(
            msm,
            8;
            ion_window=0,
            normalize=false,
            smoothing=0,
        )

        @test rawintensities(trace) == [7.0, 0.0, 0.0]
        @test intensityunit(trace) === nothing
        @test extras(trace)["source"] == "unit"
        @test extras(trace)["trace_type"] == "alkane_molecular_ion"
        @test extras(trace)["molecular_ion"] == 114

        traceattrs = attrs(first(trace))
        @test traceattrs.carbon == 8
        @test traceattrs.molecularion == 114
        @test traceattrs.centermz == [114.0]
        @test traceattrs.lowermz == [100.0]
        @test traceattrs.uppermz == [128.0]
        @test traceattrs.mzunit === nothing
        @test !traceattrs.normalized

        normalized = JuChrom.alkane_molecular_ion_trace(
            msm,
            8;
            ion_window=0,
            normalize=true,
            smoothing=0,
        )
        @test rawintensities(normalized) == [1.0, 0.0, 0.0]
        @test intensityunit(normalized) === nothing

        unitful = MassScanMatrix(
            [1.0, 2.0, 3.0],
            [100.0, 114.0, 128.0]u"Th",
            rawcounts,
        )
        unitful_trace = JuChrom.alkane_molecular_ion_trace(
            unitful,
            8;
            ion_window=0,
            normalize=false,
            smoothing=0,
        )
        @test attrs(first(unitful_trace)).centermz == [114.0]u"Th"
        @test attrs(first(unitful_trace)).mzunit == u"Th"
        @test extras(unitful_trace)["mz_unit"] == u"Th"

        unitful_kth = MassScanMatrix(
            [1.0, 2.0, 3.0],
            [0.100, 0.114, 0.128]u"kTh",
            rawcounts,
        )
        unitful_kth_trace = JuChrom.alkane_molecular_ion_trace(
            unitful_kth,
            8;
            ion_window=0,
            normalize=false,
            smoothing=0,
        )
        @test rawintensities(unitful_kth_trace) == [7.0, 0.0, 0.0]
        @test attrs(first(unitful_kth_trace)).centermz == [0.114]u"kTh"
        @test attrs(first(unitful_kth_trace)).mzunit == u"kTh"

        noninteger = MassScanMatrix([1.0], [114.2], reshape([1.0], 1, 1))
        @test_throws ArgumentError JuChrom.alkane_molecular_ion_trace(noninteger, 8)
        incompatible_mz_unit = MassScanMatrix(
            [1.0],
            [114.0]u"s",
            reshape([1.0], 1, 1),
        )
        @test_throws ArgumentError JuChrom.alkane_molecular_ion_trace(
            incompatible_mz_unit,
            8,
        )
        @test_throws ArgumentError JuChrom.alkane_molecular_ion_trace(msm, 8; ion_window=-1)
        @test_throws ArgumentError JuChrom.alkane_molecular_ion_trace(msm, 8; step_mass=0)
        @test_throws ArgumentError JuChrom.alkane_molecular_ion_trace(msm, 8; smoothing=-1)
    end

    @testset "alkanemzpeakdistancetrace" begin
        rawcounts = [
            0.0 0.0 0.0;
            10.0 10.0 10.0;
            0.0 0.0 0.0;
        ]
        variances = ones(size(rawcounts))
        msm = MassScanMatrix(
            [0.0, 1.0, 2.0],
            [43.0, 57.0, 70.0],
            rawcounts;
            sample=(name="distance test",),
            extras=Dict("source" => "unit"),
        )
        mzretentionkwargs = (
            retentionref=:start,
            scaninterval=1.0,
            mzcount=3,
            order=:descending,
            dwellref=:middle,
            dwell=:homogeneous,
        )

        trace = JuChrom.alkanemzpeakdistancetrace(
            msm,
            8;
            variances=variances,
            mzretentionkwargs=mzretentionkwargs,
            smoothing=0,
        )

        @test trace isa ChromScanSeries
        @test rawretentions(trace) == [0.0, 1.0, 2.0]
        @test rawintensities(trace) ≈ [0.5, 1.0, 0.5]
        @test intensityunit(trace) === nothing
        @test sample(trace) == sample(msm)
        @test extras(trace)["source"] == "unit"
        @test extras(trace)["trace_type"] == "alkane_mzpeak_distance_evidence"
        @test extras(trace)["carbon"] == 8
        @test extras(trace)["scan_interval"] == 1.0
        @test extras(trace)["retention_unit"] === nothing

        traceattrs = attrs(first(trace))
        @test traceattrs.carbon == 8
        @test traceattrs.referenceattrs.label == "octane"
        @test traceattrs.selectedmz == [43.0, 57.0, 70.0]
        @test traceattrs.selectedmzindices == [1, 2, 3]
        @test traceattrs.mzunit === nothing
        @test traceattrs.scaninterval == 1.0
        @test traceattrs.retentionunit === nothing
        @test traceattrs.selectedpeakcounts == [1, 1, 1]
        @test traceattrs.peakzmin == 4.0
        @test traceattrs.transform == :inverse_quadratic
        @test traceattrs.smoothing == 0
        @test traceattrs.scaledbetweenzeroandone

        flatcounts = ones(3, 3)
        flatmsm = MassScanMatrix([0.0, 1.0, 2.0], [43.0, 57.0, 70.0], flatcounts)
        empty_trace = JuChrom.alkanemzpeakdistancetrace(
            flatmsm,
            8;
            variances=ones(size(flatcounts)),
            mzretentionkwargs=mzretentionkwargs,
            smoothing=0,
        )
        @test rawintensities(empty_trace) == [0.0, 0.0, 0.0]

        unitful = MassScanMatrix(
            [0.0, 60.0, 120.0]u"s",
            [0.043, 0.057, 0.070]u"kTh",
            rawcounts .* u"pA",
        )
        unitful_kwargs = (
            retentionref=:start,
            scaninterval=60.0u"s",
            mzcount=3,
            order=:descending,
            dwellref=:middle,
            dwell=:homogeneous,
        )
        unitful_trace = JuChrom.alkanemzpeakdistancetrace(
            unitful,
            8;
            variances=variances,
            mzretentionkwargs=unitful_kwargs,
            smoothing=0,
        )
        @test rawintensities(unitful_trace) ≈ [0.5, 1.0, 0.5]
        @test retentionunit(unitful_trace) == u"s"
        @test intensityunit(unitful_trace) === nothing
        @test attrs(first(unitful_trace)).selectedmz == [0.043, 0.057, 0.070]u"kTh"
        @test attrs(first(unitful_trace)).mzunit == u"kTh"
        @test attrs(first(unitful_trace)).retentionunit == u"s"
        @test attrs(first(unitful_trace)).scaninterval == 60.0
        @test extras(unitful_trace)["mz_unit"] == u"kTh"
        @test extras(unitful_trace)["retention_unit"] == u"s"
        @test_throws ArgumentError JuChrom.raw_retention_value(unitful, 1.0u"Th")

        high_threshold_trace = JuChrom.alkanemzpeakdistancetrace(
            msm,
            8;
            variances=variances,
            mzretentionkwargs=mzretentionkwargs,
            peakzmin=100.0,
            smoothing=0,
        )
        @test rawintensities(high_threshold_trace) == [0.0, 0.0, 0.0]
        @test attrs(first(high_threshold_trace)).selectedpeakcounts == [0, 0, 0]

        @test_throws UndefKeywordError JuChrom.alkanemzpeakdistancetrace(msm, 8)
        @test_throws DimensionMismatch JuChrom.alkanemzpeakdistancetrace(
            msm,
            8;
            variances=ones(size(rawcounts, 1), size(rawcounts, 2) + 1),
            mzretentionkwargs=mzretentionkwargs,
        )
        bad_variances = copy(variances)
        bad_variances[1, 1] = -1.0
        @test_throws ArgumentError JuChrom.alkanemzpeakdistancetrace(
            msm,
            8;
            variances=bad_variances,
            mzretentionkwargs=mzretentionkwargs,
        )
        @test_throws ArgumentError JuChrom.alkanemzpeakdistancetrace(
            msm,
            8;
            variances=variances,
            mzretentionkwargs=mzretentionkwargs,
            smoothing=-1,
        )
    end

    @testset "alkaneevidencetrace" begin
        function evidence_test_trace(values; carbon=8, retentions=[1.0, 2.0, 3.0])
            chroms = [
                ChromScan(rt, nothing, y, nothing; attrs=(carbon=carbon,))
                for (rt, y) in zip(retentions, values)
            ]
            ChromScanSeries(
                chroms;
                sample=(name="evidence test",),
                extras=Dict("source" => "unit", "carbon" => carbon),
            )
        end

        function evidence_test_trace_without_carbon_attrs(
            values;
            carbon=nothing,
            retentions=[1.0, 2.0, 3.0],
        )
            chroms = [
                ChromScan(rt, nothing, y, nothing)
                for (rt, y) in zip(retentions, values)
            ]
            trace_extras = Dict{String, Any}("source" => "unit")
            isnothing(carbon) || (trace_extras["carbon"] = carbon)
            ChromScanSeries(
                chroms;
                sample=(name="evidence test",),
                extras=trace_extras,
            )
        end

        molecularion = evidence_test_trace([1.0, 0.5, 0.25])
        match = evidence_test_trace([0.2, 0.4, 0.8])
        distance = evidence_test_trace([0.5, 0.25, 0.0])

        trace = JuChrom.alkaneevidencetrace(
            molecularion,
            match;
            mzpeakdistancetrace=distance,
        )

        @test rawintensities(trace) == [0.1, 0.05, 0.0]
        @test intensityunit(trace) === nothing
        @test sample(trace) == sample(molecularion)
        @test extras(trace)["source"] == "unit"
        @test extras(trace)["trace_type"] == "alkane_evidence_product"
        @test extras(trace)["carbon"] == 8
        @test extras(trace)["distance_evidence_included"]
        @test attrs(first(trace)).carbon == 8
        @test attrs(first(trace)).distanceevidenceincluded

        without_distance = JuChrom.alkaneevidencetrace(molecularion, match)
        @test rawintensities(without_distance) == [0.2, 0.2, 0.2]
        @test !extras(without_distance)["distance_evidence_included"]
        @test !attrs(first(without_distance)).distanceevidenceincluded

        extras_carbon_molecularion = evidence_test_trace_without_carbon_attrs(
            [1.0, 0.5, 0.25];
            carbon=Int8(9),
        )
        extras_carbon_trace = JuChrom.alkaneevidencetrace(
            extras_carbon_molecularion,
            match,
        )
        @test attrs(first(extras_carbon_trace)).carbon === 9
        @test extras(extras_carbon_trace)["carbon"] === 9

        missing_carbon_molecularion = evidence_test_trace_without_carbon_attrs(
            [1.0, 0.5, 0.25],
        )
        missing_carbon_trace = JuChrom.alkaneevidencetrace(
            missing_carbon_molecularion,
            match,
        )
        @test ismissing(attrs(first(missing_carbon_trace)).carbon)
        @test ismissing(extras(missing_carbon_trace)["carbon"])

        mismatched_length = evidence_test_trace([1.0, 2.0])
        @test_throws DimensionMismatch JuChrom.alkaneevidencetrace(
            molecularion,
            mismatched_length,
        )

        mismatched_retention = evidence_test_trace([1.0, 2.0, 3.0]; retentions=[1.0, 2.0, 4.0])
        @test_throws DimensionMismatch JuChrom.alkaneevidencetrace(
            molecularion,
            mismatched_retention,
        )
    end

    @testset "trace helper branches" begin
        @test JuChrom.spectral_power_transform_value(-4.0, 0.5, false) == -2.0
        @test JuChrom.spectral_power_transform_value(3, 2, false) == 9.0
    end
end
