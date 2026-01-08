
module TestMzchrom

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core/convert.jl (including the last function)
# ─────────────────────────────────────────────────────────────────────────────

using Test
using JuChrom
using Unitful
using Logging

using JuChrom: mscanmatrix, mzchrom,
               MassScan, MassScanSeries, ChromScanSeries,
               retentions, intensities,
               intensityunit, instrument, acquisition, user, sample, extras

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

# A: No explicit levels; unitless m/z
function _toy_series_no_levels()
    s1 = MassScan(1.0u"s", [101.2, 101.6, 102.1], [10, 20, 30])
    s2 = MassScan(2.0u"s", [100.9, 101.2, 102.6], [20, 5, 30])
    MassScanSeries([s1, s2];
        instrument=(detector="Orbitrap", manufacturer="Thermo"),
        acquisition=(mode="FullScan", method="DDA", polarity="positive"),
        user=(operator="Alice", project="QC-2025"),
        sample=(ID="sample_001", matrix="plasma", prep="SPE"),
        extras=Dict("injection_volume"=>5.0, "comment"=>"QC run"))
end

# B: Same spectra but with *unitful* m/z using Thomson (Th)
function _toy_series_unitful_mz()
    s1 = MassScan(1.0u"s", [101.2, 101.6, 102.1] .* u"Th", [10, 20, 30])
    s2 = MassScan(2.0u"s", [100.9, 101.2, 102.6] .* u"Th", [20, 5, 30])
    MassScanSeries([s1, s2];
        instrument=(detector="Orbitrap", manufacturer="Thermo"),
        acquisition=(mode="FullScan", method="DDA", polarity="positive"),
        user=(operator="Alice", project="QC-2025"),
        sample=(ID="sample_001", matrix="plasma", prep="SPE"),
        extras=Dict("injection_volume"=>5.0, "comment"=>"QC run"))
end

# ─────────────────────────────────────────────────────────────────────────────
# mzchrom(series, ...)
# ─────────────────────────────────────────────────────────────────────────────

@testset "mzchrom(series) — TIC/XIC, tol, units, warnings, metadata" begin
    # A) Unitless m/z series
    seriesA = _toy_series_no_levels()  # ret: [1,2] s; m/z: [100.9,101.2,101.6,102.1,102.6]
    s1 = scans(seriesA)[1]; s2 = scans(seriesA)[2]

    # TIC: sums across all m/z per scan
    css_tic = mzchrom(seriesA, warning=false)
    @test css_tic isa ChromScanSeries
    @test retentions(css_tic) == retentions(seriesA)  # retention preserved
    @test intensityunit(css_tic) === intensityunit(seriesA)  # intensity unit preserved
    @test intensities(css_tic) == [sum(intensities(s1)), sum(intensities(s2))]

    # XIC scalar exact match: 101.2 → [10, 5]
    css_xic_scalar = mzchrom(seriesA, 101.2, warning=false)
    @test intensities(css_xic_scalar) == [10, 5]

    # XIC vector exact match (no cross-scan nearest): [101.2, 102.1] → [40, 5]
    css_xic_vec = mzchrom(seriesA, [101.2, 102.1], warning=false)
    @test intensities(css_xic_vec) == [10 + 30, 5 + 0]

    # tol behavior (numeric):
    # - too tight (no match) → zeros
    css_tol_tight = mzchrom(seriesA, 101.201; tol=1e-5, warning=false)
    @test intensities(css_tol_tight) == [0, 0]
    # - loose enough (nearest within tol) → equals exact match
    css_tol_ok = mzchrom(seriesA, 101.201; tol=1e-2, warning=false)
    @test intensities(css_tol_ok) == [10, 5]

    # warning path when no match within tol
    logger = TestLogger()
    with_logger(logger) do
        css_warn = mzchrom(seriesA, 999.0; tol=1e-6, warning=true)
        @test intensities(css_warn) == [0, 0]
    end

    @test any(l -> l.level == Logging.Warn &&
                occursin("No m/z match", string(l.message)),
            logger.logs)

    # unitful tol not allowed when series m/z are unitless
    @test_throws ArgumentError mzchrom(seriesA, 101.2; tol=1e-3u"Th")

    # B) Unitful m/z series (Th)
    seriesU = _toy_series_unitful_mz()  # m/z carry u"Th"
    # scalar selection as Quantity
    css_u_scalar = mzchrom(seriesU, 101.2u"Th")
    @test intensities(css_u_scalar) == [10, 5]
    # numeric tol is interpreted in Th
    css_u_numtol = mzchrom(seriesU, 101.201; tol=1e-2, warning=false)
    @test intensities(css_u_numtol) == [10, 5]
    # unitful tol converted to Th and honored
    css_u_qu_tol = mzchrom(seriesU, 101.201u"Th"; tol=1e-3u"Th", warning=false)
    @test intensities(css_u_qu_tol) == [10, 5]

    # incompatible unit for tol should throw
    @test_throws ArgumentError mzchrom(seriesU, 101.2u"Th"; tol=1e-3u"s")

    # metadata propagation
    css_meta = mzchrom(seriesA)
    @test instrument(css_meta) == instrument(seriesA)
    @test acquisition(css_meta) == acquisition(seriesA)
    @test user(css_meta) == user(seriesA)
    @test sample(css_meta) == sample(seriesA)
    @test extras(css_meta) == extras(seriesA)

    # C) Intensity units preserved (positive control)
    # Build a tiny series with unitful intensities to verify propagation
    sI1 = MassScan(1.0u"s", [100.0], [10.0]u"pA")
    sI2 = MassScan(2.0u"s", [100.0], [20.0]u"pA")
    seriesI = MassScanSeries([sI1, sI2])
    cssI = mzchrom(seriesI, 100.0)
    @test intensityunit(cssI) == u"pA"
    @test intensities(cssI) == [10.0, 20.0]u"pA"
end

@testset "mzchrom(msm) — TIC/XIC, by=:mz|:index, tol, dedup, warnings, metadata" begin
    # Build matrices from the same toy series used for the series tests
    seriesA = _toy_series_no_levels()  # unitless m/z
    msmA = mscanmatrix(seriesA)  # dense; grid: [100.9, 101.2, 101.6, 102.1, 102.6]
    I = intensities(msmA)
    rts = retentions(msmA)

    # Sanity: TIC equals per-row sums
    css_tic = mzchrom(msmA; warning=false)
    @test css_tic isa ChromScanSeries
    @test retentions(css_tic) == rts
    @test intensityunit(css_tic) === intensityunit(msmA)
    @test intensities(css_tic) == vec(sum(I, dims=2))

    # XIC by m/z (scalar, exact)
    css_mz_scalar = mzchrom(msmA, 101.2; by=:mz, warning=false)
    @test intensities(css_mz_scalar) == [I[1,2], I[2,2]]  # [10, 5]

    # XIC by m/z (vector, exact)
    css_mz_vec = mzchrom(msmA, [101.2, 102.1]; by=:mz, warning=false)
    @test intensities(css_mz_vec) == [I[1,2] + I[1,4], I[2,2] + I[2,4]]  # [40, 5]

    # tol behavior (numeric): too tight → zeros; loose enough → nearest matches
    css_tol_tight = mzchrom(msmA, 101.201; by=:mz, tol=1e-5, warning=false)
    @test intensities(css_tol_tight) == [0, 0]
    css_tol_ok = mzchrom(msmA, 101.201; by=:mz, tol=1e-2, warning=false)
    @test intensities(css_tol_ok) == [I[1,2], I[2,2]]  # [10, 5]

    # WARNING path (no match within tol) — capture silently
    logger = TestLogger()
    with_logger(logger) do
        css_warn = mzchrom(msmA, 999.0; by=:mz, tol=1e-6, warning=true)
        @test intensities(css_warn) == [0, 0]
    end
    @test any(l -> l.level == Logging.Warn && occursin("No m/z match", string(l.message)), 
              logger.logs)

    # Unitless m/z grid: unitful tol must throw
    @test_throws ArgumentError mzchrom(msmA, 101.2; by=:mz, tol=1e-3u"Th")

    # Matrix-specific: by=:index
    # Scalar index selection equals the same column as m/z selection above
    css_idx_scalar = mzchrom(msmA, 2; by=:index, warning=false)  # 1-based
    @test intensities(css_idx_scalar) == intensities(css_mz_scalar)

    # Vector of indices (exact); order/duplicates/out-of-bounds are handled:
    # - duplicates should be deduped
    # - out-of-bounds resolved to 0 then filtered out
    css_idx_vec = mzchrom(msmA, [2, 2, 999, -1]; by=:index, warning=false)
    @test intensities(css_idx_vec) == intensities(css_idx_scalar)

    # Sorting of resolved index vector should not affect sums
    css_idx_unsorted = mzchrom(msmA, [4, 2]; by=:index, warning=false)
    css_idx_sorted = mzchrom(msmA, [2, 4]; by=:index, warning=false)
    @test intensities(css_idx_unsorted) == intensities(css_idx_sorted)

    # Matrix-specific: dedup for by=:mz with near-equal targets mapping to same bin
    # Both selections resolve to column 2 (101.2) within tol; ensure not double-counted
    css_mz_dedup = mzchrom(msmA, [101.2001, 101.2002]; by=:mz, tol=1e-2, warning=false)
    @test intensities(css_mz_dedup) == [I[1,2], I[2,2]]  # not doubled

    # Unitful m/z grid matrix
    seriesU = _toy_series_unitful_mz()  # m/z with u"Th"
    msmU = mscanmatrix(seriesU)
    IU = intensities(msmU)

    # Scalar selection as Quantity
    cssU_q = mzchrom(msmU, 101.2u"Th"; by=:mz, warning=false)
    @test intensities(cssU_q) == [IU[1,2], IU[2,2]]  # [10, 5]

    # Numeric tol interpreted in Th
    cssU_tol_num = mzchrom(msmU, 101.201; by=:mz, tol=1e-2, warning=false)
    @test intensities(cssU_tol_num) == [IU[1,2], IU[2,2]]

    # Unitful tol converted to Th
    cssU_tol_qu = mzchrom(msmU, 101.201u"Th"; by=:mz, tol=1e-3u"Th", warning=false)
    @test intensities(cssU_tol_qu) == [IU[1,2], IU[2,2]]

    # Incompatible unit for tol should throw (catch as ArgumentError per API)
    @test_throws ArgumentError mzchrom(msmU, 101.2u"Th"; by=:mz, tol=1e-3u"s")

    # Metadata propagation from matrix to ChromScanSeries
    css_meta = mzchrom(msmA; warning=false)
    @test instrument(css_meta) == instrument(msmA)
    @test acquisition(css_meta) == acquisition(msmA)
    @test user(css_meta) == user(msmA)
    @test sample(css_meta) == sample(msmA)
    @test extras(css_meta) == extras(msmA)
end

end # module ConvertTests
