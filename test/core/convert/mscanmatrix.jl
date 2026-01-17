module TestMscanmatrix

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core/convert.jl (including the last function)
# ─────────────────────────────────────────────────────────────────────────────

using Test
using JuChrom
using Unitful
using SparseArrays

using JuChrom: mscanmatrix, SPARSE,
               MassScan, MassScanSeries, MassScanMatrix,
               retentions, intensities, mzvalues, mzunit,
               intensityunit, retentionunit, instrument, acquisition, user, sample, extras,
               level

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

# B: Explicit different levels
function _toy_series_levels()
    s1 = MassScan(1.0u"s", [101.2, 101.6, 102.1], [10, 20, 30]; level=1)
    s2 = MassScan(2.0u"s", [100.9, 101.2, 102.6], [20, 5, 30]; level=2)
    MassScanSeries([s1, s2];
        instrument=(detector="Orbitrap", manufacturer="Thermo"),
        acquisition=(mode="FullScan", method="DDA", polarity="positive"),
        user=(operator="Alice", project="QC-2025"),
        sample=(ID="sample_001", matrix="plasma", prep="SPE"),
        extras=Dict("injection_volume"=>5.0, "comment"=>"QC run"))
end

struct DummyVal
    val::Float64
end

Base.float(x::DummyVal) = x.val

# ─────────────────────────────────────────────────────────────────────────────
# mscanmatrix(::MassScanSeries, [format]; target_level, threshold)
# ─────────────────────────────────────────────────────────────────────────────

@testset "mscanmatrix — dense/sparse, threshold, levels, metadata" begin
    # A) Default: no explicit levels → both scans at level=1; union m/z across scans.
    seriesA = _toy_series_no_levels()

    msmA_dense = mscanmatrix(seriesA)  # default dense
    @test msmA_dense isa MassScanMatrix
    @test size(intensities(msmA_dense)) == (2, 5)
    @test retentions(msmA_dense) == [1.0, 2.0]u"s"
    @test retentionunit(msmA_dense) == u"s"
    @test mzvalues(msmA_dense) ≈ [100.9, 101.2, 101.6, 102.1, 102.6]
    @test mzunit(msmA_dense) === nothing
    @test intensities(msmA_dense) ≈ [0 10 20 30 0; 20 5 0 0 30]
    @test intensityunit(msmA_dense) === intensityunit(seriesA)

    msmA_sparse = mscanmatrix(seriesA, SPARSE)
    @test intensities(msmA_sparse) isa SparseMatrixCSC
    @test Array(intensities(msmA_sparse)) == intensities(msmA_dense)

    # Threshold strictly ‘>’
    msmA_thr  = mscanmatrix(seriesA; threshold=7.5)
    @test intensities(msmA_thr)  ≈ [0 10 20 30 0; 20 0 0 0 30]
    msmA_thr2 = mscanmatrix(seriesA; threshold=10.0)
    @test intensities(msmA_thr2) ≈ [0 0 20 30 0; 20 0 0 0 30]  # 10 is NOT > 10

    # Metadata propagation
    @test instrument(msmA_dense) == instrument(seriesA)
    @test acquisition(msmA_dense) == acquisition(seriesA)
    @test user(msmA_dense) == user(seriesA)
    @test sample(msmA_dense) == sample(seriesA)
    @test extras(msmA_dense) == extras(seriesA)

    # B) Per-level behavior: explicit levels → 1×3 matrices for each level.
    seriesB = _toy_series_levels()

    msmB_lvl1 = mscanmatrix(seriesB; target_level=1)
    @test size(intensities(msmB_lvl1)) == (1, 3)
    @test mzvalues(msmB_lvl1) ≈ [101.2, 101.6, 102.1]
    @test intensities(msmB_lvl1) ≈ [10 20 30]

    msmB_lvl2 = mscanmatrix(seriesB; target_level=2)
    @test level(msmB_lvl2) == 2
    @test retentions(msmB_lvl2) ≈ [2.0]u"s"
    @test mzvalues(msmB_lvl2) ≈ [100.9, 101.2, 102.6]
    @test Array(intensities(msmB_lvl2)) ≈ [20 5 30]

    # Error if no scans at requested level
    @test_throws ArgumentError mscanmatrix(seriesB; target_level=3)

    # Element type follows input intensity type (Int)
    sI1 = MassScan(1.0u"s", [101.0],        Int[7];   level=1)
    sI2 = MassScan(2.0u"s", [101.0, 102.0], Int[3,9]; level=1)
    series_int = MassScanSeries([sI1, sI2])
    msm_int = mscanmatrix(series_int; target_level=1)
    @test eltype(intensities(msm_int)) == Int
    @test intensities(msm_int) == [7 0; 3 9]
end

end # module ConvertTests
