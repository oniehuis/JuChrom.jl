module TestMzRetention

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core./scan_timing/mz_retention.jl
# ─────────────────────────────────────────────────────────────────────────────

using Test
using JuChrom
using Unitful

@testset "mzretention branches" begin
    # invalid symbols
    @test_throws ArgumentError mzretention(1.0; mzindex=1, mzcount=1, retentionref=:foo)
    @test_throws ArgumentError mzretention(1.0; mzindex=1, mzcount=1, dwellref=:foo)
    @test_throws ArgumentError mzretention(1.0; mzindex=1, mzcount=1, dwell=:foo)
    @test_throws ArgumentError mzretention(1.0; mzindex=1, mzcount=1, order=:foo)

    # mzindex resolution and validation
    @test_throws ArgumentError mzretention(1.0; mzvalue=5, mzvalues=[1, 2, 3], mzcount=3)
    @test_throws ArgumentError mzretention(1.0; mzvalue=5)
    @test_throws ArgumentError mzretention(1.0; mzindex=0, mzcount=1)
    @test_throws ArgumentError mzretention(1.0; mzindex=1)
    @test_throws ArgumentError mzretention(1.0; mzindex=3, mzcount=2)

    # mzcount inference
    r_inf1 = mzretention(0.0; mzvalue=20, mzvalues=[10, 20, 30],
        dwell=:homogeneous, scaninterval=3.0, retentionref=:start, dwellref=:start)
    @test r_inf1 == 1.0

    r_inf2 = mzretention(0.0; mzindex=2, dwell=:heterogeneous,
        dwellretentions=[1.0, 1.0, 1.0], retentionref=:start, dwellref=:start)
    @test r_inf2 == 1.0

    # homogeneous dwell
    r = mzretention(10.0; mzindex=2, mzcount=4, dwell=:homogeneous,
        scaninterval=8.0, retentionref=:start, dwellref=:middle)
    @test r == 13.0

    r2 = mzretention(10.0; mzindex=1, mzcount=4, dwell=:homogeneous,
        dwellretention=2.0, scaninterval=8.0, retentionref=:start, dwellref=:start)
    @test r2 == 10.0

    @test_throws ArgumentError mzretention(10.0; mzindex=1, mzcount=4, dwell=:homogeneous,
        dwellretention=2.0, scaninterval=7.0)
    @test_throws ArgumentError mzretention(10.0; mzindex=1, mzcount=4, dwell=:homogeneous,
        scaninterval=-1.0)

    # heterogeneous dwell
    dw = [1.0, 2.0, 3.0]
    r3 = mzretention(0.0; mzindex=2, mzcount=3, dwell=:heterogeneous,
        dwellretentions=dw, retentionref=:start, dwellref=:end, order=:ascending)
    @test r3 == 3.0

    @test_throws ArgumentError mzretention(0.0; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwellretentions=dw, scaninterval=10.0)
    r4 = mzretention(0.0; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwellretentions=dw, scaninterval=10.0, validate_span=false)
    @test r4 == 0.5
    @test_throws ArgumentError mzretention(0.0; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwellretentions=dw, scaninterval=-1.0)
    @test_throws ArgumentError mzretention(0.0; mzindex=1, mzcount=4, dwell=:heterogeneous,
        dwellretentions=dw)
    @test_throws ArgumentError mzretention(0.0; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwellretentions=[1.0, 0.0, 1.0])

    # rtol/atol handling for span validation
    rtol_ok = mzretention(0.0; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwellretentions=dw, scaninterval=6.03, rtol=1e-2)
    @test rtol_ok == 0.5

    dw_u = [1.0, 1.0, 1.0]u"s"
    atol_ok = mzretention(0.0u"s"; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwellretentions=dw_u, scaninterval=3.0005u"s", atol=1e-3u"s")
    @test atol_ok ≈ 0.5u"s"

    # order mapping
    r5 = mzretention(0.0; mzindex=1, mzcount=4, dwell=:homogeneous, dwellretention=1.0,
        order=:descending, retentionref=:start, dwellref=:start)
    @test r5 == 3.0

    # retention_ref handling
    r6 = mzretention(10.0; mzindex=1, mzcount=2, dwell=:homogeneous, dwellretention=2.0,
        retentionref=:middle, dwellref=:start)
    @test r6 == 8.0
    r7 = mzretention(10.0; mzindex=1, mzcount=2, dwell=:homogeneous, dwellretention=2.0,
        retentionref=:end, dwellref=:start)
    @test r7 == 6.0

    # unit preservation
    r8 = mzretention(1.0u"minute"; mzindex=2, mzcount=2, dwell=:homogeneous,
        scaninterval=20u"s", retentionref=:start, dwellref=:end)
    @test r8 ≈ (1.0u"minute" + 20u"s")
    @test Unitful.unit(r8) == u"minute"
end

end # module TestMzRetention
