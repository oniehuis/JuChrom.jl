module TestMzRetention

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core./scan_timing/mz_retention.jl
# ─────────────────────────────────────────────────────────────────────────────

using Test
using JuChrom
using Unitful

@testset "mzretention branches" begin
    # invalid symbols
    @test_throws ArgumentError mzretention(1.0; mzindex=1, mzcount=1, retention_ref=:foo)
    @test_throws ArgumentError mzretention(1.0; mzindex=1, mzcount=1, dwell_ref=:foo)
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
        dwell=:homogeneous, scan_interval=3.0, retention_ref=:start, dwell_ref=:start)
    @test r_inf1 == 1.0

    r_inf2 = mzretention(0.0; mzindex=2, dwell=:heterogeneous,
        dwell_retentions=[1.0, 1.0, 1.0], retention_ref=:start, dwell_ref=:start)
    @test r_inf2 == 1.0

    # homogeneous dwell
    r = mzretention(10.0; mzindex=2, mzcount=4, dwell=:homogeneous,
        scan_interval=8.0, retention_ref=:start, dwell_ref=:middle)
    @test r == 13.0

    r2 = mzretention(10.0; mzindex=1, mzcount=4, dwell=:homogeneous,
        dwell_retention=2.0, scan_interval=8.0, retention_ref=:start, dwell_ref=:start)
    @test r2 == 10.0

    @test_throws ArgumentError mzretention(10.0; mzindex=1, mzcount=4, dwell=:homogeneous,
        dwell_retention=2.0, scan_interval=7.0)
    @test_throws ArgumentError mzretention(10.0; mzindex=1, mzcount=4, dwell=:homogeneous,
        scan_interval=-1.0)

    # heterogeneous dwell
    dw = [1.0, 2.0, 3.0]
    r3 = mzretention(0.0; mzindex=2, mzcount=3, dwell=:heterogeneous,
        dwell_retentions=dw, retention_ref=:start, dwell_ref=:end, order=:ascending)
    @test r3 == 3.0

    @test_throws ArgumentError mzretention(0.0; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwell_retentions=dw, scan_interval=10.0)
    r4 = mzretention(0.0; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwell_retentions=dw, scan_interval=10.0, validate_span=false)
    @test r4 == 0.5
    @test_throws ArgumentError mzretention(0.0; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwell_retentions=dw, scan_interval=-1.0)
    @test_throws ArgumentError mzretention(0.0; mzindex=1, mzcount=4, dwell=:heterogeneous,
        dwell_retentions=dw)
    @test_throws ArgumentError mzretention(0.0; mzindex=1, mzcount=3, dwell=:heterogeneous,
        dwell_retentions=[1.0, 0.0, 1.0])

    # order mapping
    r5 = mzretention(0.0; mzindex=1, mzcount=4, dwell=:homogeneous, dwell_retention=1.0,
        order=:descending, retention_ref=:start, dwell_ref=:start)
    @test r5 == 3.0

    # retention_ref handling
    r6 = mzretention(10.0; mzindex=1, mzcount=2, dwell=:homogeneous, dwell_retention=2.0,
        retention_ref=:middle, dwell_ref=:start)
    @test r6 == 8.0
    r7 = mzretention(10.0; mzindex=1, mzcount=2, dwell=:homogeneous, dwell_retention=2.0,
        retention_ref=:end, dwell_ref=:start)
    @test r7 == 6.0

    # unit preservation
    r8 = mzretention(1.0u"minute"; mzindex=2, mzcount=2, dwell=:homogeneous,
        scan_interval=20u"s", retention_ref=:start, dwell_ref=:end)
    @test r8 ≈ (1.0u"minute" + 20u"s")
    @test Unitful.unit(r8) == u"minute"
end

end # module TestMzRetention
