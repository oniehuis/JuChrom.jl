module TestTrim

using Test
using Unitful
using JuChrom

# ── Helpers ──────────────────────────────────────────────────────────────────

_make_cs(rt, int) = ChromScan(rt, int)
_make_css(scans::Vector{<:ChromScan}; meta...) = ChromScanSeries(scans; meta...)

_make_ms(mz::AbstractVector, ints::AbstractVector; rt=1.0u"s", lvl::Int=1) =
    MassScan(rt, mz, ints; level=lvl)
_make_mss(scans::Vector{<:MassScan}; meta...) = MassScanSeries(scans; meta...)

# ─────────────────────────────────────────────────────────────────────────────
# levelscans
# ─────────────────────────────────────────────────────────────────────────────

@testset "levelscans(series::AbstractMassScanSeries, target_level::Integer=1)" begin
    s1 = _make_ms([100.0], [10.0]; rt=1.0u"s",  lvl=1)
    s2 = _make_ms([150.0], [15.0]; rt=2.0u"s",  lvl=2)
    s3 = _make_ms([200.0], [30.0]; rt=3.0u"s",  lvl=1)
    s4 = _make_ms([250.0], [35.0]; rt=4.0u"s",  lvl=2)
    mss = _make_mss([s1, s2, s3, s4])

    ms1 = levelscans(mss, 1)
    @test scancount(ms1) == 2
    @test all(level.(scans(ms1)) .== 1)

    ms2 = levelscans(mss, 2)
    @test scancount(ms2) == 2
    @test all(level.(scans(ms2)) .== 2)

    @test_throws ArgumentError levelscans(mss, 3)
end

# ─────────────────────────────────────────────────────────────────────────────
# retentiontrim / retentiontrim!
# ─────────────────────────────────────────────────────────────────────────────

@testset "retentiontrim(series; start/stop with units)" begin
    css = _make_css([
        _make_cs(0.5u"s", 1.0),
        _make_cs(1.0u"s", 2.0),
        _make_cs(1.5u"s", 3.0),
        _make_cs(2.0u"s", 4.0),
        _make_cs(2.5u"s", 5.0),
    ])

    t = retentiontrim(css; start=1.0u"s", stop=2.0u"s")
    @test scancount(t) == 3
    @test retentions(t) == [1.0, 1.5, 2.0]u"s"
    @test JuChrom.intensities(t) == [2.0, 3.0, 4.0]
end

@testset "retentiontrim(series; start/stop unitless)" begin
    css = _make_css([
        _make_cs(0.5, 1.0),
        _make_cs(1.0, 2.0),
        _make_cs(1.5, 3.0),
        _make_cs(2.0, 4.0),
        _make_cs(2.5, 5.0),
    ])

    t = retentiontrim(css; start=1.0, stop=2.0)
    @test scancount(t) == 3
    @test retentions(t) == [1.0, 1.5, 2.0]
    @test JuChrom.intensities(t) == [2.0, 3.0, 4.0]
    @test scancount(css) == 5
end

@testset "retentiontrim!(series; start/stop with units)" begin
    css = _make_css([
        _make_cs(0.5u"s", 1.0),
        _make_cs(1.0u"s", 2.0),
        _make_cs(1.5u"s", 3.0),
        _make_cs(2.0u"s", 4.0),
        _make_cs(2.5u"s", 5.0),
    ])

    retentiontrim!(css; start=1.0u"s", stop=2.0u"s")
    @test scancount(css) == 3
    @test retentions(css) == [1.0, 1.5, 2.0]u"s"
end

#= @testset "retentiontrim(msmatrix; start/stop with units)" begin
    ret = [0.5, 1.0, 1.5, 2.0, 2.5]u"s"
    mz = [100.0, 200.0]
    ints = [1.0 2.0;
            3.0 4.0;
            5.0 6.0;
            7.0 8.0;
            9.0 10.0]
    msm = MassScanMatrix(
        ret,
        mz,
        ints;
        level=2,
        instrument=(model="X",),
        acquisition=(mode="scan",),
        user=(name="tester",),
        sample=(id="S1",),
        extras=Dict("note" => "ok"),
    )

    trimmed = retentiontrim(msm; start=1.0u"s", stop=2.0u"s")
    @test retentions(trimmed) == [1.0, 1.5, 2.0]u"s"
    @test rawintensities(trimmed) == ints[2:4, :]
    @test retentionunit(trimmed) == u"s"
    @test mzvalues(trimmed) == mz
    @test mzunit(trimmed) === nothing
    @test intensityunit(trimmed) === nothing
    @test level(trimmed) == 2
    @test extras(trimmed) == Dict("note" => "ok")
    @test retentions(msm) == ret
end =#

#= @testset "retentiontrim(msmatrix; start/stop unitless)" begin
    ret = [0.5, 1.0, 1.5, 2.0, 2.5]
    mz = [100.0, 200.0]
    ints = [1.0 2.0;
            3.0 4.0;
            5.0 6.0;
            7.0 8.0;
            9.0 10.0]
    msm = MassScanMatrix(
        ret,
        mz,
        ints;
        level=1,
        extras=Dict("note" => "ok"),
    )

    trimmed = retentiontrim(msm; start=1.0, stop=2.0)
    @test retentions(trimmed) == [1.0, 1.5, 2.0]
    @test rawintensities(trimmed) == ints[2:4, :]
    @test retentionunit(trimmed) === nothing
    @test mzvalues(trimmed) == mz
    @test mzunit(trimmed) === nothing
    @test intensityunit(trimmed) === nothing
    @test level(trimmed) == 1
    @test extras(trimmed) == Dict("note" => "ok")
    @test retentions(msm) == ret
    @test_throws ArgumentError retentiontrim(msm; start=2.0, stop=1.0)
    @test_throws ArgumentError retentiontrim(msm; start=10.0, stop=11.0)
end =#


# ─────────────────────────────────────────────────────────────────────────────
# indextrim / indextrim!
# ─────────────────────────────────────────────────────────────────────────────

@testset "indextrim(series; start/stop)" begin
    css = _make_css([
        _make_cs(1.0u"s", 1.0),
        _make_cs(2.0u"s", 2.0),
        _make_cs(3.0u"s", 3.0),
        _make_cs(4.0u"s", 4.0),
        _make_cs(5.0u"s", 5.0),
    ])

    t = indextrim(css; start=2, stop=4)
    @test scancount(t) == 3
    @test retentions(t) == [2.0, 3.0, 4.0]u"s"

    @test_throws ArgumentError indextrim(css; start=4, stop=3)
    @test_throws BoundsError  indextrim(css; start=0, stop=2)
    @test_throws BoundsError  indextrim(css; start=2, stop=10)
end

@testset "indextrim!(series; start/stop)" begin
    css = _make_css([
        _make_cs(1.0u"s", 1.0),
        _make_cs(2.0u"s", 2.0),
        _make_cs(3.0u"s", 3.0),
        _make_cs(4.0u"s", 4.0),
        _make_cs(5.0u"s", 5.0),
    ])

    indextrim!(css; start=2, stop=4)
    @test scancount(css) == 3
    @test retentions(css) == [2.0, 3.0, 4.0]u"s"

    css2 = _make_css([
        _make_cs(1.0u"s", 1.0),
        _make_cs(2.0u"s", 2.0),
        _make_cs(3.0u"s", 3.0),
    ])
    @test_throws ArgumentError indextrim!(css2; start=3, stop=2)
    @test_throws BoundsError  indextrim!(css2; start=0, stop=2)
    @test_throws BoundsError  indextrim!(css2; start=2, stop=10)
end

end  # module TestTrim
