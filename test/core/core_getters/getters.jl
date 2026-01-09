module TestGetters

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core./getters.jl
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

using Test
using Unitful
using Unitful: s, pA, nA

# ── Internal Project Imports ─────────────────────────────────────────────────

using JuChrom

# ── Constants ────────────────────────────────────────────────────────────────

const TIME = 1.0s
const TIME2 = 2.0s
const TIMES = [1.0s, 2.0s, 3.0s]

const INT = 100.0
const INT2 = 200.0

const MZS = [200.0, 201.0, 202.0]
const MZS2 = [100.0, 101.0, 102.0]

const INTS = [1.0, 2.0, 3.0]
const INTS2 = [4.0, 5.0, 6.0]
const INTS3 = [7.0, 8.0, 9.0]
const IM = [1.0 2.0 3.0;
            4.0 5.0 6.0;
            7.0 8.0 9.0]

const METADATA_TUPLE = (column="1", flowrate=0.1)
const METADATA_TUPLE2 = (column="2", flowrate=0.2)

const CHROMSCAN = ChromScan(TIME, INT; attrs=METADATA_TUPLE)
const CHROMSCAN2 = ChromScan(TIME2, INT2; attrs=METADATA_TUPLE2)
const CHROMSCAN_UNITFUL_INT = ChromScan(1.5s, 50.0u"pA")
const CHROMSCAN_UNITLESS = ChromScan(3.0, 20.0)

const LEVEL = 1
const LEVEL2 = 2
const LEVELS = [LEVEL, LEVEL2]

const MASSSCAN = MassScan(TIME, MZS, INTS; level=LEVEL, attrs=METADATA_TUPLE)
const MASSSCAN2 = MassScan(TIME2, MZS2, INTS2; level=LEVEL2, attrs=METADATA_TUPLE2)

const INSTRUMENT = (model="TestModel", vendor="TestVendor")
const ACQUISITION = (method="TestMethod", date="2025-05-29")
const USER = (name="Mr./Ms. Test", id=42)
const SAMPLE = (id="1", description="Test sample")
const METADATA_DICT = Dict("extra data 1" => "1", "extra data 2" => 2)

const CHROMSCANS = [CHROMSCAN, CHROMSCAN2]
const MASSSCANS = [MASSSCAN, MASSSCAN2]

# UPDATED: use concrete series constructors
const CSS = ChromScanSeries(
    CHROMSCANS;
    instrument=INSTRUMENT,
    acquisition=ACQUISITION,
    user=USER,
    sample=SAMPLE,
    extras=METADATA_DICT)

const MSS = MassScanSeries(
    MASSSCANS;
    instrument=INSTRUMENT,
    acquisition=ACQUISITION,
    user=USER,
    sample=SAMPLE,
    extras=METADATA_DICT)

# UPDATED: use concrete matrix constructor
const SM = MassScanMatrix(TIMES, MZS, IM;
    level=LEVEL,
    instrument=INSTRUMENT,
    acquisition=ACQUISITION,
    user=USER,
    sample=SAMPLE,
    extras=METADATA_DICT)
const SM_UNITLESS = MassScanMatrix(
    [1.0, 2.0], MZS, IM[1:2, :];
    level=LEVEL,
    instrument=INSTRUMENT,
    acquisition=ACQUISITION,
    user=USER,
    sample=SAMPLE,
    extras=METADATA_DICT)

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests
# ─────────────────────────────────────────────────────────────────────────────

# ── acquisition ──────────────────────────────────────────────────────────────

@testset "acquisition(series::AbstractScanSeries)" begin
    # Test that `acquisition` works for both concrete types
    @test acquisition(CSS) == ACQUISITION
    @test acquisition(MSS) == ACQUISITION

    # Error on invalid `acquisition` type for getter
    @test_throws MethodError acquisition("invalid type")
    @test_throws MethodError acquisition([1, 2, 3])
    @test_throws MethodError acquisition(1.0)
    @test_throws MethodError acquisition(nothing)
end

@testset "acquisition(scanmatrix::AbstractScanMatrix)" begin
    # Test getter: acquisition
    @test acquisition(SM) == ACQUISITION

    # Error on invalid `acquisition` type for getter
    @test_throws MethodError acquisition("invalid type")
    @test_throws MethodError acquisition([1, 2, 3])
    @test_throws MethodError acquisition(1.0)
    @test_throws MethodError acquisition(nothing)
end

# ── attrs ────────────────────────────────────────────────────────────────────

@testset "attrs(scan::AbstractScan)" begin
    # Test that `attrs` works for both concrete types
    @test attrs(CHROMSCAN) == METADATA_TUPLE
    @test attrs(MASSSCAN) == METADATA_TUPLE

    # Error on invalid `attrs` type for getter
    @test_throws MethodError attrs("invalid type")
    @test_throws MethodError attrs([1, 2, 3])
    @test_throws MethodError attrs(1.0)
    @test_throws MethodError attrs(nothing)
end

# ── extras ───────────────────────────────────────────────────────────────────

@testset "extras(series::AbstractScanSeries)" begin
    # Test that `extras` works for both concrete types
    @test extras(CSS) == METADATA_DICT
    @test extras(MSS) == METADATA_DICT

    # Error on invalid `extras` type for getter
    @test_throws MethodError extras("invalid type")
    @test_throws MethodError extras([1, 2, 3])
    @test_throws MethodError extras(1.0)
    @test_throws MethodError extras(nothing)
end

@testset "extras(scanmatrix::AbstractScanMatrix)" begin
    # Test that `extras` works for scanmatrix
    @test extras(SM) == METADATA_DICT

    # Error on invalid `extras` type for getter
    @test_throws MethodError extras("invalid type")
    @test_throws MethodError extras([1, 2, 3])
    @test_throws MethodError extras(1.0)
    @test_throws MethodError extras(nothing)
end

# ── intensity ────────────────────────────────────────────────────────────────

@testset "intensity(scan::AbstractChromScan)" begin
    # Test that `intensity` getter works
    @test intensity(CHROMSCAN) == INT

    # Error on invalid `intensity` type for getter
    @test_throws MethodError intensity("invalid type")
    @test_throws MethodError intensity([1, 2, 3])
    @test_throws MethodError intensity(1.0)
    @test_throws MethodError intensity(nothing)
end

# ── rawintensity ─────────────────────────────────────────────────────────────

@testset "rawintensity(scan::AbstractChromScan)" begin
    # Unitless path returns raw numeric and errors if a unit is requested
    @test rawintensity(CHROMSCAN) == INT
    @test_throws ArgumentError rawintensity(CHROMSCAN; unit=u"pA")

    # Unitful path strips units and handles conversion
    @test rawintensity(CHROMSCAN_UNITFUL_INT) == 50.0
    @test rawintensity(CHROMSCAN_UNITFUL_INT; unit=u"fA") == 50_000.0

    # Error on invalid `rawintensity` type for getter
    @test_throws MethodError rawintensity("invalid type")
    @test_throws MethodError rawintensity([1, 2, 3])
    @test_throws MethodError rawintensity(1.0)
    @test_throws MethodError rawintensity(nothing)
end

# ── intensities ──────────────────────────────────────────────────────────────

@testset "intensities(scan::AbstractMassScan)" begin
    # Test that `intensities` getter works
    @test intensities(MASSSCAN) == INTS

    # Error on invalid `intensities` type for getter
    @test_throws MethodError intensities("invalid type")
    @test_throws MethodError intensities([1, 2, 3])
    @test_throws MethodError intensities(1.0)
    @test_throws MethodError intensities(nothing)
end

# ── rawintensities ───────────────────────────────────────────────────────────

@testset "rawintensities(scan::AbstractMassScan)" begin
    m_unitful = MassScan(1.0u"s", [100.0, 200.0]u"Th", [10.0, 20.0]u"pA")
    @test rawintensities(m_unitful) == [10.0, 20.0]
    @test rawintensities(m_unitful; unit=u"nA") == [0.01, 0.02]
    @test rawintensities(MASSSCAN) == INTS
end

@testset "intensities(series::AbstractScanSeries, scanindex::Integer)" begin
    # Test that `intensities` works
    @test intensities(MSS, 1) == INTS
    @test intensities(MSS, 2) == INTS2

    # Error on invalid `intensities` value
    @test_throws BoundsError intensities(MSS, 0)  # Invalid index
    @test_throws BoundsError intensities(MSS, 3)  # Out of bounds index

    # Error on invalid `intensities` type
    @test_throws MethodError intensities("invalid type", 1)
    @test_throws MethodError intensities([1, 2, 3], 1)
    @test_throws MethodError intensities(1.0, 1)
    @test_throws MethodError intensities(nothing, 1)
    @test_throws MethodError intensities(MSS, "invalid type")
    @test_throws MethodError intensities(MSS, [1, 2, 3])
    @test_throws MethodError intensities(MSS, 1.0)
    @test_throws MethodError intensities(MSS, nothing)
end

@testset "intensities(series::AbstractChromScanSeries)" begin
    # Test that `intensities` works
    @test intensities(CSS) == [INT, INT2]

    # Error on invalid `intensities` type for getter
    @test_throws MethodError intensities("invalid type")
    @test_throws MethodError intensities([1, 2, 3])
    @test_throws MethodError intensities(1.0)
    @test_throws MethodError intensities(nothing)
end

@testset "intensities(scanmatrix::AbstractScanMatrix)" begin
    # Test getter: intensities
    @test intensities(SM) == IM

    # Error on invalid `intensities` type for getter
    @test_throws MethodError intensities("invalid type")
    @test_throws MethodError intensities([1, 2, 3])
    @test_throws MethodError intensities(1.0)
    @test_throws MethodError intensities(nothing)
end

# ── ioncount → mzcount ───────────────────────────────────────────────────────

@testset "mzcount(scan::AbstractMassScan)" begin
    # Test that `mzcount` getter works
    @test mzcount(MASSSCAN) == length(MZS)

    # Error on invalid `mzcount` type for getter
    @test_throws MethodError mzcount("invalid type")
    @test_throws MethodError mzcount([1, 2, 3])
    @test_throws MethodError mzcount(1.0)
    @test_throws MethodError mzcount(nothing)
end

# ── ions → mzvalues ─────────────────────────────────────────────────────────

@testset "mzvalues(scan::AbstractMassScan)" begin
    # Test that `mzvalues` getter works
    @test mzvalues(MASSSCAN) == MZS

    # Error on invalid `mzvalues` type for getter
    @test_throws MethodError mzvalues("invalid type")
    @test_throws MethodError mzvalues([1, 2, 3])
    @test_throws MethodError mzvalues(1.0)
    @test_throws MethodError mzvalues(nothing)
end

# ── rawmzvalues ──────────────────────────────────────────────────────────────

@testset "rawmzvalues(scan::AbstractMassScan)" begin
    m_unitful = MassScan(1.0u"s", [100.0, 200.0]u"Th", [10.0, 20.0]u"pA")
    @test rawmzvalues(m_unitful) == [100.0, 200.0]
    @test rawmzvalues(m_unitful; unit=u"Th") == [100.0, 200.0]
    @test rawmzvalues(MASSSCAN) == MZS
end

@testset "mzvalues(series::AbstractScanSeries, scanindex::Integer)" begin
    # Test that `mzvalues` works
    @test mzvalues(MSS, 1) == MZS
    @test mzvalues(MSS, 2) == MZS2

    # Error on invalid `mzvalues` value
    @test_throws BoundsError mzvalues(MSS, 0)  # Invalid index
    @test_throws BoundsError mzvalues(MSS, 3)  # Out of bounds index

    # Error on invalid `mzvalues` type
    @test_throws MethodError mzvalues("invalid type", 1)
    @test_throws MethodError mzvalues([1, 2, 3], 1)
    @test_throws MethodError mzvalues(1.0, 1)
    @test_throws MethodError mzvalues(nothing, 1)
    @test_throws MethodError mzvalues(MSS, "invalid type")
    @test_throws MethodError mzvalues(MSS, [1, 2, 3])
    @test_throws MethodError mzvalues(MSS, 1.0)
    @test_throws MethodError mzvalues(MSS, nothing)
end

@testset "mzvalues(scanmatrix::AbstractScanMatrix)" begin
    # Test getter: mzvalues
    @test mzvalues(SM) == MZS

    # Error on invalid `mzvalues` type for getter
    @test_throws MethodError mzvalues("invalid type")
    @test_throws MethodError mzvalues([1, 2, 3])
    @test_throws MethodError mzvalues(1.0)
    @test_throws MethodError mzvalues(nothing)
end

# ── instrument ───────────────────────────────────────────────────────────────

@testset "instrument(series::AbstractScanSeries)" begin
    # Test getter: instrument
    @test instrument(CSS) == INSTRUMENT

    # Error on invalid `instrument` type for getter
    @test_throws MethodError instrument("invalid type")
    @test_throws MethodError instrument([1, 2, 3])
    @test_throws MethodError instrument(1.0)
    @test_throws MethodError instrument(nothing)
end

@testset "instrument(scanmatrix::AbstractScanMatrix)" begin
    # Test getter: instrument
    @test instrument(SM) == INSTRUMENT

    # Error on invalid `instrument` type for getter
    @test_throws MethodError instrument("invalid type")
    @test_throws MethodError instrument([1, 2, 3])
    @test_throws MethodError instrument(1.0)
    @test_throws MethodError instrument(nothing)
end

# ── level ────────────────────────────────────────────────────────────────────

@testset "level(scan::AbstractMassScan)" begin
    # Test getter: level
    @test level(MASSSCAN) == 1

    # Error on invalid `level` type for getter
    @test_throws MethodError level("invalid type")
    @test_throws MethodError level([1, 2, 3])
    @test_throws MethodError level(1.0)
    @test_throws MethodError level(nothing)
end

@testset "level(scanmatrix::AbstractScanMatrix)" begin
    # Test that `level` works for both concrete types
    @test level(SM) == 1

    # Error on invalid `level` type for getter
    @test_throws MethodError level("invalid type")
    @test_throws MethodError level([1, 2, 3])
    @test_throws MethodError level(1.0)
    @test_throws MethodError level(nothing)
end

# ── levels ───────────────────────────────────────────────────────────────────

@testset "levels(series::AbstractMassScanSeries)" begin
    # Test getter: instrument
    @test levels(MSS) == LEVELS

    # Error on invalid `instrument` type for getter
    @test_throws MethodError levels("invalid type")
    @test_throws MethodError levels([1, 2, 3])
    @test_throws MethodError levels(1.0)
    @test_throws MethodError levels(nothing)
end

# ── rawretention ──────────────────────────────────────────────────────────────────────────

@testset "rawretention(scan::AbstractScan)" begin
    # Unitful scan strips units and supports conversion
    @test rawretention(CHROMSCAN) == 1.0
    @test rawretention(MASSSCAN; unit=u"ms") == 1_000.0

    # Unitless scan returns raw numeric and errors if a unit is requested
    @test rawretention(CHROMSCAN_UNITLESS) == 3.0
    @test_throws ArgumentError rawretention(CHROMSCAN_UNITLESS; unit=u"s")

    # Error on invalid `rawretention` type for getter
    @test_throws MethodError rawretention("invalid type")
    @test_throws MethodError rawretention([1, 2, 3])
    @test_throws MethodError rawretention(1.0)
    @test_throws MethodError rawretention(nothing)
end

# ── rawretentions ─────────────────────────────────────────────────────────────────────────

@testset "rawretentions(scanmatrix::AbstractMassScanMatrix)" begin
    # Unitful matrix strips units and supports conversion
    @test rawretentions(SM) == [1.0, 2.0, 3.0]
    @test rawretentions(SM; unit=u"ms") == [1_000.0, 2_000.0, 3_000.0]

    # Unitless matrix returns raw numeric and errors if a unit is requested
    @test rawretentions(SM_UNITLESS) == [1.0, 2.0]
    @test_throws ArgumentError rawretentions(SM_UNITLESS; unit=u"s")

    # Error on invalid `rawretentions` type for getter
    @test_throws MethodError rawretentions("invalid type")
    @test_throws MethodError rawretentions([1, 2, 3])
    @test_throws MethodError rawretentions(1.0)
    @test_throws MethodError rawretentions(nothing)
end

@testset "rawretentions(series::AbstractScanSeries)" begin
    # Unitful series strips units and supports conversion
    @test rawretentions(CSS) == [1.0, 2.0]
    @test rawretentions(MSS; unit=u"ms") == [1_000.0, 2_000.0]

    # Unitless series returns raw numeric and errors if a unit is requested
    unitless_series = ChromScanSeries([CHROMSCAN_UNITLESS, ChromScan(4.0, 30.0)])
    @test rawretentions(unitless_series) == [3.0, 4.0]
    @test_throws ArgumentError rawretentions(unitless_series; unit=u"s")

    # Error on invalid `rawretentions` type for getter
    @test_throws MethodError rawretentions("invalid type")
    @test_throws MethodError rawretentions([1, 2, 3])
    @test_throws MethodError rawretentions(1.0)
    @test_throws MethodError rawretentions(nothing)
end

# ── retention ─────────────────────────────────────────────────────────────────────────────

@testset "retention(scan::AbstractScan)" begin
    # Test that `retention` works for both concrete types
    @test retention(CHROMSCAN) == TIME
    @test retention(MASSSCAN) == TIME

    # Error on invalid `retention` type for getter
    @test_throws MethodError retention("invalid type")
    @test_throws MethodError retention([1, 2, 3])
    @test_throws MethodError retention(1.0)
    @test_throws MethodError retention(nothing)
end

# ── retentions ────────────────────────────────────────────────────────────────────────────

@testset "retentions(scanmatrix::AbstractScanMatrix)" begin
    # Test that `retentions` works for scanmatrix
    @test retentions(SM) == TIMES

    # Error on invalid `retentions` type for getter
    @test_throws MethodError retentions("invalid type")
    @test_throws MethodError retentions([1, 2, 3])
    @test_throws MethodError retentions(1.0)
    @test_throws MethodError retentions(nothing)
end

@testset "retentions(series::AbstractScanSeries)" begin
    # Test that `retentions` works
    @test retentions(CSS) == [TIME, TIME2]
    @test retentions(MSS) == [TIME, TIME2]

    # Error on invalid `retentions` type for getter
    @test_throws MethodError retentions("invalid type")
    @test_throws MethodError retentions([1, 2, 3])
    @test_throws MethodError retentions(1.0)
    @test_throws MethodError retentions(nothing)
end

# ── sample ───────────────────────────────────────────────────────────────────

@testset "sample(series::AbstractScanSeries)" begin
    # Test that `sample` works for both concrete types
    @test sample(CSS) == SAMPLE
    @test sample(MSS) == SAMPLE

    # Error on invalid `sample` type for getter
    @test_throws MethodError sample("invalid type")
    @test_throws MethodError sample([1, 2, 3])
    @test_throws MethodError sample(1.0)
    @test_throws MethodError sample(nothing)
end

@testset "sample(scanmatrix::AbstractScanMatrix)" begin
    # Test that `sample` works for scanmatrix
    @test sample(SM) == SAMPLE

    # Error on invalid `sample` type for getter
    @test_throws MethodError sample("invalid type")
    @test_throws MethodError sample([1, 2, 3])
    @test_throws MethodError sample(1.0)
    @test_throws MethodError sample(nothing)
end

# ── scan ─────────────────────────────────────────────────────────────────────

@testset "scan(series::AbstractScanSeries, scanindex::Integer)" begin
    # Test that `scan` works for both concrete types
    @test scan(CSS, 1) == CHROMSCAN
    @test scan(MSS, 1) == MASSSCAN

    # Error on invalid `scan` value
    @test_throws BoundsError scan(CSS, 0)  # Invalid index
    @test_throws BoundsError scan(CSS, 3)  # Out of bounds index

    # Error on invalid `scan` type
    @test_throws MethodError scan("invalid type", 1)
    @test_throws MethodError scan([1, 2, 3], 1)
    @test_throws MethodError scan(1.0, 1)
    @test_throws MethodError scan(nothing, 1)
    @test_throws MethodError scan(CSS, "invalid type")
    @test_throws MethodError scan(CSS, [1, 2, 3])
    @test_throws MethodError scan(CSS, 1.0)
    @test_throws MethodError scan(CSS, nothing)
end

# ── scancount ────────────────────────────────────────────────────────────────

@testset "scancount(series::AbstractScanSeries)" begin
    # Test that `scancount` works for both concrete types
    @test scancount(CSS) == length(CHROMSCANS)
    @test scancount(MSS) == length(MASSSCANS)

    # Error on invalid `scancount` type for getter
    @test_throws MethodError scancount("invalid type")
    @test_throws MethodError scancount([1, 2, 3])
    @test_throws MethodError scancount(1.0)
    @test_throws MethodError scancount(nothing)
end

# ── scans ────────────────────────────────────────────────────────────────────

@testset "scans(series::AbstractScanSeries)" begin
    # Test that `scans` works for both concrete types
    @test scans(CSS) == CHROMSCANS
    @test scans(MSS) == MASSSCANS

    # Error on invalid `scans` type for getter
    @test_throws MethodError scans("invalid type")
    @test_throws MethodError scans([1, 2, 3])
    @test_throws MethodError scans(1.0)
    @test_throws MethodError scans(nothing)
end

# ── uniquemzvalues ────────────────────────────────────────────────────────

@testset "uniquemzvalues(series::AbstractMassScanSeries, target_level::Integer=1)" begin
    @test uniquemzvalues(MSS, 1) == sort(unique(MZS))
    @test uniquemzvalues(MSS, 2) == sort(unique(MZS2))
    @test_throws ArgumentError uniquemzvalues(MSS, 3)  # no scans at this level

    # Type errors (invalid series or level)
    @test_throws MethodError uniquemzvalues("invalid type", 1)
    @test_throws MethodError uniquemzvalues([1, 2, 3], 1)
    @test_throws MethodError uniquemzvalues(1.0, 1)
    @test_throws MethodError uniquemzvalues(nothing, 1)
    @test_throws MethodError uniquemzvalues(MSS, "invalid type")
    @test_throws MethodError uniquemzvalues(MSS, 1.01)
    @test_throws MethodError uniquemzvalues(MSS, nothing)
end

# ── user ─────────────────────────────────────────────────────────────────────

@testset "user(series::AbstractScanSeries)" begin
    # Test that `user` works for both concrete types
    @test user(CSS) == USER
    @test user(MSS) == USER

    # Error on invalid `user` type for getter
    @test_throws MethodError user("invalid type")
    @test_throws MethodError user([1, 2, 3])
    @test_throws MethodError user(1.0)
    @test_throws MethodError user(nothing)
end

@testset "user(scanmatrix::AbstractScanMatrix)" begin
    # Test that `user` works for scanmatrix
    @test user(SM) == USER

    # Error on invalid `user` type for getter
    @test_throws MethodError user("invalid type")
    @test_throws MethodError user([1, 2, 3])
    @test_throws MethodError user(1.0)
    @test_throws MethodError user(nothing)
end

end
