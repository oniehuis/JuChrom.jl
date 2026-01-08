module TestMassHunterMSLoader

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for MassHunterMSLoader.jl
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

using SHA
using Test

# ── Internal Project Imports ─────────────────────────────────────────────────

using JuChrom
using JuChrom.InputOutput
using JuChrom.MassHunterMSLoader
import JuChrom.MassHunterMSLoader: MassHunterMSOptions, MassHunterMSLoaderSpec, 
    MassHunterMSv1, MSScanBinV1, MSPeakBinV1, magicnumber, mspeakdata, msscandata, 
    readfile, read_scan_data_ms, read_scan_data_tic

# ── Helper functions and constants for unit tests ────────────────────────────

sha1hex(v::AbstractVector) = bytes2hex(sha1(reinterpret(UInt8, v)))

const TESTPATH = joinpath(JuChrom.agilent, "C7-C40_MassHunterMS.D")
const TESTDATA = joinpath(TESTPATH, "AcqData")
const MSSCAN_FILE = joinpath(TESTDATA, "MSScan.bin")
const MSPEAK_FILE = joinpath(TESTDATA, "MSPeak.bin")
const SCANDATA = msscandata(MSScanBinV1(), MSSCAN_FILE)
const TEST_TMPDIR = mkpath(joinpath(@__DIR__, "tmp"))

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests
# ─────────────────────────────────────────────────────────────────────────────

# ── Loader Configuration Structures and Loader Spec Constructors ─────────────

@testset "MassHunterMSOptions" begin
    options = MassHunterMSOptions(:ms, nothing, Nothing, Nothing, Nothing)
    @test isa(options, MassHunterMSOptions)
    @test options.mode == :ms
    @test options.level === nothing
    @test options.scantimetype === Nothing
    @test options.iontype === Nothing
    @test options.intensitytype === Nothing

    custom = MassHunterMSOptions(:tic, 2, Float32, Float32, Float64)
    @test custom.mode == :tic
    @test custom.level == 2
    @test custom.scantimetype == Float32
    @test custom.iontype == Float32
    @test custom.intensitytype == Float64
end

@testset "MassHunterMSLoaderSpec" begin
    options = MassHunterMSOptions(:tic, 1, Float32, Float32, Float64)
    direct = MassHunterMSLoaderSpec{MassHunterMSv1}("dummy_path", options, 
        Val(MassHunterMSv1))
    @test isa(direct, MassHunterMSLoaderSpec)
    @test direct.path == "dummy_path"
    @test direct.options === options

    kw = MassHunterMSLoaderSpec{MassHunterMSv1}("kw_path"; mode=:ms, level=1, 
        scantimetype=Float32, iontype=Float32, intensitytype=Float64)
    @test isa(kw, MassHunterMSLoaderSpec)
    @test kw.path == "kw_path"
    @test kw.options.mode == :ms
    @test kw.options.level == 1
    @test kw.options.scantimetype == Float32
    @test kw.options.iontype == Float32
    @test kw.options.intensitytype == Float64

    @test_throws ArgumentError MassHunterMSLoaderSpec{MassHunterMSv1}("dummy"; mode=:bad)
end

@testset "MassHunterMS" begin
    spec = MassHunterMS(TESTPATH; mode=:tic, level=1)
    @test isa(spec, MassHunterMSLoaderSpec)
    @test spec.path == TESTPATH
    @test spec.options.mode == :tic
    @test spec.options.level == 1

    @test_throws ArgumentError MassHunterMS(TESTPATH; mode=:unknown)
end

# ── Data Loading Interface ───────────────────────────────────────────────────

@testset "load" begin
    spec = MassHunterMS(TESTPATH; mode=:ms)
    mss = load(spec)
    @test isa(mss, MassScanSeries)
    @test scancount(mss) == 2405

    tic_spec = MassHunterMS(TESTPATH; mode=:tic)
    css = load(tic_spec)
    @test isa(css, ChromScanSeries)
    @test scancount(css) == 2405

    file_spec = MassHunterMS(MSSCAN_FILE; mode=:ms)
    @test_throws ArgumentError load(file_spec)

    mktempdir(TEST_TMPDIR) do tmp
        spec = MassHunterMS(tmp; mode=:ms)
        @test_throws MissingFolderError load(spec)
    end

    mktempdir(TEST_TMPDIR) do tmp
        dfolder = joinpath(tmp, "Broken.D")
        acq = joinpath(dfolder, "AcqData")
        mkpath(acq)
        open(joinpath(acq, "MSScan.bin"), "w") do io
            write(io, zeros(UInt8, 32))
        end
        open(joinpath(acq, "MSPeak.bin"), "w") do io
            write(io, zeros(UInt8, 32))
        end
        spec = MassHunterMS(dfolder)
        @test_throws FileCorruptionError load(spec)
    end

    missing_level_spec = MassHunterMS(TESTPATH; mode=:ms, level=2)
    @test_throws FileCorruptionError load(missing_level_spec)
end

# ── File Reader Interface ───────────────────────────────────────────────────

@testset "readfile (:ms)" begin
    options = MassHunterMSOptions(:ms, nothing, Nothing, Nothing, Nothing)
    mss = readfile(MassHunterMSv1, TESTDATA, options)

    @test isa(mss, MassScanSeries)
    @test scancount(mss) == 2405
    @test unit(first(retentions(mss))) == u"s"

    st = ustrip.(retentions(mss))
    @test sha1hex(st) == "867faa5d1ac9f73bbae354a71eca501b5c4daa1a"

    first_scan = first(scans(mss))
    @test level(first_scan) == 1
    @test eltype(mzvalues(first_scan)) == Float64
    @test eltype(intensities(first_scan)) == Float64
    @test sha1hex(mzvalues(first_scan)) == "ff19749a4dab34a8b6af801096a99a9c7acb6dd0"
    @test sha1hex(intensities(first_scan)) == "57c82ff5c18ce9f2d17bfd4b9e6ceeb792d64dc5"
    @test getproperty(attrs(first_scan), :ScanID) == 1

    last_scan = last(scans(mss))
    @test sha1hex(mzvalues(last_scan)) == "052afb9e54521cd128312325a28fc5e08c11827f"
    @test sha1hex(intensities(last_scan)) == "a9804cb39164f94d380e131b4f56906584d6d8ec"

    @test instrument(mss) == NamedTuple()
    @test acquisition(mss) == NamedTuple()
    @test user(mss) == NamedTuple()
    @test sample(mss) == NamedTuple()
    @test extras(mss) == Dict()
end

@testset "readfile (:tic)" begin
    options = MassHunterMSOptions(:tic, nothing, Nothing, Nothing, Nothing)
    css = readfile(MassHunterMSv1, TESTDATA, options)

    @test isa(css, ChromScanSeries)
    @test scancount(css) == 2405
    @test unit(first(retentions(css))) == u"s"

    st = ustrip.(retentions(css))
    @test sha1hex(st) == "867faa5d1ac9f73bbae354a71eca501b5c4daa1a"
    @test sha1hex(intensities(css)) == "17da30bffc6a3bfb79ff764b2d5f8963f3c791be"

    @test instrument(css) == NamedTuple()
    @test acquisition(css) == NamedTuple()
    @test user(css) == NamedTuple()
    @test sample(css) == NamedTuple()
    @test extras(css) == Dict()
end

# ── Binary Helpers and Scan Data Extraction ──────────────────────────────────

@testset "magicnumber" begin
    @test magicnumber(MSSCAN_FILE) == 257
    @test magicnumber(MSPEAK_FILE) == 259
end

@testset "msscandata" begin
    @test length(SCANDATA) == 2405

    first_scan = first(SCANDATA)
    @test first_scan.ScanID == 1
    @test first_scan.PointCount == 650
    @test first_scan.ByteCount == 10400
    @test isapprox(first_scan.ScanTime, 3.1990166666666666; atol=eps(Float64))
    @test first_scan.MSLevel == 1

    last_scan = last(SCANDATA)
    @test last_scan.ScanID == 2405
    @test last_scan.PointCount == 617
    @test last_scan.ByteCount == 9872
    @test isapprox(last_scan.ScanTime, 31.650783333333333; atol=eps(Float64))
end

@testset "read_scan_data_tic" begin
    options = MassHunterMSOptions(:tic, nothing, Float32, Nothing, Float32)
    scantimes, tic_values, attrs_vec = read_scan_data_tic(MassHunterMSv1, SCANDATA, options)

    @test length(scantimes) == length(SCANDATA)
    @test length(tic_values) == length(SCANDATA)
    @test length(attrs_vec) == length(SCANDATA)
    @test eltype(scantimes) == Float32
    @test eltype(tic_values) == Float32
    @test scantimes[1] ≈ Float32(3.1990166666666666 * 60)
    @test tic_values[1] ≈ Float32(672293.0)
    @test getproperty(first(attrs_vec), :ScanID) == 1
end

@testset "read_scan_data_ms" begin
    options = MassHunterMSOptions(:ms, nothing, Float32, Float32, Float32)
    mss = read_scan_data_ms(MassHunterMSv1, TESTDATA, SCANDATA, options)
    @test isa(mss, MassScanSeries)
    @test scancount(mss) == 2405
    @test unit(first(retentions(mss))) == u"s"
    @test eltype(mzvalues(first(scans(mss)))) == Float32
    @test eltype(intensities(first(scans(mss)))) == Float32

    level_options = MassHunterMSOptions(:ms, 2, Nothing, Nothing, Nothing)
    @test_throws FileCorruptionError read_scan_data_ms(MassHunterMSv1, TESTDATA, SCANDATA, 
        level_options)
end

@testset "mspeakdata" begin
    options = MassHunterMSOptions(:ms, nothing, Nothing, Nothing, Nothing)
    scans_vec = mspeakdata(MSPeakBinV1(), MSPEAK_FILE, SCANDATA, options)
    @test length(scans_vec) == 2405
    first_scan = first(scans_vec)
    @test ustrip(retention(first_scan)) ≈ 3.1990166666666666 * 60
    @test length(mzvalues(first_scan)) == 650
end

end  # module TestMassHunterMSLoader
