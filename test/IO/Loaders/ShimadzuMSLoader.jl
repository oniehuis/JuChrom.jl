module TestShimadzuMSLoader

# ─────────────────────────────────────────────────────────────────────────────
# Comprehensive unit tests for ShimadzuMSLoader.jl
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

using SHA
using Test
using Unitful
using PyCall

# ── Internal Project Imports ─────────────────────────────────────────────────

using JuChrom
using JuChrom.InputOutput
# Ensure the PyCall extension is loaded so ShimadzuMSLoader is available.
let ext = Base.get_extension(JuChrom, :PyCallExtension)
    if ext === nothing && isdefined(Base, :retry_load_extensions)
        Base.retry_load_extensions()
        ext = Base.get_extension(JuChrom, :PyCallExtension)
    end
    if ext === nothing && !isdefined(JuChrom, :ShimadzuMSLoader)
        @eval JuChrom include(joinpath(dirname(pathof(JuChrom)), "IO", "Loaders", "ShimadzuMSLoader.jl"))
    end
    isdefined(JuChrom, :ShimadzuMSLoader) ||
        error("JuChrom ShimadzuMSLoader could not be loaded; ensure PyCall is available.")
end
using JuChrom.ShimadzuMSLoader
import JuChrom.ShimadzuMSLoader: ShimadzuMSOptions, readbytestring, extractdata, 
    build_mass_scans, readfile, load_mass_spectra, load_tic

# ── Constants ────────────────────────────────────────────────────────────────

const TESTFILE = joinpath(JuChrom.shimadzu, "AR190311.qgd")
const TEST_TMPDIR = mkpath(joinpath(@__DIR__, "tmp"))

const SHA_MS_STREAM = "0f7e746742d659c1223de0c4801b3ecb2aa033e5"
const SHA_RT_STREAM = "809654af6f90111ac91c441a769c223e488ac4cd"
const SHA_TIC_STREAM = "9dc2f897230be6995e7e38530687c6b88d91a4ac"
const SHA_RETENTIONS = "809654af6f90111ac91c441a769c223e488ac4cd"
const SHA_COUNTS = "9935af1ee6051d6f55d3c609c5c962d21113a774"
const SHA_MZS = "c5666f2868e1f520ef838dbad9a99393b3141742"
const SHA_INTS = "55a74cf9b3f90aff0de0c030c009e32f655d72ff"
const SHA_FIRST_MZ = "e045d62d82ccb03106c6a38dea2344551d81c4b8"
const SHA_FIRST_INT = "307239925a4c872757a90069eb9d647beaee619a"
const SHA_LAST_MZ = "90dcaea4d4b520017ba627d0342779fc9ed28897"
const SHA_LAST_INT = "fe3b3dc71a4d2607ffcb491dbbc1d008391dd62a"
const SHA_TIC_VALUES = "9dc2f897230be6995e7e38530687c6b88d91a4ac"

const SCAN_COUNT = 9_210
const TOTAL_POINTS = 5_212_860

# ── Helpers ──────────────────────────────────────────────────────────────────

sha1hex(data::AbstractVector{UInt8}) = bytes2hex(sha1(Vector{UInt8}(data)))
sha1hex(v::AbstractVector{<:Number}) = bytes2hex(sha1(Vector{UInt8}(reinterpret(UInt8, v))))

# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────

@testset "ShimadzuMSOptions" begin
    opts = ShimadzuMSOptions(:ms)
    @test opts.mode == :ms
end

@testset "ShimadzuMSLoaderSpec" begin
    spec1 = ShimadzuMSLoaderSpec{ShimadzuMSv1}(TESTFILE, ShimadzuMSOptions(:ms), 
        Val(ShimadzuMSv1))
    @test spec1.path == TESTFILE
    @test spec1.options.mode == :ms

    spec2 = ShimadzuMSLoaderSpec{ShimadzuMSv1}(TESTFILE; mode=:tic)
    @test spec2.options.mode == :tic

    @test_throws ArgumentError ShimadzuMSLoaderSpec{ShimadzuMSv1}(TESTFILE; mode=:foo)
end

@testset "ShimadzuMS constructor" begin
    spec = ShimadzuMS(TESTFILE; mode=:ms)
    @test isa(spec, ShimadzuMSLoaderSpec)
    @test spec.options.mode == :ms
    @test_throws ArgumentError ShimadzuMS(TESTFILE; mode=:unknown)
end

@testset "load" begin
    ms_series = load(ShimadzuMS(TESTFILE))
    @test isa(ms_series, MassScanSeries)
    @test length(scans(ms_series)) == SCAN_COUNT

    tic_series = load(ShimadzuMS(TESTFILE, mode=:tic))
    @test isa(tic_series, ChromScanSeries)
    @test length(scans(tic_series)) == SCAN_COUNT

    dir_spec = ShimadzuMS(dirname(TESTFILE))
    @test_throws MissingFileError load(dir_spec)

    missing = joinpath(dirname(TESTFILE), "missing_shimadzu.qgd")
    isfile(missing) && rm(missing)
    spec_missing = ShimadzuMS(missing)
    @test_throws MissingFileError load(spec_missing)

    mktemp(TEST_TMPDIR) do path, io
        write(io, zeros(UInt8, 64))
        close(io)
        spec = ShimadzuMS(path)
        @test_throws FileCorruptionError load(spec)
    end
end

@testset "readbytestring" begin
    ms_data = readbytestring(TESTFILE, ["GCMS Raw Data", "MS Raw Data"])
    @test length(ms_data) == 21_605_186
    @test sha1hex(ms_data) == SHA_MS_STREAM

    rt_data = readbytestring(TESTFILE, ["GCMS Raw Data", "Retention Time"])
    @test length(rt_data) == 36_840
    @test sha1hex(rt_data) == SHA_RT_STREAM

    tic_data = readbytestring(TESTFILE, ["GCMS Raw Data", "TIC Data"])
    @test length(tic_data) == 73_680
    @test sha1hex(tic_data) == SHA_TIC_STREAM

    @test_throws FileFormatError readbytestring(TESTFILE, ["GCMS Raw Data", "DoesNotExist"])
end

@testset "extractdata" begin
    bytes = readbytestring(TESTFILE, ["GCMS Raw Data", "MS Raw Data"])
    retentions_unitfree, counts, mzs, ints = extractdata(bytes)

    @test length(retentions_unitfree) == SCAN_COUNT
    @test length(counts) == SCAN_COUNT
    @test sum(counts) == TOTAL_POINTS
    @test length(mzs) == TOTAL_POINTS
    @test length(ints) == TOTAL_POINTS

    @test sha1hex(retentions_unitfree) == SHA_RETENTIONS
    @test sha1hex(counts) == SHA_COUNTS
    @test sha1hex(mzs) == SHA_MZS
    @test sha1hex(ints) == SHA_INTS

    badbytes = Vector{UInt8}(bytes)
    badbytes[21] = 0x05  # Force bytecount > 4
    @test_throws FileFormatError extractdata(badbytes)
end

@testset "build_mass_scans" begin
    bytes = readbytestring(TESTFILE, ["GCMS Raw Data", "MS Raw Data"])
    retentions_unitfree, counts, mzs, ints = extractdata(bytes)
    scans = build_mass_scans(retentions_unitfree, counts, mzs, ints)

    @test length(scans) == SCAN_COUNT
    @test retentionunit(first(scans)) == u"ms"
    @test intensityunit(first(scans)) === nothing
    @test sha1hex(mzvalues(first(scans))) == SHA_FIRST_MZ
    @test sha1hex(intensities(first(scans))) == SHA_FIRST_INT
    @test sha1hex(mzvalues(last(scans))) == SHA_LAST_MZ
    @test sha1hex(intensities(last(scans))) == SHA_LAST_INT

    @test_throws DimensionMismatch build_mass_scans(retentions_unitfree[1:10], 
        counts[1:9], mzs, ints)
    @test_throws DimensionMismatch build_mass_scans(retentions_unitfree[1:10], 
        counts[1:10], mzs[1:10], ints[1:9])
end

@testset "readfile" begin
    mss = readfile(TESTFILE, :ms)
    @test isa(mss, MassScanSeries)
    @test length(scans(mss)) == SCAN_COUNT
    @test sha1hex(ustrip.(retentions(mss))) == SHA_RETENTIONS
    @test sha1hex(mzvalues(first(scans(mss)))) == SHA_FIRST_MZ
    @test sha1hex(intensities(last(scans(mss)))) == SHA_LAST_INT
    @test mss.instrument == NamedTuple()
    @test mss.acquisition == NamedTuple()
    @test mss.user == NamedTuple()
    @test mss.sample == NamedTuple()

    css = readfile(TESTFILE, :tic)
    @test isa(css, ChromScanSeries)
    @test length(scans(css)) == SCAN_COUNT
    @test sha1hex(ustrip.(retentions(css))) == SHA_RETENTIONS
    @test sha1hex(rawintensity.(scans(css))) == SHA_TIC_VALUES
    @test css.instrument == NamedTuple()
    @test css.acquisition == NamedTuple()
    @test css.user == NamedTuple()
    @test css.sample == NamedTuple()

    @test_throws ArgumentError readfile(TESTFILE, :unknown)

    mktemp(TEST_TMPDIR) do path, io
        write(io, zeros(UInt8, 32))
        close(io)
        @test_throws FileFormatError readfile(path, :ms)
    end
end

@testset "load_mass_spectra / load_tic" begin
    mss = load_mass_spectra(TESTFILE)
    @test isa(mss, MassScanSeries)
    @test length(scans(mss)) == SCAN_COUNT

    css = load_tic(TESTFILE)
    @test isa(css, ChromScanSeries)
    @test length(scans(css)) == SCAN_COUNT

    mktemp(TEST_TMPDIR) do path, io
        write(io, zeros(UInt8, 16))
        close(io)
        @test_throws FileCorruptionError load_mass_spectra(path)
        @test_throws FileCorruptionError load_tic(path)
    end
end

end  # module TestShimadzuMSLoader
