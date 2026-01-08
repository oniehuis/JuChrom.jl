module TestChemStationMSLoader

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ChemStationMSLoader.jl
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

using Pkg.Artifacts
using SHA
using StringEncodings
using Test

# ── Internal Project Imports ─────────────────────────────────────────────────

using JuChrom
using JuChrom.InputOutput
using JuChrom.ChemStationMSLoader
import JuChrom.ChemStationMSLoader: build_mass_scans, ChemStationMSOptions, 
    ensure_bytes_available, read_bytes_ltoh, read_len_prefixed_string_ltoh_at, 
    read_len_prefixed_string_ltoh, read_metadata, read_scalar_ltoh, read_scalar_ntoh, 
    read_scalar_ntoh_at, read_scan_data_ms, read_scan_data_tic, readfile, validate_header

# ── Helper function to for unit tests ────────────────────────────────────────

function extract_test_data(file::String)
    open(file, "r") do f
        validate_header(f, file)
        metadata = read_metadata(f)

        scancount = read_scalar_ntoh_at(f, 278, UInt32)
        
        scan_offset = 2 * read_scalar_ntoh_at(f, 266, UInt16) - 2
        
        seek(f, scan_offset)
        length_units = read_scalar_ntoh(f, UInt16)
        chunk_size_bytes = 2 * length_units
        
        seek(f, scan_offset)
        scan_chunk1 = Vector{UInt8}(undef, chunk_size_bytes)
        read!(f, scan_chunk1)
        scan_chunk2 = Vector{UInt8}(undef, chunk_size_bytes)
        read!(f, scan_chunk2)
        
        return (
            metadata = metadata,
            scancount = scancount,
            scan_offset = scan_offset,
            first_scan_chunk = scan_chunk1,
            second_scan_chunk = scan_chunk2,
        )
    end
end

sha1hex(v::AbstractVector) = bytes2hex(sha1(reinterpret(UInt8, v)))

function make_tic_chunk(scan_time::Int32, tic_value::Int32; payload_bytes::Int=8)
    payload_bytes ≥ 0 || throw(ArgumentError("payload_bytes must be non-negative"))
    chunk_size_bytes = payload_bytes + 8
    length_units = UInt16((chunk_size_bytes + 2) ÷ 2)
    io = IOBuffer()
    write(io, Base.hton(length_units))
    write(io, Base.hton(scan_time))
    if payload_bytes > 0
        write(io, zeros(UInt8, payload_bytes))
    end
    write(io, Base.hton(tic_value))
    take!(io)
end

function make_header_buf(version::String, type::String)
    return IOBuffer(UInt8[
        length(version),
        encode(version, ENCODING)...,
        0x00, 0x00,
        length(type),
        encode(type, ENCODING)...,
    ])
end

# ── Constants ────────────────────────────────────────────────────────────────

const TESTFILE = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
const ENCODING = "Windows-1252"
const testdata = extract_test_data(TESTFILE)
const TEST_TMPDIR = mkpath(joinpath(@__DIR__, "tmp"))

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests
# ─────────────────────────────────────────────────────────────────────────────

# ── Loader Configuration Structures and Loader Spec Constructors ─────────────

@testset "ChemStationMSOptions" begin
    options = ChemStationMSOptions(:ms)
    @test isa(options, ChemStationMSOptions)
    @test options.mode == :ms
end

@testset "ChemStationMSLoaderSpec" begin
    # Direct constructor
    spec1 = ChemStationMSLoaderSpec{ChemStationMSv2}("dummy_path",
        ChemStationMSOptions(:ms), Val(ChemStationMSv2))
    @test isa(spec1, ChemStationMSLoaderSpec)
    @test spec1.path == "dummy_path"
    @test isa(spec1.options, ChemStationMSOptions)
    @test spec1.options.mode == :ms

    # Keyword constructor
    spec2 = ChemStationMSLoaderSpec{ChemStationMSv2}("dummy_path"; mode=:ms)
    @test isa(spec2, ChemStationMSLoaderSpec)
    @test spec2.path == "dummy_path"
    @test isa(spec2.options, ChemStationMSOptions)
    @test spec2.options.mode == :ms
end

@testset "ChemStationMS" begin
    spec = ChemStationMS("dummy_path"; mode=:ms)
    @test isa(spec, ChemStationMSLoaderSpec)
    @test spec.path == "dummy_path"
    @test isa(spec.options, ChemStationMSOptions)
    @test spec.options.mode == :ms
end

# ── Data Loading Interface ───────────────────────────────────────────────────

@testset "load" begin
    spec = ChemStationMS(TESTFILE; mode=:ms)
    mss = load(spec)
    @test isa(mss, MassScanSeries)

    testdir = dirname(TESTFILE)
    spec = ChemStationMS(testdir; mode=:ms)
    mss = load(spec)
    @test isa(mss, MassScanSeries)

    falsefile = joinpath(testdir, "nonexistent.ms")
    spec = ChemStationMS(falsefile; mode=:ms)
    @test_throws MissingFileError load(spec)

    # Corrupted file
    mktempdir(TEST_TMPDIR) do tmp
        dst = joinpath(tmp, "data.ms")
        open(dst, "w") do io
            write(io, zeros(UInt8, 32))
        end

        req = ChemStationMSLoaderSpec{ChemStationMSv2}(dst)
        @test_throws FileCorruptionError load(req)
    end
end

# ─── File Reader Interface ───────────────────────────────────────────────────

@testset "readfile" begin
    # Extact MS data and check its contents
    mss = readfile(TESTFILE, :ms)

    @test isa(mss, MassScanSeries)
    @test eltype(scans(mss)) <: MassScan

    @test unit(first(retentions(mss))) == u"ms"
    st = ustrip.(retentions(mss))
    @test sha1hex(st) == "7b40412299cd7951433d01d240319298214a7e47"

    first_scan = first(scans(mss))
    @test sha1hex(mzvalues(first_scan)) == "c20378b8a7e4ca28900e6c0c25016e76a32bb19e"
    @test sha1hex(intensities(first_scan)) == "313c1f50ef0ba7882d10154cd4ff6e0190474b60"
    @test level(first_scan) == 1
    @test attrs(first_scan) == NamedTuple()

    last_scan = last(scans(mss))
    @test sha1hex(mzvalues(last_scan)) == "568b4386013fbe8f20440c2d4fd916b1a572fc0f"
    @test sha1hex(intensities(last_scan)) == "14d4d25bcaa7e5cbd89cf645f67af72af295d4e5"
    @test level(first_scan) == 1
    @test attrs(first_scan) == NamedTuple()

    @test hasproperty(mss, :instrument)
    @test hasproperty(mss, :acquisition)
    @test hasproperty(mss, :user)
    @test hasproperty(mss, :sample)
    @test extras(mss) == Dict()

    # Extact TIC data and check its contents
    css = readfile(TESTFILE, :tic)

    @test isa(css, ChromScanSeries)
    @test eltype(scans(css)) <: ChromScan

    @test unit(first(retentions(css))) == u"ms"
    st = ustrip.(retentions(css))
    @test sha1hex(st) == "7b40412299cd7951433d01d240319298214a7e47"

    @test sha1hex(intensities(css)
        ) == "bd194aa7c0d3d284767f7cef67ac1ece7880d51f"

    @test hasproperty(css, :instrument)
    @test hasproperty(css, :acquisition)
    @test hasproperty(css, :user)
    @test hasproperty(css, :sample)
    @test extras(css) == Dict()
end

# ─── File Header Validation and Metadata Extraction ──────────────────────────

@testset "validate_header" begin
    version = "2"
    type1 = "GC / MS Data File"
    type2 = "GC / MS DATA FILE"
    wrong_version = "3"
    wrong_type = "Unknown Type"

    @test validate_header(make_header_buf(version, type1), 
        "file")
    @test validate_header(make_header_buf(version, type2), 
        "file")

    @test_throws FileFormatError validate_header(make_header_buf(wrong_version, type1), 
        "file")
    @test_throws FileFormatError validate_header(make_header_buf(version, wrong_type), 
        "file")
end

@testset "read_metadata" begin
    md = testdata.metadata

    @test md.instrument.type == "MSD 5977B"
    @test md.instrument.inlet == "unknown"
    @test md.instrument.method == "VM_STD_Auto_SPME.M"

    @test md.acquisition.datetime == "19 Jun 24  05:06 pm"
    @test md.acquisition.sequence == 0
    @test md.acquisition.vial == 23
    @test md.acquisition.replicate == 0

    @test md.user.operator == "unknown"

    @test md.sample.sample == "C7-C40_vial4_dil4_June19_2024"
    @test md.sample.description == "Alkanstandard vial 4 4x verdünnt"
end

# ─── Scan Data Reading ───────────────────────────────────────────────────────

@testset "read_scan_data_ms" begin
    chunk_io = IOBuffer(testdata.first_scan_chunk)
    scantimes, counts, mzs, ints = read_scan_data_ms(chunk_io, 1)
    
    @test eof(chunk_io)

    @test isa(scantimes, Vector{Float32})
    @test isa(counts, Vector{Int32})
    @test isa(mzs, Vector{Float32})
    @test isa(ints, Vector{Int32})

    @test scantimes ≈ Float32[191941.0f0]
    @test counts == Int32[650]
    @test length(mzs) == sum(counts)
    @test length(ints) == sum(counts)
    @test sha1hex(mzs) == "1ea3db1b07b421abbccb556e68c09c1fe6f7b9bc"
    @test sha1hex(ints) == "468a256a60a5f776c413562e56dd69f5eab49337"
    
    seekstart(chunk_io)
    @test_throws FileCorruptionError read_scan_data_ms(chunk_io, 2)
end

@testset "read_scan_data_tic" begin
    scan_time = Int32(191941)
    tic_value = Int32(8_103_078)
    chunk = make_tic_chunk(scan_time, tic_value; payload_bytes=12)
    chunk_io = IOBuffer(chunk)
    scantimes, tic = read_scan_data_tic(chunk_io, 1)
    
    @test eof(chunk_io)

    @test scantimes == Float32[scan_time]
    @test tic == Float32[tic_value]
    
    seekstart(chunk_io)
    @test_throws FileCorruptionError read_scan_data_tic(chunk_io, 2)
end

# ─── Mass Scan Construction ──────────────────────────────────────────────────

@testset "build_mass_scans" begin
    # Sample data
    scantimes = [1.0f0, 2.0f0] * u"ms"
    counts = [3, 2]
    mzs = Float32[100.0, 50.0, 75.0, 200.0, 150.0]  # unsorted within scans
    ints = Int32[10, 20, 30, 40, 50]

    scans = build_mass_scans(ustrip.(scantimes), counts, mzs, ints)

    @test length(scans) == 2

    # Scan data checks
    @test retention.(scans) ≈ scantimes
    @test mzvalues(scans[1]) == sort(mzs[1:3])
    @test intensities(scans[1]) == ints[1:3][sortperm(mzs[1:3])]
    @test mzvalues(scans[2]) == sort(mzs[4:5])
    @test intensities(scans[2]) == ints[4:5][sortperm(mzs[4:5])]

    # Type checks
    @test eltype(mzvalues(first(scans))) == Float32
    @test eltype(intensities(first(scans))) == Int32

    # Edge case: empty input
    scantimes_empty = Vector{typeof(0.0f0)}()
    counts_empty = Int32[]
    mzs_empty = Float32[]
    ints_empty = Int32[]
    empty_result = build_mass_scans(scantimes_empty, counts_empty, mzs_empty, ints_empty)
    @test isempty(empty_result)

    # Edge case: mismatched lengths
    @test_throws DimensionMismatch build_mass_scans([1.0f0, 2.0f0], 
        [1], Float32[100.0], Int32[10])  # mismatched times vs. counts
    @test_throws DimensionMismatch build_mass_scans([1.0f0], [1], Float32[100.0], 
        Int32[10, 20])  # too many ints
end

# ─── Helper functions for binary data reading and byte order conversion ──────

@testset "ensure_bytes_available" begin
    buf = IOBuffer(rand(UInt8, 10))
    seek(buf, 0)
    @test ensure_bytes_available(buf, 5) === nothing
    seek(buf, 8)
    @test_throws UnexpectedEOFError ensure_bytes_available(buf, 4)
end

@testset "read_len_prefixed_string_ltoh" begin
    # Valid string "ABC": length = 3, followed by ASCII codes for 'A', 'B', 'C'
    s = "ABC"
    buf = IOBuffer([UInt8(length(s)); codeunits(s)...])
    @test read_len_prefixed_string_ltoh(buf) == s

    # Empty string: length = 0 means read no bytes, expect ""
    buf_empty = IOBuffer(UInt8[0])
    @test read_len_prefixed_string_ltoh(buf_empty) == ""

    # Invalidly large string length
    s = "ABC"
    buf = IOBuffer([UInt8(length(s)+1); codeunits(s)...])
    @test_throws UnexpectedEOFError read_len_prefixed_string_ltoh(buf) == s
end

@testset "read_len_prefixed_string_ltoh_at" begin
    # Valid string "ABC": length = 3, followed by ASCII codes for 'A', 'B', 'C'
    s = "ABC"
    buf = IOBuffer([UInt8(length(s)); codeunits(s)...])
    @test read_len_prefixed_string_ltoh_at(buf, 0) == s

    # Empty string
    s = ""
    buf = IOBuffer([UInt8(length(s)); codeunits(s)...])
    @test read_len_prefixed_string_ltoh_at(buf, 0) == s

    # Try reading from an invalid position (e.g., beyond buffer)
    @test_throws UnexpectedEOFError read_len_prefixed_string_ltoh_at(buf, 100)
end

@testset "read_scalar_ltoh" begin
    buf = IOBuffer()
    val = Int32(0x12345678)
    write(buf, reinterpret(UInt8, [val]))  # write raw bytes in little-endian order
    seekstart(buf)

    res = read_scalar_ltoh(buf, Int32)
    @test res == val

    # Test reading with insufficient bytes
    buf2 = IOBuffer(UInt8[0x01])
    @test_throws UnexpectedEOFError read_scalar_ltoh(buf2, Int32)
end

@testset "read_bytes_ltoh" begin
    data = UInt8[10, 20, 30, 40, 50]
    buf = IOBuffer(data)
    bytes = read_bytes_ltoh(buf, 3)
    @test bytes == data[1:3]

    # Request more bytes than available
    buf2 = IOBuffer(UInt8[1,2])
    @test_throws UnexpectedEOFError read_bytes_ltoh(buf2, 5)
end

@testset "read_scalar_ntoh" begin
    buf = IOBuffer()
    val = UInt16(0x1234)
    # Write in network byte order (big endian)
    be_bytes = ntoh(val) # or write bytes manually in big endian
    write(buf, [0x12, 0x34])
    seekstart(buf)

    res = read_scalar_ntoh(buf, UInt16)
    @test res == val

    buf2 = IOBuffer(UInt8[0x01])
    @test_throws UnexpectedEOFError read_scalar_ntoh(buf2, UInt16)
end

@testset "read_scalar_ntoh_at" begin
    val = UInt32(0x12345678)
    bytes = UInt8[0x12, 0x34, 0x56, 0x78]  # big endian
    buf = IOBuffer(bytes)

    res = read_scalar_ntoh_at(buf, 0, UInt32)
    @test res == val

    # Invalid position beyond buffer length
    @test_throws UnexpectedEOFError read_scalar_ntoh_at(buf, 10, UInt32)
end

end
