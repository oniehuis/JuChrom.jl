module TestAgilentFIDLoader

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for AgilentFIDLoader.jl
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

using SHA
using Logging
using StringEncodings
using Test
using Unitful

# ── Internal Project Imports ─────────────────────────────────────────────────

using JuChrom
using JuChrom.InputOutput
using JuChrom.AgilentFIDLoader
import JuChrom.AgilentFIDLoader: AgilentFIDLoaderSpec, AgilentFIDv179, readfile,
    validate_header, read_metadata, read_scalar_ntoh_at, read_vector_ntoh_at,
    read_len_prefixed_string_ltoh_at, read_len_prefixed_string_ltoh,
    read_scalar_ltoh, read_vector_ltoh_at, ensure_bytes_available

# ── Constants ────────────────────────────────────────────────────────────────

const TESTDIR = joinpath(JuChrom.agilent, "ZK_ONUBE_Mix1_11.D")
const TESTFILE = joinpath(TESTDIR, "FID1A.ch")

sha1hex(v) = bytes2hex(sha1(reinterpret(UInt8, v)))

function make_fid_buffer(; scancount::UInt32=UInt32(1), intensityunit::String="pA",
        scalingfactor::Float64=1.0)
    buf = zeros(UInt8, 7000)

    version = "179"
    header = [UInt8(length(version)); codeunits(version)...]
    copyto!(buf, 1, header, 1, length(header))

    type = "GC data file"
    type_bytes = encode(type, AgilentFIDLoader.ENCODING_2_bytes)
    buf[347 + 1] = UInt8(length(type))
    copyto!(buf, 347 + 2, type_bytes, 1, length(type_bytes))

    sc_bytes = reinterpret(UInt8, [hton(scancount)])
    copyto!(buf, 278 + 1, sc_bytes, 1, length(sc_bytes))

    iu_bytes = encode(intensityunit, AgilentFIDLoader.ENCODING_2_bytes)
    buf[4172 + 1] = UInt8(length(intensityunit))
    copyto!(buf, 4172 + 2, iu_bytes, 1, length(iu_bytes))

    sf_bits = hton(reinterpret(UInt64, scalingfactor))
    sf_bytes = reinterpret(UInt8, [sf_bits])
    copyto!(buf, 4732 + 1, sf_bytes, 1, length(sf_bytes))

    buf
end

struct BadBytes end

struct BadDecodeBuffer <: IO
    buf::IOBuffer
end

Base.position(f::BadDecodeBuffer) = position(f.buf)
Base.seek(f::BadDecodeBuffer, pos::Integer) = seek(f.buf, pos)
Base.seekend(f::BadDecodeBuffer) = seekend(f.buf)
Base.read(f::BadDecodeBuffer, ::Type{UInt8}) = read(f.buf, UInt8)
Base.read(f::BadDecodeBuffer, n::Integer) = read(f.buf, n)

JuChrom.AgilentFIDLoader.read_bytes_ltoh(::BadDecodeBuffer, ::Integer) = BadBytes()

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests
# ─────────────────────────────────────────────────────────────────────────────

@testset "AgilentFIDLoaderSpec" begin
    spec = AgilentFIDLoaderSpec{AgilentFIDv179}("dummy_path"; unit=nothing)
    @test isa(spec, AgilentFIDLoaderSpec)
    @test spec.path == "dummy_path"
    @test spec.options.unit === nothing

    spec2 = AgilentFID("dummy_path"; unit=u"pA")
    @test isa(spec2, AgilentFIDLoaderSpec)
    @test spec2.options.unit == u"pA"

    @test_throws ArgumentError AgilentFID("dummy_path"; unit=u"mA")
end

@testset "validate_header" begin
    # Happy path using real file
    open(TESTFILE, "r") do io
        @test validate_header(io, TESTFILE)
    end

    # Tampered version
    buf = IOBuffer(read(TESTFILE); read=true, write=true)
    seek(buf, 0)
    write(buf, UInt8(1))  # length already 3, overwrite to 1 => wrong version
    seekstart(buf)
    @test_throws FileFormatError validate_header(buf, "file")

    # Tampered type (overwrite with shorter length to stay within buffer)
    buf2 = IOBuffer(read(TESTFILE); read=true, write=true)
    seek(buf2, 347)
    write(buf2, UInt8(4))
    write(buf2, encode("WRNG", AgilentFIDLoader.ENCODING_2_bytes))
    seekstart(buf2)
    @test_throws FileFormatError validate_header(buf2, "file")
end

@testset "read_metadata" begin
    open(TESTFILE, "r") do io
        validate_header(io, TESTFILE)
        md = read_metadata(io)
        @test md.instrument.type == "GC"
        @test md.instrument.inlet == "7890"
        @test md.instrument.method == "GC-EAD.M"
        @test md.instrument.software == "Mustang ChemStation"
        @test md.instrument.signalunit == "Front Signal"
        @test md.acquisition.datetime == "14 Jan 20  02:40 pm"
        @test md.user.operator == "unknown"
        @test md.sample.sample == "ZK_ONUBE_Mix1_11"
        @test md.sample.description == "unknown"
    end
end

@testset "load/readfile" begin
    # Silence expected warnings from negative intensities in the fixture
    series = with_logger(SimpleLogger(IOBuffer(), Logging.Error)) do
        load(AgilentFID(TESTFILE))
    end
    @test isa(series, ChromScanSeries)
    @test length(series) == 4151
    @test unit(first(retentions(series))) == u"ms"
    @test sha1hex(collect(ustrip.(retentions(series)))) == "a8dfb8d4f5d69a6f5370bc1a2dbafa60dcac6d8d"
    @test sha1hex(collect(intensity.(series))) == "4cc493d8b068d4858603537795840749e29dfbdd"
    @test series.instrument.type == "GC"
    @test series.sample.sample == "ZK_ONUBE_Mix1_11"

    # Directory-based load
    series2 = with_logger(SimpleLogger(IOBuffer(), Logging.Error)) do
        load(AgilentFID(TESTDIR))
    end
    @test length(series2) == length(series)

    missing = joinpath(TESTDIR, "missing_FID1A.ch")
    isfile(missing) && rm(missing)
    @test_throws MissingFileError load(AgilentFID(missing))

    mktemp() do path, io
        close(io)
        err = try
            load(AgilentFID(path))
            nothing
        catch e
            e
        end
        @test err isa FileCorruptionError
        @test occursin("Failed to read FID data at", sprint(showerror, err))
    end
end

@testset "read_len_prefixed_string_ltoh_at" begin
    open(TESTFILE, "r") do io
        @test read_len_prefixed_string_ltoh_at(io, 4172, 2) == "pA"
    end
end

@testset "read_len_prefixed_string_ltoh" begin
    buf = IOBuffer([UInt8(1), 0x41])
    @test read_len_prefixed_string_ltoh(buf, 1) == "A"

    buf_bad = IOBuffer(UInt8[1, 0x01, 0x02, 0x03])
    err_bad = try
        read_len_prefixed_string_ltoh(buf_bad, 3)
        nothing
    catch e
        e
    end
    @test err_bad isa FileDecodingError
    @test occursin("Unsupported byte size", sprint(showerror, err_bad))

    bad_decode = BadDecodeBuffer(IOBuffer([UInt8(1), 0xFF]))
    err_decode = try
        read_len_prefixed_string_ltoh(bad_decode, 1)
        nothing
    catch e
        e
    end
    @test err_decode isa FileDecodingError
    @test occursin("Could not decode 1 bytes", sprint(showerror, err_decode))
end

@testset "read_scalar_ntoh_at" begin
    open(TESTFILE, "r") do io
        @test read_scalar_ntoh_at(io, 278, UInt32) == 4151
    end
end

@testset "readfile guards and unit handling" begin
    mktemp() do path, io
        buf = make_fid_buffer(scancount=UInt32(0))
        write(io, buf)
        close(io)
        err = try
            readfile(path, u"pA")
            nothing
        catch e
            e
        end
        @test err isa FileCorruptionError
        @test occursin("Suspicious scan count", sprint(showerror, err))
    end

    mktemp() do path, io
        buf = make_fid_buffer(scancount=UInt32(1_000_001))
        write(io, buf)
        close(io)
        err = try
            readfile(path, u"pA")
            nothing
        catch e
            e
        end
        @test err isa FileCorruptionError
        @test occursin("Suspicious scan count", sprint(showerror, err))
    end

    mktemp() do path, io
        buf = make_fid_buffer(scancount=UInt32(1), intensityunit="pA",
            scalingfactor=-1.0)
        write(io, buf)
        close(io)
        err = try
            readfile(path, u"pA")
            nothing
        catch e
            e
        end
        @test err isa FileCorruptionError
        @test occursin("Negative scaling factor", sprint(showerror, err))
    end

    mktemp() do path, io
        buf = make_fid_buffer(scancount=UInt32(1), intensityunit="pA",
            scalingfactor=1.0)
        write(io, buf)
        close(io)
        err = try
            readfile(path, u"nA")
            nothing
        catch e
            e
        end
        @test err isa FileCorruptionError
        @test occursin("Unsupported signal unit", sprint(showerror, err))
    end

    mktempdir() do dir
        err = try
            readfile(dir, u"pA")
            nothing
        catch e
            e
        end
        @test err isa FileCorruptionError
        @test occursin("Failed to open/read Agilent FID file", sprint(showerror, err))
    end
end

@testset "ensure_bytes_available" begin
    buf = IOBuffer(UInt8[0x01])
    @test_throws UnexpectedEOFError ensure_bytes_available(buf, 4)
end

@testset "read_vector_ntoh_at" begin
    open(TESTFILE, "r") do io
        vals = read_vector_ntoh_at(io, 282, Float32, 2) # start/stop
        @test vals[1] ≈ 0.437f0
        @test vals[2] ≈ 830000.4375f0
    end
end

@testset "read_vector_ltoh_at" begin
    open(TESTFILE, "r") do io
        ints = read_vector_ltoh_at(io, 6144, Float64, 2)
        @test ints[1] == 0.0
        @test ints[2] == 5.0
    end
end

end
