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
    read_len_prefixed_string_ltoh_at, read_scalar_ltoh, read_vector_ltoh_at

# ── Constants ────────────────────────────────────────────────────────────────

const TESTDIR = joinpath(JuChrom.agilent, "ZK_ONUBE_Mix1_11.D")
const TESTFILE = joinpath(TESTDIR, "FID1A.ch")

sha1hex(v) = bytes2hex(sha1(reinterpret(UInt8, v)))

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
end

@testset "read_len_prefixed_string_ltoh_at" begin
    open(TESTFILE, "r") do io
        @test read_len_prefixed_string_ltoh_at(io, 4172, 2) == "pA"
    end
end

@testset "read_scalar_ntoh_at" begin
    open(TESTFILE, "r") do io
        @test read_scalar_ntoh_at(io, 278, UInt32) == 4151
    end
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
