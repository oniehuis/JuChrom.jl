module TestMassHunterMSLoader

using Test
using SHA
using JuChrom
using JuChrom.InputOutput
using JuChrom.MassHunterMSLoader
import JuChrom.MassHunterMSLoader: MassHunterMSOptions, MassHunterMSLoaderSpec,
    MassHunterMSv1, MSScanBinV1, MSPeakBinV1, magicnumber, mspeakdata, msscandata,
    readfile, read_scan_data_ms, read_scan_data_tic, scandatum_attrs,
    MSSCAN_ATTR_EXCLUDES, MSSCAN_ATTR_EXCLUDES_TIC

@testset "MassHunterMSOptions" begin
    options = MassHunterMSOptions(:ms, nothing, Nothing, Nothing, Nothing)
    @test isa(options, MassHunterMSOptions)
    @test options.mode ≡ :ms
    @test options.level ≡ nothing
    @test options.scantimetype ≡ Nothing
    @test options.iontype ≡ Nothing
    @test options.intensitytype ≡ Nothing

    custom = MassHunterMSOptions(:tic, 2, Float32, Float32, Float64)
    @test custom.mode ≡ :tic
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
    @test direct.options ≡ options

    kw = MassHunterMSLoaderSpec{MassHunterMSv1}("kw_path"; mode=:ms, level=1, 
        scantimetype=Float32, iontype=Float32, intensitytype=Float64)
    @test isa(kw, MassHunterMSLoaderSpec)
    @test kw.path == "kw_path"
    @test kw.options.mode ≡ :ms
    @test kw.options.level == 1
    @test kw.options.scantimetype == Float32
    @test kw.options.iontype == Float32
    @test kw.options.intensitytype == Float64

    @test_throws ArgumentError MassHunterMSLoaderSpec{MassHunterMSv1}("dummy"; mode=:bad)
end

@testset "MassHunterMS" begin
    spec = MassHunterMS("example.D"; mode=:tic, level=1)
    @test isa(spec, MassHunterMSLoaderSpec)
    @test spec.path == "example.D"
    @test spec.options.mode ≡ :tic
    @test spec.options.level == 1

    @test_throws ArgumentError MassHunterMS("example.D"; mode=:unknown)
end

# Self-contained schema tests also run before the example-data artifact is refreshed.
include("MassHunterMSSchema.jl")

@testset "Loader entry points and metadata" begin
    mktempdir() do dir
        acq = synthetic_masshunter(dir)
        scanfile = joinpath(acq, "MSScan.bin")
        peakfile = joinpath(acq, "MSPeak.bin")
        data = msscandata(MSScanBinV1(), scanfile)
        @test magicnumber(scanfile) == 257
        @test magicnumber(peakfile) == 259
        @test length(data) == 2
        @test first(data).ScanTime == 0.5
        @test last(data).ScanTime == 1.0
        @test first(data).PointCount == 3
        @test first(data).ByteCount == 48

        for mode in (:ms, :tic)
            options = MassHunterMSOptions(mode, nothing, Nothing, Nothing, Nothing)
            series = readfile(MassHunterMSv1, acq, options)
            @test series isa (mode == :ms ? MassScanSeries : ChromScanSeries)
            @test scancount(series) == 2
            @test unit(first(retentions(series))) == u"minute"
            @test ustrip.(retentions(series)) == [0.5, 1.0]
            @test instrument(series) == NamedTuple()
            @test acquisition(series) == NamedTuple()
            @test user(series) == NamedTuple()
            @test sample(series) == NamedTuple()
            @test extras(series) == Dict()
        end

        attrs_ms = scandatum_attrs(first(data))
        attrs_tic = scandatum_attrs(first(data); include_tic=false)
        @test haskey(attrs_ms, :TIC)
        @test !haskey(attrs_tic, :TIC)
        @test all(name -> name ∉ MSSCAN_ATTR_EXCLUDES, keys(attrs_ms))
        @test all(name -> name ∉ MSSCAN_ATTR_EXCLUDES_TIC, keys(attrs_tic))
        @test attrs_ms.ScanID == attrs_tic.ScanID == 1
        @test attrs_ms.MinX == 40
        @test attrs_ms.MinY == 2

        options_tic = MassHunterMSOptions(:tic, nothing, Float32, Nothing, Float32)
        times, tics, scan_attrs = read_scan_data_tic(MassHunterMSv1, data, options_tic)
        @test times == Float32[0.5, 1.0]
        @test eltype(times) == Float32
        @test tics == Float32[10, 20]
        @test eltype(tics) == Float32
        @test length(scan_attrs) == 2
        @test first(scan_attrs).ScanID == 1

        options_ms = MassHunterMSOptions(:ms, nothing, Float32, Float32, Float32)
        mss = read_scan_data_ms(MassHunterMSv1, acq, data, options_ms)
        @test unit(first(retentions(mss))) == u"minute"
        @test ustrip(retention(first(scans(mss)))) === 0.5f0
        @test eltype(mzvalues(first(scans(mss)))) == Float32
        @test eltype(intensities(first(scans(mss)))) == Float32
        decoded = mspeakdata(MSPeakBinV1(), peakfile, data, options_ms)
        @test length(decoded) == 2
        @test retention(first(decoded)) == 0.5u"minute"
        @test mzvalues(first(decoded)) == [40, 50, 60]
        @test intensities(first(decoded)) == [2, 3, 5]
        @test_throws ArgumentError load(MassHunterMS(dir; level=3))
        @test_throws ArgumentError load(MassHunterMS(scanfile))
    end
end

@testset "Invalid paths and signatures" begin
    mktempdir() do dir
        @test_throws MissingFolderError load(MassHunterMS(dir))
        @test_throws MissingFolderError load(MassHunterMS(joinpath(dir, "Missing.D")))
        acq = synthetic_masshunter(dir)
        write(joinpath(acq, "MSScan.bin"), zeros(UInt8, 32))
        @test_throws FileFormatError load(MassHunterMS(dir))
    end
end

@testset "Bundled MassHunter acquisition with its original schema" begin
    example = joinpath(JuChrom.agilent, "C7-C40_MassHunterMS.D")
    schemafile = joinpath(example, "AcqData", "MSScan.xsd")
    if !isfile(schemafile)
        # The v0.0.1 artifact omitted the schema; do not invent a replacement for it.
        @test_skip isfile(schemafile)
    else
        mss = load(MassHunterMS(example))
        tic = load(MassHunterMS(example; mode=:tic))
        @test scancount(mss) == scancount(tic) == 2405
        @test unit(first(retentions(mss))) == u"minute"
        @test retentions(mss) == retentions(tic)
        @test ustrip(first(retentions(mss))) == 3.1990166666666666
        @test ustrip(last(retentions(mss))) == 31.650783333333333
        sha1hex(v) = bytes2hex(sha1(reinterpret(UInt8, v)))
        @test sha1hex(mzvalues(first(scans(mss)))) == "ff19749a4dab34a8b6af801096a99a9c7acb6dd0"
        @test sha1hex(intensities(first(scans(mss)))) == "57c82ff5c18ce9f2d17bfd4b9e6ceeb792d64dc5"
        @test sha1hex(mzvalues(last(scans(mss)))) == "052afb9e54521cd128312325a28fc5e08c11827f"
        @test sha1hex(intensities(last(scans(mss)))) == "a9804cb39164f94d380e131b4f56906584d6d8ec"
        @test sha1hex(intensities(tic)) == "17da30bffc6a3bfb79ff764b2d5f8963f3c791be"
    end
end

end  # module TestMassHunterMSLoader
