# Small synthetic files, independent of the private standards and example-data artifact.
# This file is included inside TestMassHunterMSLoader.
const MH = JuChrom.MassHunterMSLoader
const SCHEMA_PREFIX_FIELDS = [
    "ScanID" => "int", "ScanMethodID" => "int", "TimeSegmentID" => "int",
    "ScanTime" => "double", "MSLevel" => "int", "ScanType" => "int",
    "TIC" => "double", "BasePeakMZ" => "double", "BasePeakValue" => "double",
    "CycleNumber" => "int", "Status" => "int", "IonMode" => "int", "IonPolarity" => "int"]
const SCHEMA_EXTRA_FIELDS = [
    "Fragmentor" => "double", "CollisionEnergy" => "double", "MzOfInterest" => "double",
    "SamplingPeriod" => "double", "MeasuredMassRangeMin" => "double",
    "MeasuredMassRangeMax" => "double", "Threshold" => "double",
    "IsFragmentorDynamic" => "int", "IsCollisionEnergyDynamic" => "int"]
const SCHEMA_SPECTRUM_FIELDS = [
    "SpectrumFormatID" => "int", "SpectrumOffset" => "long", "ByteCount" => "int",
    "PointCount" => "int", "MinX" => "double", "MaxX" => "double",
    "MinY" => "double", "MaxY" => "double"]
const FIXTURE_TYPES = Dict("int" => Int32, "long" => Int64, "double" => Float64)

writele(io, x::Integer) = write(io, htol(x))
writele(io, x::Float32) = writele(io, reinterpret(UInt32, x))
writele(io, x::Float64) = writele(io, reinterpret(UInt64, x))

function synthetic_masshunter(dir; wide=false, unbounded=false,
    start=160, storage=(Float64, Float64), spectrum_fields=SCHEMA_SPECTRUM_FIELDS)
    acq = mkpath(joinpath(dir, "AcqData"))
    extras = wide ? SCHEMA_EXTRA_FIELDS : ["SamplingPeriod" => "double"]
    spectrum = copy(spectrum_fields)
    fields = [SCHEMA_PREFIX_FIELDS; extras; spectrum]
    element(p) = "<x:element name=\"$(first(p))\" type=\"x:$(last(p))\"/>"
    prefix = join(element.(SCHEMA_PREFIX_FIELDS))
    extra = join(element.(extras))
    spec = join(element.(spectrum))
    occurs = unbounded ? " maxOccurs=\"unbounded\"" : ""
    xml = """
    <x:schema xmlns:x="http://www.w3.org/2001/XMLSchema" xmlns:v="urn:masshunter:test"
              xmlns="urn:masshunter:test" targetNamespace="urn:masshunter:test">
      <x:complexType name="ScanRecordType"><x:sequence>
        $prefix $extra
        <!-- <x:element name="Fragmentor" type="x:double"/>
             <x:element name="CollisionEnergy" type="x:double"/> -->
        <x:element name="SpectrumParamValues" type="v:SpectrumParamsType"$occurs/>
      </x:sequence></x:complexType>
      <x:complexType name="SpectrumParamsType"><x:sequence>$spec</x:sequence></x:complexType>
    </x:schema>
    """
    write(joinpath(acq, "MSScan.xsd"), xml)
    write(joinpath(acq, "MSTS.xml"), """
        <TimeSegments><TimeSegment TimeSegmentID="1"><NumOfScans>2</NumOfScans>
        </TimeSegment></TimeSegments>""")
    records = Dict{String,Real}[]
    open(joinpath(acq, "MSPeak.bin"), "w") do io
        writele(io, UInt16(259))
        write(io, zeros(UInt8, 66))
        for i in 1:2
            offset = position(io)
            T = storage[i]
            for x in T[40, 50, 60, 2i, 3i, 5i]
                writele(io, x)
            end
            record = Dict{String,Real}(name => 0 for (name, _) in fields)
            merge!(record, Dict("ScanID" => i, "ScanTime" => 0.5i,
                "MSLevel" => i, "TIC" => 10i, "BasePeakMZ" => 60,
                "BasePeakValue" => 5i, "SpectrumFormatID" => 1,
                "SpectrumOffset" => offset, "ByteCount" => 6sizeof(T), "PointCount" => 3,
                "MinX" => 40, "MaxX" => 60, "MinY" => 2i, "MaxY" => 5i))
            wide && (record["Fragmentor"] = 71)
            push!(records, record)
        end
    end
    open(joinpath(acq, "MSScan.bin"), "w") do io
        write(io, zeros(UInt8, start))
        seekstart(io)
        writele(io, UInt16(257))
        seek(io, 88)
        writele(io, UInt32(start))
        seek(io, start)
        for record in records, (name, xtype) in fields
            writele(io, FIXTURE_TYPES[xtype](record[name]))
        end
    end
    acq
end

function change_schema(acq, old, new)
    file = joinpath(acq, "MSScan.xsd")
    write(file, replace(read(file, String), old => new))
end
function change_scan_field(acq, name, value; index=1)
    layout = MH.scanlayout(joinpath(acq, "MSScan.xsd"))
    names = fieldnames(layout.recordtype)
    n = findfirst(==(Symbol(name)), names)
    offset = sum(MH.packedsize(fieldtype(layout.recordtype, i)) for i in 1:n-1; init=0)
    file = joinpath(acq, "MSScan.bin")
    open(file, "r+") do io
        seek(io, 88)
        start = MH.readle(io, UInt32)
        seek(io, start + (index-1)*layout.recordsize + offset)
        writele(io, convert(fieldtype(layout.recordtype, n), value))
    end
end

@testset "Schema-driven MassHunter records" begin
    for wide in (false, true), T in (Float32, Float64)
        mktempdir() do dir
            acq = synthetic_masshunter(dir; wide, storage=(T, T))
            layout = MH.scanlayout(joinpath(acq, "MSScan.xsd"))
            @test layout.recordsize == (wide ? 184 : 128)
            data = msscandata(MSScanBinV1(), joinpath(acq, "MSScan.bin"))
            @test data[1].MinX == 40
            @test data[1].MaxX == 60
            @test data[1].MinY == 2
            @test data[1].MaxY == 5
            @test hasproperty(data[1], :Fragmentor) == wide
            wide && @test data[1].Fragmentor == 71
            mss = load(MassHunterMS(dir))
            css = load(MassHunterMS(dir; mode=:tic))
            @test scancount(mss) == scancount(css) == 2
            @test retentions(mss) == retentions(css) == [0.5, 1.0]u"minute"
            @test unit(first(retentions(mss))) == u"minute"
            @test unit(first(retentions(css))) == u"minute"
            @test ustrip.(retentions(mss)) == [s.ScanTime for s in data]
            @test intensities(css) == [10, 20]
            for (i, scan) in enumerate(scans(mss))
                @test eltype(mzvalues(scan)) == T
                @test eltype(intensities(scan)) == T
                @test mzvalues(scan) == [40, 50, 60]
                @test intensities(scan) == [2i, 3i, 5i]
                @test attrs(scan).TIC == sum(intensities(scan))
                @test attrs(scan).BasePeakValue == maximum(intensities(scan))
                @test attrs(scan).MinX == 40
                @test attrs(scan).MinY == 2i
                @test level(scan) == i
            end
            @test !haskey(attrs(first(scans(css))), :TIC)
            for mode in (:ms, :tic)
                selected = load(MassHunterMS(dir; mode, level=2, scantimetype=Float32,
                    iontype=Float32, intensitytype=Float32))
                @test scancount(selected) == 1
                scan = only(scans(selected))
                @test ustrip(retention(scan)) === 1f0
                @test unit(retention(scan)) == u"minute"
                @test attrs(scan).ScanID == 2
                if mode == :ms
                    @test eltype(mzvalues(scan)) == Float32
                    @test eltype(intensities(scan)) == Float32
                else
                    @test only(intensities(selected)) === 20f0
                end
            end
        end
    end
    @testset "Schema is required for every layout" begin
        for wide in (false, true)
            mktempdir() do dir
                acq = synthetic_masshunter(dir; wide)
                schemafile = joinpath(acq, "MSScan.xsd")
                rm(schemafile)
                for mode in (:ms, :tic)
                    err = try
                        load(MassHunterMS(dir; mode))
                    catch e
                        e
                    end
                    @test err isa MissingFileError
                    @test err.path == schemafile
                end
            end
        end
    end
    @testset "32-bit header offset and mixed precision" begin
        mktempdir() do dir
            synthetic_masshunter(dir; start=65536, storage=(Float32, Float64))
            mss = load(MassHunterMS(dir))
            @test scancount(mss) == 2
            @test eltype(mzvalues(first(scans(mss)))) == Float64
            @test eltype(mzvalues(last(scans(mss)))) == Float64
        end
    end
    @testset "Documented vendor repetitions and namespaces" begin
        mktempdir() do dir
            # The original 184-byte example puts the Y extrema before the X extrema.
            spectrum_fields = [SCHEMA_SPECTRUM_FIELDS[1:4];
                SCHEMA_SPECTRUM_FIELDS[7:8]; SCHEMA_SPECTRUM_FIELDS[5:6]]
            acq = synthetic_masshunter(dir; wide=true, unbounded=true, spectrum_fields)
            layout = MH.scanlayout(joinpath(acq, "MSScan.xsd"))
            @test layout.recordsize == 184
            @test fieldnames(layout.recordtype)[end-3:end] == (:MinY, :MaxY, :MinX, :MaxX)
            for mode in (:ms, :tic)
                series = load(MassHunterMS(dir; mode))
                @test scancount(series) == 2
                for (i, scan) in enumerate(scans(series))
                    @test attrs(scan).MinX == 40
                    @test attrs(scan).MaxX == 60
                    @test attrs(scan).MinY == 2i
                    @test attrs(scan).MaxY == 5i
                end
            end
        end
        mktempdir() do dir
            acq = synthetic_masshunter(dir; unbounded=true)
            change_schema(acq, "type=\"v:SpectrumParamsType\"", "type=\"SpectrumParamsType\"")
            change_schema(acq, "name=\"Status\" type=\"x:int\"", "name=\"Status\" type=\"x:int\" minOccurs=\"0\"")
            @test scancount(load(MassHunterMS(dir))) == 2
            # An independent count rules out treating a genuinely variable record as fixed.
            write(joinpath(acq, "MSTS.xml"), "<TimeSegments><TimeSegment><NumOfScans>1</NumOfScans></TimeSegment></TimeSegments>")
            @test_throws FileFormatError load(MassHunterMS(dir))
            rm(joinpath(acq, "MSTS.xml"))
            @test_throws FileFormatError load(MassHunterMS(dir))
        end
    end
    @testset "Fixed primitive repetitions and deeper nesting" begin
        mktempdir() do dir
            acq = synthetic_masshunter(dir)
            change_schema(acq, "<x:element name=\"ScanMethodID\" type=\"x:int\"/><x:element name=\"TimeSegmentID\" type=\"x:int\"/>",
                "<x:element name=\"MethodAndSegment\" type=\"x:int\" minOccurs=\"2\" maxOccurs=\"2\"/>")
            change_schema(acq, "<x:element name=\"SpectrumParamValues\" type=\"v:SpectrumParamsType\"/>",
                "<x:element name=\"Wrapper\" type=\"v:WrapperType\"/>")
            change_schema(acq, "</x:schema>", "<x:complexType name=\"WrapperType\"><x:sequence><x:element name=\"Params\" type=\"v:SpectrumParamsType\"/></x:sequence></x:complexType></x:schema>")
            @test MH.scanlayout(joinpath(acq, "MSScan.xsd")).recordsize == 128
            @test attrs(first(scans(load(MassHunterMS(dir))))).MethodAndSegment == (0, 0)
        end
    end
end

@testset "Schema and index failures" begin
    for (old, new) in (
        ("x:sequence", "x:choice"),
        ("type=\"x:int\"", "type=\"x:string\""),
        ("type=\"x:int\"", "type=\"missing:int\""),
        ("xmlns:v=\"urn:masshunter:test\"", "xmlns:v=\"urn:wrong\""),
        ("name=\"ScanID\"", "name=\"ScanTime\""),
        ("name=\"ScanID\"", "name=\"ScanID\" maxOccurs=\"unbounded\""),
        ("name=\"ScanID\"", "name=\"ScanID\" minOccurs=\"0\" maxOccurs=\"2\""),
        ("name=\"ScanID\"", "name=\"ScanID\" nillable=\"true\""),
        ("name=\"ScanID\"", "ref=\"v:ScanID\""),
        ("type=\"v:SpectrumParamsType\"", "type=\"v:ScanRecordType\""),
        ("type=\"v:SpectrumParamsType\"", "type=\"v:SpectrumParamsType\" minOccurs=\"2\" maxOccurs=\"2\""),
        ("type=\"v:SpectrumParamsType\"", "type=\"v:SpectrumParamsType\" minOccurs=\"0\""),
        ("name=\"ScanRecordType\"", "name=\"OtherRecordType\""),
        ("</x:schema>", "<x:include schemaLocation=\"external.xsd\"/></x:schema>"),
        ("</x:schema>", "<x:complexType name=\"ScanRecordType\"/></x:schema>"),
        ("</x:schema>", "</broken>"))
        mktempdir() do dir
            acq = synthetic_masshunter(dir)
            change_schema(acq, old, new)
            @test_throws FileFormatError load(MassHunterMS(dir))
        end
    end
    for offset in (0, 88, 91, 100_000, typemax(UInt32))
        mktempdir() do dir
            acq = synthetic_masshunter(dir)
            open(joinpath(acq, "MSScan.bin"), "r+") do io
                seek(io, 88)
                writele(io, UInt32(offset))
            end
            E = offset < 92 ? FileCorruptionError : UnexpectedEOFError
            @test_throws E load(MassHunterMS(dir))
        end
    end
    for length in (0, 1, 90, 160+127, 160+128)
        mktempdir() do dir
            acq = synthetic_masshunter(dir)
            open(io -> truncate(io, length), joinpath(acq, "MSScan.bin"), "r+")
            E = length == 160+128 ? FileCorruptionError : UnexpectedEOFError
            @test_throws E load(MassHunterMS(dir))
        end
    end
    mktempdir() do dir
        acq = synthetic_masshunter(dir)
        rm(joinpath(acq, "MSScan.bin"))
        @test_throws MissingFileError load(MassHunterMS(dir))
    end
end

@testset "Spectrum errors are independent of scan-index parsing" begin
    for (field, value, E) in (
        ("SpectrumOffset", -1, FileCorruptionError),
        ("SpectrumOffset", 0, FileCorruptionError),
        ("SpectrumOffset", typemax(Int64), UnexpectedEOFError),
        ("SpectrumOffset", 150, UnexpectedEOFError),
        ("ByteCount", -1, FileCorruptionError),
        ("ByteCount", typemax(Int32), UnexpectedEOFError),
        ("PointCount", -1, FileCorruptionError),
        ("PointCount", 0, FileCorruptionError),
        ("PointCount", 5, FileDecodingError),
        ("ByteCount", 12, FileDecodingError),
        ("SpectrumFormatID", 99, FileDecodingError))
        mktempdir() do dir
            acq = synthetic_masshunter(dir)
            change_scan_field(acq, field, value)
            @test scancount(load(MassHunterMS(dir; mode=:tic))) == 2
            @test_throws E load(MassHunterMS(dir))
        end
    end
    mktempdir() do dir
        acq = synthetic_masshunter(dir)
        change_scan_field(acq, "SpectrumFormatID", 99; index=2)
        @test scancount(load(MassHunterMS(dir; level=1))) == 1
        @test_throws FileDecodingError load(MassHunterMS(dir; level=2))
        rm(joinpath(acq, "MSPeak.bin"))
        @test scancount(load(MassHunterMS(dir; mode=:tic))) == 2
        @test_throws MissingFileError load(MassHunterMS(dir))
    end
    mktempdir() do dir
        acq = synthetic_masshunter(dir)
        file = joinpath(acq, "MSPeak.bin")
        open(io -> truncate(io, filesize(file)-1), file, "r+")
        @test_throws UnexpectedEOFError load(MassHunterMS(dir))
    end
    mktempdir() do dir
        acq = synthetic_masshunter(dir)
        change_scan_field(acq, "PointCount", 0)
        change_scan_field(acq, "ByteCount", 0)
        @test_throws FileDecodingError load(MassHunterMS(dir))
        @test scancount(load(MassHunterMS(dir; mode=:tic))) == 2
    end
end

@testset "Binary read failures retain file and scan context" begin
    mktempdir() do dir
        acq = synthetic_masshunter(dir)
        file = joinpath(acq, "MSScan.bin")
        layout = MH.scanlayout(joinpath(acq, "MSScan.xsd"))
        bytes = read(file)[161:end]
        for nbytes in (0, layout.recordsize - 1, layout.recordsize, length(bytes) - 1)
            err = try
                MH.read_scan_records(IOBuffer(bytes[1:nbytes]), layout.recordtype, 2, file)
            catch e
                e
            end
            @test err isa UnexpectedEOFError
            @test sprint(showerror, err) == "Truncated scan record in $file"
        end
        # Closed streams are programming/I/O errors, not truncated acquisitions.
        io = IOBuffer(bytes)
        close(io)
        @test_throws ArgumentError MH.read_scan_records(io, layout.recordtype, 2, file)
    end
    for T in (Float32, Float64)
        io = IOBuffer()
        for value in T[40, 50, 60, 2, 3, 5]
            writele(io, value)
        end
        bytes = take!(io)
        file = "MSPeak.bin"
        # Cut within the m/z array, between arrays, and within the intensity array.
        for nbytes in (0, 3sizeof(T) - 1, 3sizeof(T), length(bytes) - 1)
            err = try
                MH.read_spectrum(IOBuffer(bytes[1:nbytes]), T, 3, 7, file)
            catch e
                e
            end
            @test err isa UnexpectedEOFError
            @test sprint(showerror, err) == "Truncated spectrum at scan record 7 in $file"
        end
        io = IOBuffer(bytes)
        close(io)
        @test_throws ArgumentError MH.read_spectrum(io, T, 3, 7, file)
    end
end

@testset "Spectrum blocks at the file boundary" begin
    block = (SpectrumOffset=Int64(68), ByteCount=Int32(48), PointCount=Int32(3))
    @test MH.validate_spectrum_bounds(block, 116, 1) === nothing
    @test_throws UnexpectedEOFError MH.validate_spectrum_bounds(block, 115, 1)
    emptyblock = (SpectrumOffset=Int64(116), ByteCount=Int32(0), PointCount=Int32(0))
    # Empty blocks have valid bounds even though they cannot become a MassScan.
    @test MH.validate_spectrum_bounds(emptyblock, 116, 2) === nothing
end

@testset "Primitive endianness and required metadata" begin
    # Fixed byte sequences exercise signed values and IEEE floating-point bit patterns.
    for (T, bytes, value) in (
        (Int8, "fe", -2), (UInt8, "fe", 254),
        (Int16, "feff", -2), (UInt16, "3412", 0x1234),
        (Int32, "feffffff", -2), (UInt32, "78563412", 0x12345678),
        (Int64, "feffffffffffffff", -2), (UInt64, "8877665544332211", 0x1122334455667788),
        (Float32, "0000c03f", 1.5), (Float64, "000000000000f83f", 1.5))
        @test MH.readle(IOBuffer(hex2bytes(bytes)), T) === T(value)
    end
    for name in ("ScanTime", "MSLevel", "SpectrumOffset", "TIC")
        mktempdir() do dir
            acq = synthetic_masshunter(dir)
            change_schema(acq, "name=\"$name\"", "name=\"OtherValue\"")
            if name == "SpectrumOffset"
                @test scancount(load(MassHunterMS(dir; mode=:tic))) == 2
                @test_throws FileFormatError load(MassHunterMS(dir))
            elseif name == "TIC"
                @test scancount(load(MassHunterMS(dir))) == 2
                @test !haskey(attrs(first(scans(load(MassHunterMS(dir))))), :TIC)
                @test_throws FileFormatError load(MassHunterMS(dir; mode=:tic))
            else
                @test_throws FileFormatError load(MassHunterMS(dir))
            end
        end
    end
    for (name, value, mode) in (("ScanTime", NaN, :ms), ("MSLevel", 0, :tic), ("TIC", Inf, :tic))
        mktempdir() do dir
            acq = synthetic_masshunter(dir)
            change_scan_field(acq, name, value)
            @test_throws FileCorruptionError load(MassHunterMS(dir; mode))
        end
    end
    mktempdir() do dir
        acq = synthetic_masshunter(dir)
        change_schema(acq, "<x:sequence>", "<x:sequence>unexpected text")
        @test_throws FileFormatError load(MassHunterMS(dir))
    end
    mktempdir() do dir
        acq = synthetic_masshunter(dir)
        file = joinpath(acq, "MSScan.xsd")
        write(file, "<!DOCTYPE x:schema [<!ENTITY field 'unused'>]>" * read(file, String))
        @test_throws FileFormatError load(MassHunterMS(dir))
    end
end

@testset "Independent scan-count verification" begin
    mktempdir() do dir
        acq = synthetic_masshunter(dir; unbounded=true)
        file = joinpath(acq, "MSTS.xml")
        write(file, """
            <t:TimeSegments xmlns:t="urn:segments"><!-- <t:NumOfScans>100</t:NumOfScans> -->
            <t:TimeSegment><t:NumOfScans>1</t:NumOfScans></t:TimeSegment>
            <t:TimeSegment><t:NumOfScans>1</t:NumOfScans></t:TimeSegment></t:TimeSegments>""")
        @test scancount(load(MassHunterMS(dir))) == 2
        for (xml, E) in (
            ("<TimeSegments/>", FileFormatError),
            ("<TimeSegments><TimeSegment/></TimeSegments>", FileFormatError),
            ("<TimeSegments><TimeSegment><NumOfScans>-1</NumOfScans></TimeSegment></TimeSegments>", FileCorruptionError),
            ("<TimeSegments><TimeSegment><NumOfScans>3</NumOfScans></TimeSegment></TimeSegments>", FileCorruptionError))
            write(file, xml)
            @test_throws E load(MassHunterMS(dir))
        end
    end
    mktempdir() do dir
        acq = synthetic_masshunter(dir)
        rm(joinpath(acq, "MSTS.xml"))
        @test scancount(load(MassHunterMS(dir))) == 2
        @test scancount(load(MassHunterMS(dir; mode=:tic))) == 2
    end
    mktempdir() do dir
        acq = synthetic_masshunter(dir; unbounded=true)
        # Repetitions other than the verified spectrum block are never collapsed to one.
        change_schema(acq, "name=\"SpectrumParamValues\"", "name=\"OtherParams\"")
        @test_throws FileFormatError load(MassHunterMS(dir))
    end
    mktempdir() do dir
        acq = synthetic_masshunter(dir; unbounded=true)
        change_schema(acq, "name=\"MinX\"", "name=\"OtherMinimum\"")
        @test_throws FileFormatError load(MassHunterMS(dir))
    end
end
