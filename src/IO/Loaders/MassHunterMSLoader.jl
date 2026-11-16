module MassHunterMSLoader

# ─────────────────────────────────────────────────────────────────────────────
# Imports, Dependencies, and Public Interface
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

import FileIO: load
import EzXML
using Unitful

# ── Internal Project Imports ─────────────────────────────────────────────────
import ..JuChrom: ChromScan, MassScan, ChromScanSeries, MassScanSeries
import ..InputOutput: FileCorruptionError, FileDecodingError, FileFormatError, 
    FileIOError, MissingFileError, MissingFolderError, UnexpectedEOFError

# ── Public API ───────────────────────────────────────────────────────────────

export MassHunterMS, MassHunterMSv1, MassHunterMSLoaderSpec, load

# ─────────────────────────────────────────────────────────────────────────────
# Loader Specification and Entry Points for Agilent MassHunter MS Files
# ─────────────────────────────────────────────────────────────────────────────

# ── Constants and Format Tags ────────────────────────────────────────────────

const SUPPORTED_MODES = (:ms, :tic)

abstract type MassHunterMSFormat end
struct MassHunterMSv1 <: MassHunterMSFormat end

# ── Loader Configuration Structures ──────────────────────────────────────────

struct MassHunterMSOptions{
    T1<:Union{AbstractFloat, Nothing}, 
    T2<:Union{AbstractFloat, Nothing}, 
    T3<:Union{AbstractFloat, Nothing}}

    mode::Symbol  # :ms or :tic
    level::Union{Int, Nothing}  # MSLevel (Nothing = all)
    scantimetype::Type{T1}
    iontype::Type{T2}
    intensitytype::Type{T3}
end

struct MassHunterMSLoaderSpec{
    F<:MassHunterMSFormat}
    path::String
    options::MassHunterMSOptions
    _format::Val{F}
end

# ── Loader Spec Constructors ─────────────────────────────────────────────────

"""
    MassHunterMSLoaderSpec{F}(path::String; mode::Symbol=:ms, level=nothing,
                              scantimetype=Nothing, iontype=Nothing, intensitytype=Nothing)

Internal structure representing a configured load request for Agilent MassHunter MS data.
The `path` must point to a `.D` directory containing an `AcqData` folder with
`MSScan.bin` and `MSScan.xsd`. `mode=:ms` additionally requires `MSPeak.bin`.
In either mode, `MSTS.xml` is required when `MSScan.xsd` declares `SpectrumParamValues`
within `ScanRecordType` with `maxOccurs="unbounded"`: its scan count is used to verify
one spectrum block per scan. For other supported schemas, `MSTS.xml` is optional;
its scan count is checked if present.

The `mode` keyword selects full mass spectral data (`:ms`) or the total ion
chromatogram only (`:tic`). The `level` keyword filters by MS level, while `scantimetype`, 
`iontype`, and `intensitytype` allow conversion of scan times (minutes), m/z values, and 
intensities, respectively.

Type specification allows conversion from storage types that are wider than the effective
precision. Some MassHunter files store values in a wider floating-point type than their
effective precision (for example, Float64 storage of Float32-precision data). Converting
to an appropriate type can remove misleading numerical artifacts and reduce memory use.

Normally, this is constructed indirectly via `MassHunterMS(...)`.

See also [`MassHunterMS`](@ref MassHunterMS).
"""
function MassHunterMSLoaderSpec{F}(
    path::String;
    mode::Symbol=:ms,
    level::Union{Int, Nothing}=nothing,
    scantimetype::Type{T1}=Nothing,
    iontype::Type{T2}=Nothing,
    intensitytype::Type{T3}=Nothing
    ) where {
        F,
        T1<:Union{AbstractFloat, Nothing},
        T2<:Union{AbstractFloat, Nothing}, 
        T3<:Union{AbstractFloat, Nothing}
    }
    
    mode ∈ SUPPORTED_MODES || throw(ArgumentError(
        "Unsupported mode: $mode. Supported modes: $(SUPPORTED_MODES)"))
    MassHunterMSLoaderSpec{F}(path, MassHunterMSOptions(mode, level,
        scantimetype, iontype, intensitytype), Val{F}())
end

"""
    MassHunterMS(
        path::String;
        mode::Symbol=:ms,
        level=nothing,
        scantimetype=Nothing,
        iontype=Nothing,
        intensitytype=Nothing
    ) -> MassHunterMSLoaderSpec

Construct a loader specification for Agilent MassHunter GC/MS data.
The `path` must point to a `.D` directory containing an `AcqData` folder with
`MSScan.bin` and `MSScan.xsd`. `mode=:ms` additionally requires `MSPeak.bin`.
In either mode, `MSTS.xml` is required when `MSScan.xsd` declares `SpectrumParamValues`
within `ScanRecordType` with `maxOccurs="unbounded"`: its scan count is used to verify
one spectrum block per scan. For other supported schemas, `MSTS.xml` is optional;
its scan count is checked if present.

The `mode` keyword selects full mass spectral data (`:ms`) or the total ion
chromatogram only (`:tic`). The `level` keyword filters by MS level, while
`scantimetype`, `iontype`, and `intensitytype` allow conversion of scan times (minutes), 
m/z values, and intensities, respectively.

Type specification allows conversion from storage types that are wider than the effective
precision. Some MassHunter files store values in a wider floating-point type than their
effective precision (for example, Float64 storage of Float32-precision data). Converting
to an appropriate type can remove misleading numerical artifacts and reduce memory use.

Returns a [`MassHunterMSLoaderSpec`](@ref MassHunterMSLoaderSpec) object used for deferred 
data loading via `load(...)`.

See also [`MassHunterMSLoaderSpec`](@ref MassHunterMSLoaderSpec).
"""
function MassHunterMS(
    path::String;
    mode::Symbol=:ms,
    level::Union{Int, Nothing}=nothing,
    scantimetype::Type{T1}=Nothing,
    iontype::Type{T2}=Nothing,
    intensitytype::Type{T3}=Nothing
    ) where {T1<:Union{AbstractFloat, Nothing}, T2<:Union{AbstractFloat, Nothing}, 
        T3<:Union{AbstractFloat, Nothing}}

    mode ∈ SUPPORTED_MODES || throw(ArgumentError(
        "Unsupported mode: $mode. Supported modes: $(SUPPORTED_MODES)"))
    options = MassHunterMSOptions(mode, level, scantimetype, iontype, intensitytype)
    MassHunterMSLoaderSpec{MassHunterMSv1}(path, options, Val{MassHunterMSv1}())
end

# ── Data Loading Interface ───────────────────────────────────────────────────

"""
    load(req::MassHunterMSLoaderSpec{MassHunterMSv1}) -> AbstractScanSeries

Loads and parses mass spectrometry data from Agilent MassHunter (version 1) GC/MS files.
The request must be a `MassHunterMSLoaderSpec`, typically created via `MassHunterMS(...)`.
Returns an `AbstractScanSeries` subtype containing either a vector of 
[`MassScan`](@ref JuChrom.MassScan) objects for `mode=:ms` or a vector of 
[`ChromScan`](@ref JuChrom.ChromScan) objects for `mode=:tic`. Scan metadata is retained
by field name; sample, user, acquisition, and instrument metadata are currently empty.

See also
[`MassHunterMS`](@ref MassHunterMS),
[`MassHunterMSLoaderSpec`](@ref MassHunterMSLoaderSpec).
"""
function load(req::MassHunterMSLoaderSpec{T}) where {T<:MassHunterMSv1}
    dfolder = req.path
    if isfile(dfolder)
        throw(ArgumentError(
            "Expected path to an Agilent .D folder, but received path to a file: $dfolder"))
    elseif isdir(dfolder)
        path = joinpath(dfolder, "AcqData")
        isdir(path) || throw(MissingFolderError("Cannot find AcqData folder: $path"))
    else
        throw(MissingFolderError(
            "Path does not exist or is not a valid directory: $dfolder"))
    end
    readfile(T, path, req.options)
end

# ─────────────────────────────────────────────────────────────────────────────
# Schema and binary readers for Agilent MassHunter MS files
# ─────────────────────────────────────────────────────────────────────────────

function magicnumber(file::AbstractString)
    isfile(file) || throw(MissingFileError(file))
    filesize(file) >= 2 || throw(UnexpectedEOFError("Truncated file signature in $file"))
    open(file, "r") do io
        return readle(io, UInt16)
    end
end

struct MSScanBin end
struct MSScanBinV1 end
struct MSPeakBin end
struct MSPeakBinV1 end

const MSSCAN_ATTR_EXCLUDES = Set((:ScanTime, :MSLevel, :ByteCount, :PointCount))
const MSSCAN_ATTR_EXCLUDES_TIC = union(MSSCAN_ATTR_EXCLUDES, Set((:TIC,)))

# MSScan.xsd describes packed binary fields, not Julia's in-memory struct layout.
const XSD_NAMESPACE = "http://www.w3.org/2001/XMLSchema"
const XSD_PRIMITIVES = Dict{String,DataType}(
    "byte" => Int8, "unsignedByte" => UInt8,
    "short" => Int16, "unsignedShort" => UInt16,
    "int" => Int32, "unsignedInt" => UInt32,
    "long" => Int64, "unsignedLong" => UInt64,
    "float" => Float32, "double" => Float64)

struct MSScanLayout
    recordtype::DataType
    recordsize::Int
    single_spectrum_repetition::Bool
end

schemaerror(message) = throw(FileFormatError("Unsupported MSScan.xsd: $message"))
xmlattr(node, name, default="") = haskey(node, name) ? node[name] : default
isxsd(node, name) = EzXML.hasnamespace(node) &&
    EzXML.namespace(node) == XSD_NAMESPACE && EzXML.nodename(node) == name

function read_metadata_xml(file)
    isfile(file) || throw(MissingFileError(file))
    doc = try
        EzXML.readxml(file)
    catch e
        e isa EzXML.XMLError || rethrow()
        throw(FileFormatError("Malformed XML in $file: $(sprint(showerror, e))"))
    end
    EzXML.hasdtd(doc) && throw(FileFormatError("DTD declarations are unsupported in $file"))
    doc
end

# Comments and annotations describe the format but contribute no binary fields.
function schema_children(node)
    children = EzXML.Node[]
    for child in EzXML.eachnode(node)
        EzXML.iscomment(child) && continue
        EzXML.istext(child) && isempty(strip(EzXML.nodecontent(child))) && continue
        EzXML.iselement(child) || schemaerror("unexpected content at $(EzXML.nodepath(node))")
        EzXML.hasnamespace(child) && EzXML.namespace(child) == XSD_NAMESPACE ||
            schemaerror("foreign namespace at $(EzXML.nodepath(child))")
        isxsd(child, "annotation") || push!(children, child)
    end
    children
end

function schema_attributes(node, allowed)
    for attr in EzXML.eachattribute(node)
        !EzXML.hasnamespace(attr) && EzXML.nodename(attr) in allowed ||
            schemaerror("attribute $(EzXML.nodename(attr)) at $(EzXML.nodepath(node))")
    end
end

function schema_qname(node, value)
    parts = split(strip(value), ':')
    length(parts) in (1, 2) && all(!isempty, parts) || schemaerror("invalid QName '$value'")
    prefix, name = length(parts) == 1 ? ("", only(parts)) : (parts[1], parts[2])
    namespaces = Dict(EzXML.namespaces(node))
    isempty(prefix) || haskey(namespaces, prefix) || schemaerror("unbound prefix in '$value'")
    get(namespaces, prefix, ""), String(name)
end

function occurrence(node, attr, default)
    value = xmlattr(node, attr, default)
    count = tryparse(Int, value)
    count !== nothing && 0 <= count <= 1024 ||
        schemaerror("$attr='$value' at $(EzXML.nodepath(node)); expected a bounded count ≤ 1024")
    count
end

function scanlayout(file::AbstractString)
    doc = read_metadata_xml(file)
    root = EzXML.root(doc)
    isxsd(root, "schema") || schemaerror("root must be an XML Schema schema element")
    schema_attributes(root, ("targetNamespace", "elementFormDefault", "attributeFormDefault",
        "version", "id"))
    target = xmlattr(root, "targetNamespace")
    types = Dict{String,EzXML.Node}()
    for child in schema_children(root)
        if isxsd(child, "complexType")
            name = xmlattr(child, "name")
            isempty(name) && schemaerror("unnamed top-level complexType")
            haskey(types, name) && schemaerror("duplicate complexType '$name'")
            types[name] = child
        elseif !isxsd(child, "element")
            schemaerror("top-level $(EzXML.nodename(child)); imports/includes and custom simple types are unsupported")
        end
    end
    haskey(types, "ScanRecordType") || schemaerror("missing ScanRecordType")
    fields = Pair{Symbol,Type}[]
    repeated = Ref(false)
    expand_scan_type!(fields, types["ScanRecordType"], types, target, String[], repeated)
    isempty(fields) && schemaerror("ScanRecordType is empty")
    names = Tuple(first.(fields))
    length(unique(names)) == length(names) || schemaerror("duplicate field names after expanding nested types")
    recordtype = NamedTuple{names,Tuple{last.(fields)...}}
    recordsize = sum(packedsize, last.(fields))
    MSScanLayout(recordtype, recordsize, repeated[])
end

function expand_scan_type!(fields, node, types, target, stack, repeated)
    name = xmlattr(node, "name")
    name in stack && schemaerror("recursive complexType '$name'")
    length(stack) < 32 || schemaerror("complex types nested more than 32 levels")
    schema_attributes(node, ("name", "id"))
    children = schema_children(node)
    length(children) == 1 && isxsd(only(children), "sequence") ||
        schemaerror("complexType '$name' must contain one sequence (no choices, extensions or attributes)")
    seq = only(children)
    schema_attributes(seq, ("id", "minOccurs", "maxOccurs"))
    occurrence(seq, "minOccurs", "1") == occurrence(seq, "maxOccurs", "1") == 1 ||
        schemaerror("repeated or optional sequences in '$name'")
    nextstack = [stack; name]
    for elem in schema_children(seq)
        isxsd(elem, "element") || schemaerror("$(EzXML.nodename(elem)) in sequence '$name'")
        schema_attributes(elem, ("name", "type", "minOccurs", "maxOccurs", "id"))
        field = xmlattr(elem, "name")
        isempty(field) && schemaerror("element without a name in '$name'")
        isempty(schema_children(elem)) || schemaerror("inline types for '$field'; use a named complex type")
        uri, typename = schema_qname(elem, xmlattr(elem, "type"))
        minimum = occurrence(elem, "minOccurs", "1")
        unbounded = xmlattr(elem, "maxOccurs", "1") == "unbounded"
        maximum = unbounded ? 1 : occurrence(elem, "maxOccurs", "1")
        if uri == XSD_NAMESPACE
            unbounded && schemaerror("unbounded primitive field '$field'")
            haskey(XSD_PRIMITIVES, typename) || schemaerror("primitive type '$typename' for '$field'")
            # In the validated MassHunter layouts, active minOccurs=0 scalars occupy
            # fixed-width slots. This is a binary convention, not a general XSD rule.
            (maximum == 1 && minimum in (0, 1)) || (minimum == maximum && maximum > 1) ||
                schemaerror("variable or zero repetition of '$field'")
            T = XSD_PRIMITIVES[typename]
            push!(fields, Symbol(field) => (maximum == 1 ? T : NTuple{maximum,T}))
        else
            uri == target && haskey(types, typename) || schemaerror("unknown complex type '$typename' in namespace '$uri'")
            minimum == maximum == 1 || schemaerror("optional/repeated complex field '$field'")
            nested = Pair{Symbol,Type}[]
            if unbounded
                # MSTS.xml's independent scan count must match the index size computed
                # with one SpectrumParamsType per scan; divisibility alone is insufficient.
                name == "ScanRecordType" && field == "SpectrumParamValues" &&
                    typename == "SpectrumParamsType" && !repeated[] ||
                    schemaerror("unbounded complex field '$field'")
                repeated[] = true
            end
            expand_scan_type!(nested, types[typename], types, target, nextstack, repeated)
            if unbounded
                expected = [:SpectrumFormatID => Int32, :SpectrumOffset => Int64,
                    :ByteCount => Int32, :PointCount => Int32, :MinX => Float64,
                    :MaxX => Float64, :MinY => Float64, :MaxY => Float64]
                # The schema determines storage order; the same named fields may be reordered.
                length(nested) == length(expected) && Set(nested) == Set(expected) ||
                    schemaerror("unbounded SpectrumParamsType has an unverified layout")
            end
            append!(fields, nested)
        end
        length(fields) <= 256 || schemaerror("more than 256 expanded fields")
    end
end

packedsize(::Type{T}) where {T<:Number} = sizeof(T)
packedsize(::Type{NTuple{N,T}}) where {N,T} = N * packedsize(T)
readle(io, ::Type{T}) where {T<:Integer} = ltoh(read(io, T))
readle(io, ::Type{Float32}) = reinterpret(Float32, readle(io, UInt32))
readle(io, ::Type{Float64}) = reinterpret(Float64, readle(io, UInt64))
readle(io, ::Type{NTuple{N,T}}) where {N,T} = ntuple(_ -> readle(io, T), Val(N))

function read_scan_records(io, ::Type{T}, count) where {T}
    records = Vector{T}(undef, count)
    for i in eachindex(records)
        values = ntuple(j -> readle(io, fieldtype(T, j)), Val(fieldcount(T)))
        records[i] = T(values)
    end
    records
end

function metadata_scan_count(path; required=false)
    file = joinpath(path, "MSTS.xml")
    if !isfile(file)
        required && schemaerror("unbounded SpectrumParamValues requires MSTS.xml scan counts to verify a single spectrum per record")
        return nothing
    end
    doc = read_metadata_xml(file)
    root = EzXML.root(doc)
    EzXML.nodename(root) == "TimeSegments" ||
        throw(FileFormatError("Expected TimeSegments in $file"))
    ns = root.namespace
    total = 0
    segments = 0
    for segment in EzXML.eachelement(root)
        EzXML.nodename(segment) == "TimeSegment" && segment.namespace == ns || continue
        counts = [n for n in EzXML.eachelement(segment)
            if EzXML.nodename(n) == "NumOfScans" && n.namespace == ns]
        length(counts) == 1 || throw(FileFormatError("Expected one NumOfScans per TimeSegment in $file"))
        count = tryparse(Int, strip(EzXML.nodecontent(only(counts))))
        count !== nothing && 0 <= count <= typemax(Int) - total ||
            throw(FileCorruptionError("Invalid NumOfScans in $file"))
        total += count
        segments += 1
    end
    segments > 0 || throw(FileFormatError("No TimeSegment scan counts in $file"))
    total
end

function require_scan_fields(::Type{T}, fields) where {T}
    for (name, kind) in fields
        name in fieldnames(T) && fieldtype(T, name) <: kind ||
            throw(FileFormatError("MSScan layout requires scalar $name of type $kind"))
    end
end

function msscandata(::MSScanBinV1, file::AbstractString)
    magicnumber(file) == 257 || throw(FileFormatError("$file has an unsupported MSScan signature (expected 257)"))
    path = dirname(file)
    schemafile = joinpath(path, "MSScan.xsd")
    layout = scanlayout(schemafile)
    require_scan_fields(layout.recordtype, (:ScanTime => AbstractFloat, :MSLevel => Integer))
    expectedcount = metadata_scan_count(path; required=layout.single_spectrum_repetition)
    open(file, "r") do io
        endpos = filesize(file)
        endpos >= 92 || throw(UnexpectedEOFError("Truncated MSScan.bin header in $file (need at least 92 bytes)"))
        seek(io, 88)
        startpos = Int64(readle(io, UInt32))
        startpos >= 92 || throw(FileCorruptionError("MSScan.bin first-record offset $startpos overlaps the header"))
        startpos <= endpos || throw(UnexpectedEOFError("MSScan.bin first-record offset $startpos exceeds file size $endpos"))
        payload = endpos - startpos
        size = layout.recordsize
        size > 0 || schemaerror("zero-sized ScanRecordType")
        if layout.single_spectrum_repetition && Int128(expectedcount) * size < payload
            schemaerror("payload exceeds one SpectrumParamsType per scan according to MSTS.xml; multiple spectra/variable records are unsupported")
        end
        if payload % size != 0
            throw(UnexpectedEOFError("MSScan.bin payload ($payload bytes) is not a whole number of schema records ($size bytes); truncated record or inconsistent schema"))
        end
        count = div(payload, size)
        if expectedcount !== nothing && count != expectedcount
            throw(FileCorruptionError("MSScan.bin contains $count schema records but MSTS.xml declares $expectedcount scans"))
        end
        seek(io, startpos)
        records = try
            read_scan_records(io, layout.recordtype, count)
        catch e
            e isa EOFError || rethrow()
            throw(UnexpectedEOFError("Truncated scan record in $file"))
        end
        position(io) == endpos || throw(FileCorruptionError("MSScan.bin record boundary mismatch"))
        for (i, s) in enumerate(records)
            isfinite(s.ScanTime) && s.MSLevel >= 1 ||
                throw(FileCorruptionError("Invalid ScanTime or MSLevel at MSScan.bin record $i"))
        end
        records
    end
end

@inline function scandatum_attrs(scandatum::NamedTuple; include_tic::Bool=true)
    excludes = include_tic ? MSSCAN_ATTR_EXCLUDES : MSSCAN_ATTR_EXCLUDES_TIC
    fields = filter(name -> name ∉ excludes, keys(scandatum))
    NamedTuple{fields}(map(name -> getproperty(scandatum, name), fields))
end

function read_scan_data_tic(
    ::Type{T},
    scandata::AbstractVector,
    options::MassHunterMSOptions
    ) where {T<:MassHunterMSv1}

    require_scan_fields(eltype(scandata), (:TIC => AbstractFloat,))
    retentions_unitfree = Vector{Float64}()
    intensities_unitfree = Vector{Float64}()
    attrs = NamedTuple[]
    for scandatum in scandata
        (options.level ≡ nothing || scandatum.MSLevel == options.level) || continue
        push!(retentions_unitfree, scandatum.ScanTime)
        isfinite(scandatum.TIC) || throw(FileCorruptionError("Non-finite TIC in MSScan.bin"))
        push!(intensities_unitfree, scandatum.TIC)
        push!(attrs, scandatum_attrs(scandatum; include_tic=false))
    end

    tretention = options.scantimetype
    tintensity = options.intensitytype
    retentions_converted = tretention ≡ Nothing ? retentions_unitfree :
        convert(Vector{tretention}, retentions_unitfree)
    intensities_converted = tintensity ≡ Nothing ? intensities_unitfree :
        convert(Vector{tintensity}, intensities_unitfree)

    retentions_converted, intensities_converted, attrs
end

scantimes(scandata::AbstractVector) = [s.ScanTime for s in scandata]

function read_scan_data_ms(
    ::Type{T},
    path::AbstractString,
    scandata::AbstractVector,
    options::MassHunterMSOptions
    ) where {T<:MassHunterMSv1}
    
    mspeak_file = joinpath(path, "MSPeak.bin")
    scans = mspeakdata(MSPeakBinV1(), mspeak_file, scandata, options)
    
    MassScanSeries(
        scans,
        instrument=NamedTuple(), 
        acquisition=NamedTuple(),
        user=NamedTuple(),
        sample=NamedTuple()
    )
end

# Index parsing deliberately does not inspect spectral encodings. TIC-only loads
# can use schemas referencing encodings that the MSPeak reader cannot decode.
function validate_spectrum_bounds(s, endpos, i)
    s.PointCount >= 0 && s.ByteCount >= 0 ||
        throw(FileCorruptionError("Negative PointCount or ByteCount at scan record $i"))
    s.SpectrumOffset >= 68 ||
        throw(FileCorruptionError("SpectrumOffset $(s.SpectrumOffset) at scan record $i overlaps the MSPeak.bin header"))
    s.SpectrumOffset <= endpos && s.ByteCount <= endpos - s.SpectrumOffset ||
        throw(UnexpectedEOFError("Spectrum block at scan record $i exceeds MSPeak.bin size $endpos"))
    (s.PointCount == 0) == (s.ByteCount == 0) ||
        throw(FileCorruptionError("Inconsistent zero PointCount/ByteCount at scan record $i"))
end

function spectrum_storage_type(s, i)
    s.SpectrumFormatID == 1 || throw(FileDecodingError(
        "Unsupported SpectrumFormatID $(s.SpectrumFormatID) at scan record $i; only MSPeak format 1 paired floating-point arrays are supported (try mode=:tic)"))
    s.PointCount > 0 || throw(FileDecodingError(
        "Empty spectrum at scan record $i cannot be represented as a MassScan (try mode=:tic)"))
    s.ByteCount % s.PointCount == 0 || throw(FileDecodingError(
        "Unsupported spectrum packing at scan record $i: ByteCount=$(s.ByteCount), PointCount=$(s.PointCount)"))
    bpp = div(s.ByteCount, s.PointCount)
    bpp in (8, 16) || throw(FileDecodingError(
        "Unsupported spectrum packing at scan record $i: $bpp bytes per point; expected 8 or 16"))
    bpp == 8 ? Float32 : Float64
end

function read_spectrum_values(io, ::Type{T}, count) where {T}
    U = T == Float32 ? UInt32 : UInt64
    bits = read!(io, Vector{U}(undef, count))
    map!(ltoh, bits, bits)
    collect(reinterpret(T, bits))
end

function mspeakdata(
    ::MSPeakBinV1,
    file::AbstractString, 
    scandata::AbstractVector,
    options::MassHunterMSOptions,
    )
    
    require_scan_fields(eltype(scandata), (:SpectrumFormatID => Integer,
        :SpectrumOffset => Integer, :ByteCount => Integer, :PointCount => Integer))
    magicnumber(file) == 259 || throw(FileFormatError("Unsupported MSPeak.bin signature in $file (expected 259)"))
    endpos = filesize(file)
    endpos >= 68 || throw(UnexpectedEOFError("Truncated MSPeak.bin header in $file"))
    for (i, s) in enumerate(scandata)
        validate_spectrum_bounds(s, endpos, i)
    end
    selected = findall(s -> options.level ≡ nothing || s.MSLevel == options.level, scandata)
    level_desc = options.level ≡ nothing ? "any level" : "level $(options.level)"
    isempty(selected) && throw(ArgumentError("No MassHunter scans at $level_desc"))
    storagetypes = [spectrum_storage_type(scandata[i], i) for i in selected]
    # A series requires one concrete scan type; promote mixed storage widths losslessly.
    defaulttype = foldl(promote_type, storagetypes)
    tretention = options.scantimetype
    tmz = options.iontype ≡ Nothing ? defaulttype : options.iontype
    tintensity = options.intensitytype ≡ Nothing ? defaulttype : options.intensitytype
    retentions = tretention ≡ Nothing ? scantimes(scandata) :
        convert(Vector{tretention}, scantimes(scandata))
    scans = MassScan[]
    open(file, "r") do io
        for (i, storagetype) in zip(selected, storagetypes)
            s = scandata[i]
            seek(io, s.SpectrumOffset)
            mzs, ints = try
                (read_spectrum_values(io, storagetype, s.PointCount),
                    read_spectrum_values(io, storagetype, s.PointCount))
            catch e
                e isa EOFError || rethrow()
                throw(UnexpectedEOFError("Truncated spectrum at scan record $i in $file"))
            end
            all(isfinite, mzs) && all(>(0), mzs) && all(>(0), diff(mzs)) &&
                all(isfinite, ints) || throw(FileCorruptionError(
                    "Invalid m/z or intensity values at scan record $i in $file"))
            mzs_converted = convert(Vector{tmz}, mzs)
            ints_converted = convert(Vector{tintensity}, ints)
            push!(scans, MassScan(retentions[i], u"minute", mzs_converted, nothing,
                ints_converted, nothing; level=s.MSLevel, attrs=scandatum_attrs(s)))
        end
    end
    Vector{typeof(first(scans))}(scans)
end

# ─── File Reader Interface ───────────────────────────────────────────────────

function readfile(
    ::Type{T},
    path::AbstractString, 
    options::MassHunterMSOptions
    ) where {T<:MassHunterMSv1}

    mode = options.mode
    mode ∈ (:ms, :tic) || throw(ArgumentError(
        "Unsupported mode: $mode. Supported modes: (:ms, :tic)"))
    msscan_file = joinpath(path, "MSScan.bin")
    scandata = msscandata(MSScanBinV1(), msscan_file)
    if mode ≡ :ms
        return read_scan_data_ms(T, path, scandata, options)
    elseif mode ≡ :tic
        retentions_unitfree, intensities_unitfree, scan_attrs =
            read_scan_data_tic(T, scandata, options)
        retentionunit = u"minute"
        intensityunit = nothing

        scans = [
            ChromScan(
                retention, retentionunit, intensity, intensityunit;
                attrs=attrs
            ) for (retention, intensity, attrs) in
                zip(retentions_unitfree, intensities_unitfree, scan_attrs)
        ]
        return ChromScanSeries(
            scans,
            instrument=NamedTuple(),
            acquisition=NamedTuple(),
            user=NamedTuple(),
            sample=NamedTuple()
        )
    end
end

end  # module
