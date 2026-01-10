module MassHunterMSLoader

# ─────────────────────────────────────────────────────────────────────────────
# Imports, Dependencies, and Public Interface
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

import FileIO: load
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
    MassHunterMSLoaderSpec{F}(path::String; mode::Symbol=:ms)

Internal structure representing a configured load request for Agilent MassHunter MS data.
This type encapsulates:
- the source path,
- selected mode (`:ms` or `:tic`), and
- a type tag for format versioning.

Normally, this is constructed indirectly via `MassHunterMS(...)`.
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
    MassHunterMS(path::String; mode::Symbol=:ms) -> MassHunterMSLoaderSpec

Construct a loader specification for Agilent MassHunter GC/MS data. The `path` can be
either a directory (containing `data.ms`) or a direct path to the file. 

The keyword argument `mode` determines the type of data to load:
- `:ms`  → full mass spectral data (default)
- `:tic` → total ion chromatogram only

Returns a `MassHunterMSLoaderSpec` object used for deferred data loading via `load(...)`.
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
Takes a `MassHunterMSLoaderSpec` created via `MassHunterMS(...)`.

Returns a `AbstractScanSeries` subtype object containing either:
- a vector of `MassScan` objects (for `mode=:ms`), or
- a vector of `ChromScan` objects (for `mode=:tic`),

along with parsed metadata (sample, user, acquisition, instrument).
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
    return try
        readfile(T, path, req.options)
    catch e
        throw(FileCorruptionError("Failed to read MS data at $path: $e"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Functions for extracting data from Agilent ChemStation MS binary files 
# (e.g., 'data.ms')
# ─────────────────────────────────────────────────────────────────────────────

function magicnumber(file::AbstractString)
    open(file, "r") do io
        return ltoh(read(io, UInt16))
    end
end

struct MSScanBin end
struct MSScanBinV1 end
struct MSPeakBin end
struct MSPeakBinV1 end

struct MSScanDataV1
    ScanID::Int32
    ScanMethodID::Int32
    TimeSegmentID::Int32
    ScanTime::Float64
    MSLevel::Int32
    ScanType::Int32
    TIC::Float64
    BasePeakMZ::Float64
    BasePeakValue::Float64
    CycleNumber::Int32
    Status::Int32
    IonMode::Int32
    IonPolarity::Int32
    Fragmentor::Float64
    CollisionEnergy::Float64
    MzOfInterest::Float64
    SamplingPeriod::Float64
    MeasuredMassRangeMin::Float64
    MeasuredMassRangeMax::Float64
    Threshold::Float64
    IsFragmentorDynamic::Int32
    IsCollisionEnergyDynamic::Int32
    SpectrumFormatID::Int32
    SpectrumOffset::Int64
    ByteCount::Int32
    PointCount::Int32
    MinY::Float64
    MaxY::Float64
    MinX::Float64
    MaxX::Float64
end

const MSSCAN_FIELDNAMES = fieldnames(MSScanDataV1)
const MSSCAN_ATTR_EXCLUDES = Set((:ScanTime, :MSLevel, :ByteCount, :PointCount))
const MSSCAN_ATTR_EXCLUDES_TIC = union(MSSCAN_ATTR_EXCLUDES, Set((:TIC,)))
const MSSCAN_ATTR_FIELDS = Tuple(filter(name -> name ∉ MSSCAN_ATTR_EXCLUDES, MSSCAN_FIELDNAMES))
const MSSCAN_ATTR_FIELDS_TIC = Tuple(filter(name -> name ∉ MSSCAN_ATTR_EXCLUDES_TIC, MSSCAN_FIELDNAMES))

@inline function scandatum_attrs(scandatum::MSScanDataV1; include_tic::Bool=true)
    fields = include_tic ? MSSCAN_ATTR_FIELDS : MSSCAN_ATTR_FIELDS_TIC
    values = ntuple(i -> getfield(scandatum, fields[i]), length(fields))
    NamedTuple{fields}(values)
end

for type in (MSScanDataV1, )
    expr = Expr(:call, type, map(ft -> Expr(:call, :ltoh, Expr(:call, :read, :io, ft)), 
        fieldtypes(type))...)
    @eval Base.read(io::IO, ::Type{$type}) = $expr
end

function msscandata(::MSScanBinV1, file::AbstractString)
    size = sum(sizeof(fieldtype) for fieldtype in fieldtypes(MSScanDataV1))
    open(file, "r") do io
        seekend(io)
        endpos = position(io)
        seek(io, 88)
        startpos = ltoh(read(io, UInt16))
        seek(io, startpos)
        (endpos - startpos) % size == 0 || throw(
            FileCorruptionError("MSScan.bin has unexpected data content"))
        scancount = div(endpos - startpos, size)
        data = Vector{MSScanDataV1}(undef, scancount)
        for i in 1:scancount
            data[i] = read(io, MSScanDataV1)
        end
        data
    end
end

function read_scan_data_tic(
    ::Type{T},
    scandata::Vector{MSScanDataV1}, 
    options::MassHunterMSOptions
    ) where {T<:MassHunterMSv1}

    retentions_unitfree = Vector{Float64}()
    intensities_unitfree = Vector{Float64}()
    attrs = NamedTuple[]
    for scandatum in scandata
        (options.level === nothing || scandatum.MSLevel == options.level) || continue
        push!(retentions_unitfree, scandatum.ScanTime)
        push!(intensities_unitfree, scandatum.TIC)
        push!(attrs, scandatum_attrs(scandatum; include_tic=false))
    end

    tretention = options.scantimetype
    tintensity = options.intensitytype
    retentions_converted = tretention === Nothing ? retentions_unitfree :
        convert(Vector{tretention}, retentions_unitfree)
    intensities_converted = tintensity === Nothing ? intensities_unitfree :
        convert(Vector{tintensity}, intensities_unitfree)

    retentions_converted .* 60, intensities_converted, attrs
end

scantimes(scandata::Vector{MSScanDataV1}) = [s.ScanTime for s in scandata]

function read_scan_data_ms(
    ::Type{T},
    path::AbstractString,
    scandata::Vector{MSScanDataV1}, 
    options::MassHunterMSOptions
    ) where {T<:MassHunterMSv1}
    
    mspeak_file = joinpath(path, "MSPeak.bin")
    isfile(mspeak_file) || throw(
        MissingFileError("Cannot find MSPeak.bin file at: $mspeak_file"))
    signature = magicnumber(mspeak_file)
    signature == 259 || throw(
        FileFormatError("\"$mspeak_file\" has unexpected file signature: \"$signature\""))

    scans = mspeakdata(MSPeakBinV1(), mspeak_file, scandata, options)
    
    MassScanSeries(
        scans,
        instrument=NamedTuple(), 
        acquisition=NamedTuple(),
        user=NamedTuple(),
        sample=NamedTuple()
    )
end

function mspeakdata(
    ::MSPeakBinV1,
    file::AbstractString, 
    scandata::Vector{MSScanDataV1},
    options::MassHunterMSOptions,
    )
    
    tretention = options.scantimetype
    tmz = options.iontype
    tintensity = options.intensitytype

    retentions = tretention === Nothing ? scantimes(scandata) :
        convert(Vector{tretention}, scantimes(scandata))
    retentions_sec = retentions .* 60
    retentionunit = u"s"
    mzunit = nothing
    intensityunit = nothing
    scans = Vector{Any}()
    open(file, "r") do f
        for (i, scandatum) in enumerate(scandata)

            (options.level === nothing || scandatum.MSLevel == options.level) || continue

            seek(f, scandatum.SpectrumOffset)
            storagetype = Float32
            if scandatum.PointCount > 0
                bpp = scandatum.ByteCount / scandatum.PointCount
                bpp == 8 || bpp == 16 || throw(
                    FileDecodingError("\"$file\" has unexpected ByteCount per point: \"$bpp\""))
                storagetype = bpp == 8 ? Float32 : Float64
            end
            mzs = ltoh.(read!(f, Vector{storagetype}(undef, scandatum.PointCount)))
            ints = ltoh.(read!(f, Vector{storagetype}(undef, scandatum.PointCount)))
            mzs_converted = tmz === Nothing ? mzs : convert(Vector{tmz}, mzs)
            ints_converted = tintensity === Nothing ? ints : convert(Vector{tintensity}, ints)
            attrs = scandatum_attrs(scandatum; include_tic=true)
            scan = MassScan(
                retentions_sec[i], retentionunit,
                mzs_converted, mzunit, 
                ints_converted, intensityunit;
                level=scandatum.MSLevel,
                attrs=attrs)
            push!(scans, scan)
        end
    end
    level_desc = options.level === nothing ? "any level" :
        "level $(options.level)"
    isempty(scans) && throw(
        FileCorruptionError("MSPeak.bin did not contain any scans at $level_desc"))
    S = typeof(first(scans))
    Vector{S}(scans)
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
    try
        msscan_file = joinpath(path, "MSScan.bin")
        signature = magicnumber(msscan_file)
        signature == 257 || throw(
            FileFormatError("\"$msscan_file\" has unexpected file signature: \"$signature\""))
        scandata = msscandata(MSScanBinV1(), msscan_file)
        if mode == :ms
            return read_scan_data_ms(T, path, scandata, options)
        elseif mode == :tic
            retentions_unitfree, intensities_unitfree, scan_attrs = 
                read_scan_data_tic(T, scandata, options)
            retentionunit = u"s"
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
    catch e
        throw(FileCorruptionError("Failed to open/read MassHunterMS data \"$path\": $e"))
    end
end

end  # module
