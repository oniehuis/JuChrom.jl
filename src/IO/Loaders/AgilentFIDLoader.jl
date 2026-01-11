module AgilentFIDLoader

# ─────────────────────────────────────────────────────────────────────────────
# Imports, Dependencies, and Public Interface
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

import FileIO: load
using StringEncodings
using Unitful

# ── Internal Project Imports ─────────────────────────────────────────────────
import ..JuChrom: ChromScan, ChromScanSeries
import ..InputOutput: FileCorruptionError, FileDecodingError, FileFormatError, 
    FileIOError, MissingFileError, UnexpectedEOFError

# ── Public API ───────────────────────────────────────────────────────────────

export AgilentFID, AgilentFIDLoaderSpec, load

# ─────────────────────────────────────────────────────────────────────────────
# Loader Specification and Entry Points for Agilent FID Files
# ─────────────────────────────────────────────────────────────────────────────

# ── Constants and Format Tags ────────────────────────────────────────────────

# const AFIDSample01 = artifact"AFIDSample01"
const ENCODING_1_byte = "Windows-1252"
const ENCODING_2_bytes = "UTF-16LE"
const SUPPORTED_UNITS = (nothing, u"pA")

abstract type AgilentFIDFormat end
struct AgilentFIDv179 <: AgilentFIDFormat end

# ── Loader Configuration Structures ──────────────────────────────────────────

struct ChemStationMSOptions{U<:Union{Nothing, Unitful.Units}}
    unit::U  # :arbitrary or :pA
end

struct AgilentFIDLoaderSpec{F<:AgilentFIDFormat}
    path::String
    options::ChemStationMSOptions
    _format::Val{F}
end

# ── Loader Spec Constructors ─────────────────────────────────────────────────

"""
    AgilentFIDLoaderSpec{F}(path::String; unit=nothing)

Internal structure representing a configured load request for Agilent FID data. The
`path` can be either a directory (containing `FID1A.ch`) or a direct path to the file.
The `unit` keyword sets the signal units (`nothing` for unitless or `u"pA"` for pA).

Normally, this is constructed indirectly via `AgilentFID(...)`.

See also [`AgilentFID`](@ref AgilentFID).
"""
function AgilentFIDLoaderSpec{F}(path::String; unit::U=nothing) where {F, U<:Union{Nothing, Unitful.Units}}
    unit ∈ SUPPORTED_UNITS || throw(ArgumentError(
        "Unsupported mode: $unit. Supported modes: $(SUPPORTED_UNITS)"))
    AgilentFIDLoaderSpec{F}(path, ChemStationMSOptions(unit), Val{F}())
end

"""
    AgilentFID(path::String) -> AgilentFIDLoaderSpec

Construct a loader specification for Agilent FID data. The `path` can be either a directory
(containing `FID1A.ch`) or a direct path to the file. The `unit` keyword sets the signal
units (`nothing` for unitless or `u"pA"` for pA). This convenience constructor currently
uses the `AgilentFIDv179` reader, but it may select other reader versions automatically in
future releases.

Returns a `AgilentFIDLoaderSpec` object used for deferred data loading via `load(...)`.

See also [`AgilentFIDLoaderSpec`](@ref AgilentFIDLoaderSpec).
"""
function AgilentFID(path::String; unit::Union{Nothing, Unitful.Units}=nothing)
    unit ∈ SUPPORTED_UNITS || throw(ArgumentError(
        "Unsupported mode: $unit. Supported modes: $(SUPPORTED_UNITS)"))
    AgilentFIDLoaderSpec{AgilentFIDv179}(path; unit=unit)
end

# ── Data Loading Interface ───────────────────────────────────────────────────

"""
    load(req::AgilentFIDLoaderSpec{AgilentFIDv179}) -> ChromScanSeries

Loads and parses FID data from Agilent FID (version 179) data files. The request must be
an `AgilentFIDLoaderSpec`, typically created via `AgilentFID(...)`. Returns a
[`ChromScanSeries`](@ref JuChrom.ChromScanSeries) containing a vector of
[`ChromScan`](@ref JuChrom.ChromScan) objects, along with parsed metadata (sample, user,
acquisition, instrument).

See also 
[`AgilentFID`](@ref AgilentFID),
[`AgilentFIDLoaderSpec`](@ref AgilentFIDLoaderSpec).
"""
function load(req::AgilentFIDLoaderSpec{AgilentFIDv179})
    path = req.path
    if isdir(path)
        datafile = joinpath(path, "FID1A.ch")
        isfile(datafile) || throw(MissingFileError(joinpath(path, "FID1A.ch")))
    elseif isfile(path)
        datafile = path
    else
        throw(MissingFileError(
            "Path does not exist or is not a valid file/directory: $path"))
    end
    return try
        readfile(datafile, req.options.unit)
    catch e
        throw(FileCorruptionError("Failed to read FID data at $datafile: $e"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Functions for extracting data from Agilent FID binary files (e.g., 'FID1A.ch')
# ─────────────────────────────────────────────────────────────────────────────

function readfile(file::AbstractString, unit::Union{Nothing, Unitful.Units})
    try
        open(file, "r") do f
            validate_header(f, file)

            metadata = read_metadata(f)

            scancount = read_scalar_ntoh_at(f, 278, UInt32)
            if scancount == 0 || scancount > 10^6
                throw(FileCorruptionError("Suspicious scan count: $scancount"))
            end
        
            # Start and stop time
            starttime, stoptime = read_vector_ntoh_at(f, 282, Float32, 2)
            retentions_unitfree = LinRange(starttime, stoptime, scancount)
            retentionunit = u"ms"

            # Signal unit
            intensityunit = read_len_prefixed_string_ltoh_at(f, 4172, 2)
            intensityunit == "pA" || throw(FileFormatError(
                "Unsupported intensity unit \"$intensityunit\" in file: $file"))

            if unit == u"pA"
                scalingfactor = read_scalar_ntoh_at(f, 4732, Float64)
                scalingfactor < 0 && 
                    throw(FileDecodingError("Negative scaling factor found in file: $file"))
            elseif unit === nothing
                scalingfactor = 1
            else
                throw(ArgumentError("Unsupported signal unit: $unit. Supported units: $(SUPPORTED_UNITS)"))
            end

            # Signals
            raw_intensities = read_vector_ltoh_at(f, 6144, Float64, scancount)
            any(raw_intensities .< 0) && @warn "Negative intensity values found in file: $file"

            scans = Vector{ChromScan{eltype(retentions_unitfree), typeof(retentionunit),
                eltype(raw_intensities), typeof(unit), @NamedTuple{}}}(undef, scancount)
            @inbounds for i in 1:scancount
                intensity = raw_intensities[i] * scalingfactor
                scans[i] = ChromScan(retentions_unitfree[i], retentionunit, intensity, 
                    unit)
            end
            ChromScanSeries(
                scans;
                instrument=metadata.instrument, 
                acquisition=metadata.acquisition,
                user=metadata.user,
                sample=metadata.sample
            )
        end
    catch e
        throw(FileCorruptionError("Failed to open/read Agilent FID file \"$file\": " * 
            "$e"))
    end
end

function validate_header(f::IO, file::AbstractString)
    ver = read_len_prefixed_string_ltoh(f, 1)
    ver == "179" || throw(FileFormatError("Unsupported version \"$ver\" in \"$file\""))
    typ = read_len_prefixed_string_ltoh_at(f, 347, 2)
    lowercase(typ) == "gc data file" || throw(FileFormatError(
        "Unsupported type \"$typ\" in \"$file\""))
    true
end

function read_metadata(f::IO)
    posmap = Dict(
        "sample"      =>  858,
        "description" => 1369,
        "operator"    => 1880,
        "datetime"    => 2391,
        "inlet"       => 2492,
        "type"        => 2533,
        "method"      => 2574,
        "software"    => 3089,
        "signal"      => 4213
    )


    rawmeta = Dict{String, Any}()
    for (key, pos) in pairs(posmap)
        str = read_len_prefixed_string_ltoh_at(f, pos, 2)
        rawmeta[key] = str ≠ "" ? str : "unknown"
    end

    instrument = (
        type        = rawmeta["type"], #
        inlet       = rawmeta["inlet"],  #
        method      = rawmeta["method"], #
        software    = rawmeta["software"], #
        signalunit  = rawmeta["signal"], #
    )

    acquisition = (
        datetime    = rawmeta["datetime"], # 
    )
    
    user = (
        operator    = rawmeta["operator"], # 
    )

    sample = (
        sample      = rawmeta["sample"], # 
        description = rawmeta["description"] # 
    )

    (instrument = instrument, acquisition = acquisition, user = user, sample = sample)
end

# ─────────────────────────────────────────────────────────────────────────────
# Helper functions for binary data reading and byte order conversion
# ─────────────────────────────────────────────────────────────────────────────

function ensure_bytes_available(f::IO, nbytes::Integer)
    curr_pos = position(f)
    seekend(f)
    end_pos = position(f)
    seek(f, curr_pos)
    bytes_left = end_pos - curr_pos
    if bytes_left < nbytes
        throw(UnexpectedEOFError(
            "Unexpected end of file: need $nbytes bytes, only $bytes_left available"))
    end
end

function read_scalar_ntoh(f::IO, ::Type{T}) where {T<:Union{Integer, AbstractFloat}}
    nbytes = sizeof(T)
    ensure_bytes_available(f, nbytes)
    value = read(f, T)
    ntoh(value)
end

read_scalar_ntoh_at(f::IO, pos::Integer, T::Type) = (seek(f, pos); read_scalar_ntoh(f, T))

function read_vector_ntoh_at(f::IO, pos::Integer, ::Type{T}, n::Integer) where {T<:Real}
    seek(f, pos)
    nbytes = sizeof(T) * n
    ensure_bytes_available(f, nbytes)
    value = read!(f, Vector{T}(undef, n))
    ntoh.(value)
end

function read_vector_ltoh_at(f::IO, pos::Integer, ::Type{T}, n::Integer) where {T<:Real}
    seek(f, pos)
    nbytes = sizeof(T) * n
    ensure_bytes_available(f, nbytes)
    value = read!(f, Vector{T}(undef, n))
    ltoh.(value)
end

function read_len_prefixed_string_ltoh(f::IO, bytes::Int=1)
    len = read_scalar_ltoh(f, UInt8)
    len == 0 && return ""

    raw = read_bytes_ltoh(f, len * bytes)
    ENCODING = bytes == 1 ? ENCODING_1_byte : bytes == 2 ? ENCODING_2_bytes : "unknown"

    try
        bytes == 1 || bytes == 2 ||
            throw(ArgumentError("Unsupported byte size for string encoding: $bytes"))
        decode(raw, ENCODING)
    catch e
        throw(FileDecodingError("Could not decode $len bytes as $ENCODING: $e"))
    end
end

function read_len_prefixed_string_ltoh_at(f::IO, pos::Integer, bytes::Int)
    seek(f, pos)
    read_len_prefixed_string_ltoh(f, bytes)
end

function read_bytes_ltoh(f::IO, nbytes::Integer)
    ensure_bytes_available(f, nbytes)
    value = read(f, nbytes)
    ltoh(value)
end

function read_scalar_ltoh(f::IO, ::Type{T}) where {T<:Union{Integer, AbstractFloat}}
    nbytes = sizeof(T)
    ensure_bytes_available(f, nbytes)
    value = read(f, T)
    ltoh(value)
end

end  # module
