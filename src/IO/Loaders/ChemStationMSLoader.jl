module ChemStationMSLoader

# ─────────────────────────────────────────────────────────────────────────────
# Imports, Dependencies, and Public Interface
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

import FileIO: load
using StringEncodings
using Unitful

# ── Internal Project Imports ─────────────────────────────────────────────────
import ..JuChrom: ChromScan, MassScan, ChromScanSeries, MassScanSeries
import ..InputOutput: FileCorruptionError, FileDecodingError, FileFormatError, 
    FileIOError, MissingFileError, UnexpectedEOFError

# ── Public API ───────────────────────────────────────────────────────────────

export ChemStationMS, ChemStationMSv2, ChemStationMSLoaderSpec, load

# ─────────────────────────────────────────────────────────────────────────────
# Loader Specification and Entry Points for Agilent ChemStation MS Files
# ─────────────────────────────────────────────────────────────────────────────

# ── Constants and Format Tags ────────────────────────────────────────────────

# const CSMSample01 = artifact"CSMSample01"
const ENCODING = "Windows-1252"
const MIN_CHUNK_SIZE = 16
const SUPPORTED_MODES = (:ms, :tic)

abstract type ChemStationMSFormat end
struct ChemStationMSv2 <: ChemStationMSFormat end

# ── Loader Configuration Structures ──────────────────────────────────────────

struct ChemStationMSOptions
    mode::Symbol  # :ms or :tic
end

struct ChemStationMSLoaderSpec{F<:ChemStationMSFormat}
    path::String
    options::ChemStationMSOptions
    _format::Val{F}
end

# ── Loader Spec Constructors ─────────────────────────────────────────────────

"""
    ChemStationMSLoaderSpec{F}(path::String; mode::Symbol=:ms)

Internal structure representing a configured load request for ChemStation MS data. The
`path` can be either a directory (containing `data.ms`) or a direct path to the file. The
`mode` keyword selects full mass spectral data (`:ms`) or the total ion chromatogram only
(`:tic`).

Normally, this is constructed indirectly via `ChemStationMS(...)`.

See also [`ChemStationMS`](@ref ChemStationMS).
"""
function ChemStationMSLoaderSpec{F}(path::String; mode::Symbol=:ms) where {F}
    mode ∈ SUPPORTED_MODES || throw(ArgumentError(
        "Unsupported mode: $mode. Supported modes: $(SUPPORTED_MODES)"))
    ChemStationMSLoaderSpec{F}(path, ChemStationMSOptions(mode), Val{F}())
end

"""
    ChemStationMS(path::String; mode::Symbol=:ms) -> ChemStationMSLoaderSpec

Construct a loader specification for Agilent ChemStation GC/MS data. The `path` can be
either a directory (containing `data.ms`) or a direct path to the file. The `mode` keyword
selects full mass spectral data (`:ms`) or the total ion chromatogram only (`:tic`). This
convenience constructor currently uses the `ChemStationMSv2` reader, but it may select other 
reader versions automatically in future releases.

Returns a `ChemStationMSLoaderSpec` object used for deferred data loading via `load(...)`.

See also 
[`ChemStationMSLoaderSpec`](@ref ChemStationMSLoaderSpec).
"""
function ChemStationMS(path::String; mode::Symbol=:ms)
    mode ∈ SUPPORTED_MODES || throw(ArgumentError(
        "Unsupported mode: $mode. Supported modes: $(SUPPORTED_MODES)"))
    ChemStationMSLoaderSpec{ChemStationMSv2}(path, mode=mode)
end

# ── Data Loading Interface ───────────────────────────────────────────────────

"""
    load(req::ChemStationMSLoaderSpec{ChemStationMSv2}) -> AbstractScanSeries

Loads and parses mass spectrometry data from Agilent ChemStation (version 2) GC/MS files.
The request must be a `ChemStationMSLoaderSpec`, typically created via `ChemStationMS(...)`.
Returns an `AbstractScanSeries` subtype containing either a vector of
[`MassScan`](@ref JuChrom.MassScan) objects for `mode=:ms` or a vector of
[`ChromScan`](@ref JuChrom.ChromScan) objects for `mode=:tic`, along with parsed metadata
(sample, user, acquisition, instrument).

See also 
[`ChemStationMS`](@ref ChemStationMS),
[`ChemStationMSLoaderSpec`](@ref ChemStationMSLoaderSpec).
"""
function load(req::ChemStationMSLoaderSpec{ChemStationMSv2})
    path = req.path
    if isdir(path)
        datafile = joinpath(path, "data.ms")
        isfile(datafile) || throw(MissingFileError(joinpath(path, "data.ms")))
    elseif isfile(path)
        datafile = path
    else
        throw(MissingFileError(
            "Path does not exist or is not a valid file/directory: $path"))
    end
    return try
        readfile(datafile, req.options.mode)
    catch e
        throw(FileCorruptionError("Failed to read MS data at $datafile: $e"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Functions for extracting data from Agilent ChemStation MS binary files 
# (e.g., 'data.ms')
# ─────────────────────────────────────────────────────────────────────────────

# ─── File Reader Interface ───────────────────────────────────────────────────

function readfile(file::AbstractString, mode::Symbol)
    mode ∈ (:ms, :tic) || throw(ArgumentError("Unsupported mode: $mode"))
    try
        open(file, "r") do f
            validate_header(f, file)

            metadata = read_metadata(f)

            scancount = read_scalar_ntoh_at(f, 278, UInt32)

            if scancount == 0 || scancount > 10^6
                throw(FileCorruptionError("Suspicious scan count: $scancount"))
            end

            scan_offset = 2 * read_scalar_ntoh_at(f, 266, UInt16) - 2
            seek(f, scan_offset)

            if mode == :ms
                retentions_unitfree, counts, mzs, intensities_unitfree = 
                    read_scan_data_ms(f, scancount)
                scans = build_mass_scans(retentions_unitfree, counts, mzs, 
                    intensities_unitfree)
                return MassScanSeries(
                    scans,
                    instrument=metadata.instrument, 
                    acquisition=metadata.acquisition,
                    user=metadata.user,
                    sample=metadata.sample
                )

            elseif mode == :tic
                retentions_unitfree, intensities_unitfree = read_scan_data_tic(f, scancount)
                retentionunit = u"ms"
                intensityunit = nothing
                scans = [ChromScan(retention_unitfree, retentionunit, intensity_unitfree, 
                    intensityunit) for (retention_unitfree, 
                    intensity_unitfree) in zip(retentions_unitfree, intensities_unitfree)]
                return ChromScanSeries(
                    scans,
                    instrument=metadata.instrument, 
                    acquisition=metadata.acquisition,
                    user=metadata.user,
                    sample=metadata.sample
                )
            end
        end
    catch e
        throw(FileCorruptionError("Failed to open/read ChemStationMS file \"$file\": " * 
            "$e"))
    end
end

# ─── File Header Validation and Metadata Extraction ──────────────────────────

function validate_header(f::IO, file::AbstractString)
    ver = read_len_prefixed_string_ltoh(f)
    ver == "2" || throw(FileFormatError("Unsupported version \"$ver\" in \"$file\""))
    seek(f, 4)
    typ = read_len_prefixed_string_ltoh(f)
    lowercase(typ) == "gc / ms data file" || throw(FileFormatError(
        "Unsupported type \"$typ\" in \"$file\""))
    true
end

function read_metadata(f::IO)
    posmap = Dict(
        "sample"      =>  24,
        "description" =>  86,
        "operator"    => 148,
        "datetime"    => 178,
        "type"        => 208,
        "inlet"       => 218,
        "method"      => 228,
    )

    rawmeta = Dict{String, Any}()
    for (key, pos) in pairs(posmap)
        str = read_len_prefixed_string_ltoh_at(f, pos)
        rawmeta[key] = str ≠ "" ? str : "unknown"
    end

    rawmeta["sequence"]  = Int(read_scalar_ntoh_at(f, 252, Int16))
    rawmeta["vial"]      = Int(read_scalar_ntoh_at(f, 254, Int16))
    rawmeta["replicate"] = Int(read_scalar_ntoh_at(f, 256, Int16))

    instrument = (
        type        = rawmeta["type"],
        inlet       = rawmeta["inlet"], 
        method      = rawmeta["method"],
    )

    acquisition = (
        datetime    = rawmeta["datetime"],
        sequence    = rawmeta["sequence"], 
        vial        = rawmeta["vial"],
        replicate   = rawmeta["replicate"],
    )
    
    user = (
        operator    = rawmeta["operator"],
    )

    sample = (
        sample      = rawmeta["sample"],
        description = rawmeta["description"]
    )

    (instrument = instrument, acquisition = acquisition, user = user, sample = sample)
end

# ─── Scan Data Reading ───────────────────────────────────────────────────────

function read_scan_data_ms(f::IO, scancount::Integer)

    retentions_unitfree = Vector{Float32}(undef, scancount)
    counts = Vector{Int32}(undef, scancount)
    mzs = Float32[]
    ints = Int32[]

    for i in 1:scancount
        try
            # Read chunk length (in 2-byte units)
            length_units = read_scalar_ntoh(f, UInt16)
            chunk_size_bytes = 2 * length_units - 2  # subtract the 2 bytes just read

            # Basic corruption check
            chunk_size_bytes < MIN_CHUNK_SIZE && throw(FileCorruptionError(
                "Scan $i chunk too small ($chunk_size_bytes bytes) – possibly corrupt"))

            # Check enough bytes available
            ensure_bytes_available(f, chunk_size_bytes)

            # Read the chunk into a buffer and wrap in IOBuffer
            chunk = Vector{UInt8}(undef, chunk_size_bytes)
            read!(f, chunk)
            io = IOBuffer(chunk)

            # Read scan time (4 bytes Int32)
            retentions_unitfree[i] = ntoh(read(io, Int32))

            # Read number of points
            total_data_length = ntoh(read(io, UInt16))  # includes 6 header bytes
            n = fld(total_data_length - 6, 2)
            n < 0 && throw(FileCorruptionError(
                "Negative number of points ($n) in scan $i – corrupt header"))
            n > 2048 && @warn "Unusually high point count ($n) in scan $i – " * 
                "continuing anyway"
            counts[i] = n

            # Skip 10 bytes of additional header
            seek(io, position(io) + 10)

            # Read and process the packed data: 2n UInt16 values
            raw_data = Vector{UInt16}(undef, 2n)
            read!(io, raw_data)
            raw_data = ntoh.(raw_data)

            # Sanity check
            @assert length(raw_data) == 2n "Corrupt scan $i: expected 2n values, " * 
                "got $(length(raw_data))"

            # Split into m/z and intensity raw data
            mz_raw = raw_data[1:2:end]
            int_raw = raw_data[2:2:end]

            # Scale m/z values
            scan_mzs = Float32.(mz_raw) ./ 20
            any(isnan, scan_mzs) && throw(FileCorruptionError(
                "NaN detected in scan $i: corrupt m/z values"))

            # Unpack intensity: exponent in bits 14–15, mantissa in lower 14 bits
            exponents = int_raw .>> 14
            mantissas = int_raw .& 0x3FFF
            scan_ints = Int32.(mantissas .* (8 .^ exponents))

            # Append scan data to overall arrays
            append!(mzs, scan_mzs)
            append!(ints, scan_ints)
        catch e
            throw(FileCorruptionError("Error reading MS scan $i of $scancount: $e"))
        end
    end

    retentions_unitfree, counts, mzs, ints
end

function read_scan_data_tic(
    f::IO, 
    scancount::Integer
    )::Tuple{Vector{Float32}, Vector{Float32}}
    
    retentions_unitfree = Vector{Float32}(undef, scancount)
    tic = Vector{Float32}(undef, scancount)

    for i in 1:scancount
        try
            # Read chunk length (in 2-byte units)
            length_units = read_scalar_ntoh(f, UInt16)
            chunk_size_bytes = 2 * length_units - 2  # subtract the 2 bytes just read
            
            # Basic corruption check
            chunk_size_bytes < MIN_CHUNK_SIZE && throw(FileCorruptionError(
                "Scan $i chunk too small ($chunk_size_bytes bytes) – possibly corrupt"))

            # Check enough bytes available
            ensure_bytes_available(f, chunk_size_bytes)

            # Read scan time (4 bytes Int32)
            retentions_unitfree[i] = read_scalar_ntoh(f, Int32)

            # Skip the rest of the chunk except the last 4 bytes (TIC value)
            skip_bytes = chunk_size_bytes - 8  # 4 bytes for scan time and TIC value each
            seek(f, position(f) + skip_bytes)

            # Read TIC intensity (4 bytes Int32)
            tic[i] = read_scalar_ntoh(f, Int32) 
        catch e
            throw(FileCorruptionError("Error reading TIC scan $i of $scancount: $e"))
        end
    end

    retentions_unitfree, tic
end

# ─── Mass Scan Construction ──────────────────────────────────────────────────

function build_mass_scans(
    retentions_unitfree::AbstractVector{<:Real},
    counts::AbstractVector{<:Integer},
    mzs::AbstractVector{<:Real},
    ints::AbstractVector{<:Real})

    length(retentions_unitfree) == length(counts) || throw(DimensionMismatch(
        "Retentions and counts must have the same length: " * 
        "$(length(retentions_unitfree)) vs $(length(counts))"))
    length(mzs) == sum(counts) || throw(DimensionMismatch(
        "Total m/z count $(length(mzs)) does not match sum of counts $(sum(counts))"))
    length(ints) == sum(counts) || throw(DimensionMismatch(
        "Total intensity count $(length(ints)) does not match " * 
        "sum of counts $(sum(counts))"))
    
    retentionunit = u"ms"
    intensityunit = nothing
    level = 1
    scans = Vector{MassScan{Float32, typeof(retentionunit), Vector{Float32}, Nothing, 
        Vector{Int32}, typeof(intensityunit), typeof(level), @NamedTuple{}}}(undef, 
        length(counts))

    offset = 0
    for i in eachindex(counts)
        n = counts[i]
        x = @view mzs[offset+1:offset+n]
        y = @view ints[offset+1:offset+n]
        p = sortperm(x)
        scans[i] = MassScan(retentions_unitfree[i], retentionunit, x[p], nothing, y[p], 
            intensityunit)
        offset += n
    end

    scans
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

function read_len_prefixed_string_ltoh(f::IO)
    len = read_scalar_ltoh(f, UInt8)
    len == 0 && return ""
    bytesting = read_bytes_ltoh(f, len)
    return try
        decode(bytesting, ENCODING)
    catch e
        throw(FileDecodingError("Could not decode $len bytes as $ENCODING: $e"))
    end
end

function read_len_prefixed_string_ltoh_at(f::IO, pos::Integer)
    seek(f, pos)
    val = read_len_prefixed_string_ltoh(f)
end

function read_scalar_ltoh(f::IO, ::Type{T}) where {T<:Union{Integer, AbstractFloat}}
    nbytes = sizeof(T)
    ensure_bytes_available(f, nbytes)
    value = read(f, T)
    ltoh(value)
end

function read_bytes_ltoh(f::IO, nbytes::Integer)
    ensure_bytes_available(f, nbytes)
    value = read(f, nbytes)
    ltoh(value)
end

function read_scalar_ntoh(f::IO, ::Type{T}) where {T<:Union{Integer, AbstractFloat}}
    nbytes = sizeof(T)
    ensure_bytes_available(f, nbytes)
    value = read(f, T)
    ntoh(value)
end

read_scalar_ntoh_at(f::IO, pos::Integer, T::Type) = (seek(f, pos); read_scalar_ntoh(f, T))

end  # module
