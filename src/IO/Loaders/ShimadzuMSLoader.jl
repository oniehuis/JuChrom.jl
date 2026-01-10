module ShimadzuMSLoader

# ─────────────────────────────────────────────────────────────────────────────
# Imports, Dependencies, and Public Interface
# ─────────────────────────────────────────────────────────────────────────────

# ── External Dependencies ────────────────────────────────────────────────────

import FileIO: load
using Unitful
using PyCall

const olefile = PyNULL()  # shared PyCall handle for the olefile module

function __init__()
    copy!(olefile, pyimport_conda("olefile", "olefile"))
end

# ── Internal Project Imports ─────────────────────────────────────────────────
import JuChrom: ChromScan, MassScan, ChromScanSeries, MassScanSeries
import JuChrom.InputOutput: FileCorruptionError, FileFormatError, MissingFileError

# ── Public API ───────────────────────────────────────────────────────────────

export ShimadzuMS, ShimadzuMSv1, ShimadzuMSLoaderSpec, load

# ─────────────────────────────────────────────────────────────────────────────
# Loader Specification and Entry Points for Shimadzu MS Files
# ─────────────────────────────────────────────────────────────────────────────

# ── Constants and Format Tags ────────────────────────────────────────────────

const SUPPORTED_MODES = (:ms, :tic)  # :ms → full spectra, :tic → chromatogram only

abstract type ShimadzuMSFormat end
struct ShimadzuMSv1 <: ShimadzuMSFormat end

# ── Loader Configuration Structures ──────────────────────────────────────────

struct ShimadzuMSOptions
    mode::Symbol  # :ms or :tic
end

"""
    ShimadzuMSLoaderSpec{F}

Loader specification for Shimadzu GC/MS files. Stores the file path, selected mode,
and a format tag used by `load`.
"""
struct ShimadzuMSLoaderSpec{F<:ShimadzuMSFormat}
    path::String
    options::ShimadzuMSOptions
    _format::Val{F}
end

# ── Loader Spec Constructors ─────────────────────────────────────────────────

"""
    ShimadzuMSLoaderSpec{F}(path::String; mode::Symbol=:ms)

Internal structure representing a configured load request for Shimadzu MS data.
This type encapsulates:
- the source path,
- selected mode (`:ms` or `:tic`), and
- a type tag for format versioning.

Normally, this is constructed indirectly via `ShimadzuMS(...)`.
"""
function ShimadzuMSLoaderSpec{F}(path::String; mode::Symbol=:ms) where {F}
    mode ∈ SUPPORTED_MODES || throw(ArgumentError(
        "Unsupported mode: $mode. Supported modes: $(SUPPORTED_MODES)"))
    ShimadzuMSLoaderSpec{F}(path, ShimadzuMSOptions(mode), Val{F}())
end

"""
    ShimadzuMS(path::String; mode::Symbol=:ms) -> ShimadzuMSLoaderSpec

Construct a loader specification for Agilent Shimadzu GC/MS data. The `path` must be
an explicit path to the `.qgd`/`.ms` file. 

The keyword argument `mode` determines the type of data to load:
- `:ms`  → full mass spectral data (default)
- `:tic` → total ion chromatogram only

Returns a `ShimadzuMSLoaderSpec` object used for deferred data loading via `load(...)`.
"""
function ShimadzuMS(path::String; mode::Symbol=:ms)
    mode ∈ SUPPORTED_MODES || throw(ArgumentError(
        "Unsupported mode: $mode. Supported modes: $(SUPPORTED_MODES)"))
    ShimadzuMSLoaderSpec{ShimadzuMSv1}(path, mode=mode)
end

# ── Data Loading Interface ───────────────────────────────────────────────────

"""
    load(req::ShimadzuMSLoaderSpec{ShimadzuMSv1}) -> AbstractScanSeries

Loads and parses mass spectrometry data from Agilent Shimadzu (version 1) GC/MS files.
Takes a `ShimadzuMSLoaderSpec` created via `ShimadzuMS(...)`.

Returns a `AbstractScanSeries` subtype object containing either:
- a vector of `MassScan` objects (for `mode=:ms`), or
- a vector of `ChromScan` objects (for `mode=:tic`),

along with parsed metadata (sample, user, acquisition, instrument).
"""
function load(req::ShimadzuMSLoaderSpec{ShimadzuMSv1})
    path = req.path
    # Shimadzu loader requires explicit file paths (no directory discovery)
    isfile(path) || throw(MissingFileError(
        "ShimadzuMS expects a direct file path, but got missing/non-file path $path"))
    try
        readfile(path, req.options.mode)
    catch e
        isa(e, FileCorruptionError) && rethrow()
        throw(FileCorruptionError("Failed to read Shimadzu MS data at $path: $e"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Functions for extracting data from Shimadzu MS files
# ─────────────────────────────────────────────────────────────────────────────

function readbytestring(file::AbstractString, streamname::Vector{String})
    # Open the OLE container and pull a raw byte buffer for the requested stream
    ole::PyObject = olefile.OleFileIO(file)
    msrawdata::String = try
        ole.exists(streamname)::Bool || throw(
            FileFormatError("file does not contain streamname $streamname"))
        msrawdatastream::PyObject = ole.openstream(streamname)
        msrawdatastream.read()
    finally
        ole.close()
    end
    # Transcode the data into a byte string
    transcode(UInt8, msrawdata)::Base.CodeUnits{UInt8, String}
end

function extractdata(bytesting)
    # Shimadzu encodes each scan as a 32-byte header followed by packed m/z/int data
    i = 1  # byte offset within the overall stream
    scantimes = Vector{Int32}()
    pointcounts = Vector{Int16}()
    mzvec = Vector{Float64}()
    intsvec = Vector{UInt32}()
    while i < length(bytesting)
        push!(scantimes, ltoh.(reinterpret(Int32, bytesting[i+4:i+7]))[1])
        bytecount = ltoh.(reinterpret(Int16, bytesting[i+20:i+21]))[1]
        if !(1 ≤ bytecount ≤ 4)
            throw(FileFormatError(
                "Unexpected Shimadzu MS byte count $bytecount (expected 1–4)"))
        end
        pointcount = ltoh.(reinterpret(Int16, bytesting[i+22:i+23]))[1]
        push!(pointcounts, pointcount)
        start = i + 32
        stepwidth = 2 + bytecount
        stop = start + stepwidth * (pointcount - 1)
        scanmzs = Vector{Float64}(undef, pointcount)
        scanints = Vector{UInt32}(undef, pointcount)
        for (s, z) in enumerate(start:stepwidth:stop)
            # m/z values are scaled Int16, intensities are variable-width UInts
            scanmzs[s] = ltoh.(reinterpret(Int16, bytesting[z:z+1]))[1] / 20
            intbytes = zeros(UInt8, 4) 
            intbytes[1:bytecount] = bytesting[z+2:z+1+bytecount]
            int = ltoh.(reinterpret(UInt32, intbytes))[1]
            scanints[s] = int
        end
        append!(mzvec, scanmzs)
        append!(intsvec, scanints)
        i += 32 + (2 + bytecount) * pointcount
    end
    scantimes, pointcounts, mzvec, intsvec
end

function build_mass_scans(
    retentions_unitfree::AbstractVector{<:Real},
    counts::AbstractVector{<:Integer},
    mzs::AbstractVector{<:Real},
    ints::AbstractVector{<:Real})

    # Defensive length checks catch inconsistencies introduced during parsing
    length(retentions_unitfree) == length(counts) || throw(DimensionMismatch(
        "Retentions and counts must have the same length: " * 
        "$(length(retentions_unitfree)) vs $(length(counts))"))
    length(mzs) == sum(counts) || throw(DimensionMismatch(
        "Total m/z count $(length(mzs)) does not match sum of counts $(sum(counts))"))
    length(ints) == sum(counts) || throw(DimensionMismatch(
        "Total intensity count $(length(ints)) does not match " * 
        "sum of counts $(sum(counts))"))
    
    retention_unit = u"ms"
    intensity_unit = nothing
    level = 1
    scans = Vector{MassScan{Int32, typeof(retention_unit), Vector{Float64}, Nothing, 
        Vector{UInt32}, typeof(intensity_unit), typeof(level), @NamedTuple{}}}(undef, 
        length(counts))

    offset = 0
    for i in eachindex(counts)
        n = counts[i]
        x = @view mzs[offset+1:offset+n]
        y = @view ints[offset+1:offset+n]
        p = sortperm(x)
        scans[i] = MassScan(retentions_unitfree[i], retention_unit, x[p], nothing, y[p], 
            intensity_unit)
        offset += n
    end

    scans
end

# ─── File Reader Interface ───────────────────────────────────────────────────

function readfile(file::AbstractString, mode::Symbol)
    mode ∈ (:ms, :tic) || throw(ArgumentError("Unsupported mode: $mode"))
    # All supported Shimadzu formats are stored as OLE compound documents
    olefile.isOleFile(file) || throw(FileFormatError("not an OLE file: $file"))
    # Delegate to the appropriate stream decoder
    if mode == :ms
        return load_mass_spectra(file)
    elseif mode == :tic
        return load_tic(file)
    else
        throw(ArgumentError("Unsupported mode: $mode"))
    end
end

function load_mass_spectra(file::AbstractString)
    try
        # Read the packed MS Raw Data stream and decode into scan-level arrays
        streamname = ["GCMS Raw Data", "MS Raw Data"]
        bytestring = readbytestring(file, streamname)
        retentions_unitfree, counts, mzs, ints = extractdata(bytestring)
        scans = build_mass_scans(retentions_unitfree, counts, mzs, ints)
        MassScanSeries(
            scans,
            instrument=NamedTuple(),
            acquisition=NamedTuple(),
            user=NamedTuple(),
            sample=NamedTuple()
        )
    catch e
        throw(FileCorruptionError(
            "Failed to read Shimadzu MS spectral data from \"$file\": $e"))
    end
end

function load_tic(file::AbstractString)
    try
        # Shimadzu stores TIC retention times and intensities in separate streams
        streamname = ["GCMS Raw Data", "Retention Time"]
        bytestring = readbytestring(file, streamname)
        if length(bytestring) % 4 != 0
            throw(FileFormatError(
                "Shimadzu retention time stream size is not divisible by 4 bytes"))
        end
        retentions_unitfree = ltoh.(reinterpret(Int32, bytestring))
        retention_unit = u"ms"

        streamname = ["GCMS Raw Data", "TIC Data"]
        bytestring = readbytestring(file, streamname)
        if length(bytestring) % 8 != 0
            throw(FileFormatError(
                "Shimadzu TIC data stream size is not divisible by 8 bytes"))
        end
        intensities_unitfree = ltoh.(reinterpret(Int64, bytestring))
        intensity_unit = nothing

        # Pair retention time/intensity into individual ChromScan entries
        scans = [ChromScan(retention_unitfree, retention_unit, intensity_unitfree,
            intensity_unit) for (retention_unitfree,
            intensity_unitfree) in zip(retentions_unitfree, intensities_unitfree)]

        ChromScanSeries(
            scans,
            instrument=NamedTuple(),
            acquisition=NamedTuple(),
            user=NamedTuple(),
            sample=NamedTuple()
        )
    catch e
        throw(FileCorruptionError(
            "Failed to read Shimadzu MS TIC data from \"$file\": $e"))
    end
end

end  # module
