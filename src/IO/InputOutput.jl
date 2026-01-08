module InputOutput

using Logging
import FileIO: load

export as, @as
export DEFAULT_VERBOSITY
export FileCorruptionError
export FileDecodingError
export FileFormatError
export FileIOError
export MissingFileError
export MissingFolderError
export set_verbosity
export UnexpectedEOFError

# export load
# export LoaderSpec

# Default global verbosity setting for all loaders using this base
const DEFAULT_VERBOSITY = Ref{LogLevel}(Logging.Info)

# """
#     set_verbosity(level::LogLevel)

# Set global verbosity level for MyIO logging.

# Example:
#     set_verbosity(Logging.Debug)
# """
function set_verbosity(level::LogLevel)
    DEFAULT_VERBOSITY[] = level
    global_logger(SimpleLogger(stderr, level))
end

# """
#     LoaderSpec{T}

# A wrapper struct to specify the loader format type T and associated file path.

# Used to disambiguate between loaders when calling load.
# """
# struct LoaderSpec{T}
#     path::String
# end


# """
#     load(req::LoaderSpec)

# Generic interface to load data from the given LoaderSpec.

# Specific loaders should extend this method for their types.
# """
# function load(req::LoaderSpec)
#     throw(MethodError(load, (req,)))
# end


# """
#     as(::Type{T}, path::String) where T

# Construct a LoaderSpec instance specifying the loader format type T and path.

# Example:
#     req = as(AgilentFID, "/path/to/data")
#     data = load(req)
# """
as(::Type{T}, path::String) where T = LoaderSpec{T}(path)

# """
#     @as T "path"

# Macro version of as for creating a LoaderSpec instance.

# Example:
#     data = load(@as AgilentFID "/path/to/data")
# """
macro as(T, str)
    return :(LoaderSpec{$T}($str))
end

# """
#     FileIOError <: Exception

# Abstract base class for errors occurring during file input/output operations,
# such as malformed content, missing files, or incompatible formats.
# """
abstract type FileIOError <: Exception end

# """
#     FileFormatError(msg)

# Thrown when the format of a file is invalid or unsupported.
# """
struct FileFormatError <: FileIOError
    msg::String
end
Base.showerror(io::IO, e::FileFormatError) = print(io, "FileFormatError: ", e.msg)

# """
#     MissingFileError(path::AbstractString)

# Thrown when a required file is missing from the specified path.
# """
struct MissingFileError <: FileIOError
    path::String
end
Base.showerror(io::IO, e::MissingFileError) = 
    print(io, "MissingFileError: File not found at \"", e.path, "\"")

# """
#     MissingFolderError(path::AbstractString)

# Thrown when a required folder is missing from the specified path.
# """
struct MissingFolderError <: FileIOError
    path::String
end
Base.showerror(io::IO, e::MissingFolderError) = 
    print(io, "MissingFolderError: Folder not found at \"", e.path, "\"")

# """
#     FileCorruptionError(detail::AbstractString)

# Thrown when file content is corrupt or inconsistent.
# """
struct FileCorruptionError <: FileIOError
    msg::String
end
Base.showerror(io::IO, e::FileCorruptionError) = print(io, "FileCorruptionError: ", e.msg)

# """
#     FileDecodingError(msg)

# Thrown when byte sequences in a file cannot be decoded with the expected encoding.
# """
struct FileDecodingError <: FileIOError
    msg::String
end
Base.showerror(io::IO, e::FileDecodingError) = print(io, "FileDecodingError: ", e.msg)


# """
#     UnexpectedEOFError(msg::String)

# Thrown when an unexpected end-of-file (EOF) condition is encountered during file reading 
# operations, indicating that the data is incomplete or truncated.
# """
struct UnexpectedEOFError <: Exception
    msg::String
end
Base.showerror(io::IO, e::UnexpectedEOFError) = print(io, e.msg)

end # module
