module InputOutputTests

# ──────────────────────────────────────────────────────────────────────────────────────────
# Unit tests for ./IO/InputOutput.jl
# ──────────────────────────────────────────────────────────────────────────────────────────

using Test
using JuChrom
using Logging
using JuChrom.AgilentFIDLoader: AgilentFID
using JuChrom.InputOutput: FileDecodingError, FileFormatError, FileCorruptionError, 
    FileIOError

# ──────────────────────────────────────────────────────────────────────────────────────────
# Unit tests
# ──────────────────────────────────────────────────────────────────────────────────────────

# ── FileIOError ───────────────────────────────────────────────────────────────────────────

@testset "FileIOError" begin
    fd = FileDecodingError("oops")
    ff = FileFormatError("bad")
    fi = FileCorruptionError("corrupt")
    fe = FileIOError  # abstract

    @test fd isa FileDecodingError
    @test fd isa FileIOError
    @test ff isa FileFormatError
    @test ff isa FileIOError
    @test fi isa FileCorruptionError
    @test fi isa FileIOError
end

# ── set_verbosity ─────────────────────────────────────────────────────────────────────────

@testset "set_verbosity" begin
    # Set to Debug and confirm the Ref is updated
    set_verbosity(Logging.Debug)
    @test JuChrom.InputOutput.DEFAULT_VERBOSITY[] == Logging.Debug

    # Set back to Info and confirm update
    set_verbosity(Logging.Info)
    @test JuChrom.InputOutput.DEFAULT_VERBOSITY[] == Logging.Info

    # Invalid type should error (non-LogLevel)
    @test_throws MethodError set_verbosity("not a level")
end

end  # module InputOutputTests
