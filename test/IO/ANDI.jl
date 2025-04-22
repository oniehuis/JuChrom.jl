using JuChrom
using Test

@testset "ANDI" begin
    @test Nothing == ANDI().scantimetype
    @test Float32 == ANDI(; scantimetype=Float32).scantimetype
    @test Nothing == ANDI().iontype
    @test Float32 == ANDI(; iontype=Float32).iontype
    @test Nothing == ANDI().intensitytype
    @test Float32 == ANDI(; intensitytype=Float32).intensitytype
    @test isa(ANDI(), FileFormat)
end

@testset "importdata ANDI()" begin
    file = joinpath(JuChrom.andi, "C7-C40_13_Nov_1.CDF")
    chrom = importdata(file, ANDI())
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float64, ions: Float32, intensities: Float32}\n",
        "5221 scans; scan time range: 191.942 s - 3898.719 s\n",
        "5248 ions; range: m/z 29.0 - 562.9\n",
        "intensity range: 0.0 - 1.051136e6\n",
        "metadata: 2 entries")

    chrom = importdata(file, ANDI(; scantimetype=Float32))
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float32, ions: Float32, intensities: Float32}\n",
        "5221 scans; scan time range: 191.942f0 s - 3898.719f0 s\n",
        "5248 ions; range: m/z 29.0 - 562.9\n",
        "intensity range: 0.0 - 1.051136e6\n",
        "metadata: 2 entries")

    chrom = importdata(file, ANDI(; iontype=Float64))
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float64, ions: Float64, intensities: Float32}\n",
        "5221 scans; scan time range: 191.942 s - 3898.719 s\n",
        "5248 ions; range: m/z 29.0 - 562.9000244140625\n",
        "intensity range: 0.0 - 1.051136e6\n",
        "metadata: 2 entries")

    chrom = importdata(file, ANDI(; intensitytype=Float64))
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float64, ions: Float32, intensities: Float64}\n",
        "5221 scans; scan time range: 191.942 s - 3898.719 s\n",
        "5248 ions; range: m/z 29.0 - 562.9\n",
        "intensity range: 0.0 - 1.051136e6\n",
        "metadata: 2 entries")

    chrom = importdata(file, ANDI(; scantimetype=Float32, iontype=Float64, 
        intensitytype=Float64))
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float32, ions: Float64, intensities: Float64}\n",
        "5221 scans; scan time range: 191.942f0 s - 3898.719f0 s\n",
        "5248 ions; range: m/z 29.0 - 562.9000244140625\n",
        "intensity range: 0.0 - 1.051136e6\n",
        "metadata: 2 entries")

    # ERROR: IOError: unsupported source for ANDI data import 
    file = JuChrom.andi
    @test_throws JuChrom.InputOutput.IOError importdata(file, ANDI())
end
