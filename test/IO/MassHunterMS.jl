using JuChrom
using Test

@testset "MassHunterMS" begin
    @test Nothing == MassHunterMS().scantimetype
    @test Float32 == MassHunterMS(; scantimetype=Float32).scantimetype
    @test Nothing == MassHunterMS().iontype
    @test Float32 == MassHunterMS(; iontype=Float32).iontype
    @test Nothing == MassHunterMS().intensitytype
    @test Float32 == MassHunterMS(; intensitytype=Float32).intensitytype
    @test isa(MassHunterMS(), FileFormat)
end

@testset "importdata MassHunterMS()" begin
    dfolder = joinpath(JuChrom.agilent, "C7-C40_MassHunterMS.D")
    chrom = importdata(dfolder, MassHunterMS())
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float64, ions: Float64, intensities: Float64}\n",
        "2405 scans; scan time range: 191.941 s - 1899.047 s\n",
        "50275 ions; range: m/z 29.020000457763672 - 562.8900146484375\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 2 entries")

    chrom = importdata(dfolder, MassHunterMS(; scantimetype=Float32))
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float32, ions: Float64, intensities: Float64}\n",
        "2405 scans; scan time range: 191.941f0 s - 1899.047f0 s\n",
        "50275 ions; range: m/z 29.020000457763672 - 562.8900146484375\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 2 entries")

    chrom = importdata(dfolder, MassHunterMS(; iontype=Float32))
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float64, ions: Float32, intensities: Float64}\n",
        "2405 scans; scan time range: 191.941 s - 1899.047 s\n",
        "50275 ions; range: m/z 29.02 - 562.89\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 2 entries")

    chrom = importdata(dfolder, MassHunterMS(; intensitytype=Float32))
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float64, ions: Float64, intensities: Float32}\n",
        "2405 scans; scan time range: 191.941 s - 1899.047 s\n",
        "50275 ions; range: m/z 29.020000457763672 - 562.8900146484375\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 2 entries")

    chrom = importdata(dfolder, MassHunterMS(; scantimetype=Float32, iontype=Float32, 
        intensitytype=Float32))
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float32, ions: Float32, intensities: Float32}\n",
        "2405 scans; scan time range: 191.941f0 s - 1899.047f0 s\n",
        "50275 ions; range: m/z 29.02 - 562.89\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 2 entries")

    # ERROR: IOError: unsupported source for MassHunterMS data import 
    dfolder = joinpath(JuChrom.agilent, "C40_MassHunterMS", "imaginary")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, MassHunterMS())

    # IOError: not an Agilent .D folder
    dfolder = joinpath(JuChrom.agilent)
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, MassHunterMS())

    # ERROR: IOError: cannot find folder AcqData in Agilent .D folder
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, MassHunterMS())

    # ERROR: IOError: cannot find "MSPeak.bin" file in path
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, MassHunterMS())

    # ERROR: IOError: cannot find "MSScan.bin" file in path
    dfolder = joinpath(JuChrom.agilent, "ZK_ONUBE_Mix1_11.D")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, MassHunterMS())
end
