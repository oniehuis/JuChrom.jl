using JuChrom
using Test

@testset "MassHunterMS" begin
    @test Float32 == MassHunterMS().scantimetype
    @test Float64 == MassHunterMS(; scantimetype=Float64).scantimetype
    @test Float32 == MassHunterMS().iontype
    @test Float64 == MassHunterMS(; iontype=Float64).iontype
    @test Float32 == MassHunterMS().intensitytype
    @test Float64 == MassHunterMS(; intensitytype=Float64).intensitytype
    @test isa(MassHunterMS(), FileFormat)
end

@testset "importdata MassHunterMS()" begin
    dfolder = joinpath(JuChrom.agilent, "C7-C40_MassHunterMS.D")
    gcms = importdata(dfolder, MassHunterMS())
    io = IOBuffer()
    Base.show(io, gcms)
    @test String(take!(io)) == string(
        "JuChrom.GCMS {scan times: Float32, ions: Float32, intensities: Float32}\n",
        "2405 scans; scan time range: 3.1990166f0 minute - 31.650784f0 minute\n",
        "50275 ions; range: m/z 29.02 - 562.89\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 3 entries")

    gcms = importdata(dfolder, MassHunterMS(; scantimetype=Float64))
    io = IOBuffer()
    Base.show(io, gcms)
    @test String(take!(io)) == string(
        "JuChrom.GCMS {scan times: Float64, ions: Float32, intensities: Float32}\n",
        "2405 scans; scan time range: 3.1990166666666666 minute - 31.650783333333333 ", 
            "minute\n",
        "50275 ions; range: m/z 29.02 - 562.89\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 3 entries")

    gcms = importdata(dfolder, MassHunterMS(; iontype=Float64))
    io = IOBuffer()
    Base.show(io, gcms)
    @test String(take!(io)) == string(
        "JuChrom.GCMS {scan times: Float32, ions: Float64, intensities: Float32}\n",
        "2405 scans; scan time range: 3.1990166f0 minute - 31.650784f0 minute\n",
        "50275 ions; range: m/z 29.020000457763672 - 562.8900146484375\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 3 entries")

    gcms = importdata(dfolder, MassHunterMS(; intensitytype=Float64))
    io = IOBuffer()
    Base.show(io, gcms)
    @test String(take!(io)) == string(
        "JuChrom.GCMS {scan times: Float32, ions: Float32, intensities: Float64}\n",
        "2405 scans; scan time range: 3.1990166f0 minute - 31.650784f0 minute\n",
        "50275 ions; range: m/z 29.02 - 562.89\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 3 entries")

    gcms = importdata(dfolder, MassHunterMS(; scantimetype=Float64, iontype=Float64, 
        intensitytype=Float64))
    io = IOBuffer()
    Base.show(io, gcms)
    @test String(take!(io)) == string(
        "JuChrom.GCMS {scan times: Float64, ions: Float64, intensities: Float64}\n",
        "2405 scans; scan time range: 3.1990166666666666 minute - 31.650783333333333 ", 
        "minute\n",
        "50275 ions; range: m/z 29.020000457763672 - 562.8900146484375\n",
        "intensity range: 0.0 - 1.1872475e6\n",
        "metadata: 3 entries")

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
