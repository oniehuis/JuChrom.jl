using JuChrom
using Test

@testset "ChemStationMS" begin
    @test "data.ms" == JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS().datafilename
    @test "datasim.ms" == JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS(; datafilename="datasim.ms").datafilename

    # Check supertype
    @test isa(JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS(), FileFormat)
end

@testset "importdata ChemStationMS()" begin
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    gcms = importdata(dfolder, JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS())
    io = IOBuffer()
    Base.show(io, gcms)
    @test String(take!(io)) == string(
        "JuChrom.GCMS {scan times: Float32, ions: Float32, intensities: Int64}\n",
        "2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms\n",
        "5176 ions; range: m/z 29.0 - 562.9\n",
        "intensity range: 0 - 1186816\n",
        "metadata: 0 entries")

    dfile = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
    gcms = importdata(dfile, JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS())
    io = IOBuffer()
    Base.show(io, gcms)
    @test String(take!(io)) == string(
        "JuChrom.GCMS {scan times: Float32, ions: Float32, intensities: Int64}\n",
        "2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms\n",
        "5176 ions; range: m/z 29.0 - 562.9\n",
        "intensity range: 0 - 1186816\n",
        "metadata: 0 entries")
end

