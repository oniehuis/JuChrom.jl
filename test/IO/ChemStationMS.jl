using JuChrom
using Test

@testset "ChemStationMS" begin
    @test "data.ms" == JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS().datafilename
    @test "datasim.ms" == JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS(; 
        datafilename="datasim.ms").datafilename
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

    # IOError: unsupported source for ChemStationMS data import, because file does not 
    # exist
    dfile = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data2.ms")
    @test_throws JuChrom.InputOutput.IOError importdata(dfile, 
        JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS())
    
    # IOError: unsupported source for ChemStationMS data import, because specified D. 
    # folder does not exist
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "imaginary")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, 
        JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS())

    # IOError: not an Agilent .D folder, because folder exists, but does not end with .D
    dfolder = JuChrom.agilent
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, 
        JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS())

    # target file not in D. folder
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, 
        JuChrom.InputOutput.ChemStationMSReaders.ChemStationMS(datafilename="datasim.ms"))
end

