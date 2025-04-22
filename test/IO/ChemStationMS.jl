using JuChrom
using Test

@testset "ChemStationMS" begin
    @test "data.ms" == ChemStationMS().datafilename
    @test "datasim.ms" == ChemStationMS(; datafilename="datasim.ms").datafilename
    @test isa(ChemStationMS(), FileFormat)
end

@testset "importdata ChemStationMS()" begin
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    chrom = importdata(dfolder, ChemStationMS())
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float32, ions: Float32, intensities: Int32}\n",
        "2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms\n",
        "5176 ions; range: m/z 29.0 - 562.9\n",
        "intensity range: 0 - 1186816\n",
        "metadata: 10 entries")

    dfile = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data.ms")
    chrom = importdata(dfile, ChemStationMS())
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Float32, ions: Float32, intensities: Int32}\n",
        "2405 scans; scan time range: 191941.0f0 ms - 1.899047f6 ms\n",
        "5176 ions; range: m/z 29.0 - 562.9\n",
        "intensity range: 0 - 1186816\n",
        "metadata: 10 entries")

    # IOError: unsupported source for ChemStationMS data import, because file does not 
    # exist
    dfile = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "data2.ms")
    @test_throws JuChrom.InputOutput.IOError importdata(dfile, ChemStationMS())
    
    # IOError: unsupported source for ChemStationMS data import, because specified D. 
    # folder does not exist
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D", "imaginary")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, ChemStationMS())

    # IOError: not an Agilent .D folder, because folder exists, but does not end with .D
    dfolder = JuChrom.agilent
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, ChemStationMS())

    # target file not in D. folder
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, ChemStationMS(
        datafilename="datasim.ms"))
end
