using JuChrom
using Test

@testset "AgilentFID" begin
    @test "FID1A.ch" == AgilentFID().datafilename
    @test "FID2A.ch" == AgilentFID(; datafilename="FID2A.ch").datafilename
    @test isa(AgilentFID(), FileFormat)
end

@testset "importdata AgilentFID()" begin
    dfolder = joinpath(JuChrom.agilent, "ZK_ONUBE_Mix1_11.D")
    fid = importdata(dfolder, AgilentFID())
    io = IOBuffer()
    Base.show(io, fid)
    @test String(take!(io)) == string(
        "JuChrom.Chrom {scan times: Float32, intensities: Float64}\n",
        "4151 scans; scan time range: 0.437f0 ms - 830000.44f0 ms\n",
        "intensity range: 0.0 - 1.0738316309895832e6\n",
        "metadata: 10 entries")

    dfile = joinpath(JuChrom.agilent, "ZK_ONUBE_Mix1_11.D", "FID1A.ch")
    fid = importdata(dfile, AgilentFID())
    io = IOBuffer()
    Base.show(io, fid)
    @test String(take!(io)) == string(
        "JuChrom.Chrom {scan times: Float32, intensities: Float64}\n",
        "4151 scans; scan time range: 0.437f0 ms - 830000.44f0 ms\n",
        "intensity range: 0.0 - 1.0738316309895832e6\n",
        "metadata: 10 entries")

    # IOError: unsupported source for AgilentFID data import, because file does not 
    # exist
    dfile = joinpath(JuChrom.agilent, "ZK_ONUBE_Mix1_11.D", "FID2A.ch")
    @test_throws JuChrom.InputOutput.IOError importdata(dfile, AgilentFID())
    
    # IOError: unsupported source for AgilentFID data import, because specified D. 
    # folder does not exist
    dfolder = joinpath(JuChrom.agilent, "ZK_ONUBE_Mix1_11.D", "imaginary")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, AgilentFID())

    # IOError: not an Agilent .D folder, because folder exists, but does not end with .D
    dfolder = JuChrom.agilent
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, AgilentFID())

    # target file not in D. folder
    dfolder = joinpath(JuChrom.agilent, "ZK_ONUBE_Mix1_11.D")
    @test_throws JuChrom.InputOutput.IOError importdata(dfolder, AgilentFID(
        datafilename="FID2A.ch"))
end
