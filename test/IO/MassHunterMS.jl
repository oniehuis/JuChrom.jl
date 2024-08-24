using JuChrom
using Test

# @testset "MassHunterMS" begin
#     @test "data.ms" == MassHunterMS().datafilename
#     @test "datasim.ms" == MassHunterMS(; datafilename="datasim.ms").datafilename
#     @test isa(MassHunterMS(), FileFormat)
# end

@testset "importdata MassHunterMS()" begin
    dfolder = joinpath(JuChrom.agilent, "C7-C40_MassHunterMS.D")
    gcms = importdata(dfolder, MassHunterMS())
    io = IOBuffer()
    Base.show(io, gcms)
    @test String(take!(io)) == string(
        "JuChrom.GCMS {scan times: Float64, ions: Float64, intensities: Int64}\n",
        "2405 scans; scan time range: 3.199 minute - 31.651 minute\n",
        "5174 ions; range: m/z 29.0 - 562.9\n",
        "intensity range: 0 - 1187248\n",
        "metadata: 0 entries")

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
