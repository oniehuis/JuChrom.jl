using JuChrom
using Test

@testset "ShimadzuMS" begin
    @test isa(ShimadzuMS(), FileFormat)
end

@testset "importdata ShimadzuMS()" begin
    file = joinpath(JuChrom.shimadzu, "AR190311.qgd")
    chrom = importdata(file, ShimadzuMS())
    io = IOBuffer()
    Base.show(io, chrom)
    @test String(take!(io)) == string(
        "JuChrom.ChromMS {scan times: Int32, ions: Float64, intensities: UInt32}\n",
        "9210 scans; scan time range: 183000 ms - 2024800 ms\n",
        "3620 ions; range: m/z 34.9 - 600.4\n",
        "intensity range: 0 - 924954\n",
        "metadata: 0 entries")

    # ERROR: IOError: unsupported source for ANDI data import
    @test_throws JuChrom.InputOutput.IOError importdata(JuChrom.shimadzu, ShimadzuMS())
    file = joinpath(JuChrom.andi, "C7-C40_13_Nov_1.CDF")
    @test_throws JuChrom.InputOutput.IOError importdata(file, ShimadzuMS())
end
