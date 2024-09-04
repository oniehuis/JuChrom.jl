using JuChrom
using Test
using DelimitedFiles

@testset "DelimitedTextWriter" begin
    rootpath = joinpath(splitpath(pathof(JuChrom))[begin:end-2])
    pathname = joinpath(rootpath, "tmp", "delimtest")
    chrom = Chrom((1:2)u"minute", [123, 224])
    exportdata(chrom, pathname, DelimitedText(), overwrite=true)
    data = readdlm(string(pathname, ".tsv"), '\t')
    @test data == Any["scan time [minute]" "intensity"; 1 123; 2 224]
end