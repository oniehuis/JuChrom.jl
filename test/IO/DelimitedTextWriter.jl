using JuChrom
using Test
using DelimitedFiles

@testset "DelimitedTextWriter AbstractChrom" begin
    rootpath = joinpath(splitpath(pathof(JuChrom))[begin:end-2])
    tmppath = joinpath(rootpath, "tmp")
    mkpath(tmppath)
    file = joinpath(tmppath, "./delimtest")
    chrom = Chrom((1:2)u"minute", [123, 224])
    exportdata(chrom, file, DelimitedText(), overwrite=true)
    data = readdlm(string(file, ".tsv"), '\t')
    @test data == Any["scan time [minute]" "intensity"; 1 123; 2 224]

    exportdata(chrom, file, DelimitedText(delim=","), overwrite=true)
    data = readdlm(string(file, ".csv"), ',')
    @test data == Any["scan time [minute]" "intensity"; 1 123; 2 224]

    exportdata(chrom, file, DelimitedText(delim=";"), overwrite=true)
    data = readdlm(string(file, ".csv"), ';')
    @test data == Any["scan time [minute]" "intensity"; 1 123; 2 224]

    exportdata(chrom, file, DelimitedText(delim="*"), overwrite=true)
    data = readdlm(string(file, ".txt"), '*')
    @test data == Any["scan time [minute]" "intensity"; 1 123; 2 224]
end

@testset "DelimitedTextWriter AbstractChromMS" begin
    rootpath = joinpath(splitpath(pathof(JuChrom))[begin:end-2])
    tmppath = joinpath(rootpath, "tmp")
    mkpath(tmppath)
    file = joinpath(tmppath, "./delimtest")
    chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    exportdata(chrom, file, DelimitedText(), overwrite=true)
    data = readdlm(string(file, ".tsv"), '\t')
    @test data == Any["scan times [s]" "m/z 85" "m/z 100"; 1 0 12; 2 34 956; 3 23 1]

    exportdata(chrom, file, DelimitedText(delim=","), overwrite=true)
    data = readdlm(string(file, ".csv"), ',')
    @test data == Any["scan times [s]" "m/z 85" "m/z 100"; 1 0 12; 2 34 956; 3 23 1]

    exportdata(chrom, file, DelimitedText(delim=";"), overwrite=true)
    data = readdlm(string(file, ".csv"), ';')
    @test data == Any["scan times [s]" "m/z 85" "m/z 100"; 1 0 12; 2 34 956; 3 23 1]

    exportdata(chrom, file, DelimitedText(delim="*"), overwrite=true)
    data = readdlm(string(file, ".txt"), '*')
    @test data == Any["scan times [s]" "m/z 85" "m/z 100"; 1 0 12; 2 34 956; 3 23 1]
end