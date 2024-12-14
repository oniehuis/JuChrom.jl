using Test
using Unitful
using XLSX
using JuChrom

@testset "Excel AbstractChrom" begin
    rootpath = joinpath(splitpath(pathof(JuChrom))[begin:end-2])
    tmppath = joinpath(rootpath, "tmp")
    mkpath(tmppath)
    file = joinpath(tmppath, "./exceltest")
    chrom = Chrom((1:2)u"minute", [123, 224])
    exportdata(chrom, file, Excel(), overwrite=true)
    
    xf = XLSX.readxlsx(string(file, ".xlsx"))
    sheet = xf[first(XLSX.sheetnames(xf))]
    data = sheet["A1:B3"]
    @test data == Any["scan time [minute]" "intensity"; 1 123; 2 224]

    exportdata(chrom, file, Excel(sheetname="TIC"), overwrite=true)
    data = XLSX.readdata(string(file, ".xlsx"), "TIC", "A1:B3")
    @test data == Any["scan time [minute]" "intensity"; 1 123; 2 224]
end


@testset "Excel AbstractChromMS" begin
    rootpath = joinpath(splitpath(pathof(JuChrom))[begin:end-2])
    tmppath = joinpath(rootpath, "tmp")
    mkpath(tmppath)
    file = joinpath(tmppath, "./exceltest")
    chrom = ChromMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    exportdata(chrom, file, Excel(), overwrite=true)
    
    xf = XLSX.readxlsx(string(file, ".xlsx"))
    sheet = xf[first(XLSX.sheetnames(xf))]
    data = sheet["A1:C4"]
    @test data == Any["scan times [s]" "m/z 85" "m/z 100"; 1 0 12; 2 34 956; 3 23 1]

    exportdata(chrom, file, Excel(sheetname="GCMS"), overwrite=true)
    data = XLSX.readdata(string(file, ".xlsx"), "GCMS", "A1:C4")
    @test data == Any["scan times [s]" "m/z 85" "m/z 100"; 1 0 12; 2 34 956; 3 23 1]
end