using Test
using Unitful
using XLSX
using JuChrom

@testset "Excel" begin
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