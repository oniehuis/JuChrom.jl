using JuChrom
using Test
using Unitful: 𝐓

@testset "scantimes" begin
    # FID
    @test [1, 2, 3]u"s" == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test [1, 2, 3] == scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(FID([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true)

    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓 , nothing}}} == typeof(scantimes(FID(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓 , nothing}}} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23])))
    @test Vector{Quantity{Rational{Int64}, 𝐓 , Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓 , nothing}}} == typeof(scantimes(FID(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(FID(Int[1, 2, 3]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(FID(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), ustripped=true))

    # GCMS
    @test [1, 2, 3]u"s" == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test [1, 2, 3] == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute", ustripped=true)

    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓 , nothing}}} == typeof(scantimes(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓 , nothing}}} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1])))
    @test Vector{Quantity{Rational{Int64}, 𝐓 , Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓 , nothing}}} == typeof(scantimes(GCMS(Int64[1, 2, 3]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", 
        [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(GCMS(Int[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))
    @test Vector{Float64} == typeof(scantimes(GCMS(Float64[1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true))

    # TIC
    @test [1, 2, 3]u"s" == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute")
    @test [1, 2, 3] == scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), ustripped=true)
    @test [1/60, 1/30, 1/20] ≈ scantimes(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        timeunit=u"minute", ustripped=true)

    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),
        ), 𝐓 , nothing}}} == typeof(scantimes(TIC(Int64[1, 2, 3]u"s", [12, 956, 23])))
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1
        ),), 𝐓 , nothing}}} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23])))
    @test Vector{Quantity{Rational{Int64}, 𝐓 , Unitful.FreeUnits{ (Unitful.Unit{:Minute, 
        𝐓}(0, 1//1),), 𝐓 , nothing}}} == typeof(scantimes(TIC(Int64[1, 2, 3]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), timeunit=u"minute"))
    @test Vector{Int} == typeof(scantimes(TIC(Int[1, 2, 3]u"s", [12, 956, 23]), 
        ustripped=true))
    @test Vector{Float64} == typeof(scantimes(TIC(Float64[1.0, 2.0, 3.0]u"s", 
        [12, 956, 23]), ustripped=true))
end

@testset "intensities" begin
    # FID
    @test [12, 956, 23] == intensities(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test Vector{Int} == typeof(intensities(FID([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Vector{Float64} == typeof(intensities(FID([1, 2, 3]u"s", 
        Float64[12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(intensities(FID([1, 2, 3]u"s", 10:5:20)))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(FID([1, 2, 3]u"s", 10.0:5:20.0)))

    # GCMS
    @test [0 12; 34 956; 23 1] == intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]))
    @test reshape([0, 956, 1], (:,1)) == intensities(GCMS([1, 2, 3]u"s", [85],
        reshape([0, 956, 1], (:,1))))
    @test Matrix{Int64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Int64[0 12; 34 956; 23 1])))
    @test Matrix{Float64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1])))
    @test Matrix{Int64} == typeof(intensities(GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1],
         (:,1)))))

    # TIC
    @test [12, 956, 23] == intensities(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test Vector{Int} == typeof(intensities(TIC([1, 2, 3]u"s", Int[12, 956, 23])))
    @test Vector{Float64} == typeof(intensities(TIC([1, 2, 3]u"s", 
        Float64[12.0, 956.0, 23.0])))
    @test StepRange{Int64, Int64} == typeof(intensities(TIC([1, 2, 3]u"s", 10:5:20)))
    @test StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, 
        Int64} == typeof(intensities(TIC([1, 2, 3]u"s", 10.0:5:20.0)))
end

@testset "metadata" begin
    # FID
    @test Dict() == metadata(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test Dict(:id => 1, :name => "sample") == metadata(FID([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 1, :name => "sample")))
    @test Dict{Any, Any} == typeof(metadata(FID([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 1, :name => "sample"))))
    
    # GCMS
    @test Dict() == metadata(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test Dict(:id => 1, :name => "sample") == metadata(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 1, :name => "sample")))
    @test Dict{Any, Any} == typeof(metadata(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 1, :name => "sample"))))   

    # TIC
    @test Dict() == metadata(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test Dict(:id => 1, :name => "sample") == metadata(TIC([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 1, :name => "sample")))
    @test Dict{Any, Any} == typeof(metadata(TIC([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 1, :name => "sample"))))
end
