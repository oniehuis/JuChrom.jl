using JuChrom
using Test
using Unitful: 𝐓

@testset "gcms" begin
    @test isa(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), AbstractChromatogram)
    @test isa(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), AbstractGCMS)

    # # Preserve the provided argument types
    # @test (typeof(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1])) ==
    #     GCMS{Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓 , nothing}}}, Vector{Int64}, Matrix{Int64}, Nothing})
    # @test (typeof(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1])) ==
    #     GCMS{Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓 , nothing}}}, Vector{Int64}, Matrix{Int64}, Nothing})
    # @test (typeof(GCMS([1, 2, 3]u"s", [85.0, 100.0], [0 12; 34 956; 23 1])) ==
    #     GCMS{Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓 , nothing}}}, Vector{Float64}, Matrix{Int64}, Nothing})
    # @test (typeof(GCMS([1, 2, 3]u"s", [85, 100], [0.0 12.0; 34.0 956.0; 23.0 1.0])) ==
    #     GCMS{Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓 , nothing}}}, Vector{Int64}, Matrix{Float64}, Nothing})
    # @test (typeof(GCMS([1, 2, 3]u"s", [85, 100], [0.0 12.0; 34.0 956.0; 23.0 1.0], nothing)) ==
    #     GCMS{Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓 , nothing}}}, Vector{Int64}, Matrix{Float64}, Nothing})
    # @test (typeof(GCMS([1, 2, 3]u"s", [85, 100], [0.0 12.0; 34.0 956.0; 23.0 1.0], "sample name")) ==
    #     GCMS{Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓 , nothing}}}, Vector{Int64}, Matrix{Float64}, String})
    # @test (typeof(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1], SubString("sample name"))) ==
    #     GCMS{Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓 , nothing}}}, Vector{Int64}, Matrix{Int64},  SubString{String}})

    # # Scan time vector must contain quantities in time units
    # @test_throws MethodError GCMS([1, 2, 3], [85, 100], [0 12; 34 956; 23 1])
    # @test_throws MethodError GCMS([1, 2, 3]u"m", [85, 100], [0 12; 34 956; 23 1])

    # # Ion vector must have no a unit
    # @test_throws MethodError GCMS([1, 2, 3], [85, 100]u"s", [0 12; 34 956; 23 1])
    
    # # Intensity matrix must have no a unit    
    # @test_throws MethodError GCMS([1, 2, 3], [85, 100], [0 12; 34 956; 23 1]u"s")

    # # Intensity matrix must have no a unit    
    # @test_throws MethodError GCMS([1, 2, 3], [85, 100], [0 12; 34 956; 23 1]u"s")

    # # Source cannot something else than an AbstractString object or nothing
    # @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], [0.0 12.0; 34.0 956.0; 23.0 1.0], :sample_name)
    # @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], [0.0 12.0; 34.0 956.0; 23.0 1.0], 1)
    # @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], [0.0 12.0; 34.0 956.0; 23.0 1.0], 1.0)
    # @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], [0.0 12.0; 34.0 956.0; 23.0 1.0], 'S')

    # # Too many scan times in comparison to intensity matrix
    # @test_throws DimensionMismatch GCMS([1, 2, 3, 4]u"s", [85, 100], [0 12; 34 956; 23 1])
    # # Too few scan times in comparison to intensity matrix
    # @test_throws DimensionMismatch GCMS([1, 2]u"s", [85, 100], [0 12; 34 956; 23 1])
    # # Too many ions in comparison to intensity matrix
    # @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85, 100, 200], [0 12; 34 956; 23 1])
    # # Too few ions in comparison to intensity matrix
    # @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85], [0 12; 34 956; 23 1])
    # # Too many intensities in comparison to intensity matrix
    # @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85, 100], [0 12 43; 34 956 65; 23 1 73])
    # # Too few intensities in comparison to intensity matrix
    # @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85, 100, 200], [0 12; 34 956; 23 1])
    # # Too many intensities in comparison to scan times
    # @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1; 40 98])
    # # Too few intensities in comparison to scan times
    # @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956])
    # # Scan times not in ascending order
    # @test_throws ArgumentError GCMS([2, 1, 3]u"s", [85, 100], [0 12; 34 956; 23 1])
    # # Ions not in ascending order
    # @test_throws ArgumentError GCMS([1, 2, 3]u"s", [100, 85], [0 12; 34 956; 23 1])
    # # Negative intensity value
    # @test_throws ArgumentError GCMS([1, 2, 3]u"s", [85, 100], [-1 12; 34 956; 23 1])

    # # Testing getters
    # @test [85, 100] == ions(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))

    # # @test [0 12; 34 956; 23 1] == intensities(gcms)

    # @test [1, 2, 3]u"s" == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    # @test [1000, 2000, 3000]u"ms" == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"ms")
    # @test [1//60, 1//30, 1//20]u"minute" == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")
    # @test [ 1//3600, 1//1800, 1//1200]u"hr" == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"hr")

    # @test [1000.0, 2000.0, 3000.0]u"ms" == scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"ms")
    # @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")
    # @test [1/3600, 1/1800, 1/1200]u"hr" ≈ scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"hr")
    # @test [1, 2, 3] == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), ustripped=true)
    # @test [1.0, 2.0, 3.0] == scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), ustripped=true)
    # @test [1000, 2000, 3000] == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"ms", ustripped=true)
    # @test [1000.0, 2000.0, 3000.0] == scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"ms", ustripped=true)

    # @test (typeof(scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))) ==
    #     Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓 , nothing}}})
    # @test (typeof(scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]))) ==
    #     Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 𝐓 , nothing}}})
    # @test (typeof(scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")) ==
    #     Vector{Quantity{Rational{Int64}, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓 , nothing}}})
    # @test (typeof(scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"minute")) ==
    #     Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓 , nothing}}})
    # @test (typeof(scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"hr")) ==
    #     Vector{Quantity{Rational{Int64}, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Hour, 𝐓}(0, 1//1),), 𝐓 , nothing}}})
    # @test (typeof(scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"hr")) ==
    #     Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Hour, 𝐓}(0, 1//1),), 𝐓 , nothing}}})

    # @test (typeof(scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), timeunit=u"ms")) ==
    #     Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1000),), 𝐓 , nothing}}})


    # @test isnothing(source(gcms))
    # @test dimension(eltype(scantimes(gcms))) == 𝐓
    # @test unit(eltype(scantimes(gcms))) == u"s"

end