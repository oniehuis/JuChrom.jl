using JuChrom
using Test
using Unitful: 𝐓

@testset "GCMS" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(GCMS(Int64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(GCMS(Float64[1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).scantimes)

    # Ions
    @test [85, 100] == GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]).ions
    @test [85] == GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1], (:,1))).ions
    @test Vector{Int64} == typeof(GCMS([1, 2, 3]u"s", Int64[85, 100], 
        [0 12; 34 956; 23 1]).ions)
    @test Vector{Float64} == typeof(GCMS([1, 2, 3]u"s", Float64[85, 100], 
        [0 12; 34 956; 23 1]).ions)
    @test Vector{Int64} == typeof(GCMS([1, 2, 3]u"s", Int64[85], reshape([0, 956, 1], 
        (:,1))).ions)

    # Intensities
    @test [0 12; 34 956; 23 1] == GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).intensities
    @test reshape([0, 956, 1], (:,1)) == GCMS([1, 2, 3]u"s", [85],
        reshape([0, 956, 1], (:,1))).intensities
    @test Matrix{Int64} == typeof(GCMS([1, 2, 3]u"s", [85, 100], 
        Int64[0 12; 34 956; 23 1]).intensities)
    @test Matrix{Float64} == typeof(GCMS([1, 2, 3]u"s", [85, 100], 
        Float64[0 12; 34 956; 23 1]).intensities)
    @test Matrix{Int64} == typeof(GCMS([1, 2, 3]u"s", [85], reshape([0, 956, 1],
         (:,1))).intensities)

    # Metadata
    @test Dict{Any, Any} == typeof(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]).metadata)
    @test Dict() == GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]).metadata
    @test Dict(:id => 1, :name => "sample") == GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 1, :name => "sample")).metadata

    # Check the associated supertypes
    @test isa(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), AbstractChromatogram)
    @test isa(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), AbstractGCMS)

    # Scan time vector accepts only time values
    @test_throws MethodError GCMS([1, 2, 3], [85, 100], [0 12; 34 956; 23 1])
    @test_throws MethodError GCMS([1, 2, 3]u"m", [85, 100], [0 12; 34 956; 23 1])

    # Ion vector must not accept a unit
    @test_throws MethodError GCMS([1, 2, 3], [85, 100]u"s", [0 12; 34 956; 23 1])
    
    # Intensity matrix must not accept a unit
    @test_throws MethodError GCMS([1, 2, 3], [85, 100], [0 12; 34 956; 23 1]u"s")

    # Metadata cannot be anything other than dictionary
    @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], :sample_name)
    @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], 1)
    @test_throws MethodError GCMS([1, 2, 3]u"s", [85, 100], 
        [0.0 12.0; 34.0 956.0; 23.0 1.0], 'S')

    # Check size compatibility of intensities, scantimes, and ions
    # Too many scan times in comparison to intensity matrix
    @test_throws DimensionMismatch GCMS([1, 2, 3, 4]u"s", [85, 100], [0 12; 34 956; 23 1])
    # Too few scan times in comparison to intensity matrix
    @test_throws DimensionMismatch GCMS([1, 2]u"s", [85, 100], [0 12; 34 956; 23 1])
    # Too many ions in comparison to intensity matrix
    @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85, 100, 200], 
        [0 12; 34 956; 23 1])
    # Too few ions in comparison to intensity matrix
    @test_throws DimensionMismatch GCMS([1, 2, 3]u"s", [85], [0 12; 34 956; 23 1])

    # Scan times must be in ascending order
    @test_throws ArgumentError GCMS([2, 1, 3]u"s", [85, 100], [0 12; 34 956; 23 1])

    # Ions must be in ascending order
    @test_throws ArgumentError GCMS([1, 2, 3]u"s", [100, 85], [0 12; 34 956; 23 1])

    # Intensity values cannot be less than zero
    @test_throws ArgumentError GCMS([1, 2, 3]u"s", [85, 100], [-1 12; 34 956; 23 1])

    # GCMS is broadcastable
    @test [1, 2]u"s" == ((gcms, i) -> gcms.scantimes[i]).(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), [1, 2])

    # Testing getters
    # Scantimes
    @test [85, 100] == ions(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))

    @test [1, 2, 3]u"s" == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test [1000, 2000, 3000]u"ms" == scantimes(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"ms")
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test [1/3600, 1/1800, 1/1200]u"hr" ≈ scantimes(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"hr")

    @test [1000.0, 2000.0, 3000.0]u"ms" == scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"ms")
    @test [1/60, 1/30, 1/20]u"minute" ≈ scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"minute")
    @test [1/3600, 1/1800, 1/1200]u"hr" ≈ scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"hr")
    @test [1, 2, 3] == scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        ustripped=true)
    @test [1.0, 2.0, 3.0] == scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), ustripped=true)
    @test [1000, 2000, 3000] == scantimes(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"ms", ustripped=true)
    @test [1000.0, 2000.0, 3000.0] == scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), timeunit=u"ms", ustripped=true)

    @test (typeof(scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))) ==
        Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 1//1),), 
        𝐓 , nothing}}})
    @test (typeof(scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]))) ==
        Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}})
    @test (typeof(scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute")) == Vector{Quantity{Rational{Int64}, 𝐓 , Unitful.FreeUnits{
        (Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓 , nothing}}})
    @test (typeof(scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"minute")) == Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{
        (Unitful.Unit{:Minute, 𝐓}(0, 1//1),), 𝐓 , nothing}}})
    @test (typeof(scantimes(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"hr")) == Vector{Quantity{Rational{Int64}, 𝐓 , Unitful.FreeUnits{
        (Unitful.Unit{:Hour, 𝐓}(0, 1//1),), 𝐓 , nothing}}})
    @test (typeof(scantimes(GCMS([1.0, 2.0, 3.0]u"s", [85, 100], [0 12; 34 956; 23 1]), 
        timeunit=u"hr")) == Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{
        (Unitful.Unit{:Hour, 𝐓}(0, 1//1),), 𝐓 , nothing}}})

    # Ions
    @test Int64[85, 100] == ions(GCMS([1, 2, 3]u"s", Int64[85, 100], [0 12; 34 956; 23 1]))
    @test Float64[85, 100] == ions(GCMS([1, 2, 3]u"s", Float64[85.0, 100.0], 
        [0 12; 34 956; 23 1]))
    @test_throws MethodError ions(TIC([1, 2, 3]u"s", [0, 34, 956]))

    # Intensities
    @test Int64[12, 956, 23] == intensities(FID([1, 2, 3]u"s",
        Int64[12, 956, 23]))
    @test Float64[12.0, 956.0, 23.0] == intensities(FID([1, 2, 3]u"s",
        Float64[12.0, 956.0, 23.0]))
    @test Int64[12, 956, 23] == intensities(TIC([1, 2, 3]u"s",
        Int64[12, 956, 23]))
    @test Float64[12.0, 956.0, 23.0] == intensities(TIC([1, 2, 3]u"s",
        Float64[12.0, 956.0, 23.0]))
    @test Int64[0 12; 34 956; 23 1] == intensities(GCMS([1, 2, 3]u"s", [85, 100],
        Int64[0 12; 34 956; 23 1]))
    @test Float64[0 12; 34 956; 23 1] == intensities(GCMS([1, 2, 3]u"s", [85, 100],
        Float64[0 12; 34 956; 23 1]))
        
    @test Dict() == metadata(FID([1, 2, 3]u"s", [12, 956, 23]))
    @test Dict(:id => 1, :name => "sample") == metadata(FID([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 1, :name => "sample")))
    @test Dict() == metadata(GCMS([1, 2, 3]u"s", [85, 100], [0 12; 34 956; 23 1]))
    @test Dict(:id => 1, :name => "sample") == metadata(GCMS([1, 2, 3]u"s", [85, 100], 
        [0 12; 34 956; 23 1], Dict(:id => 1, :name => "sample")))
    @test Dict() == metadata(TIC([1, 2, 3]u"s", [12, 956, 23]))
    @test Dict(:id => 1, :name => "sample") == metadata(TIC([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 1, :name => "sample")))

end

@testset "FID" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == FID([1, 2, 3]u"s", [12, 956, 23]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(FID(Int64[1, 2, 3]u"s", [12, 956, 23]).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(FID(Float64[1, 2, 3]u"s", 
            [12, 956, 23]).scantimes)

    # Intensities
    @test [12, 956, 23] == FID([1, 2, 3]u"s", [12, 956, 23]).intensities
    @test Vector{Int64} == typeof(FID([1, 2, 3]u"s", Int64[12, 956, 23]).intensities)
    @test Vector{Float64} == typeof(FID([1, 2, 3]u"s",  Float64[12, 956, 23]).intensities)

    # Metadata
    @test Dict{Any, Any} == typeof(FID([1, 2, 3]u"s", [12, 956, 23]).metadata)
    @test Dict() == FID([1, 2, 3]u"s", [12, 956, 23]).metadata
    @test Dict(:id => 1, :name => "sample") == FID([1, 2, 3]u"s", [12, 956, 23], 
        Dict(:id => 1, :name => "sample")).metadata

    # Check the associated supertypes
    @test isa(FID([1, 2, 3]u"s", [12, 956, 23]), AbstractChromatogram)
    @test isa(FID([1, 2, 3]u"s", [12, 956, 23]), AbstractGC)
    @test isa(FID([1, 2, 3]u"s", [12, 956, 23]), AbstractFID)

    # Scan time vector accepts only time values
    @test_throws MethodError FID([1, 2, 3], [12, 956, 23])
    @test_throws MethodError FID([1, 2, 3]u"m", [12, 956, 23])
    
    # Intensity vector must not accept a unit
    @test_throws MethodError FID([1, 2, 3], [12, 956, 23]u"s")

    # Metadata cannot be anything other than a dictionary
    @test_throws MethodError FID([1, 2, 3]u"s", [12, 956, 23], :sample_name)
    @test_throws MethodError FID([1, 2, 3]u"s", [12, 956, 23], 1)
    @test_throws MethodError FID([1, 2, 3]u"s", [12, 956, 23], 'S')

    # Number of scantimes and intensities must be identical
    @test_throws DimensionMismatch FID([1, 2, 3, 4]u"s", [12, 956, 23])
    @test_throws DimensionMismatch FID([1, 2]u"s", [12, 956, 23])

    # Scan times must be in ascending order
    @test_throws ArgumentError FID([2, 1, 3]u"s", [12, 956, 23])

    # Intensity values cannot be less than zero
    @test_throws ArgumentError FID([1, 2, 3]u"s", [-1, 956, 23])

    # FID is broadcastable
    @test [1, 2]u"s" == ((fid, i) -> fid.scantimes[i]).(FID([1, 2, 3]u"s", [12, 956, 23]), 
        [1, 2])

end

@testset "TIC" begin
    # Verify object construction and field types depending on constructor arguments
    # Scantimes
    @test [1, 2, 3]u"s" == TIC([1, 2, 3]u"s", [12, 956, 23]).scantimes
    @test Vector{Quantity{Int64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(TIC(Int64[1, 2, 3]u"s", [12, 956, 23]).scantimes)
    @test Vector{Quantity{Float64, 𝐓 , Unitful.FreeUnits{(Unitful.Unit{:Second, 𝐓}(0, 
        1//1),), 𝐓 , nothing}}} == typeof(TIC(Float64[1, 2, 3]u"s", 
            [12, 956, 23]).scantimes)

    # Intensities
    @test [12, 956, 23] == TIC([1, 2, 3]u"s", [12, 956, 23]).intensities
    @test Vector{Int64} == typeof(TIC([1, 2, 3]u"s", Int64[12, 956, 23]).intensities)
    @test Vector{Float64} == typeof(TIC([1, 2, 3]u"s",  Float64[12, 956, 23]).intensities)

     # Metadata
     @test Dict{Any, Any} == typeof(TIC([1, 2, 3]u"s", [12, 956, 23]).metadata)
     @test Dict() == TIC([1, 2, 3]u"s", [12, 956, 23]).metadata
     @test Dict(:id => 1, :name => "sample") == TIC([1, 2, 3]u"s", [12, 956, 23], 
         Dict(:id => 1, :name => "sample")).metadata

    # Check the associated supertypes
    @test isa(TIC([1, 2, 3]u"s", [12, 956, 23]), AbstractChromatogram)
    @test isa(TIC([1, 2, 3]u"s", [12, 956, 23]), AbstractGC)
    @test isa(TIC([1, 2, 3]u"s", [12, 956, 23]), AbstractTIC)

    # Scan time vector accepts only time values
    @test_throws MethodError TIC([1, 2, 3], [12, 956, 23])
    @test_throws MethodError TIC([1, 2, 3]u"m", [12, 956, 23])
    
    # Intensity vector must not accept a unit
    @test_throws MethodError TIC([1, 2, 3], [12, 956, 23]u"s")

    # Metadata cannot be anything other than a dictionary.
    @test_throws MethodError TIC([1, 2, 3]u"s", [12, 956, 23], :sample_name)
    @test_throws MethodError TIC([1, 2, 3]u"s", [12, 956, 23], 1)
    @test_throws MethodError TIC([1, 2, 3]u"s", [12, 956, 23], 'S')

    # Number of scantimes and intensities must be identical
    @test_throws DimensionMismatch TIC([1, 2, 3, 4]u"s", [12, 956, 23])
    @test_throws DimensionMismatch TIC([1, 2]u"s", [12, 956, 23])

    # Scan times must be in ascending order
    @test_throws ArgumentError TIC([2, 1, 3]u"s", [12, 956, 23])

    # Intensity values cannot be less than zero
    @test_throws ArgumentError TIC([1, 2, 3]u"s", [-1, 956, 23])

    # FID is broadcastable
    @test [1, 2]u"s" == ((tic, i) -> tic.scantimes[i]).(TIC([1, 2, 3]u"s", [12, 956, 23]), 
        [1, 2])
end

@testset "showGCMS" begin
    io = IOBuffer()
    show(io, GCMS(Int64[1, 2, 3]u"s", Int64[85, 100], Int64[0 12; 34 956; 23 1], Dict(:id => 1, 
        :name => "sample")))
    @test String(take!(io)) == string(
        "GCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n", 
        "3 scans; time range: 1 s - 3 s\n",
        "2 ions; range: m/z 85 - 100\n",
        "intensity range: 0 - 956\n", 
        "metadata: :id, :name")
end

@testset "showFID" begin
    io = IOBuffer()
    show(io, FID(Int64[1, 2, 3]u"s", Int64[12, 956, 23], Dict(:id => 1, :name => "sample")))
    @test String(take!(io)) == string(
        "FID {scantimes: Int64, intensities: Int64}\n",
        "3 scans; time range: 1 s - 3 s\n",
        "intensity range: 12 - 956\n",
        "metadata: :id, :name")
end

@testset "showTIC" begin
    io = IOBuffer()
    show(io, TIC(Int64[1, 2, 3]u"s", Int64[12, 956, 23], Dict(:id => 1, :name => "sample")))
    @test String(take!(io)) == string(
        "TIC {scantimes: Int64, intensities: Int64}\n",
        "3 scans; time range: 1 s - 3 s\n",
        "intensity range: 12 - 956\n",
        "metadata: :id, :name")
end