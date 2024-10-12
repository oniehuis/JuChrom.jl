using JuChrom
using Test
using Unitful: 𝐓


############################################################################################
# Base.iterate(LM::LocalMaxima, currentindex=LM.startindex)
############################################################################################
@testset "LocalMaxima" begin
    maxima = []
    for lm in JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7])
        push!(maxima, lm)
    end
    @test [3:3, 5:7] == maxima

    empty!(maxima)
    for lm in JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, stopindex=4)
        push!(maxima, lm)
    end
    @test [3:3] == maxima
end


############################################################################################
# JuChrom.LocalMaxima(values::AbstractVector{<:Real}; 
# startindex::Integer=firstindex(values), stopindex::Integer=lastindex(values))
############################################################################################
@testset "LocalMaxima" begin
    @test [4, 3, 5, 3, 6, 6, 6, 4, 7] == JuChrom.LocalMaxima(
        [4, 3, 5, 3, 6, 6, 6, 4, 7]).values
    @test 1 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7]).startindex
    @test 9 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7]).stopindex
    @test 9 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7]).length
    @test [4, 3, 5, 3, 6, 6, 6, 4, 7] == JuChrom.LocalMaxima(
        [4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, stopindex=8).values
    @test 2 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, 
        stopindex=8).startindex
    @test 8 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, 
        stopindex=8).stopindex
    @test 9 == JuChrom.LocalMaxima([4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=2, 
        stopindex=8).length

    @test Vector{Int} == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7]).values)
    @test Int == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7]).startindex)
    @test UInt8 == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0x1).startindex)
    @test UInt8 == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=0x9).stopindex)
    @test Int == typeof(JuChrom.LocalMaxima(Int[4, 3, 5, 3, 6, 6, 6, 4, 7]).length)

    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=10)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=0)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=10)
    @test_throws ArgumentError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=1, stopindex=2)
    @test_throws ArgumentError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=8)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0, stopindex=1)
    @test_throws BoundsError JuChrom.LocalMaxima(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=9, stopindex=10)
 end


############################################################################################
# JuChrom.nextlocalmaximum(values::AbstractVector{<:Real}; 
# startindex::Integer=firstindex(values), stopindex::Integer=lastindex(values))
############################################################################################
@testset "nextlocalmaximum" begin
    @test 3:3 == JuChrom.nextlocalmaximum(Int[4, 3, 5, 3, 6, 6, 6, 4, 7])
    @test 3:3 == JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7])
    @test 5:7 == JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], startindex=4)
    @test nothing === JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=4, stopindex=7)
    @test 3:3 === JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=2, stopindex=7)

    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=10)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=0)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        stopindex=10)
    @test_throws ArgumentError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=1, stopindex=2)
    @test_throws ArgumentError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=8)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=0, stopindex=1)
    @test_throws BoundsError JuChrom.nextlocalmaximum(Float64[4, 3, 5, 3, 6, 6, 6, 4, 7], 
        startindex=9, stopindex=10)
 end
 

############################################################################################
# JuChrom.lsfit(xs, ys)
############################################################################################
@testset "lsfit" begin
    # Basic linear relationship
    xs1 = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys1 = [2.0, 3.0, 4.0, 5.0, 6.0]
    intercept1, slope1 = JuChrom.Deconvolution.lsfit(xs1, ys1)
    @test intercept1 ≈ 1.0 atol=1e-10
    @test slope1 ≈ 1.0 atol=1e-10

    # Negative slope
    xs2 = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys2 = [5.0, 4.0, 3.0, 2.0, 1.0]
    intercept2, slope2 = JuChrom.Deconvolution.lsfit(xs2, ys2)
    @test intercept2 ≈ 6.0 atol=1e-10
    @test slope2 ≈ -1.0 atol=1e-10

    # Zero slope (horizontal line)
    xs3 = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys3 = [3.0, 3.0, 3.0, 3.0, 3.0]
    intercept3, slope3 = JuChrom.Deconvolution.lsfit(xs3, ys3)
    @test intercept3 ≈ 3.0 atol=1e-10
    @test slope3 ≈ 0.0 atol=1e-10

    # Insufficient elements in vectors (should trigger an error)
    xs4 = [1.0]
    ys4 = [2.0]
    @test_throws ArgumentError JuChrom.Deconvolution.lsfit(xs4, ys4)

    # Vectors of different lengths (should throw an error)
    xs5 = [1.0, 2.0, 3.0]
    ys5 = [1.0, 2.0]
    @test_throws ArgumentError JuChrom.Deconvolution.lsfit(xs5, ys5)
end


############################################################################################
# JuChrom.stddev(chrom::AbstractChromMS, ionindices; windowsize::Integer=13, 
#   threshold::Real=0)
############################################################################################
@testset "JuChrom.stddev(chrom, ionindices)" begin
    # Ensure the function works on basic, well-defined input
    chrom = ChromMS((1:5)u"s", [84.8, 85.2, 100.9], [1 2 4; 2 5 1; 3 1 2; 1 3 3; 3 5 4])
    ionindices = [1, 2]
    result = JuChrom.stddev(chrom, ionindices, windowsize=5, threshold=0)
    @test result ≈ [0.7071067811865475, 1.1547005383792517]

    # Ensure the function works on real input
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    chrom = binions(importdata(dfolder, ChemStationMS()))
    ionindices = [2]
    result = JuChrom.stddev(chrom, ionindices, windowsize=13, threshold=0)
    @test result ≈ [1.1470786693528088, 1.2792042981336627, 1.6431676725154982, 
        1.3764944032233704, 0.7844645405527362, 1.1766968108291043, 0.4082482904638631]

    # Test with no transition
    chrom = ChromMS((1:5)u"s", [84.8, 85.2], [1 2; 1 2; 1 2; 1 2; 1 2])
    ionindices = [1, 2]
    result = JuChrom.stddev(chrom, ionindices, windowsize=5, threshold=0)
    @test result == Float64[]

    # Test with window size larger than scan count
    chrom = ChromMS((1:5)u"s", [84.8, 85.2, 100.9], [1 2 4; 2 5 1; 3 1 2; 1 3 3; 3 5 4])
    ionindices = [1, 2]
    @test_throws ArgumentError JuChrom.stddev(chrom, ionindices, windowsize=6, threshold=0)

    # Test with window size samller 5 scans
    chrom = ChromMS((1:5)u"s", [84.8, 85.2, 100.9], [1 2 4; 2 5 1; 3 1 2; 1 3 3; 3 5 4])
    ionindices = [1, 2]
    @test_throws ArgumentError JuChrom.stddev(chrom, ionindices, windowsize=4, threshold=0)

    # Test with threshold cutoff
    chrom = ChromMS((1:5)u"s", [84.8, 85.2, 100.9], [1 2 4; 2 5 1; 3 1 2; 1 3 3; 3 5 4])
    ionindices = [2]
    result1 = JuChrom.stddev(chrom, ionindices, windowsize=5, threshold=0)
    @test result1 ≈ [1.1547005383792517]
    result2 = JuChrom.stddev(chrom, ionindices, windowsize=5, threshold=1)
    @test result2 == Float64[]

    # Test with less transitions than half the window size
    chrom = ChromMS((1:13)u"s", [1, 2], [1 1; 2 2; 3 3; 1 1; 1 1; 3 3; 3 3; 1 1; 1 1; 3 3; 
        3 3; 1 1; 1 1])
    ionindices = [1]
    result1 = JuChrom.stddev(chrom, ionindices, windowsize=13, threshold=0)
    @test result2 == Float64[]

    # Test with three consecutive scans having an intensity above (or below) the mean ints.
    chrom = ChromMS((1:13)u"s", [1, 2], [1 1; 2 2; 3 3; 1 1; 1 1; 1 1; 3 3; 1 1; 3 3; 1 1; 
        3 3; 1 1; 3 3])
    ionindices = [1]
    result1 = JuChrom.stddev(chrom, ionindices, windowsize=13, threshold=0)
    @test result2 == Float64[]

    # Edge case: empty input for ionindices
    chrom = ChromMS((1:5)u"s", [84.8, 85.2, 100.9], [1 2 4; 2 5 1; 3 1 2; 1 3 3; 3 5 4])
    ionindices = []
    result1 = JuChrom.stddev(chrom, ionindices, windowsize=5, threshold=0)
    @test result2 == Float64[]
end


############################################################################################
# stddev(chrom::AbstractChromMS; windowsize::Integer=13, threshold::Real=0, 
#   nthreads::Integer)
############################################################################################
@testset "JuChrom.stddev(chrom)" begin

    # Ensure the function operates correctly with real input, supporting multi-threading
    if Threads.nthreads > 1
        dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
        chrom = binions(importdata(dfolder, ChemStationMS()))
        σ, n = JuChrom.stddev(chrom, windowsize=13, threshold=0)
        @test σ ≈ 1.989116630064713
        @test n == 9874
    end

    # Ensure the function works correctly with real input on a single thread
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    chrom = binions(importdata(dfolder, ChemStationMS()))
    σ, n = JuChrom.stddev(chrom, windowsize=13, threshold=0, nthreads=1)
    @test σ ≈ 1.989116630064713
    @test n == 9874

    # Ensure the function throws an error if the number of threads is less than one
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    chrom = binions(importdata(dfolder, ChemStationMS()))
    @test_throws ArgumentError JuChrom.stddev(chrom, windowsize=13, threshold=0, nthreads=0)

    # Ensure the function throws an error if the number of threads exceeds the maximum 
    # available
    dfolder = joinpath(JuChrom.agilent, "C7-C40_ChemStationMS.D")
    chrom = binions(importdata(dfolder, ChemStationMS()))
    @test_throws ArgumentError JuChrom.stddev(chrom, windowsize=13, threshold=0, 
        nthreads=Threads.nthreads()+1)

    # Ensure the function works on basic, well-defined input
    chrom = ChromMS((1:5)u"s", [84.8, 85.2, 100.9], [1 2 4; 2 5 1; 3 1 2; 1 3 3; 3 5 4])
    σ, n = JuChrom.stddev(chrom, windowsize=5, threshold=0)
    @test σ ≈ 1.3801577659941269
    @test n == 2

    # Test with no transition
    chrom = ChromMS((1:5)u"s", [84.8, 85.2], [1 2; 1 2; 1 2; 1 2; 1 2])
    σ, n = JuChrom.stddev(chrom, windowsize=5, threshold=0)
    @test σ ≈ 1.3801577659941269
    @test n == 2

    # Test with window size larger than scan count
    chrom = ChromMS((1:5)u"s", [84.8, 85.2, 100.9], [1 2 4; 2 5 1; 3 1 2; 1 3 3; 3 5 4])
    @test_throws TaskFailedException JuChrom.stddev(chrom, windowsize=6, threshold=0)

    # Test with window size samller 5 scans
    chrom = ChromMS((1:5)u"s", [84.8, 85.2, 100.9], [1 2 4; 2 5 1; 3 1 2; 1 3 3; 3 5 4])
    @test_throws TaskFailedException JuChrom.stddev(chrom, windowsize=4, threshold=0)

    # Test with threshold cutoff
    chrom = ChromMS((1:5)u"s", [84.8, 85.2, 100.9], [1 2 4; 2 5 1; 3 1 2; 1 3 3; 3 5 4])
    σ, n = JuChrom.stddev(chrom, windowsize=5, threshold=0)
    @test σ ≈ 1.3801577659941269
    @test n == 2
    σ, n = JuChrom.stddev(chrom, windowsize=5, threshold=1)
    @test isnothing(σ)
    @test n == 0

    # Test with less transitions than half the window size
    chrom = ChromMS((1:13)u"s", [1, 2], [1 1; 2 2; 3 3; 1 1; 1 1; 3 3; 3 3; 1 1; 1 1; 3 3; 
        3 3; 1 1; 1 1])
    σ, n = JuChrom.stddev(chrom, windowsize=13, threshold=0)
    @test isnothing(σ)
    @test n == 0

    # Test with three consecutive scans having an intensity above (or below) the mean ints.
    chrom = ChromMS((1:13)u"s", [1, 2], [1 1; 2 2; 3 3; 1 1; 1 1; 1 1; 3 3; 1 1; 3 3; 1 1; 
        3 3; 1 1; 3 3])
    result1 = JuChrom.stddev(chrom, windowsize=13, threshold=0)
    @test isnothing(σ)
    @test n == 0
end