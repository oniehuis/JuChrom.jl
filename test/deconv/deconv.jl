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
    if Threads.nthreads() > 1
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
    @test isnothing(σ)
    @test n == 0

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


# ############################################################################################
# # DeconvolutionWindow
# ############################################################################################
# @testset "DeconvolutionWindow" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     @test isa(dw.chromatogram, AbstractChromMS)
#     @test isa(dw.δtᵢ, Function)
#     @test isa(dw.ionindex, Int)
#     @test dw.ionindex == 1
#     @test isa(dw.maxindexrange, UnitRange{Int})
#     @test dw.maxindexrange == 4:4
#     @test isa(dw.leftindex, Int)
#     @test dw.leftindex == 1
#     @test isa(dw.rightindex, Int)
#     @test dw.rightindex == 7
#     @test isa(dw.leftlowindex, Int)
#     @test dw.leftlowindex == 1
#     @test isa(dw.rightlowindex, Int)
#     @test dw.rightlowindex == 7
# end


# ############################################################################################
# # Direction()
# ############################################################################################
# @testset "Reverse()" begin
#     @test isa(JuChrom.Deconvolution.Direction, Any)
# end


# ############################################################################################
# # Forward()
# ############################################################################################
# @testset "Forward()" begin
#     @test isa(JuChrom.Deconvolution.Forward(), JuChrom.Deconvolution.Direction)
# end


# ############################################################################################
# # Reverse()
# ############################################################################################
# @testset "Reverse()" begin
#     @test isa(JuChrom.Deconvolution.Reverse(), JuChrom.Deconvolution.Direction)
# end


# ############################################################################################
# # onestep()
# ############################################################################################
# @testset "onestep()" begin
#     @test JuChrom.Deconvolution.onestep(JuChrom.Deconvolution.Forward()) == 1
#     @test JuChrom.Deconvolution.onestep(JuChrom.Deconvolution.Reverse()) == -1
# end


# ############################################################################################
# # chromatogram(dw)
# ############################################################################################
# @testset "chromatogram(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.chromatogram(dw) == chrom
# end


# ############################################################################################
# # δtᵢ(dw)
# ############################################################################################
# @testset "δtᵢ(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.δtᵢ(dw) == δtᵢ
# end


# ############################################################################################
# # ionindex(dw)
# ############################################################################################
# @testset "ionindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     iᵢ = 1
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, iᵢ, 4:4, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.ionindex(dw) == iᵢ
# end


# ############################################################################################
# # maxindexrange(dw)
# ############################################################################################
# @testset "maxindexrange(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 5, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     lm = 4:5
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, lm, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.maxindexrange(dw) == lm
# end


# ############################################################################################
# # maxstartindex(dw)
# ############################################################################################
# @testset "maxstartindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 5, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     lm = 4:5
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, lm, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.maxstartindex(dw) == first(lm)
# end


# ############################################################################################
# # maxstopindex(dw)
# ############################################################################################
# @testset "maxstopindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 5, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     lm = 4:5
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, lm, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.maxstopindex(dw) == last(lm)
# end


# ############################################################################################
# # maxintensity(dw)
# ############################################################################################
# @testset "maxintensity(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     lm = 4:4
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, lm, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.maxintensity(dw) == 5
# end


# ############################################################################################
# # maxsize(dw)
# ############################################################################################
# @testset "maxsize(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     lm = 4:4
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, lm, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.maxsize(dw) == 1
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 5, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     lm = 4:5
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, lm, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.maxsize(dw) == 2
# end


# ############################################################################################
# # leftindex(dw)
# ############################################################################################
# @testset "leftindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     iₗ = 1
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, iₗ, 7, 1, 7)
#     @test JuChrom.Deconvolution.leftindex(dw) == iₗ
# end


# ############################################################################################
# # rightindex(dw)
# ############################################################################################
# @testset "rightindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     iᵣ = 7
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, iᵣ, 1, 7)
#     @test JuChrom.Deconvolution.rightindex(dw) == iᵣ
# end


# ############################################################################################
# # windowsize(dw)
# ############################################################################################
# @testset "windowsize(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.windowsize(dw) == 7
# end


# ############################################################################################
# # windowindices(dw)
# ############################################################################################
# @testset "windowindices(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.windowindices(dw) == 1:7
# end


# ############################################################################################
# # currentindex(dw, ::Direction)
# ############################################################################################
# @testset "currentindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     @test JuChrom.Deconvolution.currentindex(dw, JuChrom.Deconvolution.Reverse()) == 1
#     @test JuChrom.Deconvolution.currentindex(dw, JuChrom.Deconvolution.Forward()) == 7
# end


# ############################################################################################
# # currentindex!(dw, ::Direction, i::Integer)
# ############################################################################################
# @testset "currentindex!(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     JuChrom.Deconvolution.currentindex!(dw, JuChrom.Deconvolution.Reverse(), 2)
#     @test JuChrom.Deconvolution.currentindex(dw, JuChrom.Deconvolution.Reverse()) == 2
#     JuChrom.Deconvolution.currentindex!(dw, JuChrom.Deconvolution.Forward(), 6)
#     @test JuChrom.Deconvolution.currentindex(dw, JuChrom.Deconvolution.Forward()) == 6
# end


# ############################################################################################
# # leftlowindex(dw)
# ############################################################################################
# @testset "leftlowindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     iₗₗ = 2
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, iₗₗ, 7)
#     JuChrom.Deconvolution.leftlowindex(dw) == iₗₗ
# end


# ############################################################################################
# # rightlowindex(dw)
# ############################################################################################
# @testset "leftlowindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     iᵣₗ = 6
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, iᵣₗ)
#     JuChrom.Deconvolution.rightlowindex(dw) == iᵣₗ
# end


# ############################################################################################
# # leftlowindex!(dw, i::Integer)
# ############################################################################################
# @testset "leftlowindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     JuChrom.Deconvolution.leftlowindex!(dw, 2)
#     JuChrom.Deconvolution.leftlowindex(dw) == 2
# end


# ############################################################################################
# # rightlowindex!(dw, i::Integer)
# ############################################################################################
# @testset "leftlowindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     JuChrom.Deconvolution.rightlowindex!(dw, 6)
#     JuChrom.Deconvolution.rightlowindex(dw) == 6
# end


# ############################################################################################
# # lowindex(dw::DeconvolutionWindow, ::Direction)
# ############################################################################################
# @testset "leftlowindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     iₗₗ = 2
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, iₗₗ, 7)
#     @test JuChrom.Deconvolution.lowindex(dw, JuChrom.Deconvolution.Reverse()) == iₗₗ
#     iᵣₗ = 6
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, iᵣₗ)
#     @test JuChrom.Deconvolution.lowindex(dw, JuChrom.Deconvolution.Forward()) == iᵣₗ
# end


# ############################################################################################
# # lowindex(dw::DeconvolutionWindow, ::Direction)
# ############################################################################################
# @testset "leftlowindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     JuChrom.Deconvolution.lowindex!(dw, JuChrom.Deconvolution.Reverse(), 2)
#     @test JuChrom.Deconvolution.lowindex(dw, JuChrom.Deconvolution.Reverse()) == 2
#     JuChrom.Deconvolution.lowindex!(dw, JuChrom.Deconvolution.Forward(), 6)
#     @test JuChrom.Deconvolution.lowindex(dw, JuChrom.Deconvolution.Forward()) == 6
# end


# ############################################################################################
# # lowintensity(dw::DeconvolutionWindow, direction::Direction)
# ############################################################################################
# @testset "leftlowindex(dw)" begin
#     chrom = ChromMS((1:7)u"s", [29], reshape([0, 1, 2, 5, 3, 2, 1], (:,1)))
#     δtᵢ = ionscantimeshift(chrom, LinearDescending())
#     dw = JuChrom.Deconvolution.DeconvolutionWindow(chrom, δtᵢ, 1, 4:4, 1, 7, 1, 7)
#     JuChrom.Deconvolution.lowintensity(dw, JuChrom.Deconvolution.Reverse()) == 0
#     JuChrom.Deconvolution.lowintensity(dw, JuChrom.Deconvolution.Forward()) == 1
# end
