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
 