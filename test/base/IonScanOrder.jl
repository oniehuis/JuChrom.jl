using JuChrom
using Test
using Unitful: 𝐓

############################################################################################
# LinearAscending(; start::Real=0, stop::Real=1) <: IonScanOrder
############################################################################################
@testset "LinearAscending" begin
    # Verify object construction and field types depending on constructor arguments
    @test 0 == LinearAscending().start
    @test Int == typeof(LinearAscending().start)
    @test 1 == LinearAscending().stop
    @test Int == typeof(LinearAscending().stop)
    @test 0.5 == LinearAscending(; start=0.5).start
    @test Float64 == typeof(LinearAscending(; start=0.5).start)
    @test 1 == LinearAscending(; start=0.5).stop
    @test 0 == LinearAscending(; stop=0.5).start
    @test Float64 == typeof(LinearAscending(; stop=0.5).stop)
    @test 0.5 == LinearAscending(; stop=0.5).stop
    @test 0.1 == LinearAscending(; start=0.1, stop=0.5).start
    @test 0.5 == LinearAscending(; start=0.1, stop=0.5).stop
    @test Float64 == typeof(LinearAscending(; start=0.1, stop=0.5).start)
    @test Float64 == typeof(LinearAscending(; start=0.1, stop=0.5).stop)

    # Check the associated supertype
    @test isa(LinearAscending(), IonScanOrder)

    # Check if error is thrown if start and stop do not satisfy the assumptions
    # 0 ≤ start < stop ≤ 1
    @test_throws ArgumentError LinearAscending(; start=-0.1)
    @test_throws ArgumentError LinearAscending(; start=1.1)
    @test_throws ArgumentError LinearAscending(; stop=-0.1)
    @test_throws ArgumentError LinearAscending(; stop=1.1)
    @test_throws ArgumentError LinearAscending(; start=-0.1, stop=1.1)
    @test_throws ArgumentError LinearAscending(; start=0.5, stop=0.5)
end


############################################################################################
# LinearDescending(; start::Real=0, stop::Real=1) <: IonScanOrder
############################################################################################
@testset "LinearDescending" begin
    # Verify object construction and field types depending on constructor arguments
    @test 0 == LinearDescending().start
    @test Int == typeof(LinearDescending().start)
    @test 1 == LinearDescending().stop
    @test Int == typeof(LinearDescending().stop)
    @test 0.5 == LinearDescending(; start=0.5).start
    @test Float64 == typeof(LinearDescending(; start=0.5).start)
    @test 1 == LinearDescending(; start=0.5).stop
    @test 0 == LinearDescending(; stop=0.5).start
    @test Float64 == typeof(LinearDescending(; stop=0.5).stop)
    @test 0.5 == LinearDescending(; stop=0.5).stop
    @test 0.1 == LinearDescending(; start=0.1, stop=0.5).start
    @test 0.5 == LinearDescending(; start=0.1, stop=0.5).stop
    @test Float64 == typeof(LinearDescending(; start=0.1, stop=0.5).start)
    @test Float64 == typeof(LinearDescending(; start=0.1, stop=0.5).stop)

    # Check the associated supertype
    @test isa(LinearDescending(), IonScanOrder)

    # Check if error is thrown if start and stop do not satisfy the assumptions
    # 0 ≤ start < stop ≤ 1
    @test_throws ArgumentError LinearDescending(; start=-0.1)
    @test_throws ArgumentError LinearDescending(; start=1.1)
    @test_throws ArgumentError LinearDescending(; stop=-0.1)
    @test_throws ArgumentError LinearDescending(; stop=1.1)
    @test_throws ArgumentError LinearDescending(; start=-0.1, stop=1.1)
    @test_throws ArgumentError LinearDescending(; start=0.5, stop=0.5)
end



############################################################################################
# ionscantimeshift(chrom::AbstractChromMS, ionscanorder::IonScanOrder)
############################################################################################
@testset "ionscantimeshift LinearAscending" begin
    @test -0.5u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearAscending())(1)
    @test 0u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearAscending())(2)
    @test -0.25u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearAscending(start=0.5))(1)
    @test 0u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearAscending(start=0.5))(2)
    @test -0.75u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearAscending(stop=0.5))(1)
    @test -0.5u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearAscending(stop=0.5))(2)
    @test -0.5u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearAscending(start=0.25, stop=0.75))(1)
    @test -0.25u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearAscending(start=0.25, stop=0.75))(2)
end

############################################################################################
# ionscantimeshift(chrom::AbstractChromMS, ionscanorder::IonScanOrder)
############################################################################################
@testset "ionscantimeshift LinearDescending" begin
    @test -0u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearDescending())(1)
    @test -0.5u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearDescending())(2)
    @test 0u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearDescending(start=0.5))(1)
    @test -0.25u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearDescending(start=0.5))(2)
    @test -0.5u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearDescending(stop=0.5))(1)
    @test -0.75u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearDescending(stop=0.5))(2)
    @test -0.25u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearDescending(start=0.25, stop=0.75))(1)
    @test -0.5u"s" == ionscantimeshift(ChromMS([1.0, 2.0, 3.0]u"s", [85, 100], 
        [0 12; 34 956; 23 1]), LinearDescending(start=0.25, stop=0.75))(2)
end