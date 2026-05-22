using Test
using JuChrom

function series_path_test_trace(values; retentions=collect(1.0:length(values)))
    chroms = [
        ChromScan(rt, nothing, y, nothing)
        for (rt, y) in zip(retentions, values)
    ]
    ChromScanSeries(chroms)
end

@testset "alkane series path DP" begin
    @testset "selects a smooth increasing path" begin
        traces = Dict(
            8 => series_path_test_trace([0.0, 10.0, 0.0, 7.0, 0.0, 0.0, 0.0]),
            9 => series_path_test_trace([0.0, 0.0, 0.0, 9.0, 0.0, 2.0, 0.0]),
            10 => series_path_test_trace([0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0]),
        )

        result = JuChrom.alkaneseriespath(
            traces;
            minsteps=3,
            stepreward=0.0,
            minrelativepeakheight=0.0,
            minpeakheight=0.1,
            maxcandidatespertrace=nothing,
            spacingweight=1.0,
            gapincreaseweight=0.25,
        )

        @test result.success
        @test result.failurereason === nothing
        @test result.carbonnumbers == [8, 9, 10]
        @test result.scanindices == [2, 4, 6]
        @test result.retentions == [2.0, 4.0, 6.0]
        @test result.scores == [10.0, 9.0, 8.0]
        @test result.gaps == [2.0, 2.0]
        @test result.objective == 27.0
        @test length(result.candidatesbycarbon[8]) == 2
        @test result.candidateruns == [[8, 9, 10]]
        @test only(result.runresults).success
        @test result.settings.minsteps == 3
    end

    @testset "limits candidates and keeps debugging data on failure" begin
        traces = Dict(
            8 => series_path_test_trace([0.0, 1.0, 0.0, 5.0, 0.0]),
            9 => series_path_test_trace([0.0, 4.0, 0.0, 0.0, 0.0]),
        )

        result = JuChrom.alkaneseriespath(
            traces;
            minsteps=2,
            stepreward=0.0,
            minrelativepeakheight=0.0,
            minpeakheight=0.1,
            maxcandidatespertrace=1,
        )

        @test !result.success
        @test occursin("no valid path", result.failurereason)
        @test result.path == NamedTuple[]
        @test result.scanindices == Int[]
        @test result.candidatesbycarbon[8] == [
            (carbon=8, scanindex=4, retention=4.0, score=5.0),
        ]
        @test result.candidatesbycarbon[9] == [
            (carbon=9, scanindex=2, retention=2.0, score=4.0),
        ]
        @test result.candidateruns == [[8, 9]]
        @test only(result.runresults).success == false
    end

    @testset "validates settings and traces" begin
        traces = Dict(8 => series_path_test_trace([0.0, 1.0, 0.0]))

        @test_throws ArgumentError JuChrom.alkaneseriespath(traces; minsteps=0)
        @test_throws ArgumentError JuChrom.alkaneseriespath(traces; stepreward=Inf)
        @test_throws ArgumentError JuChrom.alkaneseriespath(
            traces;
            minrelativepeakheight=-0.1,
        )
        @test_throws ArgumentError JuChrom.alkaneseriespath(traces; minpeakheight=-0.1)
        @test_throws ArgumentError JuChrom.alkaneseriespath(
            traces;
            maxcandidatespertrace=0,
        )
        @test_throws ArgumentError JuChrom.alkaneseriespath(traces; spacingweight=-1.0)
        @test_throws ArgumentError JuChrom.alkaneseriespath(
            traces;
            gapincreaseweight=-1.0,
        )
    end
end
