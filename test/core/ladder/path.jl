using Test
using JuChrom

path_contrast(scan, score; z=score, centerz=2.0, isolationz=2.0) = (
    ladderstep=0,
    leftindex=scan - 1,
    apexindex=scan,
    rightindex=scan + 1,
    scanindices=[scan - 1, scan, scan + 1],
    peakmodel=[0.5, 1.0, 0.5],
    molecularionscore=score,
    z=z,
    centerz=centerz,
    isolationz=isolationz,
    centervslowerz=isolationz,
    centervsupperz=isolationz,
)

@testset "alkaneladderpath selects an increasing scored path" begin
    contrasts = Dict(
        8 => [path_contrast(10, 3.0)],
        9 => [path_contrast(20, 1.0)],
        10 => [
            path_contrast(15, 100.0),
            path_contrast(30, 5.0),
        ],
        11 => [path_contrast(40, 0.0)],
    )

    pathinfo = alkaneladderpath(
        contrasts;
        minsteps=3,
        stepreward=0.05,
        maxmissingsteps=1,
        missingsteppenalty=2.0,
    )

    @test pathinfo.status == :success
    @test pathinfo.laddersteps == [8, 9, 10]
    @test pathinfo.scanindices == [10, 20, 30]
    @test pathinfo.score ≈ 3.0 + 1.0 + 5.0 + 3 * 0.05
    @test length(pathinfo.windows) == 3
end

@testset "alkaneladderpath can bridge one missing step with a penalty" begin
    contrasts = Dict(
        8 => [path_contrast(10, 3.0)],
        10 => [path_contrast(30, 5.0)],
    )

    pathinfo = alkaneladderpath(
        contrasts;
        minsteps=2,
        stepreward=0.05,
        maxmissingsteps=1,
        missingsteppenalty=2.0,
    )

    @test pathinfo.status == :success
    @test pathinfo.laddersteps == [8, 10]
    @test pathinfo.score ≈ 3.0 + 5.0 + 2 * 0.05 - 2.0

    failed = alkaneladderpath(
        contrasts;
        minsteps=2,
        maxmissingsteps=0,
    )
    @test failed.status == :failed
    @test failed.reason == :no_valid_path
end

@testset "alkaneladderpath ignores nonpositive molecular-ion scores" begin
    contrasts = Dict(
        8 => [
            path_contrast(10, 0.0),
            (apexindex=20, z=100.0, centerz=2.0, isolationz=2.0),
        ],
    )

    pathinfo = alkaneladderpath(contrasts; minsteps=1)

    @test pathinfo.status == :failed
    @test pathinfo.reason == :no_candidates
    @test isempty(pathinfo.path)
end

@testset "alkane ladder path validation" begin
    @test JuChrom.validate_alkane_ladder_path_settings(
        1.645,
        1.645,
        1,
        0.05,
        100,
        25.0,
        5.0,
        2.5,
        1,
        2.0,
    ) === nothing
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_path_settings(
        1.645,
        1.645,
        0,
        0.05,
        100,
        25.0,
        5.0,
        2.5,
        1,
        2.0,
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_path_settings(
        1.645,
        1.645,
        1,
        Inf,
        100,
        25.0,
        5.0,
        2.5,
        1,
        2.0,
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_path_settings(
        1.645,
        1.645,
        1,
        0.05,
        0,
        25.0,
        5.0,
        2.5,
        1,
        2.0,
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_path_settings(
        1.645,
        1.645,
        1,
        0.05,
        100,
        25.0,
        5.0,
        2.5,
        -1,
        2.0,
    )
    @test_throws ArgumentError JuChrom.validate_alkane_ladder_path_settings(
        1.645,
        1.645,
        1,
        0.05,
        100,
        25.0,
        5.0,
        2.5,
        1,
        -1.0,
    )
end
