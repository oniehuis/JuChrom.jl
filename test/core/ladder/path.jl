using Test
using JuChrom

function path_contrast(scan, score; z=score, centerz=2.0, isolationz=2.0)
    AlkaneMolecularIonContrast(
        0,
        scan - 1,
        scan,
        scan + 1,
        [scan - 1, scan, scan + 1],
        [0.5, 1.0, 0.5],
        1.0,
        114,
        Int[],
        Int[],
        Int[],
        Float64[],
        Float64[],
        Float64[],
        false,
        false,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        centerz,
        0.0,
        0.0,
        isolationz,
        isolationz,
        isolationz,
        score,
        z
    )
end

path_test_settings() = (
    centerzmin=1.645,
    isolationzmin=1.645,
    minsteps=5,
    stepreward=0.05,
    maxcandidatesperstep=100,
    spacingweight=25.0,
    gapincreaseweight=5.0,
    maxgapratio=2.5,
    maxmissingsteps=1,
    missingsteppenalty=2.0,
    msm=nothing,
    variances=nothing,
    standard=nothing,
    massspectrummatch=false,
    massspectrummatchdistanceweight=5.0,
    massspectrumvariancefloor=1.0
)

@testset "alkaneladderpath selects an increasing scored path" begin
    contrasts = Dict(
        8 => [path_contrast(10, 3.0)],
        9 => [path_contrast(20, 1.0)],
        10 => [
            path_contrast(15, 100.0),
            path_contrast(30, 5.0)
        ],
        11 => [path_contrast(40, 0.0)]
    )
    settings = merge(
        path_test_settings(),
        (
            minsteps=3,
            stepreward=0.05,
            maxmissingsteps=1,
            missingsteppenalty=2.0
        )
    )

    pathinfo = alkaneladderpath(
        contrasts;
        settings...
    )

    @test pathinfo.status ≡ :success
    @test pathinfo.laddersteps == [8, 9, 10]
    @test pathinfo.scanindices == [10, 20, 30]
    @test pathinfo.score ≈ 3.0 + 1.0 + 5.0 + 3 * 0.05
    @test length(pathinfo.windows) == 3
end

@testset "alkaneladderpath settings preserve supplied values" begin
    contrasts = Dict(8 => [path_contrast(10, 3.0)])
    settings = merge(
        path_test_settings(),
        (
            centerzmin=Float32(1.645),
            isolationzmin=Float32(1.645),
            minsteps=Int32(1),
            stepreward=Float32(0.05),
            maxcandidatesperstep=Int32(100),
            spacingweight=Float32(25.0),
            gapincreaseweight=Float32(5.0),
            maxgapratio=Float32(2.5),
            maxmissingsteps=Int32(1),
            missingsteppenalty=Float32(2.0),
            massspectrummatchdistanceweight=Float32(5.0),
            massspectrumvariancefloor=Float32(1.0)
        )
    )

    pathinfo = alkaneladderpath(contrasts; settings...)

    @test pathinfo.settings.centerzmin ≡ settings.centerzmin
    @test pathinfo.settings.isolationzmin ≡ settings.isolationzmin
    @test pathinfo.settings.minsteps ≡ settings.minsteps
    @test pathinfo.settings.stepreward ≡ settings.stepreward
    @test pathinfo.settings.maxcandidatesperstep ≡ settings.maxcandidatesperstep
    @test pathinfo.settings.spacingweight ≡ settings.spacingweight
    @test pathinfo.settings.gapincreaseweight ≡ settings.gapincreaseweight
    @test pathinfo.settings.maxgapratio ≡ settings.maxgapratio
    @test pathinfo.settings.maxmissingsteps ≡ settings.maxmissingsteps
    @test pathinfo.settings.missingsteppenalty ≡ settings.missingsteppenalty
    @test pathinfo.settings.massspectrummatchdistanceweight ≡
        settings.massspectrummatchdistanceweight
    @test pathinfo.settings.massspectrumvariancefloor ≡
        settings.massspectrumvariancefloor
end

@testset "alkaneladderpath can bridge one missing step with a penalty" begin
    contrasts = Dict(
        8 => [path_contrast(10, 3.0)],
        10 => [path_contrast(30, 5.0)]
    )
    settings = merge(
        path_test_settings(),
        (
            minsteps=2,
            stepreward=0.05,
            maxmissingsteps=1,
            missingsteppenalty=2.0
        )
    )

    pathinfo = alkaneladderpath(
        contrasts;
        settings...
    )

    @test pathinfo.status ≡ :success
    @test pathinfo.laddersteps == [8, 10]
    @test pathinfo.score ≈ 3.0 + 5.0 + 2 * 0.05 - 2.0

    failedsettings = merge(
        path_test_settings(),
        (
            minsteps=2,
            maxmissingsteps=0
        )
    )
    failed = alkaneladderpath(contrasts; failedsettings...)
    @test failed.status ≡ :failed
    @test failed.reason ≡ :no_valid_path
end

@testset "alkaneladderpath ignores nonpositive molecular-ion scores" begin
    contrasts = Dict(
        8 => [
            path_contrast(10, 0.0)
        ]
    )
    settings = merge(path_test_settings(), (minsteps=1,))

    pathinfo = alkaneladderpath(contrasts; settings...)

    @test pathinfo.status ≡ :failed
    @test pathinfo.reason ≡ :no_candidates
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
        5.0,
        1.0
    ) ≡ nothing
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
        5.0,
        1.0
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
        5.0,
        1.0
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
        5.0,
        1.0
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
        5.0,
        1.0
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
        5.0,
        1.0
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
        2.0,
        -1.0,
        1.0
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
        2.0,
        5.0,
        0.0
    )
end
