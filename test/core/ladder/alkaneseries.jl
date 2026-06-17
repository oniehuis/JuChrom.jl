using Test
using JuChrom

function test_ladder_step_result()
    apex(c, scan; source=:molecularion, cosine=NaN, required=NaN, good=true) = (
        success=true,
        ladderstep=c,
        apexscanindex=Float64(scan),
        apexretention=Float64(scan) * 10.0,
        source=source,
        mass_spectrum_cosine=Float64(cosine),
        required_cosine=Float64(required),
        good_for_calibration=good
    )

    AlkaneSeriesResult(
        nothing,
        ones(1, 1),
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        (path=NamedTuple[],),
        (
            apexes=[
                apex(10, 30; source=:molecularion, good=true),
                apex(8, 10; source=:molecularion, good=false),
                merge(apex(99, 99), (success=false,))
            ],
        ),
        (
            additions=[
                (
                    source=:gapfilled,
                    apex=apex(9, 20; source=:gapfilled, cosine=0.95, required=0.85)
                ),
                (
                    source=:rightextended,
                    apex=apex(11, 40; source=:rightextended, cosine=0.93, required=0.9)
                )
            ],
            gapfilled=[
                (
                    source=:gapfilled,
                    apex=apex(9, 20; source=:gapfilled, cosine=0.95, required=0.85)
                )
            ],
            leftextended=NamedTuple[],
            rightextended=[
                (
                    source=:rightextended,
                    apex=apex(11, 40; source=:rightextended, cosine=0.93, required=0.9)
                )
            ]
        )
    )
end

@testset "alkaneladdersteps returns merged refined ladder step view" begin
    result = test_ladder_step_result()

    steps = alkaneladdersteps(result)

    @test steps isa Vector{AlkaneLadderStep}
    @test [step.ladderstep for step in steps] == [8, 9, 10, 11]
    @test [step.source for step in steps] ==
        [:molecularion, :gapfilled, :molecularion, :rightextended]
    @test steps[2].massspectrumcosine == 0.95
    @test steps[2].requiredcosine == 0.85
    @test steps[1].goodforcalibration == false
    @test steps[2].apex.ladderstep == 9
end

@testset "alkaneladdersteps filters by evidence source" begin
    result = test_ladder_step_result()

    mi = alkaneladdersteps(result; gapfilled=false, edgeextended=false)
    gap = alkaneladdersteps(result; molecularion=false, edgeextended=false)
    edge = alkaneladdersteps(result; molecularion=false, gapfilled=false)

    @test [step.ladderstep for step in mi] == [8, 10]
    @test [step.ladderstep for step in gap] == [9]
    @test [step.ladderstep for step in edge] == [11]
end
