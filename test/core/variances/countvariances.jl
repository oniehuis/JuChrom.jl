using Test
using Unitful
using JuChrom

function _countvariance_test_matrix()
    validwindow(base, step) = base .+ [0.0, step, -step, step, -step]
    vcat(
        [0.0 10.0 10.0;
         10.0 10.0 10.0;
         10.0 10.0 10.0;
         10.0 10.0 10.0;
         10.0 10.0 10.0],
        hcat(validwindow(30.0, 3.0), validwindow(50.0, 5.0), validwindow(70.0, 7.0)),
        hcat(validwindow(35.0, 3.5), validwindow(55.0, 5.5), validwindow(75.0, 7.5)),
        hcat(validwindow(40.0, 4.0), validwindow(60.0, 6.0), validwindow(80.0, 8.0)),
    )
end

@testset "countvariances" begin
    rawcounts = _countvariance_test_matrix()

    @testset "flat result and variance formula" begin
        result = countvariances(
            rawcounts;
            windowsize=5,
            mintransitioncount=1,
            zerothresholdquantile=1.0,
            positivecountquantile=0.0,
        )

        @test result isa CountVarianceEstimate
        @test varianceunit(result) ≡ nothing
        @test result.variances isa Matrix{Float64}
        @test size(result.variances) == size(rawcounts)
        @test result.noisesigma > 0
        @test result.variancefactor == result.noisesigma^2
        @test result.zerowindowthreshold == 10.0
        @test result.positivecountthreshold == 10.0
        @test result.zerowindowcount == 1
        @test result.normalizeddeviationcount == 45
        @test result.windowsize == 5
        @test result.mintransitioncount == 1
        @test result.zerothresholdquantile == 1.0
        @test result.positivecountquantile == 0.0
        @test result.intensityfloor == 10.0
        @test result.variances ≈ result.variancefactor .* max.(rawcounts, result.intensityfloor)
    end

    @testset "explicit intensity floor" begin
        result = countvariances(
            rawcounts;
            windowsize=5,
            mintransitioncount=1,
            intensityfloor=25.0,
        )

        @test result.intensityfloor == 25.0
        @test result.variances ≈ result.variancefactor .* max.(rawcounts, 25.0)
    end

    @testset "MassScanMatrix wrapper" begin
        msm = MassScanMatrix(collect(1.0:size(rawcounts, 1)), [50.0, 60.0, 70.0], rawcounts)
        frommatrix = countvariances(rawcounts; windowsize=5, mintransitioncount=1)
        frommsm = countvariances(msm; windowsize=5, mintransitioncount=1)

        @test frommsm isa CountVarianceEstimate
        @test frommsm.variances isa Matrix{Float64}
        @test varianceunit(frommsm) ≡ nothing
        @test rawvariances(frommsm) ≈ frommatrix.variances
        @test frommsm.noisesigma == frommatrix.noisesigma
        @test frommsm.intensityfloor == frommatrix.intensityfloor
    end

    @testset "MassScanMatrix wrapper preserves variance units" begin
        msm = MassScanMatrix(
            collect(1.0:size(rawcounts, 1)),
            u"s",
            [50.0, 60.0, 70.0],
            nothing,
            rawcounts,
            u"pA"
        )
        result = countvariances(msm; windowsize=5, mintransitioncount=1)

        @test result isa CountVarianceEstimate
        @test result.variances isa Matrix{Float64}
        @test varianceunit(result) == u"pA"^2
        @test rawvariances(result) ≈
            countvariances(rawcounts; windowsize=5, mintransitioncount=1).variances
        @test variances(result) == result.variances .* u"pA"^2
    end

    @testset "unitful matrix wrapper preserves variance units" begin
        result = countvariances(
            rawcounts .* u"pA";
            windowsize=5,
            mintransitioncount=1
        )

        @test result isa CountVarianceEstimate
        @test result.variances isa Matrix{Float64}
        @test varianceunit(result) == u"pA"^2
        @test rawvariances(result; unit=u"nA"^2) ≈
            countvariances(rawcounts; windowsize=5, mintransitioncount=1).variances * 1e-6
        floored = countvariances(
            rawcounts .* u"pA";
            windowsize=5,
            mintransitioncount=1,
            intensityfloor=25.0u"pA"
        )
        @test floored.intensityfloor == 25.0
        @test_throws ArgumentError countvariances(
            rawcounts;
            windowsize=5,
            mintransitioncount=1,
            intensityfloor=25.0u"pA"
        )
        @test_throws ArgumentError countvariances(
            rawcounts .* u"pA";
            windowsize=5,
            mintransitioncount=1,
            intensityfloor=25.0u"s"
        )
    end

    @testset "count noise helpers" begin
        msm = MassScanMatrix(collect(1.0:size(rawcounts, 1)), [50.0, 60.0, 70.0], rawcounts)
        stats = JuChrom.count_noise_stats(
            rawcounts;
            windowsize=5,
            mintransitioncount=1,
            zerothresholdquantile=1.0,
        )
        msmstats = JuChrom.count_noise_stats(
            msm;
            windowsize=5,
            mintransitioncount=1,
            zerothresholdquantile=1.0,
        )

        @test JuChrom.count_noise_sigma(
            rawcounts;
            windowsize=5,
            mintransitioncount=1,
            zerothresholdquantile=1.0,
        ) == stats.noisesigma
        @test JuChrom.count_noise_sigma(
            msm;
            windowsize=5,
            mintransitioncount=1,
            zerothresholdquantile=1.0,
        ) == stats.noisesigma
        @test msmstats == stats
    end

    @testset "input validation" begin
        negativecounts = copy(rawcounts)
        negativecounts[1, 1] = -1.0
        @test_throws ArgumentError countvariances(
            negativecounts; windowsize=5, mintransitioncount=1)

        nonfinitecounts = copy(rawcounts)
        nonfinitecounts[1, 1] = Inf
        @test_throws ArgumentError countvariances(
            nonfinitecounts; windowsize=5, mintransitioncount=1)

        @test_throws ArgumentError countvariances(
            rawcounts; windowsize=1, mintransitioncount=0)
        @test_throws ArgumentError countvariances(
            rawcounts; windowsize=5, mintransitioncount=5)
        @test_throws ArgumentError countvariances(
            rawcounts; windowsize=5, mintransitioncount=1, zerothresholdquantile=-0.1)
        @test_throws ArgumentError countvariances(
            rawcounts; windowsize=5, mintransitioncount=1, positivecountquantile=1.1)
        @test_throws ArgumentError countvariances(
            rawcounts; windowsize=size(rawcounts, 1) + 1, mintransitioncount=1)
        @test_throws ArgumentError countvariances(
            rawcounts; windowsize=5, mintransitioncount=1, intensityfloor=0.0)
        @test_throws ArgumentError countvariances(
            rawcounts; windowsize=5, mintransitioncount=1, intensityfloor=Inf)
    end
end
