module TestIntensityUnits

using Test
using Unitful
using JuChrom
using SparseArrays

function dwell_fixture(storage, withvariance; mz=[57, 71, 85], ret=[1.0, 2.0]u"s")
    matrix = MassScanMatrix(ret, mz, storage([10.0 20.0 0.0; 40.0 50.0 0.0]);
        level=2, instrument=(model=["MS"],), acquisition=(mode=["SIM"],),
        user=(name=["tester"],), sample=(id=["S1"],), extras=Dict("note" => ["counts"]))
    withvariance ? VarianceMassScanMatrix(matrix, storage([1.0 4.0 0.0; 16.0 25.0 0.0])) : matrix
end

@testset "explicit dwell ($storage, variance=$withvariance)" for storage in (identity, sparse),
        withvariance in (false, true)
    matrix = dwell_fixture(storage, withvariance)
    source_intensities = copy(rawintensities(matrix))
    expected = [1000 1000 0; 4000 2500 0]
    spec = (85 => 40u"ms", 57 => 0.01u"s", 71 => 20000u"μs", 99 => 1u"s")
    for dwell in (spec, reverse(spec), collect(spec), Dict(spec), Any[spec...],
                  view(collect(spec), :))
        normalized = dwellnormalize(matrix, dwell)
        @test rawintensities(normalized) ≈ expected
        @test intensityunit(normalized) == u"s^-1"
        @test issparse(rawintensities(normalized)) == (storage === sparse)
        @test rawintensities(matrix) == source_intensities
        @test rawretentions(normalized) == rawretentions(matrix)
        @test rawretentions(normalized) !== rawretentions(matrix)
        @test rawmzvalues(normalized) == rawmzvalues(matrix)
        @test rawmzvalues(normalized) !== rawmzvalues(matrix)
        @test retentionunit(normalized) == retentionunit(matrix)
        @test mzunit(normalized) == mzunit(matrix)
        @test level(normalized) == level(matrix)
        for getter in (instrument, acquisition, user, sample, extras)
            @test getter(normalized) == getter(matrix)
            @test getter(normalized) !== getter(matrix)
        end
        if withvariance
            @test normalized isa VarianceMassScanMatrix
            @test rawvariances(normalized) ≈ [10000 10000 0; 160000 62500 0]
            @test varianceunit(normalized) == u"s^-2"
            @test issparse(rawvariances(normalized)) == (storage === sparse)
            @test rawvariances(normalized) !== rawvariances(matrix)
            @test rawvariances(matrix) == [1 4 0; 16 25 0]
        end
    end

    # Equivalent duration units yield the same numeric rates and output unit.
    for dwell in (10u"ms", 0.01u"s", 10000u"μs")
        normalized = dwellnormalize(matrix, dwell)
        @test rawintensities(normalized) ≈ source_intensities .* 100
        @test intensityunit(normalized) == u"s^-1"
        @test issparse(rawintensities(normalized)) == (storage === sparse)
        if withvariance
            @test rawvariances(normalized) ≈ rawvariances(matrix) .* 10000
            @test varianceunit(normalized) == u"s^-2"
            @test issparse(rawvariances(normalized)) == (storage === sparse)
        end
    end

    mixedprecision = dwellnormalize(matrix, (57 => 0.25f0*u"s", 71 => 500.0u"ms", 85 => 1u"s"))
    @test rawintensities(mixedprecision) == [40 40 0; 160 100 0]
    if withvariance
        @test rawvariances(mixedprecision) == [16 16 0; 256 100 0]
    end

    normalized = dwellnormalize(matrix, spec)
    rawintensities(normalized)[1, 1] = -1
    rawretentions(normalized)[1] = 0
    rawmzvalues(normalized)[1] = 1
    instrument(normalized).model[1] = "changed"
    extras(normalized)["note"][1] = "changed"
    @test rawintensities(matrix)[1, 1] == 10
    @test rawretentions(matrix)[1] == 1
    @test rawmzvalues(matrix)[1] == 57
    @test instrument(matrix).model == ["MS"]
    @test extras(matrix)["note"] == ["counts"]
    if withvariance
        rawvariances(normalized)[1, 1] = 0
        @test rawvariances(matrix)[1, 1] == 1
    end

    # Every column needs a duration even if all its intensities are zero.
    @test_throws ArgumentError dwellnormalize(matrix, (57 => 1u"s", 71 => 1u"s"))
    for bad in (0u"s", -1u"ms", NaN*u"s", Inf*u"s", -Inf*u"s",
                1u"kg", complex(1, 1)*u"s", floatmax(Float64)*u"d")
        @test_throws ArgumentError dwellnormalize(matrix, bad)
    end
    for bad in (0, 0u"s", -1u"s", NaN*u"s", Inf*u"s", 1u"kg", "1 s", nothing,
                missing, complex(1, 1)*u"s")
        # Invalid extra-ion durations are also rejected.
        @test_throws ArgumentError dwellnormalize(matrix, (spec..., 100 => bad))
    end
    for bad in ((), [], Dict(), [1u"s", 2u"s", 3u"s"], (1u"s", 2u"s"),
                (spec..., 5), (spec..., 57 => 1u"s"), (spec..., 99.0 => 2u"s"))
        @test_throws ArgumentError dwellnormalize(matrix, bad)
    end
    for badkey in (0, -1, NaN, Inf, -Inf, "100", missing, complex(100, 1), 100u"Th")
        @test_throws ArgumentError dwellnormalize(matrix, (spec..., badkey => 1u"s"))
    end
    @test_throws ArgumentError dwellnormalize(matrix, (57u"Th" => 1u"s",))
    annotated = withintensityunit(matrix, u"pA")
    @test_throws ArgumentError dwellnormalize(annotated, 1u"s")
    @test_throws ArgumentError dwellnormalize(annotated, spec)
    @test_throws ArgumentError dwellnormalize(dwellnormalize(matrix, spec), spec)
end

@testset "dwell m/z units and exact matching" for withvariance in (false, true)
    matrix = dwell_fixture(identity, withvariance; mz=[57, 71, 85]u"Th")
    expected = [1000 1000 0; 4000 2500 0]
    for spec in ((57 => 10u"ms", 71 => 20u"ms", 85 => 40u"ms"),
                 (85u"Th" => 40u"ms", 57000u"mTh" => 10u"ms", 71u"Th" => 20u"ms"),
                 Dict(57000u"mTh" => 10u"ms", 71000u"mTh" => 20u"ms", 85000u"mTh" => 40u"ms"))
        @test rawintensities(dwellnormalize(matrix, spec)) ≈ expected
        @test mzunit(dwellnormalize(matrix, spec)) == u"Th"
    end
    @test_throws ArgumentError dwellnormalize(matrix, (57u"Th" => 1u"s", 57000u"mTh" => 1u"s"))
    @test_throws ArgumentError dwellnormalize(matrix, (57u"Th" => 1u"s", 71 => 1u"s"))
    @test_throws ArgumentError dwellnormalize(matrix, (57u"kg" => 1u"s",))
    @test_throws ArgumentError dwellnormalize(matrix, (0u"Th" => 1u"s",))
    @test_throws ArgumentError dwellnormalize(matrix, (NaN*u"Th" => 1u"s",))

    f32 = MassScanMatrix([1.0f0], Float32[57.1], fill(10.0f0, 1, 1))
    narrow = withvariance ? VarianceMassScanMatrix(f32, fill(1.0f0, 1, 1)) : f32
    @test_throws ArgumentError dwellnormalize(narrow, (57.1 => 1u"s",))
    @test rawintensities(dwellnormalize(narrow, (57.1f0 => 1u"s",))) == fill(10, 1, 1)
    @test rawintensities(dwellnormalize(narrow, (Float64(57.1f0) => 1u"s",))) == fill(10, 1, 1)

    precise = MassScanMatrix([1], [nextfloat(57.0)], fill(10, 1, 1))
    @test_throws ArgumentError dwellnormalize(precise, (57.0f0 => 1u"s",))
    large = 2^53 + 1
    large_matrix = MassScanMatrix([1], [large], fill(10, 1, 1))
    @test rawintensities(dwellnormalize(large_matrix, (1.0 => 1u"s", large => 2u"s"))) == fill(5, 1, 1)
    @test_throws ArgumentError dwellnormalize(large_matrix, (Float64(large) => 2u"s",))
end

@testset "dwell normalization is independent of selection and retention coordinates" for withvariance in (false, true)
    # Single scans, non-time retention units, repeated/descending coordinates, and
    # durations longer than adjacent retention spacing do not constrain normalization.
    for ret in ([1.0], [1000.0], [1.0]u"s", [1000.0]u"Th", [0.0, 0.001]u"s",
                [1.0, 1.0], [2.0, 1.0], [1.0, 1000.0]u"s")
        parentmatrix = MassScanMatrix(ret, [57, 71], repeat([10 20], length(ret), 1))
        matrix = withvariance ? VarianceMassScanMatrix(parentmatrix, repeat([1 4], length(ret), 1)) : parentmatrix
        for dwell in (2u"s", (57 => 2u"s", 71 => 2u"s", 99 => 3u"s"))
            normalized = dwellnormalize(matrix, dwell)
            @test rawintensities(normalized) == repeat([5 10], length(ret), 1)
            @test retentions(normalized) == ret
            if withvariance
                @test rawvariances(normalized) == repeat([0.25 1.0], length(ret), 1)
            end
        end
    end
    spec = (57 => 10u"ms", 71 => 20u"ms", 85 => 40u"ms")
    full = dwell_fixture(identity, withvariance)
    reducedparent = MassScanMatrix([2.0], [71], fill(50, 1, 1))
    reduced = withvariance ? VarianceMassScanMatrix(reducedparent, fill(25, 1, 1)) : reducedparent
    @test rawintensities(dwellnormalize(reduced, spec))[1] == rawintensities(dwellnormalize(full, spec))[2, 2]
    if withvariance
        @test rawvariances(dwellnormalize(reduced, spec))[1] == rawvariances(dwellnormalize(full, spec))[2, 2]
    end
end

@testset "short dwell durations do not overflow integer variance scaling" for storage in (identity, sparse)
    matrix = MassScanMatrix([1], [57], storage(fill(2.0, 1, 1)))
    vmsm = VarianceMassScanMatrix(matrix, storage(fill(4.0, 1, 1)))
    for dwell in (1u"ps", (57 => 1u"ps",))
        normalized = dwellnormalize(vmsm, dwell)
        @test rawintensities(normalized) ≈ fill(2e12, 1, 1)
        @test rawvariances(normalized) ≈ fill(4e24, 1, 1)
        @test issparse(rawvariances(normalized)) == (storage === sparse)
    end
end

@testset "dwellnormalize removed APIs" for withvariance in (false, true)
    matrix = dwell_fixture(identity, withvariance)
    @test_throws MethodError dwellnormalize(matrix)
    @test_throws MethodError dwellnormalize(matrix; acquisition=:sequential)
    @test_throws MethodError dwellnormalize(matrix, 1u"s"; acquisition=:simultaneous)
    @test_throws MethodError dwellnormalize(matrix, 1)
    @test_throws MethodError dwellnormalize(matrix, 1, u"s")
    @test_throws MethodError dwellnormalize(matrix, [1, 2, 3], u"s")
    @test_throws ArgumentError dwellnormalize(matrix, [1, 2, 3]u"s")
end

@testset "withintensityunit MassScanMatrix" begin
    ret = [1.0, 2.0]u"s"
    mz = [100.0, 200.0, 300.0]
    ints = [10.0 20.0 30.0; 40.0 50.0 60.0]
    msm = MassScanMatrix(
        ret,
        mz,
        ints;
        level=2,
        instrument=(model="MS",),
        acquisition=(mode="scan",),
        user=(name="tester",),
        sample=(id="S1",),
        extras=Dict("note" => "raw counts"),
    )

    annotated = withintensityunit(msm, u"pA")

    @test rawintensities(annotated) == rawintensities(msm)
    @test intensityunit(annotated) == u"pA"
    @test rawretentions(annotated) == rawretentions(msm)
    @test retentionunit(annotated) == retentionunit(msm)
    @test rawmzvalues(annotated) == rawmzvalues(msm)
    @test mzunit(annotated) == mzunit(msm)
    @test level(annotated) == level(msm)
    @test instrument(annotated) == instrument(msm)
    @test acquisition(annotated) == acquisition(msm)
    @test user(annotated) == user(msm)
    @test sample(annotated) == sample(msm)
    @test extras(annotated) == extras(msm)

    extras(msm)["note"] = "changed"
    @test extras(annotated)["note"] == "raw counts"

    unitful_msm = MassScanMatrix(ret, mz, ints .* u"pA")
    @test_throws ArgumentError withintensityunit(unitful_msm, u"nA")
end

@testset "withintensityunit VarianceMassScanMatrix" begin
    ret = [1.0, 2.0]u"s"
    mz = [100.0, 200.0, 300.0]
    ints = [10.0 20.0 30.0; 40.0 50.0 60.0]
    vars = [1.0 4.0 9.0; 16.0 25.0 36.0]
    msm = MassScanMatrix(
        ret,
        mz,
        ints;
        level=2,
        instrument=(model="MS",),
        acquisition=(mode="scan",),
        user=(name="tester",),
        sample=(id="S1",),
        extras=Dict("note" => "raw counts")
    )
    vmsm = VarianceMassScanMatrix(msm, vars)

    annotated = withintensityunit(vmsm, u"pA")

    @test annotated isa VarianceMassScanMatrix
    @test rawintensities(annotated) == rawintensities(vmsm)
    @test rawvariances(annotated) == rawvariances(vmsm)
    @test intensityunit(annotated) == u"pA"
    @test varianceunit(annotated) == u"pA"^2
    @test rawretentions(annotated) == rawretentions(vmsm)
    @test retentionunit(annotated) == retentionunit(vmsm)
    @test rawmzvalues(annotated) == rawmzvalues(vmsm)
    @test mzunit(annotated) == mzunit(vmsm)
    @test level(annotated) == level(vmsm)
    @test instrument(annotated) == instrument(vmsm)
    @test acquisition(annotated) == acquisition(vmsm)
    @test user(annotated) == user(vmsm)
    @test sample(annotated) == sample(vmsm)
    @test extras(annotated) == extras(vmsm)

    extras(msm)["note"] = "changed"
    @test extras(annotated)["note"] == "raw counts"

    unitful_vmsm = VarianceMassScanMatrix(withintensityunit(msm, u"pA"), vars)
    @test_throws ArgumentError withintensityunit(unitful_vmsm, u"nA")
end

end
