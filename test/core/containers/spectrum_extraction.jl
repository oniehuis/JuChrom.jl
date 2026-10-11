module TestSpectrumExtraction
using Test
using JuChrom
using SparseArrays: sparse

@testset "Exact m/z matching and mean scan subtraction" begin
    series = MassScanSeries([
        MassScan(1.0, [40.0, 50.0, 60.0], [10.0, 2.0, 5.0]),
        MassScan(2.0, [40.0, 50.01, 70.0], [3.0, 4.0, 8.0]),
        MassScan(3.0, [40.0, 50.0, 80.0], [20.0, 1.0, 7.0]),
    ])
    before = [copy(rawintensities(s)) for s in scans(series)]
    result = extractmassspectrum(series, [2, 1, 1]; subtract=[3, 3])
    @test mzvalues(result) == [40, 50, 50.01, 60, 70, 80]
    @test intensities(result) == [0, 0, 2, 2.5, 4, 0]
    @test attrs(result).aggregation == :mean
    @test attrs(result).added_scans == [1, 2]
    @test attrs(result).subtracted_scans == [3]
    @test [rawintensities(s) for s in scans(series)] == before
    @test all(iszero, intensities(extractmassspectrum(series, Int[]; subtract=[3])))
    @test all(iszero, intensities(extractmassspectrum(series, [1]; subtract=[1])))
    @test_throws ArgumentError extractmassspectrum(series, Int[])
    @test_throws BoundsError extractmassspectrum(series, [0])
    @test_throws BoundsError extractmassspectrum(series, [1]; subtract=[4])

    # Floor only after combining every contribution, independent of source order.
    signed = MassScanSeries([MassScan(i, [50], [value]) for (i, value) in
        enumerate([-5.0, 8.0, 1.0])])
    @test intensities(extractmassspectrum(signed, [1, 2]; subtract=[3])) == [0.5]
    near = MassScanSeries([MassScan(1, Float32[57.1, nextfloat(57.1f0)], [1, 2])])
    exact = extractmassspectrum(near, [1])
    @test mzvalues(exact) == rawmzvalues(scan(near, 1))
    @test eltype(mzvalues(exact)) == Float32
end

@testset "Group means use all scans, including missing channels" begin
    series = MassScanSeries([
        MassScan(1, [40, 50, 60], [12.0, 8.0, 6.0]),
        MassScan(2, [40, 50], [4.0, 0.0]),
        MassScan(3, [40, 50, 70], [3.0, 2.0, 10.0]),
        MassScan(4, [40, 50], [1.0, 4.0]),
        MassScan(5, [40], [2.0]),
    ])
    result = extractmassspectrum(series, [1, 2]; subtract=[3, 4, 5])
    @test mzvalues(result) == [40, 50, 60, 70]
    @test intensities(result) ≈ [6, 2, 3, 0]
    @test intensities(extractmassspectrum(series, [1, 2])) == [8, 4, 3]
    # Shared indices contribute with their own group's weight, not by cancellation
    # before computing means. Duplicate indices still count only once per group.
    @test intensities(extractmassspectrum(series, [1, 2, 2]; subtract=[2, 2])) == [4, 4, 3]
    repeated = MassScanSeries([
        MassScan(1, [40, 50], [12.0, 8.0]),
        MassScan(2, [40, 50], [12.0, 8.0]),
        MassScan(3, [40, 50], [3.0, 2.0]),
        MassScan(4, [40, 50], [3.0, 2.0]),
    ])
    for (positive, negative) in (([1], [3]), ([1, 2], [3]), ([1], [3, 4]))
        @test intensities(extractmassspectrum(repeated, positive; subtract=negative)) == [9, 6]
    end
end

@testset "Integer binning, units, and sufficiently wide accumulation" begin
    series = MassScanSeries([
        MassScan(1.0u"s", [50.1, 50.2]u"Th", [1000.0, 2000.0]u"pA"),
        MassScan(2.0u"s", [50.1, 60.1]u"Th", [500.0, 100.0]u"pA"),
    ])
    binned = binmzvalues(series)
    result = extractmassspectrum(binned, [1]; subtract=[2])
    @test mzvalues(result) == [50, 60]u"Th"
    @test intensities(result) == [2500, 0]u"pA"
    @test intensityunit(result) == u"pA"
    @test mzunit(result) == u"Th"
    average = extractmassspectrum(binned, [1, 2])
    @test intensities(average) == [1750, 50]u"pA"
    large = MassScanSeries([MassScan(i, [50], Int32[typemax(Int32)]) for i in 1:3])
    @test intensities(extractmassspectrum(large, [1, 2, 3])) == [Float64(typemax(Int32))]
    @test intensities(extractmassspectrum(large, [1, 2, 3]; subtract=[1])) == [0]
end

@testset "Direct matrix and variance-matrix extraction" begin
    mzs = Float32[40, 50, 60]
    values = [12.0 -2.0 6.0; 4.0 8.0 0.0; 3.0 2.0 10.0]
    dense = MassScanMatrix([1000.0, 1100.0, 1200.0], nothing,
        mzs, u"Th", values, u"pA")
    sparse_matrix = MassScanMatrix([1000.0, 1100.0, 1200.0], nothing,
        mzs, u"Th", sparse(values), u"pA")
    variance_matrix = VarianceMassScanMatrix(dense, ones(3, 3), u"pA^2")
    for matrix in (dense, sparse_matrix, variance_matrix)
        result = extractmassspectrum(matrix, [2, 1, 1]; subtract=[3, 3])
        @test mzvalues(result) == mzs * u"Th"
        @test eltype(result.mzvalues) == Float32
        @test intensities(result) == [5, 1, 0]u"pA"
        @test attrs(result).added_scans == [1, 2]
        @test attrs(result).subtracted_scans == [3]
        @test intensities(extractmassspectrum(matrix, [1, 2])) == [8, 3, 3]u"pA"
        @test all(iszero, intensities(extractmassspectrum(matrix, Int[]; subtract=[3])))
        @test intensities(extractmassspectrum(matrix, [1, 2]; subtract=[1])) == [0, 5, 0]u"pA"
        @test rawintensities(matrix) == values
        @test_throws ArgumentError extractmassspectrum(matrix, Int[])
        @test_throws BoundsError extractmassspectrum(matrix, [0])
        @test_throws BoundsError extractmassspectrum(matrix, [1]; subtract=[4])
    end
    @test rawvariances(variance_matrix) == ones(3, 3)
    @test parent(variance_matrix) === dense
end
end
