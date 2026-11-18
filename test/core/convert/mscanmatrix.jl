module TestMscanmatrix

using Test
using JuChrom
using Unitful
using SparseArrays

function toy_series()
    MassScanSeries([
        MassScan(1.0u"s", [57, 71], Int16[10, 20]),
        MassScan(2.0u"s", [57, 85], Int16[5, 30])];
        instrument=(model="test",), acquisition=(mode="SIM",),
        user=(operator="Alice",), sample=(id="sample",), extras=Dict("note" => "test"))
end

@testset "mscanmatrix explicit grid ($format)" for format in (DENSE, SPARSE)
    series = toy_series()
    expected = Int16[10 20 0 0; 5 0 30 0]
    for grid in (57:14:99, [57, 71, 85, 99], (57, 71, 85, 99),
                 view([57, 71, 85, 99], :), Any[57, 71, 85, 99])
        before = collect(grid)
        matrix = mscanmatrix(series, grid; format)
        @test matrix isa MassScanMatrix
        @test rawmzvalues(matrix) == [57, 71, 85, 99]
        @test rawmzvalues(matrix) isa Vector{Int}
        @test rawintensities(matrix) == expected
        @test eltype(rawintensities(matrix)) === Int16
        @test retentions(matrix) == [1.0, 2.0]u"s"
        @test retentionunit(matrix) === u"s"
        @test mzunit(matrix) === nothing
        @test intensityunit(matrix) === nothing
        @test level(matrix) == 1
        @test instrument(matrix) == instrument(series)
        @test acquisition(matrix) == acquisition(series)
        @test user(matrix) == user(series)
        @test sample(matrix) == sample(series)
        @test extras(matrix) == extras(series)
        @test collect(grid) == before
        @test rawintensities(first(scans(series))) == [10, 20]
        @test rawmzvalues(first(scans(series))) == [57, 71]
        if format === DENSE
            @test rawintensities(matrix) isa Matrix{Int16}
        else
            @test rawintensities(matrix) isa SparseMatrixCSC{Int16}
            @test nnz(rawintensities(matrix)) == 4
        end
    end

    # The requested axis and the source's arrays are not aliased by the output.
    grid = [57, 71, 85, 99]
    matrix = mscanmatrix(series, grid; format)
    grid[end] = 113
    @test rawmzvalues(matrix)[end] == 99
    rawintensities(matrix)[1, 1] = 100
    rawretentions(matrix)[1] = 0
    @test rawintensities(first(scans(series)))[1] == 10
    @test rawretentions(series)[1] == 1

    # A nonuniform SIM grid is supported; missing interior/edge columns remain.
    irregular = mscanmatrix(series, (1, 57, 60, 71, 85, 113); format)
    @test rawintensities(irregular) == [0 10 0 20 0 0; 0 5 0 0 30 0]

    filtered = mscanmatrix(series, 57:14:99; format, threshold=10)
    @test rawintensities(filtered) == [0 20 0 0; 0 0 30 0]
    @test rawmzvalues(filtered) == [57, 71, 85, 99]
    @test rawintensities(mscanmatrix(series, 57:14:99; format, threshold=30)) ==
        zeros(Int16, 2, 4)

    # Selection is explicit and never changes the surviving intensity values.
    selected = mscanmatrix(series, (57, 71); format, unmatched=:drop)
    @test rawintensities(selected) == [10 20; 5 0]
    @test rawintensities(mscanmatrix(series, (57,); format, unmatched=:drop)) ==
        reshape(Int16[10, 5], 2, 1)
    absent = mscanmatrix(series, (99, 113); format, unmatched=:drop)
    @test rawmzvalues(absent) == [99, 113]
    @test size(rawintensities(absent)) == (2, 2)
    @test iszero(rawintensities(absent))
    if format === SPARSE
        @test nnz(rawintensities(absent)) == 0
    end

    # Unmatched values below/equal to the threshold must still raise an error.
    for threshold in (0, 30, Inf)
        @test_throws ArgumentError mscanmatrix(series, (57, 71); format, threshold)
    end
    zero_signal = MassScanSeries([MassScan(1u"s", [57, 71], [0, -1])])
    @test_throws ArgumentError mscanmatrix(zero_signal, (57,); format)
    @test iszero(rawintensities(mscanmatrix(zero_signal, (57, 71); format)))
    @test rawintensities(mscanmatrix(zero_signal, (57, 71); format, threshold=-2)) == [0 -1]

    # Scans at other levels do not participate in matching or the output.
    levels = MassScanSeries([
        MassScan(1u"s", [57], [10]; level=1),
        MassScan(2u"s", [1000], [20]; level=2)])
    ms1 = mscanmatrix(levels, (57, 71); format)
    @test rawintensities(ms1) == [10 0]
    @test retentions(ms1) == [1]u"s"
    ms2 = mscanmatrix(levels, (999, 1000); format, target_level=2)
    @test rawintensities(ms2) == [0 20]
    @test retentions(ms2) == [2]u"s"
    @test level(ms2) == 2
    @test_throws ArgumentError mscanmatrix(levels, (57,); format, target_level=3)

    @test_throws ArgumentError mscanmatrix(series, 57:99; format, unmatched=:warn)
    for grid in ((), [], 57:56, (71, 57), (57, 57), (57, 0), (-1, 57),
                 (NaN, 57), (57, Inf), (-Inf, 57), (57, "71"), (57 + 0im, 71),
                 (57, missing), (57, 71u"Th"), Any[57u"Th", 71])
        @test_throws ArgumentError mscanmatrix(series, grid; format, unmatched=:drop)
    end
    @test_throws ArgumentError mscanmatrix(series, [57, 71, 85]u"Th"; format)
end

@testset "mscanmatrix exact numeric matching ($format)" for format in (DENSE, SPARSE)
    grid = 20.1:0.1:200.4
    series = MassScanSeries([MassScan(1.0, [grid[1], grid[2], grid[end]], Float32[1, 2, 3])])
    matrix = mscanmatrix(series, grid; format)
    @test rawmzvalues(matrix) == collect(grid)
    @test size(rawintensities(matrix)) == (1, length(grid))
    @test rawintensities(matrix)[1, [1, 2, end]] == [1, 2, 3]
    @test sum(rawintensities(matrix)) == 6
    @test eltype(rawintensities(matrix)) === Float32
    @test retentionunit(matrix) === nothing

    f32series = MassScanSeries([MassScan(1.0f0, Float32[57.1], Float32[10])])
    @test_throws ArgumentError mscanmatrix(f32series, (57.1,); format)
    matched = mscanmatrix(f32series, (57.1f0,); format)
    @test rawintensities(matched) == fill(10, 1, 1)
    @test eltype(rawmzvalues(matched)) === Float32
    @test eltype(rawretentions(matched)) === Float32
    @test rawintensities(mscanmatrix(f32series, (Float64(57.1f0),); format)) == fill(10, 1, 1)

    # No coercion of a recorded value into the grid's narrower numeric type.
    nearby = MassScanSeries([MassScan(1, [nextfloat(57.0)], [10])])
    @test_throws ArgumentError mscanmatrix(nearby, Float32[57]; format)
    @test iszero(rawintensities(mscanmatrix(nearby, Float32[57]; format, unmatched=:drop)))
    @test rawintensities(mscanmatrix(toy_series(), (57.0, 71.0, 85.0); format)) ==
        [10 20 0; 5 0 30]

    # Heterogeneous numeric grids must not round large integers through Float64.
    large = 2^53 + 1
    large_series = MassScanSeries([MassScan(1, [large], [10])])
    large_matrix = mscanmatrix(large_series, (1.0, large); format)
    @test rawmzvalues(large_matrix)[2] === large
    @test rawintensities(large_matrix) == [0 10]
    @test_throws ArgumentError mscanmatrix(large_series, (Float64(large),); format)

    rational = MassScanSeries([MassScan(1, [115//2], [7])])
    @test rawintensities(mscanmatrix(rational, (57.5,); format)) == fill(7, 1, 1)
end

@testset "mscanmatrix unit conversion ($format)" for format in (DENSE, SPARSE)
    series = MassScanSeries([
        MassScan(1.0f0u"s", [57, 71]u"Th", Float32[10, 20]u"pA"),
        MassScan(2.0f0u"s", [57, 85]u"Th", Float32[5, 30]u"pA")])
    for grid in ((57, 71, 85, 99), [57, 71, 85, 99]u"Th",
                 (57u"Th", 71u"Th", 85u"Th", 99u"Th"),
                 (57u"Th"):(14u"Th"):(99u"Th"),
                 (57000u"mTh", 71000u"mTh", 85000u"mTh", 99000u"mTh"),
                 (57u"Th", 71000u"mTh", 85u"Th", 99000u"mTh"),
                 Any[57u"Th", 71u"Th", 85u"Th", 99u"Th"])
        matrix = mscanmatrix(series, grid; format)
        @test mzvalues(matrix) == [57, 71, 85, 99]u"Th"
        @test rawintensities(matrix) == [10 20 0 0; 5 0 30 0]
        @test intensities(matrix) == Float32[10 20 0 0; 5 0 30 0]u"pA"
        @test mzunit(matrix) === u"Th"
        @test intensityunit(matrix) === u"pA"
        @test retentionunit(matrix) === u"s"
        @test eltype(rawintensities(matrix)) === Float32
        @test eltype(rawretentions(matrix)) === Float32
    end
    dropped = mscanmatrix(series, (57000u"mTh", 99000u"mTh"); format, unmatched=:drop)
    @test intensities(dropped) == Float32[10 0; 5 0]u"pA"
    @test_throws ArgumentError mscanmatrix(series, (57000u"mTh",); format, threshold=Inf)

    for grid in ((57u"s", 71u"s"), (57u"Th", 71u"s"),
                 (57u"Th", 71), Any[57, 71u"Th"], (57u"Th", 57000u"mTh"),
                 (71000u"mTh", 57u"Th"), (NaN*u"Th",), (0u"Th",),
                 (complex(57, 1)*u"Th",))
        @test_throws ArgumentError mscanmatrix(series, grid; format, unmatched=:drop)
    end

    floating_grid = (20.1u"Th"):(0.1u"Th"):(200.4u"Th")
    floating_series = MassScanSeries([
        MassScan(1u"s", [floating_grid[1], floating_grid[end]], [1, 2])])
    matrix = mscanmatrix(floating_series, floating_grid; format)
    @test mzvalues(matrix) == collect(floating_grid)
    @test sum(rawintensities(matrix)) == 3
    @test rawintensities(matrix)[1, [1, end]] == [1, 2]
end

@testset "mscanmatrix API" begin
    series = toy_series()
    @test rawintensities(mscanmatrix(series, 57:14:99)) isa Matrix
    @test_throws MethodError mscanmatrix(series)
    @test_throws MethodError mscanmatrix(series, SPARSE)
    @test_throws MethodError mscanmatrix(series, 57:14:99, SPARSE)
    @test_throws MethodError mscanmatrix(series, 57)
    @test_throws MethodError mscanmatrix(series, [57 71 85])
end

end # module
