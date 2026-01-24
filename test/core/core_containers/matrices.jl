module TestMatrices

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core/containers/matrices.jl
# ─────────────────────────────────────────────────────────────────────────────

using Test
using Unitful
using Unitful: s, pA
using JuChrom  # MassScanMatrix, AbstractMassScanMatrix

# Helpers to build small arrays
mkret(n) = collect(1.0:n)                    # Vector{Float64}
mkmzs(n) = collect(100.0:1.0:100.0+n-1)      # strictly increasing
mkints(r, c) = reshape(collect(1.0:(r*c)), r, c)

@testset "MassScanMatrix – inner constructor" begin
    ret = mkret(3)
    mzs = mkmzs(4)
    im  = mkints(length(ret), length(mzs))

    msm = MassScanMatrix{typeof(ret), typeof(u"s"),
                         typeof(mzs), Nothing,
                         typeof(im), typeof(pA),
                         Int,
                         typeof((vendor="Acme",)),
                         typeof((mode="DIA",)),
                         typeof((user="U",)),
                         typeof((sample="S",)),
                         Dict{String,Any}}(
        ret, u"s", mzs, nothing, im, pA, 2,
        (vendor="Acme",), (mode="DIA",), (user="U",), (sample="S",),
        Dict{String,Any}("k"=>"v"))

    # Stored fields
    @test msm.retentions === ret
    @test msm.retentionunit == u"s"
    @test msm.mzvalues === mzs
    @test msm.mzunit === nothing
    @test msm.intensities === im
    @test msm.intensityunit == pA
    @test msm.level === 2
    @test msm.instrument == (vendor="Acme",)
    @test msm.acquisition == (mode="DIA",)
    @test msm.user == (user="U",)
    @test msm.sample == (sample="S",)
    @test msm.extras == Dict("k"=>"v")

    # Subtyping relation: parameterization matches {R,M,I}
    @test msm isa AbstractMassScanMatrix{typeof(u"s"), Nothing, typeof(pA)}

    # Error paths enforced by inner constructor
    @test_throws ArgumentError MassScanMatrix{typeof(Float64[]), Nothing,
                                              typeof(mzs), Nothing,
                                              typeof(im), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        Float64[], nothing, mzs, nothing, im, nothing, 1,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    @test_throws ArgumentError MassScanMatrix{typeof([1.0, NaN]), Nothing,
                                              typeof(mzs), Nothing,
                                              typeof(im), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        [1.0, NaN], nothing, mzs, nothing, im, nothing, 1,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    @test_throws ArgumentError MassScanMatrix{typeof(ret), Nothing,
                                              typeof(Float64[]), Nothing,
                                              typeof(im), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        ret, nothing, Float64[], nothing, im, nothing, 1,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    @test_throws ArgumentError MassScanMatrix{typeof(ret), Nothing,
                                              typeof([100.0, 100.0]), Nothing,
                                              typeof(im[:, 1:2]), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        ret, nothing, [100.0, 100.0], nothing, im[:, 1:2], nothing, 1,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    @test_throws ArgumentError MassScanMatrix{typeof(ret), Nothing,
                                              typeof([-1.0, 2.0]), Nothing,
                                              typeof(im[:,1:2]), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        ret, nothing, [-1.0, 2.0], nothing, im[:,1:2], nothing, 1,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    @test_throws ArgumentError MassScanMatrix{typeof(ret), Nothing,
                                              typeof(mzs), Nothing,
                                              typeof(zeros(0, length(mzs))), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        ret, nothing, mzs, nothing, zeros(0, length(mzs)), nothing, 1,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    @test_throws ArgumentError MassScanMatrix{typeof(ret), Nothing,
                                              typeof(mzs), Nothing,
                                              typeof([1.0 2.0; 3.0 Inf]), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        ret, nothing, mzs[1:2], nothing, [1.0 2.0; 3.0 Inf], nothing, 1,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    # Shape mismatches
    @test_throws ArgumentError MassScanMatrix{typeof(ret[1:2]), Nothing,
                                              typeof(mzs), Nothing,
                                              typeof(im), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        ret[1:2], nothing, mzs, nothing, im, nothing, 1,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())
    @test_throws ArgumentError MassScanMatrix{typeof(ret), Nothing,
                                              typeof(mzs[1:2]), Nothing,
                                              typeof(im), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        ret, nothing, mzs[1:2], nothing, im, nothing, 1,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    # Level must be ≥ 1
    @test_throws ArgumentError MassScanMatrix{typeof(ret), Nothing,
                                              typeof(mzs), Nothing,
                                              typeof(im), Nothing,
                                              Int, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
                                              Dict{String,Any}}(
        ret, nothing, mzs, nothing, im, nothing, 0,
        NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())
end

@testset "MassScanMatrix – outer constructor (unitless arrays + explicit units)" begin
    ret = mkret(2)
    mzs = mkmzs(3)
    im  = mkints(length(ret), length(mzs))

    msm = MassScanMatrix(ret, u"s", mzs, nothing, im, pA;
                         level=3,
                         instrument=(vendor="A",),
                         acquisition=(mode="DDA",),
                         user=(u="x",),
                         sample=(s="y",),
                         extras=Dict(SubString("k",1,1)=>1))

    @test msm.retentions === ret
    @test msm.retentionunit == u"s"
    @test msm.mzvalues === mzs
    @test msm.mzunit === nothing
    @test msm.intensities === im
    @test msm.intensityunit == pA
    @test msm.level === 3
    @test msm.instrument == (vendor="A",)
    @test msm.acquisition == (mode="DDA",)
    @test msm.user == (u="x",)
    @test msm.sample == (s="y",)
    @test msm.extras isa Dict{String,Any}
    @test first(keys(msm.extras)) isa String

    # Subtyping relation holds
    @test msm isa AbstractMassScanMatrix{typeof(u"s"), Nothing, typeof(pA)}
end

@testset "MassScanMatrix – outer constructor (arrays may carry Unitful quantities)" begin
    ret = [1.0, 2.0]u"s"                          # unitful vector
    mzs = [100.0, 200.0, 300.0]                   # unitless vector
    im  = [1.0 2.0 3.0; 4.0 5.0 6.0]u"pA"         # unitful matrix

    msm = MassScanMatrix(ret, mzs, im; instrument=(id=1,))
    @test msm.retentions == [1.0, 2.0]       # stripped
    @test msm.retentionunit == u"s"
    @test msm.mzvalues == [100.0, 200.0, 300.0]
    @test msm.mzunit === nothing
    @test msm.intensities == [1.0 2.0 3.0; 4.0 5.0 6.0]
    @test msm.intensityunit == pA
    @test msm.level === 1
    @test msm.instrument == (id=1,)
    @test msm.acquisition == NamedTuple()
    @test msm.user == NamedTuple()
    @test msm.sample == NamedTuple()
    @test msm.extras == Dict{String,Any}()

    # extras with non-String keys are normalized
    s = "zz"; k = SubString(s, 1, 1)
    msm2 = MassScanMatrix(ret, mzs, im; extras=Dict(k=>:ok))
    @test msm2.extras isa Dict{String,Any}
    @test first(keys(msm2.extras)) isa String

    # Error cases via outer ctor
    @test_throws ArgumentError MassScanMatrix([1.0]u"s", [100.0, 200.0], [1.0 2.0 3.0])  # row mismatch
    @test_throws ArgumentError MassScanMatrix([1.0, 2.0], [100.0, 100.0], [1.0 2.0; 3.0 4.0]) # non-strict mz
    @test_throws ArgumentError MassScanMatrix([1.0, 2.0], [-1.0, 2.0], [1.0 2.0; 3.0 4.0])     # nonpositive mz
    @test_throws ArgumentError MassScanMatrix([1.0, 2.0], [100.0, 200.0], [1.0 2.0; 3.0 4.0]; level=0) # level
end

@testset "MassScanMatrix – subtraction" begin
    ret = mkret(2)
    mzs = mkmzs(2)
    im1 = [1.0 2.0; 3.0 4.0]
    im2 = [0.5 1.0; 1.5 2.0]

    msm1 = MassScanMatrix(ret, mzs, im1;
                          level=2,
                          instrument=(vendor="A",),
                          acquisition=(mode="DDA",),
                          user=(u="x",),
                          sample=(s="y",),
                          extras=Dict("k"=>1))
    msm2 = MassScanMatrix(ret, mzs, im2;
                          level=2,
                          instrument=(vendor="A",),
                          acquisition=(mode="DDA",),
                          user=(u="x",),
                          sample=(s="y",),
                          extras=Dict("k"=>1))

    msm_diff = msm1 - msm2
    @test intensities(msm_diff) == im1 .- im2
    @test retentions(msm_diff) == ret
    @test mzvalues(msm_diff) == mzs
    @test retentions(msm_diff) !== ret
    @test mzvalues(msm_diff) !== mzs
    @test extras(msm_diff) == extras(msm1)
    @test extras(msm_diff) !== extras(msm1)

    msm_ret_mismatch = MassScanMatrix(ret .+ 1, mzs, im1;
                                      level=2,
                                      instrument=(vendor="A",),
                                      acquisition=(mode="DDA",),
                                      user=(u="x",),
                                      sample=(s="y",),
                                      extras=Dict("k"=>1))
    @test_throws DimensionMismatch msm1 - msm_ret_mismatch

    msm_extra_mismatch = MassScanMatrix(ret, mzs, im1;
                                        level=2,
                                        instrument=(vendor="A",),
                                        acquisition=(mode="DDA",),
                                        user=(u="x",),
                                        sample=(s="y",),
                                        extras=Dict("k"=>2))
    @test_throws DimensionMismatch msm1 - msm_extra_mismatch
end

@testset "MassScanMatrix show methods" begin
    msm_plain = MassScanMatrix([1.0, 2.0], [100.0, 200.0], [1.0 2.0; 3.0 4.0])
    out_plain = sprint(show, msm_plain)
    @test occursin("MassScanMatrix with 2 scans", out_plain)
    @test occursin("├─ Retention:", out_plain)
    @test occursin("Range: 1.0 to 2.0 (unitless)", out_plain)
    @test occursin("├─ M/Z values:", out_plain)
    @test occursin("Total data points: 2", out_plain)
    @test occursin("└─ Intensity:", out_plain)
    @test occursin("MS Level: 1", out_plain)
    @test !occursin("Scan type", out_plain)

    msm_meta = MassScanMatrix([1.0, 2.0]u"s", [100.0, 200.0]u"Th",
                              [1.0 2.0; 3.0 4.0]u"pA";
                              instrument=(vendor="A",))
    out_meta = sprint(io -> show(io, MIME"text/plain"(), msm_meta))
    @test occursin("MassScanMatrix with 2 scans", out_meta)
    @test occursin("Range: 1.0 to 2.0 (s)", out_meta)
    @test occursin("Range: 100.0 to 200.0 (Th)", out_meta)
    @test occursin("Range: 1.0 to 4.0 (pA)", out_meta)
    @test occursin("├─ Intensity:", out_meta)
    @test occursin("Instrument", out_meta)
    @test occursin("vendor = A", out_meta)
    @test !occursin("Scan type", out_meta)
end

end # module
