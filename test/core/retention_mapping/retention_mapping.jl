module TestRetentionMapping

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests for ./core./retention_mapping.jl
# ─────────────────────────────────────────────────────────────────────────────

using Test
using JuChrom
import Base: broadcastable
using BSplineKit
using LinearAlgebra
using Logging
using Printf
using Random
using SparseArrays
using Statistics
using Unitful
using Unitful: AbstractQuantity
using JuChrom: tune_lambda_for_monotonic_spline

# ─────────────────────────────────────────────────────────────────────────────
# Unit tests
# ─────────────────────────────────────────────────────────────────────────────

# ── applymap(::RetentionMapper, ::Union{Real, AbstractQuantity}) ─────────────────────

@testset "applymap(rm, retention::Union{Real,AbstractQuantity})" begin

    # Helper: constant-spline mapper → always maps to rB_min (predictable)
    make_mapper(; unitful::Bool) = begin
        rA = [0.0, 0.5, 1.0]
        rB = [10.0, 20.0, 30.0]
        rA_min, rA_max = 0.0, 1.0
        rB_min, rB_max = 10.0, 30.0
        rA_norm_min, rA_norm_max = 0.0, 1.0
        rB_pred_min, rB_pred_max = rB_min, rB_max
        rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

        order = BSplineOrder(4)
        knots = collect(LinRange(0.0, 1.0, 8))
        B = BSplineBasis(order, knots)
        coefs  = zeros(length(B))  # constant spline ⇒ rB_norm ≡ 0
        spline = Spline(B, coefs)

        rA_unit = unitful ? u"minute" : nothing
        rB_unit = unitful ? u"Th"     : nothing

        JuChrom.RetentionMapper(
            rA, rA_unit,
            rA_min, rA_max,
            rA_norm_min, rA_norm_max,
            rB, rB_unit,
            rB_min, rB_max,
            rB_pred_min, rB_pred_max,
            rB_pred_norm_min, rB_pred_norm_max,
            knots, coefs, spline,
            Dict{String, Any}()
        )
    end

    # Unitless mapper: Real input required; output is raw number = rB_min
    rm0 = make_mapper(unitful=false)
    @test applymap(rm0, 0.25) ≈ rm0.rB_min  # constant spline → rB_min
    @test applymap(rm0, -1.0) ≈ rm0.rB_min  # extrapolation still returns rB_min
    @test_throws ArgumentError applymap(rm0, 0.25u"minute")  # unitless mapper rejects AbstractQuantity

    # Unitful mapper: AbstractQuantity input required; output carries rB_unit (or requested unit)
    rmU = make_mapper(unitful=true)
    y1 = applymap(rmU, 0.5u"minute")  # accepted; time unit OK
    @test y1 ≈ rmU.rB_min * u"Th"  # constant spline ⇒ rB_min with B-unit
    y2 = applymap(rmU, 30.0u"s")  # convertible time units; ustrip handles conversion
    @test y2 ≈ rmU.rB_min * u"Th"

    # Explicit unit override returns converted AbstractQuantity
    y3 = applymap(rmU, 0.5u"minute"; unit=u"Th")  # same unit → equal value
    @test y3 ≈ rmU.rB_min * u"Th"

    # Validation errors: wrong presence/absence of units
    # unitful mapper rejects unitless input
    @test_throws ArgumentError applymap(rmU, 0.5)  
end

# ── applymap(::RetentionMapper, ::MassScanMatrix) ────────────────────────────

@testset "applymap(rmap, ::MassScanMatrix)" begin
    # Minimal strictly increasing RetentionMapper (unitless; no optimization)
    rA = [0.0, 0.5, 1.0]
    rB = [10.0, 20.0, 30.0]
    rA_min, rA_max = 0.0, 1.0
    rB_min, rB_max = 10.0, 30.0
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 10))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # monotone-ish spline
    spline = Spline(B, coefs)

    rmap = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max,
        rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max,
        rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # MassScanMatrix: rows=scans, cols=m/z (so #rows == length(ret))
    ret = [0.25, 0.50, 0.75]  # 3 scans ⇒ 3 rows
    mz  = [50.0, 51.0, 52.0, 53.0]  # 4 m/z bins ⇒ 4 cols
    X   = [100.0 110.0 120.0 130.0;  # row 1 (scan at 0.25)
           200.0 210.0 220.0 230.0;  # row 2 (scan at 0.50)
           300.0 310.0 320.0 330.0]  # row 3 (scan at 0.75) -> size 3×4

    msm = MassScanMatrix(
        ret,
        nothing,  # retention unit
        mz,
        nothing,  # m/z unit
        X,
        nothing,  # intensity unit
    )

    # Apply mapping
    out = applymap(rmap, msm)

    # Expected mapped retentions (vectorized scalar path)
    new_rets = applymap.(rmap, ret)
    @test retentions(out) ≈ new_rets
    @test retentionunit(out) === nothing

    # Jacobians: one per scan/row
    J = rawderivmap.(rmap, ret; rB_unit=nothing)
    @test all(isfinite, J) && all(>(0), J)

    # Intensities scaled row-wise by 1 / Jacobian for each corresponding scan
    X_out = rawintensities(out)
    @test size(X_out) == size(X)
    @test X_out[1, :] ≈ X[1, :] ./ J[1]
    @test X_out[2, :] ≈ X[2, :] ./ J[2]
    @test X_out[3, :] ≈ X[3, :] ./ J[3]

    # m/z values and units propagated unchanged
    @test rawmzvalues(out) == rawmzvalues(msm)
    @test mzunit(out) === mzunit(msm)
    @test intensityunit(out) === intensityunit(msm)
end

# ── applymap(::RetentionMapper, ::MassScanSeries) ────────────────────────────

@testset "applymap(rmap, ::MassScanSeries)" begin

    # Minimal, valid RetentionMapper with strictly increasing spline (no optimization)
    rA = [0.0, 0.5, 1.0]
    rB = [10.0, 20.0, 30.0]
    rA_min, rA_max = 0.0, 1.0
    rB_min, rB_max = 10.0, 30.0
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 10))
    B = BSplineBasis(order, knots)
    # Monotone-ish spline: increasing coefficients → positive derivative everywhere
    coefs = collect(range(0.0, 1.0, length(B)))
    spline = Spline(B, coefs)

    rmap = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max,
        rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max,
        rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Build a tiny MassScanSeries with two scans (unitless for simplicity)
    mz = [50.0, 51.0, 52.0]
    ints1 = [100.0, 200.0, 300.0]
    ints2 = [150.0, 250.0, 350.0]

    s1 = MassScan(0.25, nothing, mz, nothing, ints1, nothing)
    s2 = MassScan(0.75, nothing, mz, nothing, ints2, nothing)
    series = MassScanSeries([s1, s2])

    # Run mapping
    out = applymap(rmap, series)

    # Expected mapped retentions and Jacobian-corrected intensities
    r_new_1 = applymap(rmap, 0.25)  # per-scan call uses same path
    r_new_2 = applymap(rmap, 0.75)

    j1 = rawderivmap(rmap, 0.25; rB_unit=nothing)  # Jacobian dB/dA (unitless case)
    j2 = rawderivmap(rmap, 0.75; rB_unit=nothing)

    # Check new retentions
    @test retention(out[1]) ≈ r_new_1
    @test retention(out[2]) ≈ r_new_2

    # Check intensity Jacobian scaling (I_new = I_old / jacobian)
    @test intensities(out[1]) ≈ (ints1 ./ j1)
    @test intensities(out[2]) ≈ (ints2 ./ j2)

    # Units of new retentions: still unitless for unitless mapper
    @test retentionunit(out[1]) === nothing
    @test retentionunit(out[2]) === nothing
end

# ── Base.broadcastable(::RetentionMapper) ─────────────────────────────────────

@testset "Base.broadcastable(::RetentionMapper)" begin
    # Minimal, manual RetentionMapper (no optimization; tiny identity-like spline)
    rA = [0.0, 0.5, 1.0]
    rB = [0.0, 0.5, 1.0]

    rA_min, rA_max = 0.0, 1.0
    rB_min, rB_max = 0.0, 1.0

    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    # Tiny cubic B-spline basis over [0,1] and a simple coefficient vector
    order = BSplineOrder(4)  # cubic
    knots = collect(LinRange(0.0, 1.0, 8))  # uniform knots for a tiny basis
    B = BSplineBasis(order, knots)
    ncoefs = length(B)  # dimension of the basis
    coefs = zeros(ncoefs)  # constant spline
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, nothing,  # rA, rA_unit
        rA_min, rA_max,
        rA_norm_min, rA_norm_max,
        rB, nothing,  # rB, rB_unit
        rB_min, rB_max,
        rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline, 
        Dict{String, Any}()
    )

    # Broadcastable wraps as Ref ⇒ scalar-like broadcasting
    b = broadcastable(rm)
    @test b isa Ref
    @test b[] === rm

    # Broadcasting with identity uses scalar semantics and returns the object itself
    @test identity.(rm) === rm

    # Broadcasting with `===` stays scalar (not an array)
    @test (rm .=== rm) === true
end

# ── derivinvmap(rm::RetentionMapper, retention::Union{<:Real, <:AbstractQuantity}) ───

@testset "derivinvmap(rm, retention::Union{Real,AbstractQuantity})" begin

    # Build a small monotone mapper (minutes → Th)
    rA = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
    rB = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]u"Th"
    rm = fitmap(rA, rB)

    # Unitful input required? (API accepts AbstractQuantity or Real; 
    # with unitful mapper, AbstractQuantity is the valid path)
    @test_throws ArgumentError derivinvmap(rm, 350.0)  # unitless input rejected for unitful rm

    # Success: returns AbstractQuantity with unit rA_unit / rB_unit
    dAB = derivinvmap(rm, 350.0u"Th")  # dA/dB at B=350 Th
    @test isfinite(ustrip(dAB)) && ustrip(dAB) > 0
    @test Unitful.unit(dAB) == u"minute"/u"Th"

    # Reciprocity with forward derivative: derivinvmap(B) ≈ 1 / derivmap(A) at matched A,B
    A = 5.0u"minute"
    Bv = applymap(rm, A)  # map A → B
    dBA = derivmap(rm, A)  # dB/dA (unit = Th/min)
    dAB2 = derivinvmap(rm, Bv)  # dA/dB (unit = min/Th)
    @test ustrip(dAB2 * dBA) ≈ 1.0 atol=1e-10

    # Unit conversion via rA_unit / rB_unit arguments (seconds for A; keep Th for B)
    dAB_s = derivinvmap(rm, 350.0u"Th"; rA_unit=u"s", rB_unit=u"Th")
    @test Unitful.unit(dAB_s) == u"s"/u"Th"
    @test ustrip(dAB_s) ≈ ustrip(u"s"/u"Th", dAB)  # same physical value in requested units

    # Extrapolation branches (silence @info logs)
    with_logger(SimpleLogger(devnull, Logging.Warn)) do
        d_lo = derivinvmap(rm, 50.0u"Th"; warn=true)
        d_hi = derivinvmap(rm, 800.0u"Th"; warn=true)
        @test isfinite(ustrip(d_lo)) && isfinite(ustrip(d_hi))
        @test Unitful.unit(d_lo) == u"minute"/u"Th"
        @test Unitful.unit(d_hi) == u"minute"/u"Th"
    end

    # Unitless mapper variant: accepts Real; returns unitless; rejects AbstractQuantity
    rA0 = ustrip.(rA)  # drop units
    rB0 = ustrip.(rB)
    rm0 = fitmap(rA0, rB0)
    @test_throws ArgumentError derivinvmap(rm0, 350.0u"Th")
    d0 = derivinvmap(rm0, 350.0)
    @test (d0 isa Real) && isfinite(d0) && d0 > 0
end

# ── derivmap(rm, retention::Union{Real,AbstractQuantity}) ────────────────────────────

@testset "derivmap(rm, retention::Union{Real,AbstractQuantity})" begin
    # Helper: simple monotone mapper (tiny spline; no optimization)
    make_mapper(; unitful::Bool) = begin
        rA = [0.0, 0.5, 1.0]
        rB = [10.0, 20.0, 30.0]
        rA_min, rA_max = 0.0, 1.0
        rB_min, rB_max = 10.0, 30.0
        rA_norm_min, rA_norm_max = 0.0, 1.0
        rB_pred_min, rB_pred_max = rB_min, rB_max
        rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

        order = BSplineOrder(4)
        knots = collect(LinRange(0.0, 1.0, 10))
        B = BSplineBasis(order, knots)
        coefs = collect(range(0.0, 1.0, length(B)))  # increasing → positive slope
        spline = Spline(B, coefs)

        rA_unit = unitful ? u"minute" : nothing
        rB_unit = unitful ? u"Th"     : nothing

        JuChrom.RetentionMapper(
            rA, rA_unit,
            rA_min, rA_max,
            rA_norm_min, rA_norm_max,
            rB, rB_unit,
            rB_min, rB_max,
            rB_pred_min, rB_pred_max,
            rB_pred_norm_min, rB_pred_norm_max,
            knots, coefs, spline,
            Dict{String, Any}()
        )
    end

    # Unitful mapper: input must be AbstractQuantity; result has unit rB_unit/rA_unit
    rmU = make_mapper(unitful=true)
    q = derivmap(rmU, 0.5u"minute")  # ok, AbstractQuantity input
    @test isfinite(ustrip(q)) && ustrip(q) > 0
    @test Unitful.unit(q) == u"Th"/u"minute"  # unit = rB_unit / rA_unit
    # Consistency with rawderivmap once units are stripped
    raw = rawderivmap(rmU, 0.5u"minute")  # unitless numeric
    @test ustrip(u"Th"/u"minute", q) ≈ raw

    # Unit validation: unitful mapper rejects unitless input
    @test_throws ArgumentError derivmap(rmU, 0.5)

    # Unitless mapper: accepts Real; returns unitless; rejects AbstractQuantity
    rm0 = make_mapper(unitful=false)
    d0 = derivmap(rm0, 0.5)
    @test (d0 isa Real) && isfinite(d0) && d0 > 0
    @test_throws ArgumentError derivmap(rm0, 0.5u"minute")
end

# ── extras(rm::AbstractRetentionMapper) ──────────────────────────────────────

@testset "extras(rm::AbstractRetentionMapper)" begin

    # Minimal, unitless RetentionMapper (no optimization)
    rA = [0.0, 0.5, 1.0]
    rB = [10.0, 20.0, 30.0]
    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 10))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # monotone-ish slope
    spline = Spline(B, coefs)

    meta = Dict{String, Any}("cal_file" => "test.cal", "instrument" => "GC-MS")

    rm = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        copy(meta),
    )

    # Extras returns the metadata Dict{String, Any}
    ex = extras(rm)
    @test ex isa Dict{String, Any}
    @test haskey(ex, "cal_file") && ex["cal_file"] == "test.cal"
    @test haskey(ex, "instrument") && ex["instrument"] == "GC-MS"

    # Metadata is mutable and mutations are reflected
    ex["analysis_date"] = "2025-10-13"
    @test haskey(extras(rm), "analysis_date") && extras(rm)["analysis_date"] == "2025-10-13"

    # Mapping is unaffected by metadata changes
    y_before = applymap(rm, 0.6)
    ex["note"] = "does not affect mapping"
    y_after  = applymap(rm, 0.6)
    @test y_after ≈ y_before atol=1e-12  # unitless mapper → unitless tolerance
end

# ── invmap(rm, retention::Union{<:Real, <:AbstractQuantity}) ─────────────────────────

@testset "invmap(rm, retention::Union{<:Real, <:AbstractQuantity})" begin

    # Synthetic monotone data (minutes → Th)
    rA = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
    rB = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]u"Th"
    rm = fitmap(rA, rB)

    # Basic success path
    xB = 350.0u"Th"
    yA = invmap(rm, xB)
    @test yA isa AbstractQuantity
    @test Unitful.unit(yA) == u"minute"
    @test isfinite(ustrip(yA))
    @test yA ≈ 5.312425320906996u"minute" atol=1e-10u"minute"  # regression target from docs

    # Unit conversion
    yA_s = invmap(rm, xB; unit=u"s")
    @test yA_s isa AbstractQuantity
    @test Unitful.unit(yA_s) == u"s"
    @test yA_s ≈ (ustrip(yA) * 60.0)u"s" atol=1e-8u"s"

    # Broadcasting over a small vector of outputs
    xsB = [250.0u"Th", 450.0u"Th"]
    ysA = invmap.(rm, xsB)
    @test all(yi -> yi isa AbstractQuantity && Unitful.unit(yi) == u"minute", ysA)
    @test isfinite(sum(ustrip, ysA))
    @test ysA ≈ [3.2283166885897936, 8.08686182716771]u"minute" atol=1e-10u"minute"

    # Extrapolation below/above domain with warn=true (just exercise the branch)
    with_logger(SimpleLogger(devnull, Logging.Warn)) do
        y_low  = invmap(rm, 50.0u"Th"; warn=true)
        y_high = invmap(rm, 800.0u"Th"; warn=true)
        @test isfinite(ustrip(y_low))
        @test isfinite(ustrip(y_high))
    end

    # Error: missing units when rm expects units
    @test_throws ArgumentError invmap(rm, 350.0)

    # Round-trip consistency with forward mapping
    tA = 5.0u"minute"
    bB = applymap(rm, tA)
    @test invmap(rm, bB) ≈ tA atol=1e-8u"minute"

    # Respect requested output units in round-trip (seconds)
    bB2 = applymap(rm, 7.5u"minute")
    @test invmap(rm, bB2; unit=u"s") ≈ (7.5 * 60)u"s" atol=1e-6u"s"
end

# ── invmapmax(rm; unit=…) ────────────────────────────────────────────────────

@testset "invmapmax(rm; unit=…)" begin

    # Minimal RetentionMapper (A unitful minutes → B unitful Th)
    rA = [1.0, 2.0, 3.0]
    rB = [100.0, 200.0, 300.0]

    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # simple increasing spline
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String,Any}(),
    )

    # Basic success path: unitful mapper returns B-domain maximum with its unit
    @test invmapmax(rm) ≈ rB_max * u"Th" atol=1e-12u"Th"

    # Unit conversion on output
    @test invmapmax(rm; unit=u"kTh") ≈ (rB_max/1000) * u"kTh" atol=1e-12u"kTh"

    # Unitless variant
    rm0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String,Any}(),
    )

    # Unitless mapper returns raw numeric; asking for a unit should error
    @test invmapmax(rm0) ≈ rB_max atol=1e-12
    @test_throws ArgumentError invmapmax(rm0; unit=u"Th")
end

# ── invmapmin(rm; unit=…) ────────────────────────────────────────────────────

@testset "invmapmin(rm; unit=…)" begin

    # Minimal RetentionMapper (A unitful minutes → B unitful Th)
    rA = [1.0, 2.0, 3.0]
    rB = [100.0, 200.0, 300.0]

    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # simple increasing spline
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Basic success path: unitful mapper returns B-domain minimum with its unit
    @test invmapmin(rm) ≈ rB_min * u"Th" atol=1e-12u"Th"

    # Unit conversion on output
    @test invmapmin(rm; unit=u"kTh") ≈ (rB_min/1000) * u"kTh" atol=1e-12u"kTh"

    # Unitless variant
    rm0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Unitless mapper returns raw numeric; asking for a unit should error
    @test invmapmin(rm0) ≈ rB_min atol=1e-12
    @test_throws ArgumentError invmapmin(rm0; unit=u"Th")
end

# ── mapmax(rm; unit=…) ──────────────────────────────────────────────────────

@testset "mapmax(rm; unit=…)" begin

    # Build two tiny mappers sharing geometry; one unitful (minute→Th), one unitless
    rA = [1.2, 2.5, 4.1]
    rB = [100.0, 200.0, 300.0]

    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    rm0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String,Any}(),
    )

    # Unitful: returns AbstractQuantity with rm.rA_unit
    @test mapmax(rm) ≈ 4.1u"minute" atol=1e-12u"minute"

    # Unitful with conversion
    @test mapmax(rm, unit=u"s") ≈ 246.0u"s" atol=1e-10u"s"

    # Unitless: returns raw numeric
    @test mapmax(rm0) ≈ 4.1 atol=1e-12

    # Unitless + unit argument → error
    @test_throws ArgumentError mapmax(rm0, unit=u"s")
end

# ── mapmin(rm; unit=…) ──────────────────────────────────────────────────────

@testset "mapmin(rm; unit=…)" begin

    # Build two tiny mappers sharing geometry; one unitful (minute→Th), one unitless
    rA = [1.2, 2.5, 4.1]
    rB = [100.0, 200.0, 300.0]

    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}(),
    )

    rm0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Unitful: returns AbstractQuantity with rm.rA_unit
    @test mapmin(rm) ≈ 1.2u"minute" atol=1e-12u"minute"

    # Unitful with conversion
    @test mapmin(rm, unit=u"s") ≈ 72.0u"s" atol=1e-10u"s"

    # Unitless: returns raw numeric
    @test mapmin(rm0) ≈ 1.2 atol=1e-12

    # Unitless + unit argument → error (per API)
    @test_throws ArgumentError mapmin(rm0, unit=u"s")
end

# ── rawapplymap(rm, retention::AbstractQuantity) ────────────────────────────────────

@testset "rawapplymap(rm, retention::AbstractQuantity)" begin
    # Minimal strictly increasing mapper (minutes → Th), no optimizer needed
    rA = [0.0, 0.5, 1.0]  # domain A (unitful in mapper)
    rB = [10.0, 20.0, 30.0]  # domain B (unitful in mapper)
    rA_min, rA_max = 0.0, 1.0
    rB_min, rB_max = 10.0, 30.0
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0
    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 10))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # monotone-ish slope
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}(),
    )

    # Basic success path: returns a unitless number
    xA = 0.4u"minute"
    y_raw = rawapplymap(rm, xA)
    @test y_raw isa Real && isfinite(y_raw)

    # Consistency with applymap after stripping units
    y_unitful = applymap(rm, xA)  # has u"Th"
    @test y_raw ≈ ustrip(y_unitful)

    # Explicit unit conversion via `unit` argument
    y_raw_kTh = rawapplymap(rm, xA; unit=u"kTh")  # numeric in kTh
    @test y_raw_kTh ≈ ustrip(u"kTh", y_unitful)

    # Validation: unitless input rejected for unitful mapper
    @test_throws ArgumentError rawapplymap(rm, 0.4)

    # Unitless mapper variant accepts AbstractQuantity? → should error per docs
    rm0 = JuChrom.RetentionMapper(
        rA, nothing, rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing, rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline, Dict{String, Any}()
    )
    @test_throws ArgumentError rawapplymap(rm0, 0.4u"minute")
end

# ── rawderivinvmap(rm, retention::Union{Real,AbstractQuantity}) ─────────────────────

@testset "rawderivinvmap(rm, retention::Union{Real,AbstractQuantity})" begin

    # Build a small, strictly increasing mapper (unitful)
    rA = [0.0, 0.5, 1.0]
    rB = [10.0, 20.0, 30.0]
    rA_min, rA_max = 0.0, 1.0
    rB_min, rB_max = 10.0, 30.0
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 10))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # monotone-ish → positive derivative
    spline = Spline(B, coefs)

    rmU = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max,
        rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max,
        rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Unit validation (matches doc/implementation): unitful rm requires AbstractQuantity input
    @test_throws ArgumentError rawderivinvmap(rmU, 20.0)

    # Finite and positive derivative for a valid AbstractQuantity input (B domain)
    dinv = rawderivinvmap(rmU, 20.0u"Th")
    @test isfinite(dinv) && dinv > 0

    # Check reciprocal relationship with rawderivmap at matched A/B points
    A = 0.6u"minute"
    Bv = applymap(rmU, A)  # map A→B (AbstractQuantity out)
    dAB = rawderivmap(rmU, A)  # dB/dA (unitless)
    dBA = rawderivinvmap(rmU, Bv)  # dA/dB (unitless)
    @test dBA ≈ 1 / dAB  # reciprocals

    # Unit parameter effects: change rA_unit/rB_unit and confirm proper rescaling/stripping.
    dBA_s = rawderivinvmap(rmU, Bv; rA_unit=u"s", rB_unit=u"Th")  # input in Th, output stripped to s/Th
    @test dBA_s ≈ dBA * 60  # A in seconds ⇒ derivative ×60 compared to minutes

    # Unitless mapper variant: must receive Real; returns unitless number.
    rm0 = JuChrom.RetentionMapper(
        rA, nothing, rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing, rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline, Dict{String, Any}()
    )
    @test_throws ArgumentError rawderivinvmap(rm0, 20.0u"Th")
    d0 = rawderivinvmap(rm0, 20.0)
    @test isfinite(d0) && d0 > 0
end

# ── rawderivmap(::RetentionMapper, ::Union{Real,AbstractQuantity}) ──────────────────

@testset "rawderivmap(rm, retention::Union{Real,AbstractQuantity})" begin
    # Helper: simple monotone mapper (no optimization; tiny spline)
    make_mapper(; unitful::Bool) = begin
        rA = [0.0, 0.5, 1.0]
        rB = [10.0, 20.0, 30.0]
        rA_min, rA_max = 0.0, 1.0
        rB_min, rB_max = 10.0, 30.0
        rA_norm_min, rA_norm_max = 0.0, 1.0
        rB_pred_min, rB_pred_max = rB_min, rB_max
        rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

        order = BSplineOrder(4)
        knots = collect(LinRange(0.0, 1.0, 10))
        B = BSplineBasis(order, knots)
        # Increasing coefficients → positive derivative everywhere (monotone-ish)
        coefs = collect(range(0.0, 1.0, length(B)))
        spline = Spline(B, coefs)

        rA_unit = unitful ? u"minute" : nothing
        rB_unit = unitful ? u"Th"     : nothing

        JuChrom.RetentionMapper(
            rA, rA_unit,
            rA_min, rA_max,
            rA_norm_min, rA_norm_max,
            rB, rB_unit,
            rB_min, rB_max,
            rB_pred_min, rB_pred_max,
            rB_pred_norm_min, rB_pred_norm_max,
            knots, coefs, spline,
            Dict{String, Any}()
        )
    end

    # Unitless mapper: must receive Real, returns finite positive number
    rm0 = make_mapper(unitful=false)
    d0 = rawderivmap(rm0, 0.5)  # unitless input accepted
    @test isfinite(d0) && d0 > 0
    @test_throws ArgumentError rawderivmap(rm0, 0.5u"minute")  # unitless mapper rejects AbstractQuantity

    # Unitful mapper: must receive AbstractQuantity; derivative finite and positive
    rmU = make_mapper(unitful=true)
    d1 = rawderivmap(rmU, 0.5u"minute")  # accepted
    @test isfinite(d1) && d1 > 0
    @test_throws ArgumentError rawderivmap(rmU, 0.5)  # unitful mapper rejects unitless input

    # Changing input unit rescales dB/dA inversely with A's unit scale
    # (minute → second multiplies A by 60 ⇒ derivative divided by 60)
    d2 = rawderivmap(rmU, 30.0u"s"; rA_unit=u"s")  # same physical input; different unit
    @test d2 ≈ d1 / 60
end

# ── rawinvmap(rm, retention::AbstractQuantity) ──────────────────────────────────────

@testset "rawinvmap(rm, retention::AbstractQuantity)" begin

    # Build a tiny, strictly increasing mapper (minutes → Th), no optimization
    rA = [0.0, 0.5, 1.0]
    rB = [10.0, 20.0, 30.0]
    rA_min, rA_max = 0.0, 1.0
    rB_min, rB_max = 10.0, 30.0
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 10))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # monotone-ish slope
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Basic success: AbstractQuantity in B-domain → raw (unitless) number in A-units
    xB = 18.0u"Th"
    y_raw = rawinvmap(rm, xB)
    @test y_raw isa Real && isfinite(y_raw)

    # Consistency with invmap after stripping units
    y_unitful = invmap(rm, xB)  # has u"minute"
    @test y_raw ≈ ustrip(y_unitful)

    # Explicit output-unit conversion (numeric in seconds)
    y_raw_s = rawinvmap(rm, xB; unit=u"s")
    @test y_raw_s ≈ ustrip(u"s", y_unitful)

    # Validation: unitless input rejected for unitful mapper
    @test_throws ArgumentError rawinvmap(rm, 18.0)

    # Unitless mapper variant rejects AbstractQuantity (per validation logic)
    rm0 = JuChrom.RetentionMapper(
        rA, nothing, rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing, rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline, Dict{String, Any}()
    )
    @test_throws ArgumentError rawinvmap(rm0, 18.0u"Th")
end

# ── rawinvmapmax(rm; unit=…) ─────────────────────────────────────────────────

@testset "rawinvmapmax(rm; unit=…)" begin

    # Minimal RetentionMapper (A unitful minutes → B unitful Th)
    rA = [1.0, 2.0, 3.0]
    rB = [100.0, 200.0, 300.0]

    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String,Any}(),
    )

    # Basic success path: returns a raw (numeric) value equal to B-max in its native unit
    @test rawinvmapmax(rm) ≈ rB_max atol=1e-12

    # Unit conversion (numeric in requested unit)
    @test rawinvmapmax(rm; unit=u"kTh") ≈ (rB_max / 1000) atol=1e-12

    # Unitless variant
    rm0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String,Any}(),
    )

    # Unitless mapper returns raw numeric; asking for a unit should error
    @test rawinvmapmax(rm0) ≈ rB_max atol=1e-12
    @test_throws ArgumentError rawinvmapmax(rm0; unit=u"Th")
end

# ── rawinvmapmin(rm; unit=…) ─────────────────────────────────────────────────

@testset "rawinvmapmin(rm; unit=…)" begin

    # Minimal RetentionMapper (A unitful minutes → B unitful Th)
    rA = [1.0, 2.0, 3.0]
    rB = [100.0, 200.0, 300.0]

    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # simple increasing spline
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String,Any}(),
    )

    # Basic success path: returns a raw (numeric) value equal to B-min in its native unit
    @test rawinvmapmin(rm) ≈ rB_min atol=1e-12

    # Unit conversion (numeric in requested unit)
    @test rawinvmapmin(rm; unit=u"kTh") ≈ (rB_min / 1000) atol=1e-12

    # Unitless variant
    rm0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String,Any}(),
    )

    # Unitless mapper returns raw numeric; asking for a unit should error
    @test rawinvmapmin(rm0) ≈ rB_min atol=1e-12
    @test_throws ArgumentError rawinvmapmin(rm0; unit=u"Th")
end

# ── rawmapmax ────────────────────────────────────────────────────────────────

@testset "rawmapmax(rm::AbstractRetentionMapper)" begin
    rA_unitful = [1.2, 2.5, 4.1]u"minute"
    rB_unitful = [10.0, 20.0, 30.0]u"minute"
    mapper_unitful = fitmap(rA_unitful, rB_unitful)

    rA_unitless = [1.2, 2.5, 4.1]
    rB_unitless = [10.0, 20.0, 30.0]
    mapper_unitless = fitmap(rA_unitless, rB_unitless)

    @test rawmapmax(mapper_unitful) ≈ 4.1
    @test rawmapmax(mapper_unitful; unit=u"s") ≈ 246.0

    @test rawmapmax(mapper_unitless) ≈ 4.1
    @test_throws ArgumentError rawmapmax(mapper_unitless; unit=u"s")
end

# ── rawmapmin ────────────────────────────────────────────────────────────────

@testset "rawmapmin(rm::AbstractRetentionMapper)" begin
    rA_unitful = [1.2, 2.5, 4.1]u"minute"
    rB_unitful = [10.0, 20.0, 30.0]u"minute"
    mapper_unitful = fitmap(rA_unitful, rB_unitful)

    rA_unitless = [1.2, 2.5, 4.1]
    rB_unitless = [10.0, 20.0, 30.0]
    mapper_unitless = fitmap(rA_unitless, rB_unitless)

    @test rawmapmin(mapper_unitful) ≈ 1.2
    @test rawmapmin(mapper_unitful; unit=u"s") ≈ 72.0

    @test rawmapmin(mapper_unitless) ≈ 1.2
    @test_throws ArgumentError rawmapmin(mapper_unitless; unit=u"s")
end

# ── rawretentions_A ──────────────────────────────────────────────────────────

@testset "rawretentions_A(rm::AbstractRetentionMapper)" begin
    rA_unitful = [1.2, 2.5, 4.1]u"minute"
    rB_unitful = [10.0, 20.0, 30.0]u"minute"
    mapper_unitful = fitmap(rA_unitful, rB_unitful)

    rA_unitless = [1.2, 2.5, 4.1]
    rB_unitless = [10.0, 20.0, 30.0]
    mapper_unitless = fitmap(rA_unitless, rB_unitless)

    @test rawretentions_A(mapper_unitful) ≈ [1.2, 2.5, 4.1]
    @test rawretentions_A(mapper_unitful; unit=u"s") ≈ [72.0, 150.0, 246.0]

    @test rawretentions_A(mapper_unitless) ≈ [1.2, 2.5, 4.1]
    @test_throws ArgumentError rawretentions_A(mapper_unitless; unit=u"s")
end

# ── rawretentions_B ──────────────────────────────────────────────────────────

@testset "rawretentions_B(rm::AbstractRetentionMapper)" begin
    rA_unitful = [1.2, 2.5, 4.1]u"minute"
    rB_unitful = [10.0, 20.0, 30.0]u"minute"
    mapper_unitful = fitmap(rA_unitful, rB_unitful)

    rA_unitless = [1.2, 2.5, 4.1]
    rB_unitless = [10.0, 20.0, 30.0]
    mapper_unitless = fitmap(rA_unitless, rB_unitless)

    @test rawretentions_B(mapper_unitful) ≈ [10.0, 20.0, 30.0]
    @test rawretentions_B(mapper_unitful; unit=u"s") ≈ [600.0, 1_200.0, 1_800.0]

    @test rawretentions_B(mapper_unitless) ≈ [10.0, 20.0, 30.0]
    @test_throws ArgumentError rawretentions_B(mapper_unitless; unit=u"s")
end

# ── retentions_A(rm; unit=…) ─────────────────────────────────────────────────

@testset "retentions_B(rm; unit=…)" begin

    # Hand-built RetentionMapper with unitful B (Th)
    rA = [1.0, 2.0, 3.0]
    rB = [100.0, 200.0, 300.0]

    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # simple increasing spline
    spline = Spline(B, coefs)

    rmB = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Returns original B-values with their unit by default
    @test retentions_B(rmB) ≈ (rB .* u"Th")

    # Unit conversion works: request same unit (no-op) and a scaled unit
    @test retentions_B(rmB; unit=u"Th")  ≈ (rB .* u"Th")
    @test retentions_B(rmB; unit=u"kTh") ≈ ((rB ./ 1000) .* u"kTh")

    # Hand-built RetentionMapper with unitless B
    rmB0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Returns raw numbers for unitless mapper
    @test retentions_B(rmB0) == rB

    # Asking for a unit on a unitless mapper errors
    @test_throws ArgumentError retentions_B(rmB0; unit=u"Th")
end

# ── retentions_B(rm; unit=…) ─────────────────────────────────────────────────

@testset "retentions_B(rm; unit=…)" begin

    # Hand-built RetentionMapper with unitful B (Th)
    rA = [1.0, 2.0, 3.0]
    rB = [100.0, 200.0, 300.0]

    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))  # simple increasing spline
    spline = Spline(B, coefs)

    rmB = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Returns original B-values with their unit by default
    @test retentions_B(rmB) ≈ (rB .* u"Th")

    # Unit conversion works: request same unit (no-op) and a scaled unit
    @test retentions_B(rmB; unit=u"Th")  ≈ (rB .* u"Th")
    @test retentions_B(rmB; unit=u"kTh") ≈ ((rB ./ 1000) .* u"kTh")

    # Hand-built RetentionMapper with unitless B
    rmB0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    # Returns raw numbers for unitless mapper
    @test retentions_B(rmB0) == rB

    # Asking for a unit on a unitless mapper errors
    @test_throws ArgumentError retentions_B(rmB0; unit=u"Th")
end

# ── retentionunit_A(rm) ──────────────────────────────────────────────────────

@testset "retentionunit_A(rm)" begin

    # Hand-built mapper: A unitful (minute), B unitless
    rA = [1.0, 2.0, 3.0]
    rB = [100.0, 200.0, 300.0]
    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0
    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))
    spline = Spline(B, coefs)

    rmA = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline, Dict{String,Any}(),
    )

    # Unitful mapper → returns the A-unit
    @test retentionunit_A(rmA) == u"minute"

    # Unitless mapper → returns nothing
    rmA0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline, Dict{String,Any}(),
    )
    @test retentionunit_A(rmA0) === nothing
end

# ── retentionunit_B(rm) ──────────────────────────────────────────────────────

@testset "retentionunit_B(rm)" begin

    # Hand-built mapper: B unitful (Th), A unitless
    rA = [1.0, 2.0, 3.0]
    rB = [100.0, 200.0, 300.0]
    rA_min, rA_max = minimum(rA), maximum(rA)
    rB_min, rB_max = minimum(rB), maximum(rB)
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0
    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = collect(range(0.0, 1.0, length(B)))
    spline = Spline(B, coefs)

    rmB = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline, Dict{String,Any}(),
    )

    # Unitful mapper → returns the B-unit
    @test retentionunit_B(rmB) == u"Th"

    # Unitless mapper → returns nothing
    rmB0 = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max, rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline, Dict{String,Any}(),
    )
    @test retentionunit_B(rmB0) === nothing
end

# ── Base.show(::MIME"text/plain", ::RetentionMapper) ─────────────────────────

@testset "Base.show(::MIME\"text/plain\", ::RetentionMapper)" begin
    # Minimal valid RetentionMapper (fast; no optimization)
    rA = [0.0, 0.5, 1.0]
    rB = [10.0, 20.0, 30.0]
    rA_min, rA_max = 0.0, 1.0
    rB_min, rB_max = 10.0, 30.0
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = zeros(length(B))
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max,
        rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max,
        rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}()
    )

    s_plain = sprint(io -> show(io, MIME"text/plain"(), rm))
    s_comp  = sprint(io -> show(io, rm))
    @test s_plain == s_comp  # must delegate to the non-MIME show
end

# ── Base.show(::MIME"text/plain", ::Vector{<:RetentionMapper}) ───────────────

@testset "Base.show(::MIME\"text/plain\", ::Vector{<:RetentionMapper})" begin

    # Helper to make a minimal, valid RetentionMapper quickly
    make_rm(rB_shift) = begin
        rA = [0.0, 0.5, 1.0]
        rB = [10.0, 20.0, 30.0] .+ rB_shift
        rA_min, rA_max = 0.0, 1.0
        rB_min, rB_max = minimum(rB), maximum(rB)
        rA_norm_min, rA_norm_max = 0.0, 1.0
        rB_pred_min, rB_pred_max = rB_min, rB_max
        rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

        order = BSplineOrder(4)
        knots = collect(LinRange(0.0, 1.0, 8))
        B = BSplineBasis(order, knots)
        coefs = zeros(length(B))  # simple constant spline
        spline = Spline(B, coefs)

        JuChrom.RetentionMapper(
            rA, nothing,
            rA_min, rA_max, rA_norm_min, rA_norm_max,
            rB, nothing,
            rB_min, rB_max, rB_pred_min, rB_pred_max,
            rB_pred_norm_min, rB_pred_norm_max,
            knots, coefs, spline,
            Dict{String, Any}("tag" => "rm$(rB_shift)")
        )
    end

    make_rm_units(rB_shift) = begin
        rA = [0.0, 0.5, 1.0]
        rB = [10.0, 20.0, 30.0] .+ rB_shift
        rA_min, rA_max = 0.0, 1.0
        rB_min, rB_max = minimum(rB), maximum(rB)
        rA_norm_min, rA_norm_max = 0.0, 1.0
        rB_pred_min, rB_pred_max = rB_min, rB_max
        rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

        order = BSplineOrder(4)
        knots = collect(LinRange(0.0, 1.0, 8))
        B = BSplineBasis(order, knots)
        coefs = zeros(length(B))
        spline = Spline(B, coefs)

        JuChrom.RetentionMapper(
            rA, u"minute",
            rA_min, rA_max, rA_norm_min, rA_norm_max,
            rB, u"Th",
            rB_min, rB_max, rB_pred_min, rB_pred_max,
            rB_pred_norm_min, rB_pred_norm_max,
            knots, coefs, spline,
            Dict{String, Any}("tag" => "rm$(rB_shift)")
        )
    end

    mappers = [make_rm(0.0), make_rm(5.0)]

    s = sprint(io -> show(io, MIME"text/plain"(), mappers))

    # Header shows count and pluralization
    @test startswith(s, "2 RetentionMapper")
    # First entry uses tree prefix and includes key facts
    @test occursin("├─ [1] ", s)
    @test occursin("pts, ", s)
    @test occursin("unitless→unitless", s)  # units are Nothing here
    @test occursin("avg residual:", s)
    # Last entry uses terminal tree prefix
    @test occursin("└─ [2] ", s)

    mappers_u = [make_rm_units(0.0), make_rm_units(5.0)]
    su = sprint(io -> show(io, MIME"text/plain"(), mappers_u))
    @test occursin("minute→Th", su)
end

# ── tune_lambda_for_monotonic_spline ────────────────────────────────────────

@testset "tune_lambda_for_monotonic_spline" begin
    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    rA_norm_min, rA_norm_max = 0.0, 1.0

    function fit_spline_with_penalty(::Real)
        coefs = zeros(length(B))
        return coefs, JuChrom.MOI.OPTIMAL
    end

    @test_throws ArgumentError tune_lambda_for_monotonic_spline(
        fit_spline_with_penalty,
        B,
        rA_norm_min,
        rA_norm_max,
        1e-3,
        1e-1,
        1e-6,
        5,
        10,
    )

    order2 = BSplineOrder(2)
    knots2 = collect(0.0:0.25:1.0)
    B2 = BSplineBasis(order2, knots2)

    function fit_spline_with_penalty_monotonic(::Real)
        coefs = collect(0.0:(length(B2)-1))
        return coefs, JuChrom.MOI.OPTIMAL
    end

    @test_logs (:warn, r"Reached maximum iterations") begin
        tune_lambda_for_monotonic_spline(
            fit_spline_with_penalty_monotonic,
            B2,
            rA_norm_min,
            rA_norm_max,
            1e-6,
            1.0,
            1e-12,
            1,
            10,
        )
    end
end

# ── Base.show(::OptimizationError) ───────────────────────────────────────────

@testset "Base.show(::OptimizationError)" begin
    err = JuChrom.OptimizationError(2, :INFEASIBLE)
    s = sprint(io -> show(io, err))
    @test occursin("OptimizationError:", s)
    @test occursin("λ = 2", s)
    @test occursin("INFEASIBLE", s)
end

# ── Base.showerror(::OptimizationError) ──────────────────────────────────────

@testset "Base.showerror(::OptimizationError)" begin
    err = JuChrom.OptimizationError(5, :NUMERICAL_ERROR)
    s = sprint(io -> showerror(io, err))
    @test occursin("λ = 5", s)
    @test occursin("NUMERICAL_ERROR", s)
    @test !occursin("OptimizationError:", s)  # concise message only
end

# ── Base.show(::RetentionMapper) ─────────────────────────────────────────────

@testset "Base.show(::RetentionMapper)" begin
    # Minimal, valid RetentionMapper with tiny cubic spline and some metadata
    rA = [0.0, 0.5, 1.0]
    rB = [10.0, 20.0, 30.0]
    rA_min, rA_max = 0.0, 1.0
    rB_min, rB_max = 10.0, 30.0
    rA_norm_min, rA_norm_max = 0.0, 1.0
    rB_pred_min, rB_pred_max = rB_min, rB_max
    rB_pred_norm_min, rB_pred_norm_max = 0.0, 1.0

    order = BSplineOrder(4)
    knots = collect(LinRange(0.0, 1.0, 8))
    B = BSplineBasis(order, knots)
    coefs = zeros(length(B))  # simple, well-defined spline
    spline = Spline(B, coefs)

    rm = JuChrom.RetentionMapper(
        rA, nothing,  # rA, rA_unit
        rA_min, rA_max,
        rA_norm_min, rA_norm_max,
        rB, nothing,  # rB, rB_unit
        rB_min, rB_max,
        rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}("source" => "unit test", "note" => "demo"),
    )

    s = sprint(io -> show(io, rm))

    # Header and domains (assert only stable, documented substrings)
    @test occursin("RetentionMapper with 3 calibration points", s)
    @test occursin("Domain A (input):", s)
    @test occursin("Domain B (output):", s)

    # Spline info block
    @test occursin("Spline:", s)
    @test occursin("Order: 4", s)  # avoid overfitting exact wording
    @test occursin("Knots:", s)
    @test occursin("Coefficients:", s)

    # Fit quality block: residual summaries appear
    @test occursin("Fit quality:", s)
    @test occursin("Min residual:", s)
    @test occursin("Avg residual:", s)
    @test occursin("Max residual:", s)

    # Metadata block prints when extras is nonempty
    @test occursin("Metadata:", s)
    @test occursin("source", s)
    @test occursin("note", s)

    # Unitful rA should include units in residual calculation path
    rm_u = JuChrom.RetentionMapper(
        rA, u"minute",
        rA_min, rA_max,
        rA_norm_min, rA_norm_max,
        rB, u"Th",
        rB_min, rB_max,
        rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}(),
    )
    su = sprint(io -> show(io, rm_u))
    @test occursin("unitless→unitless", su) == false
    @test occursin("minute", su)

    rm_complex = JuChrom.RetentionMapper(
        rA, nothing,
        rA_min, rA_max,
        rA_norm_min, rA_norm_max,
        rB, nothing,
        rB_min, rB_max,
        rB_pred_min, rB_pred_max,
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        Dict{String, Any}(
            "nested" => Dict("a" => 1, "b" => 2),
            "long" => repeat("x", 80),
        ),
    )
    scomplex = sprint(io -> show(io, rm_complex))
    @test occursin("nested", scomplex)
    @test occursin("a =", scomplex)
    @test occursin("long =", scomplex)
    @test occursin("\n", scomplex)
end

end  # module
