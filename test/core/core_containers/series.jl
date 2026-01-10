module TestScanSeriesContainers

using Test
using Unitful
using Unitful: s, minute, pA
using JuChrom  # ChromScan, MassScan, ChromScanSeries, MassScanSeries

# Helpers to build scans quickly
mkchrom(i, rt_unit, I_unit) = ChromScan(i*1.0*rt_unit, (10i)*I_unit)
mkms(rt, mzs, ints; level=1) = MassScan(rt, mzs, ints; level)

# ─────────────────────────────────────────────────────────────────────────────
# ChromScanSeries
# ─────────────────────────────────────────────────────────────────────────────

@testset "ChromScanSeries – inner constructor" begin
    scans = [ChromScan(1.0u"minute", 10.0u"pA"),
             ChromScan(2.0u"minute", 20.0u"pA")]
    S  = eltype(scans)
    T1 = typeof(scans)
    T2 = typeof((vendor="Acme",))
    T3 = typeof((method="M1",))
    T4 = typeof((user="U",))
    T5 = typeof((sample="S",))
    T6 = Dict{String,Any}

    css = ChromScanSeries{S,
                          typeof(first(scans).retention_unit),
                          typeof(first(scans).intensity_unit),
                          T1, T2, T3, T4, T5, T6}(
        scans, (vendor="Acme",), (method="M1",), (user="U",), (sample="S",),
        Dict{String,Any}("x"=>1))
    @test css.scans === scans
    @test css.instrument == (vendor="Acme",)
    @test css.acquisition == (method="M1",)
    @test css.user == (user="U",)
    @test css.sample == (sample="S",)
    @test css.extras == Dict("x"=>1)

    # Errors: empty scans, abstract eltype
    empty_vec = Vector{S}()  # CONCRETE eltype; ensures we hit the ctor's non-empty check
    @test_throws ArgumentError ChromScanSeries{S, Nothing, Nothing,
                                               typeof(empty_vec),
                                               NamedTuple, NamedTuple, NamedTuple, NamedTuple, Dict{String,Any}}(
        empty_vec, NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    abs_vec = AbstractVector{AbstractChromScan}(scans)  # abstract eltype should be rejected
    @test_throws ArgumentError ChromScanSeries{eltype(abs_vec), Nothing, Nothing,
                                               typeof(abs_vec),
                                               NamedTuple, NamedTuple, NamedTuple, NamedTuple, Dict{String,Any}}(
        abs_vec, NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())
end

@testset "ChromScanSeries – outer constructor" begin
    scans = [mkchrom(1, u"minute", u"pA"), mkchrom(2, u"minute", u"pA")]

    # Defaults & type-parameter inference (R=minute, I=pA)
    css = ChromScanSeries(scans)
    @test css.scans === scans
    @test css.instrument == NamedTuple()
    @test css.acquisition == NamedTuple()
    @test css.user == NamedTuple()
    @test css.sample == NamedTuple()
    @test css.extras isa Dict{String,Any}
    @test typeof(css).parameters[2] === typeof(first(scans).retention_unit)
    @test typeof(css).parameters[3] === typeof(first(scans).intensity_unit)

    # extras: AbstractString keys are converted to String; values widened to Any
    s = "ab"; k = SubString(s, 1, 1)
    css2 = ChromScanSeries(scans; extras=Dict(k=>1))
    @test css2.extras isa Dict{String,Any}
    @test first(keys(css2.extras)) isa String

    # Error paths on outer ctor (non-empty, concrete eltype)
    @test_throws ArgumentError ChromScanSeries(Vector{eltype(scans)}())  # empty concrete
    abs_vec = AbstractVector{AbstractChromScan}(scans)
    @test_throws ArgumentError ChromScanSeries(abs_vec)                  # abstract eltype
end

@testset "ChromScanSeries – interface and equality" begin
    scans = [mkchrom(1, u"minute", u"pA"), mkchrom(2, u"minute", u"pA")]
    css1 = ChromScanSeries(scans; instrument=(vendor="Acme",))
    css2 = ChromScanSeries(scans; instrument=(vendor="Acme",))
    css3 = ChromScanSeries(scans; instrument=(vendor="Other",))

    @test first(css1) == scans[1]
    @test last(css1) == scans[end]
    @test css1 == css2
    @test css1 != css3
end

# ─────────────────────────────────────────────────────────────────────────────
# MassScanSeries
# ─────────────────────────────────────────────────────────────────────────────

@testset "MassScanSeries – inner constructor" begin
    scans = [mkms(1.0u"s", [100.0, 150.0], [10.0, 20.0]),
             mkms(2.0u"s", [100.0, 150.0], [12.0, 22.0]; level=2)]
    S  = eltype(scans)
    T1 = typeof(scans)
    T2 = typeof((instr="MS",))
    T3 = typeof((mode="FullScan",))
    T4 = typeof((user="U",))
    T5 = typeof((sample="S",))
    T6 = Dict{String,Any}

    mss = MassScanSeries{S,
                         typeof(first(scans).retention_unit),
                         typeof(first(scans).mz_unit),
                         typeof(first(scans).intensity_unit),
                         T1, T2, T3, T4, T5, T6}(
        scans, (instr="MS",), (mode="FullScan",), (user="U",), (sample="S",),
        Dict{String,Any}("y"=>2))
    @test mss.scans === scans
    @test mss.instrument == (instr="MS",)
    @test mss.acquisition == (mode="FullScan",)
    @test mss.user == (user="U",)
    @test mss.sample == (sample="S",)
    @test mss.extras == Dict("y"=>2)

    # Errors: empty scans, abstract eltype
    empty_vec = typeof(scans)([])  # already concrete eltype
    @test_throws ArgumentError MassScanSeries{S, Nothing, Nothing, Nothing,
                                              typeof(empty_vec),
                                              NamedTuple, NamedTuple, NamedTuple, NamedTuple, Dict{String,Any}}(
        empty_vec, NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())

    abs_vec = AbstractVector{AbstractMassScan}(scans)
    @test_throws ArgumentError MassScanSeries{eltype(abs_vec), Nothing, Nothing, Nothing,
                                              typeof(abs_vec),
                                              NamedTuple, NamedTuple, NamedTuple, NamedTuple, Dict{String,Any}}(
        abs_vec, NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Dict{String,Any}())
end

@testset "MassScanSeries – outer constructor" begin
    scans = [mkms(1.0u"s", [100.0, 150.0], [10.0, 20.0]),
             mkms(2.0u"s", [101.0, 151.0], [12.0, 22.0]; level=2)]

    # Defaults & type-parameter inference (R, M, I)
    mss = MassScanSeries(scans)
    @test mss.scans === scans
    @test mss.instrument == NamedTuple()
    @test mss.acquisition == NamedTuple()
    @test mss.user == NamedTuple()
    @test mss.sample == NamedTuple()
    @test mss.extras isa Dict{String,Any}
    @test typeof(mss).parameters[2] === typeof(first(scans).retention_unit)
    @test typeof(mss).parameters[3] === typeof(first(scans).mz_unit)
    @test typeof(mss).parameters[4] === typeof(first(scans).intensity_unit)

    # extras: AbstractString keys conversion; Any value type
    s = "xyz"; k = SubString(s, 2, 3)
    mss2 = MassScanSeries(scans; extras=Dict(k=>:ok))
    @test mss2.extras isa Dict{String,Any}
    @test first(keys(mss2.extras)) isa String

    # Error paths on outer ctor (non-empty, concrete eltype)
    @test_throws ArgumentError MassScanSeries(typeof(scans)([]))
    abs_vec = AbstractVector{AbstractMassScan}(scans)
    @test_throws ArgumentError MassScanSeries(abs_vec)
end

@testset "MassScanSeries – interface and equality" begin
    scans = [mkms(1.0u"s", [100.0, 150.0], [10.0, 20.0]),
             mkms(2.0u"s", [101.0, 151.0], [12.0, 22.0]; level=2)]
    mss1 = MassScanSeries(scans; acquisition=(mode="FullScan",))
    mss2 = MassScanSeries(scans; acquisition=(mode="FullScan",))
    mss3 = MassScanSeries(scans; acquisition=(mode="SIM",))

    @test first(mss1) == scans[1]
    @test last(mss1) == scans[end]
    @test mss1 == mss2
    @test mss1 != mss3
end

@testset "ScanSeriesDisplay helpers" begin
    io = IOBuffer()
    short = "short"
    long = "x" ^ 80

    @test JuChrom.ScanSeriesDisplay.format_value_with_wrap(short, 10) == short
    @test occursin("\n", JuChrom.ScanSeriesDisplay.format_value_with_wrap(long, 10))

    JuChrom.ScanSeriesDisplay.print_wrapped_value(io, "pfx", short, "  ")
    JuChrom.ScanSeriesDisplay.print_wrapped_value(io, "pfx", "line1\nline2", "  ")

    JuChrom.ScanSeriesDisplay.print_complex_value(io, [1, 2, 3], "", "", true)
    JuChrom.ScanSeriesDisplay.print_complex_value(io, [1, 2, 3, 4], "", "", true)
    JuChrom.ScanSeriesDisplay.print_complex_value(io, Dict("a" => 1), "", "", true)
    JuChrom.ScanSeriesDisplay.print_complex_value(io, (x=1,), "", "", true)
    JuChrom.ScanSeriesDisplay.print_complex_value(io, [1 2; 3 4], "", "", true)
    JuChrom.ScanSeriesDisplay.print_complex_value(io, [1 2 3; 4 5 6; 7 8 9], "", "", true)

    empty_sections = [("Instrument", NamedTuple())]
    @test JuChrom.ScanSeriesDisplay.print_annotations(IOBuffer(), empty_sections, Dict{String, Any}()) == false

    sections = [("Instrument", (model="X",)),
                ("User", (name="A",))]
    @test JuChrom.ScanSeriesDisplay.print_annotations(IOBuffer(), sections, Dict{String, Any}()) == true

    extras = Dict("note" => "x" ^ 70, "nested" => (a=1, b=2))
    @test JuChrom.ScanSeriesDisplay.print_annotations(IOBuffer(), sections, extras) == true

    single_complex = [("Instrument", (config=(alpha=1, beta=2),))]
    @test JuChrom.ScanSeriesDisplay.print_annotations(IOBuffer(), single_complex, Dict{String, Any}()) == true

    single_long = [("Instrument", (comment="x" ^ 120,))]
    @test JuChrom.ScanSeriesDisplay.print_annotations(IOBuffer(), single_long, Dict{String, Any}()) == true

    multi_mixed = [("Instrument", (model="X", notes="x" ^ 120)),
                   ("User", (name="A",))]
    @test JuChrom.ScanSeriesDisplay.print_annotations(IOBuffer(), multi_mixed, Dict{String, Any}()) == true
end

@testset "ScanSeries show methods" begin
    scans = [mkchrom(1, u"minute", u"pA"), mkchrom(2, u"minute", u"pA")]
    css = ChromScanSeries(scans; instrument=(vendor="Acme",))
    css_out = sprint(show, css)
    @test occursin("ChromScanSeries", css_out)

    mscans = [mkms(1.0u"s", [100.0, 150.0], [10.0, 20.0]),
              mkms(2.0u"s", [101.0, 151.0], [12.0, 22.0]; level=2)]
    mss = MassScanSeries(mscans; acquisition=(mode="FullScan",))
    mss_out = sprint(show, mss)
    @test occursin("MassScanSeries", mss_out)
end

end # module
