module TestUnits

using Test
using Unitful
using JuChrom

# ─────────────────────────────────────────────────────────────────────────────
# consistentunits
# ─────────────────────────────────────────────────────────────────────────────
@testset "consistentunits" begin
    @test JuChrom.consistentunits(Int[]) === true
    @test JuChrom.consistentunits([1.0, 2.0, 3.0]) === true
    @test JuChrom.consistentunits([1, 2, 3] * u"s") === true
    @test JuChrom.consistentunits([1 * u"s", 2 * u"minute", 3 * u"s"]) === true
    @test JuChrom.consistentunits([1 * u"m", 2.0 * u"minute", 3.0 * u"kg"]) === false
    @test JuChrom.consistentunits([1 * u"s", 2.0, 3.0]) === false
end

# ─────────────────────────────────────────────────────────────────────────────
# strip_units_checked
# ─────────────────────────────────────────────────────────────────────────────
@testset "strip_units_checked" begin
    vals, u = JuChrom.strip_units_checked([1, 2]u"s", "times")
    @test vals == [1, 2]
    @test u == u"s"

    vals2, u2 = JuChrom.strip_units_checked([1.0, 2.0], "values")
    @test vals2 == [1.0, 2.0]
    @test u2 === nothing

    empty_vals, u3 = JuChrom.strip_units_checked(Int[], "empty")
    @test empty_vals == Int[]    # equality, not identity
    @test u3 === nothing

    @test_throws ArgumentError JuChrom.strip_units_checked([1u"s", 2u"m"], "mixed")
end

# ─────────────────────────────────────────────────────────────────────────────
# inverse
# ─────────────────────────────────────────────────────────────────────────────
@testset "inverse(unit)" begin
    @test JuChrom.inverse(u"m") == u"m^-1"
    @test JuChrom.inverse(u"kg") == u"kg^-1"
    @test JuChrom.inverse(u"s") == u"s^-1"
end

# ─────────────────────────────────────────────────────────────────────────────
# Core unit conversion helpers (_handle_*): vectors & scalars
# ─────────────────────────────────────────────────────────────────────────────
@testset "_handle_unitless" begin
    v = [1.0, 2.0]
    @test JuChrom._handle_unitless(v, nothing, "retentions") === v
    @test_throws ArgumentError JuChrom._handle_unitless(v, u"s", "retentions")
end

@testset "_handle_unitful_convert (vector)" begin
    vals = [1.0, 2.0]
    @test JuChrom._handle_unitful_convert(vals, u"s", nothing) == [1.0, 2.0]u"s"
    @test JuChrom._handle_unitful_convert(vals, u"s", u"ms") == [1000.0, 2000.0]u"ms"
    @test_throws Unitful.DimensionError JuChrom._handle_unitful_convert(vals, u"s", u"m")
end

@testset "_handle_unitful_strip (vector)" begin
    vals = [1.0, 2.0]
    @test JuChrom._handle_unitful_strip(vals, u"s", nothing) == [1.0, 2.0]
    @test JuChrom._handle_unitful_strip(vals, u"s", u"ms") == [1000.0, 2000.0]
    @test_throws Unitful.DimensionError JuChrom._handle_unitful_strip(vals, u"s", u"m")
end

@testset "_handle_unitful_convert_scalar" begin
    @test JuChrom._handle_unitful_convert_scalar(2.0, u"s", nothing) == 2.0u"s"
    @test JuChrom._handle_unitful_convert_scalar(2.0, u"s", u"ms") == 2000.0u"ms"
    @test_throws Unitful.DimensionError JuChrom._handle_unitful_convert_scalar(1.0, u"s", u"m")
end

@testset "_handle_unitful_strip_scalar" begin
    @test JuChrom._handle_unitful_strip_scalar(2.0, u"s", nothing) == 2.0
    @test JuChrom._handle_unitful_strip_scalar(2.0, u"s", u"ms") == 2000.0
    @test_throws Unitful.DimensionError JuChrom._handle_unitful_strip_scalar(1.0, u"s", u"m")
end

end # module TestUnits
