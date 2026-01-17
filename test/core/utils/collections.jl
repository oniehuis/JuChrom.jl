module TestCollections

using Test
using Unitful
using JuChrom

# ─────────────────────────────────────────────────────────────────────────────
# copy_with
# ─────────────────────────────────────────────────────────────────────────────
@testset "copy_with" begin
    struct MyStruct
        x::Int
        y::String
        z::Vector{Int}
    end
    obj = MyStruct(1, "a", [1,2])
    obj2 = JuChrom.copy_with(obj, (y = "b",))
    @test obj2 isa MyStruct
    @test obj2.x == 1
    @test obj2.y == "b"
    @test obj2.z == [1,2]

    # deep copy behavior on untouched fields
    obj2.z[1] = 99
    @test obj.z == [1,2]  # original not mutated

    # multiple overrides
    obj3 = JuChrom.copy_with(obj, (x=42, y="hello"))
    @test obj3.x == 42
    @test obj3.y == "hello"
    @test obj3.z == [1,2]
end

# ─────────────────────────────────────────────────────────────────────────────
# typify (dict & vector)
# ─────────────────────────────────────────────────────────────────────────────
@testset "typify" begin
    # Dict with Int values
    d = Dict("a" => 1, "b" => 2)
    td = JuChrom.typify(d)
    @test td isa Dict{String,Int}
    @test td == d

    # Dict with Float64 values
    df = Dict("a" => 1.0, "b" => 2.0)
    tdf = JuChrom.typify(df)
    @test tdf isa Dict{String,Float64}
    @test tdf == df

    # Dict with mixed numeric values → promotes to Float64
    dmix = Dict("a" => 1, "b" => 2.0)
    tdmix = JuChrom.typify(dmix)
    @test tdmix isa Dict{String,Float64}
    @test tdmix == Dict("a" => 1.0, "b" => 2.0)

    # Dict with incompatible types → falls back to Any, logs a warning
    dm = Dict("a" => 1, "b" => "x")
    @test_logs (:warn, r"Type Any promoted by typify is not concrete") match_mode = :any begin
        tdm = JuChrom.typify(dm)
        @test tdm isa Dict{String,Any}
        @test tdm == Dict("a" => 1, "b" => "x")
    end

    # Empty Dict → Union{} promotes to Any, logs a warning; result is Dict{Any,Any}()
    de = Dict{String,Int}()
    @test_logs (:warn, r"Type Union\{\} promoted by typify is not concrete") match_mode = :any begin
        tde = JuChrom.typify(de)
        @test tde isa Dict{Any,Any}
        @test isempty(tde)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# findclosest
# ─────────────────────────────────────────────────────────────────────────────
@testset "findclosest" begin
    A = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    @test JuChrom.findclosest(A, 0.0) == 3
    @test JuChrom.findclosest(A, 1.4) == 4
    @test JuChrom.findclosest(A, 1.5) == 5  # tie -> larger
    @test JuChrom.findclosest(A, -1.5) == 2 # tie -> larger
    @test JuChrom.findclosest(A, -10.0) == 1
    @test JuChrom.findclosest(A, 10.0) == 6

    # Single element vector
    @test JuChrom.findclosest([42.0], 0.0) == 1

    # Works with integer inputs
    @test JuChrom.findclosest([1, 3, 5, 7], 6) == 4  # tie -> larger
end

end  # module TestCollections
