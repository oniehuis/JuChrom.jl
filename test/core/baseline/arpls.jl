module ArplsTests

using Test
using LinearAlgebra
using Logging
using Random
using SparseArrays
using Statistics
using JuChrom
using JuChrom: arpls

function arplsreference(
    y; λ=1e5, ratio=1e-3, maxiter=10_000, variances=nothing,
    variancefloor=1e-12, peakthreshold=4.0, peakslope=1.0, zerothreshold=1e-8,
    zerofractionthreshold=0.2, zeroweight=0.01
    )

    n = length(y)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for i in 1:(n - 2)
        push!(rows, i); push!(cols, i); push!(vals, 1.0)
        push!(rows, i); push!(cols, i + 1); push!(vals, -2.0)
        push!(rows, i); push!(cols, i + 2); push!(vals, 1.0)
    end
    D = sparse(rows, cols, vals, n - 2, n)
    H = λ * (D' * D)
    w = ones(Float64, n)
    z = zeros(Float64, n)
    precisionweights = if isnothing(variances)
        ones(Float64, n)
    else
        variancefloor_eff = max(Float64(variancefloor), eps(Float64))
        p = @. 1 / max(Float64(variances), variancefloor_eff)
        p ./= mean(p)
        p
    end
    residualscales = @. 1 / sqrt(precisionweights)
    zeromask = y .< zerothreshold
    if count(zeromask) / n > zerofractionthreshold
        fill!(zeromask, false)
    end
    zerofactors = ones(Float64, n)
    zerofactors[zeromask] .= zeroweight

    function referenceweights(dw)
        dn = dw[(.!zeromask) .& (dw .< 0)]
        length(dn) < 2 && return nothing

        m = mean(dn)
        s = std(dn)
        if !(isfinite(s) && s > 0)
            return nothing
        end

        wt = @. 1 / (1 + exp(peakslope * (dw - (peakthreshold * s - m)) / s))
        wt
    end

    for _ in 1:maxiter
        fitweights = w .* precisionweights .* zerofactors
        W = Diagonal(fitweights)
        z .= cholesky(W + H) \ (fitweights .* y)
        d = y - z
        dw = isnothing(variances) ? d : d ./ residualscales
        wt = referenceweights(dw)
        isnothing(wt) && break
        norm(w - wt) / norm(w) < ratio && break
        w .= wt
    end
    z
end

@testset "arpls(::AbstractVector; ...)" begin
    Random.seed!(121)
    n = 180
    x = collect(range(0.0, 1.0; length=n))
    btrue = @. 0.2 + 0.15 * x + 0.05 * x^2
    y = copy(btrue)
    @. y += 1.4 * exp(-((x - 0.35)^2) / (2 * 0.018^2))
    @. y += 0.8 * exp(-((x - 0.72)^2) / (2 * 0.030^2))
    y .+= 0.005 .* randn(n)

    b = arpls(y; λ=1e5, ratio=1e-6, maxiter=500)
    bref = arplsreference(y; λ=1e5, ratio=1e-6, maxiter=500)
    bdefault = arpls(y; λ=1e5, maxiter=500)
    brefdefault = arplsreference(y; λ=1e5, ratio=1e-3, maxiter=500)
    boriginal = arpls(y; λ=1e5, ratio=1e-6, maxiter=500,
                      peakthreshold=2.0, peakslope=2.0)
    breforiginal = arplsreference(y; λ=1e5, ratio=1e-6, maxiter=500,
                                  peakthreshold=2.0, peakslope=2.0)

    @test length(b) == n
    @test all(isfinite, b)
    @test b ≈ bref
    @test bdefault ≈ brefdefault
    @test arpls(y; λ=1e5, ratio=1e-6, maxiter=500,
                peakthreshold=4.0, peakslope=1.0) ≈ b
    @test boriginal ≈ breforiginal
    @test norm(boriginal - b) > 1e-3
    @test norm(arpls(y; λ=1e5, ratio=1e-6, maxiter=500,
                     peakslope=0.5) - b) > 1e-3

    @test all(≥(0), arpls(-ones(12); nonnegative=true))
    @test any(<(0), arpls(-ones(12); nonnegative=false))

    variances = @. 0.2 + 0.8 * (1 + sin(2π * x))^2
    bvar = arpls(y; λ=1e5, ratio=1e-6, maxiter=500, variances=variances)
    bvarref = arplsreference(y; λ=1e5, ratio=1e-6, maxiter=500,
                             variances=variances)
    bvar_scaled = arpls(y; λ=1e5, ratio=1e-6, maxiter=500,
                        variances=1e6 .* variances)
    smallvariances = @. 1e-6 * (0.2 + (1 + sin(2π * x))^2)
    bsmall_z0 = arpls(y; λ=1e5, ratio=1e-6, maxiter=500,
                      variances=smallvariances, zerothreshold=0.0)
    bsmall_zpos = arpls(y; λ=1e5, ratio=1e-6, maxiter=500,
                        variances=smallvariances, zerothreshold=0.05)
    bsmall_floor_high = arpls(y; λ=1e5, ratio=1e-6, maxiter=500,
                              variances=smallvariances, variancefloor=1e-3)
    @test bvar ≈ bvarref
    @test bvar_scaled ≈ bvar rtol=1e-8 atol=1e-8
    @test bsmall_zpos ≈ bsmall_z0
    @test norm(bsmall_floor_high - bsmall_z0) > 1e-6
    @test arpls(y; λ=1e5, ratio=1e-6, maxiter=500,
                variances=ones(n)) ≈ b
    @test arpls(y; λ=1e5, ratio=1e-6, maxiter=500,
                variances=zeros(n)) ≈ b

    ydrop = copy(y)
    ydrop[104:112] .= 0.0
    bdrop = arpls(ydrop; λ=1e5, ratio=1e-6, maxiter=500, zeroweight=1.0)
    bdrop_zero = arpls(ydrop; λ=1e5, ratio=1e-6, maxiter=500)
    @test mean(bdrop_zero[104:112]) > mean(bdrop[104:112])

    ysparse = copy(y)
    ysparse[1:45] .= 0.0
    bsparse = arpls(ysparse; λ=1e5, ratio=1e-6, maxiter=500)
    bsparse_ref = arplsreference(ysparse; λ=1e5, ratio=1e-6, maxiter=500)
    bsparse_forcedropout = arpls(
        ysparse; λ=1e5, ratio=1e-6, maxiter=500, zerofractionthreshold=1.0)
    @test bsparse ≈ bsparse_ref
    @test norm(bsparse - bsparse_forcedropout) > 1e-3

    @test_throws ArgumentError arpls(y[1:2])
    @test_throws ArgumentError arpls([1.0, Inf, 2.0])
    @test_throws ArgumentError arpls(y; λ=0.0)
    @test_throws ArgumentError arpls(y; ratio=0.0)
    @test_throws ArgumentError arpls(y; maxiter=0)
    @test_throws ArgumentError arpls(fill(Inf, n))
    @test_throws MethodError arpls(x, y)
    @test_throws ArgumentError arpls(y; variances=ones(n - 1))
    @test_throws ArgumentError arpls(y; variances=fill(-1.0, n))
    @test_throws ArgumentError arpls(y; variances=ones(n, 1))
    @test_throws ArgumentError arpls(y; peakthreshold=0.0)
    @test_throws ArgumentError arpls(y; peakslope=0.0)
    @test_throws ArgumentError arpls(y; zerothreshold=-1e-9)
    @test_throws ArgumentError arpls(y; zerofractionthreshold=-1e-9)
    @test_throws ArgumentError arpls(y; zerofractionthreshold=1.01)
    @test_throws ArgumentError arpls(y; zeroweight=0.0)
    @test_throws ArgumentError arpls(y; variancefloor=0.0)
    @test_throws ArgumentError arpls(y; variancefloor=-1e-12)
end

@testset "arpls(::MassScanMatrix; ...)" begin
    Random.seed!(122)
    n_scans = 70
    n_mzs = 3
    steps = 0.8 .+ 0.4 .* rand(n_scans - 1)
    ret = [0.0; cumsum(steps)]
    ret .*= 4.0 / last(ret)
    mzs = [100.0, 150.0, 200.0]
    ints = zeros(n_scans, n_mzs)
    base = @. 0.4 + 0.03 * ret
    for j in 1:n_mzs
        @. ints[:, j] = base + (0.6 + 0.2j) * exp(-((ret - (0.7 + j))^2) /
                                                  (2 * (0.09 + 0.01j)^2))
    end
    ints .+= 0.002 .* randn(n_scans, n_mzs)

    msm = MassScanMatrix(ret, mzs, ints; sample=(name="arpls",), extras=Dict("k" => 2))
    bmsm = arpls(msm; λ=5e4, ratio=1e-6, maxiter=500)
    expected = reduce(hcat, [arpls(@view(ints[:, j]); λ=5e4, ratio=1e-6,
                                   maxiter=500)
                             for j in axes(ints, 2)])
    variances_ones = ones(size(ints))
    bmsm_var_ones = arpls(msm; λ=5e4, ratio=1e-6, maxiter=500,
                          variances=variances_ones)
    scaledret = ret ./ last(ret)
    variances = reduce(hcat, [@. 0.3 + 0.2j + scaledret^2
                              for j in axes(ints, 2)])
    bmsm_var = arpls(msm; λ=5e4, ratio=1e-6, maxiter=500, variances=variances)
    expected_var = reduce(hcat, [arpls(@view(ints[:, j]); λ=5e4, ratio=1e-6,
                                       maxiter=500, variances=@view(variances[:, j]))
                                 for j in axes(ints, 2)])
    bmsm_tuned = arpls(msm; λ=5e4, ratio=1e-6, maxiter=500,
                       peakthreshold=2.5, peakslope=1.5)
    expected_tuned = reduce(hcat, [arpls(@view(ints[:, j]); λ=5e4,
                                         ratio=1e-6, maxiter=500,
                                         peakthreshold=2.5, peakslope=1.5)
                                   for j in axes(ints, 2)])
    logger = TestLogger()
    with_logger(logger) do
        arpls(msm; λ=5e4, ratio=1e-12, maxiter=1)
    end
    warnmessages = [
        string(log.message) for log in logger.logs if log.level == Logging.Warn
    ]
    @test length(warnmessages) == n_mzs
    for (i, mz) in enumerate(mzs)
        @test any(
            occursin("m/z channel $i (m/z=$mz)", message) for message in warnmessages)
    end

    @test bmsm isa MassScanMatrix
    @test JuChrom.rawretentions(bmsm) == ret
    @test JuChrom.rawmzvalues(bmsm) == mzs
    @test JuChrom.rawintensities(bmsm) ≈ expected
    @test JuChrom.rawintensities(bmsm_var_ones) ≈ expected
    @test JuChrom.rawintensities(bmsm_var) ≈ expected_var
    @test JuChrom.rawintensities(bmsm_tuned) ≈ expected_tuned
    @test bmsm.sample == msm.sample
    @test bmsm.extras == msm.extras
    @test JuChrom.rawretentions(bmsm) !== JuChrom.rawretentions(msm)
    @test JuChrom.rawmzvalues(bmsm) !== JuChrom.rawmzvalues(msm)

    msmneg = MassScanMatrix(ret[1:6], mzs[1:2], -ones(6, 2))
    bmsmneg = arpls(msmneg; nonnegative=true, maxiter=20)
    @test all(≥(0), JuChrom.rawintensities(bmsmneg))
    @test_throws ArgumentError arpls(msm; variances=ones(n_scans - 1, n_mzs))
    @test_throws ArgumentError arpls(msm; variances=ones(n_scans))
end

end # module
