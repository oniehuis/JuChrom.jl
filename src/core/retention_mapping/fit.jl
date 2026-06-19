# ── fitmap ────────────────────────────────────────────────────────────────────────────────

const DEFAULT_RETENTION_MAPPING_SOLVER_TIME_LIMIT = 1.0

"""
    fitmap(retentions_A::AbstractVector{<:Union{<:Real, <:AbstractQuantity{<:Real}}},
           retentions_B::AbstractVector{<:Union{<:Real, <:AbstractQuantity{<:Real}}};
           smoothness_penalty_n::Integer=100,
           monotonicity_grid_size::Integer=10^5,
           λ::Real=3e-9,
           solver_time_limit::Union{Nothing, Real}=DEFAULT_RETENTION_MAPPING_SOLVER_TIME_LIMIT,
           extras::AbstractDict{<:AbstractString, <:Any}=Dict{String, Any}(),
           optimizer_factory=HiGHS.Optimizer)

Constructs a strictly increasing, smooth cubic B-spline that maps values from domain A 
(`retentions_A`) to domain B (`retentions_B`), subject to a monotonicity constraint and 
smoothness penalty.

The fitted spline minimizes a regularized least-squares objective: data fidelity
(squared residuals) and a smoothness penalty (second derivative norm), subject to
non-negativity of the first derivative.

The smoothing parameter `λ` controls the tradeoff between exact anchor fit and curvature.
The default `λ=3e-9` is chosen for normalized RT -> RI calibration data and favors a stable
inverse mapping over boundary-of-monotonicity interpolation.

Inputs are monotonic `retentions_A`/`retentions_B` vectors (unitful or unitless), optional
controls for smoothness (`smoothness_penalty_n`, `λ`) and monotonicity validation
(`monotonicity_grid_size`), optional metadata (`extras`), and an
optimizer constructor (`optimizer_factory`, e.g. `HiGHS.Optimizer`) used to build the
JuMP model. `solver_time_limit` sets a per-solve time limit in seconds (if supported by
the optimizer). It defaults to `$(DEFAULT_RETENTION_MAPPING_SOLVER_TIME_LIMIT)` seconds;
pass `solver_time_limit=nothing` to disable the time limit or another positive number to
override it. If the fitted spline is not strictly monotonic on the validation grid, the
function throws an `ArgumentError`; use a larger `λ` or curate the calibration points.

Returns a `RetentionMapper` containing the fitted spline, the chosen `λ`, and metadata.

See also 
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper), 
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`extras`](@ref extras(::AbstractRetentionMapper)),
[`fitmap`](@ref),
[`invmap`](@ref),
[`invmapmax`](@ref),
[`invmapmin`](@ref),
[`mapmax`](@ref),
[`mapmin`](@ref),
[`rawapplymap`](@ref),
[`rawderivinvmap`](@ref),
[`rawderivmap`](@ref),
[`rawinvmap`](@ref),
[`rawinvmapmax`](@ref),
[`rawinvmapmin`](@ref),
[`rawmapmax`](@ref),
[`rawmapmin`](@ref),
[`rawretentions_A`](@ref), 
[`rawretentions_B`](@ref), 
[`retentions_A`](@ref), 
[`retentions_B`](@ref), 
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref). 

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute";
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0];

julia> mapper = fitmap(retention_times, retention_indices, 
                       extras=Dict("instrument" => "GC-MS", "column" => "DB-5"));

julia> mapped = applymap(mapper, 5.0u"minute");
       invmap(mapper, mapped) ≈ 5.0u"minute"
true

julia> retention = invmap(mapper, 150);
       applymap(mapper, retention) ≈ 150
true
```
"""
function fitmap(
    retentions_A::AbstractVector{<:Union{<:Real, <:AbstractQuantity{<:Real}}},
    retentions_B::AbstractVector{<:Union{<:Real, <:AbstractQuantity{<:Real}}};
    smoothness_penalty_n::Integer=100,
    monotonicity_grid_size::Integer=10^5,
    λ::Real=3e-9,
    solver_time_limit::Union{Nothing, Real}=DEFAULT_RETENTION_MAPPING_SOLVER_TIME_LIMIT,
    extras::AbstractDict{<:AbstractString, <:Any}=Dict{String, Any}(),
    optimizer_factory=HiGHS.Optimizer
    )
    
    converted_extras = Dict{String, Any}(string(k) => v for (k, v) in extras)

    # Validate input data: ensure timepoints are usable for monotonic spline fitting
    length(retentions_A) > 2 || throw(ArgumentError(
        "retentions_A must have at least 3 points"))
    length(retentions_A) == length(retentions_B) || throw(ArgumentError(
        "retentions_A and retentions_B must have the same length"))
    all(diff(ustrip.(retentions_A)) .> 0) || throw(ArgumentError(
        "retentions_A must be strictly increasing"))
    all(u -> unit(u) == unit(first(retentions_A)), retentions_A) || throw(ArgumentError(
        "retentions_A must all share the same unit"))
    all(u -> unit(u) == unit(first(retentions_B)), retentions_B) || throw(ArgumentError(
        "retentions_B must all share the same unit"))
    if !isnothing(solver_time_limit)
        isfinite(solver_time_limit) && solver_time_limit > 0 || throw(ArgumentError(
            "solver_time_limit must be nothing or a finite positive number"))
    end
    isfinite(λ) && λ ≥ 0 || throw(ArgumentError(
        "λ must be finite and nonnegative"))
    monotonicity_grid_size ≥ 2 || throw(ArgumentError(
        "monotonicity_grid_size must be at least 2"))

    # Strip units for numerical processing
    if eltype(retentions_A) <: AbstractQuantity
        rA = ustrip.(retentions_A)
        rA_unit = unit(first(retentions_A))
    else
        rA = copy(retentions_A)
        rA_unit = nothing
    end
   
    if eltype(retentions_B) <: AbstractQuantity
        rB = ustrip.(retentions_B)
        rB_unit = unit(first(retentions_B))
    else
        rB = copy(retentions_B)
        rB_unit = nothing
    end

    # Get value ranges for input scaling
    rA_min, rA_max = extrema(rA)
    rB_min, rB_max = extrema(rB)

    # Normalize inputs to [0, 1] for numerical stability (force Float64)
    rA_norm = Float64.(minmax_scale.(rA, rA_min, rA_max))
    rB_norm = Float64.(minmax_scale.(rB, rB_min, rB_max))
    rA_norm_min, rA_norm_max = extrema(rA_norm)

    # Create a cubic spline basis over normalized domain
    spline_order = 4
    n_points = length(rA_norm)
    n_knots = n_points + spline_order

    function build_basis(n_knots::Int)
        knots = range(rA_norm_min, stop=rA_norm_max, length=n_knots)
        B = RecombinedBSplineBasis(BSplineBasis(BSplineOrder(spline_order), knots), Natural())
        return B, knots
    end

    B, knots = build_basis(n_knots)
    coefs_n = length(B)

    # Reduce knot count until coefs_n ≤ n_points (avoid underdetermined fit)
    while coefs_n > n_points && n_knots > spline_order + 2
        n_knots -= 1
        B, knots = build_basis(n_knots)
        coefs_n = length(B)
    end

    # Construct interpolation matrix (0th derivative)
    C = collocation_matrix(B, rA_norm, Derivative(0), SparseMatrixCSC{Float64})

    # Construct matrix for monotonicity constraint (1st derivative)
    D1 = collocation_matrix(B, LinRange(rA_norm_min, rA_norm_max, coefs_n), 
        Derivative(1), SparseMatrixCSC{Float64})
    
    # Construct matrix for smoothness penalty (2nd derivative)
    smooth_n = min(smoothness_penalty_n, 5 * coefs_n)
    smooth_n = max(smooth_n, coefs_n)
    D2 = collocation_matrix(B, LinRange(rA_norm_min, rA_norm_max, smooth_n), 
        Derivative(2), SparseMatrixCSC{Float64})

    # Build optimizer (optionally set time limit)
    optimizer = if solver_time_limit ≡ nothing
        optimizer_factory
    else
        () -> begin
            opt = optimizer_factory()
            try
                MOI.set(opt, MOI.TimeLimitSec(), Float64(solver_time_limit))
            catch e
                isa(e, MOI.UnsupportedAttribute) || rethrow()
            end
            opt
        end
    end

    # Define spline fitting procedure with L2 penalty and monotonicity constraint
    # Add a tiny ridge term for numerical stability of the QP.
    ridge = 1e-10 * sum(abs2, C)
    function fit_spline_with_penalty(λ)
        # Set up optimization model
        model = Model(optimizer)
        try
            set_silent(model)
        catch e
            isa(e, MOI.UnsupportedAttribute) || rethrow()
        end

        # Define spline coefficients as decision variables
        @variable(model, coefs[1:coefs_n])

        # Objective: squared error + λ * smoothness penalty + tiny ridge
        @objective(model, Min,
            sum((C * coefs .- rB_norm).^2) + λ * sum((D2 * coefs).^2) + ridge * sum(coefs.^2))

        # Enforce monotonicity via first derivative ≥ 0 over grid
        for i in axes(D1, 1)
            @constraint(model, sum(D1[i, j] * coefs[j] for j in 1:coefs_n) ≥ 0)
        end

        # Solve the QP
        optimize!(model)

        # Return fitted coefficients and solver status
        status = termination_status(model)
        pstatus = MOI.get(model, MOI.PrimalStatus())
        ok = status == MOI.OPTIMAL ||
             (status == MOI.TIME_LIMIT &&
              (pstatus == MOI.FEASIBLE_POINT || pstatus == MOI.NEARLY_FEASIBLE_POINT))
        ok || throw(OptimizationError(λ, status))
        value.(coefs), status
    end

    coefs, status = fit_spline_with_penalty(λ)

    # Abort if solver did not return an optimal solution
    status == MOI.OPTIMAL || throw(OptimizationError(λ, status))

    # Construct the final fitted spline
    spline = Spline(B, coefs)
    check_grid = LinRange(rA_norm_min, rA_norm_max, monotonicity_grid_size)
    is_strictly_monotonic(spline, check_grid) || throw(ArgumentError(
        "fixed λ = $(λ) does not produce a strictly monotonic spline; " *
        "increase λ or curate the calibration points"))

    # Convert predicted y-range back to original scale
    rB_pred_norm_min = spline(rA_norm_min)
    rB_pred_norm_max = spline(rA_norm_max)
    rB_pred_min = inverse_minmax_scale(rB_pred_norm_min, rB_min, rB_max)
    rB_pred_max = inverse_minmax_scale(rB_pred_norm_max, rB_min, rB_max)

    # Return a RetentionMapper object holding all necessary spline info and metadata
    RetentionMapper(
        rA, rA_unit, rA_min, rA_max, rA_norm_min, rA_norm_max,
        rB, rB_unit, rB_min, rB_max, rB_pred_min, rB_pred_max, 
        rB_pred_norm_min, rB_pred_norm_max,
        knots, coefs, spline,
        λ,
        converted_extras
    )
end

function is_strictly_monotonic(spline::Spline, x::AbstractVector{<:Real}; tol=1e-8)
    y = spline.(x)
    all(diff(y) .> tol)
end
