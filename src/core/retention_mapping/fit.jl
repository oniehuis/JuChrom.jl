
# ── fitmap ────────────────────────────────────────────────────────────────────────────────

"""
    fitmap(retentions_A::AbstractVector{<:Union{<:Real, <:AbstractQuantity{<:Real}}},
           retentions_B::AbstractVector{<:Union{<:Real, <:AbstractQuantity{<:Real}}};
           smoothness_penalty_n::Integer=100,
           monotonicity_grid_size::Integer=10^5,
           λ_max::Real=1e2,
           λ_min::Real=1e-20,
           logλ_tolerance::Real=1e-4,
           max_lambda_iters::Integer=100,
           extras::AbstractDict{<:AbstractString, <:Any}=Dict{String, Any}(),
           optimizer_factory=HiGHS.Optimizer)

Constructs a strictly increasing, smooth cubic B-spline that maps values from domain A 
(`retentions_A`) to domain B (`retentions_B`), subject to a monotonicity constraint and 
smoothness penalty.

The fitted spline minimizes a regularized least-squares objective: data fidelity
(squared residuals) and a smoothness penalty (second derivative norm), subject to
non-negativity of the first derivative.

The smoothing parameter `λ` is automatically determined via binary search: it selects the 
smallest value that yields a monotonic spline while fitting the data as closely as possible, 
thus ensuring both stability and high fidelity without introducing spurious oscillations.

Inputs are monotonic `retentions_A`/`retentions_B` vectors (unitful or unitless), optional
controls for smoothness (`smoothness_penalty_n`, `λ_min`, `λ_max`, `max_lambda_iters`) and
monotonicity enforcement (`monotonicity_grid_size`), optional metadata (`extras`), and an
optimizer constructor. `logλ_tolerance` is the stopping tolerance on `log10(λ)` used by the
binary search that tunes the smoothing parameter: when the log-range width falls below
this value, the search stops (smaller values mean tighter λ tuning at higher cost).

Returns a `RetentionMapper` containing the fitted spline and metadata.

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
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0];

julia> mapper = fitmap(retention_times, retention_indices, 
                      extras=Dict("instrument" => "GC-MS", "column" => "DB-5"));

julia> applymap(mapper, 5.0u"minute") ≈ 338.71144090340465
true

julia> invmap(mapper, 150) ≈ 1.8368240074613091u"minute"
true
```
"""
function fitmap(
    retentions_A::AbstractVector{<:Union{<:Real, <:AbstractQuantity{<:Real}}},
    retentions_B::AbstractVector{<:Union{<:Real, <:AbstractQuantity{<:Real}}};
    smoothness_penalty_n::Integer=100,
    monotonicity_grid_size::Integer=10^5,
    λ_max::Real=1e2,
    λ_min::Real=1e-20,
    logλ_tolerance::Real=1e-4,
    max_lambda_iters::Integer=100,
    extras::AbstractDict{<:AbstractString, <:Any}=Dict{String, Any}(),
    optimizer_factory=HiGHS.Optimizer,
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

    # Normalize inputs to [0, 1] for numerical stability
    rA_norm = minmax_scale.(rA, rA_min, rA_max)
    rB_norm = minmax_scale.(rB, rB_min, rB_max)
    rA_norm_min, rA_norm_max = extrema(rA_norm)

    # Create a cubic spline basis over normalized domain
    spline_order = 4
    knots = range(rA_norm_min, stop=rA_norm_max, length=length(rA_norm) + spline_order)
    B = RecombinedBSplineBasis(BSplineBasis(BSplineOrder(spline_order), knots), Natural())
    coefs_n = length(B)

    # Construct interpolation matrix (0th derivative)
    C = collocation_matrix(B, rA_norm, Derivative(0), SparseMatrixCSC{Float64})

    # Construct matrix for monotonicity constraint (1st derivative)
    D1 = collocation_matrix(B, LinRange(rA_norm_min, rA_norm_max, coefs_n), 
        Derivative(1), SparseMatrixCSC{Float64})
    
    # Construct matrix for smoothness penalty (2nd derivative)
    D2 = collocation_matrix(B, LinRange(rA_norm_min, rA_norm_max, smoothness_penalty_n), 
        Derivative(2), SparseMatrixCSC{Float64})

    # Define spline fitting procedure with L2 penalty and monotonicity constraint
    function fit_spline_with_penalty(λ)
        # Set up optimization model
        model = Model(optimizer_factory)
        try
            set_silent(model)
        catch e
            isa(e, MOI.UnsupportedAttribute) || rethrow()
        end

        # Define spline coefficients as decision variables
        @variable(model, coefs[1:coefs_n])

        # Objective: squared error + λ * smoothness penalty
        @objective(model, Min, sum((C * coefs .- rB_norm).^2) + λ * sum((D2 * coefs).^2))

        # Enforce monotonicity via first derivative ≥ 0 over grid
        for i in axes(D1, 1)
            @constraint(model, sum(D1[i, j] * coefs[j] for j in 1:coefs_n) ≥ 0)
        end

        # Solve the QP
        optimize!(model)

        # Return fitted coefficients and solver status
        status = termination_status(model)
        status == MOI.OPTIMAL || throw(OptimizationError(λ, status))
        value.(coefs), status
    end

    # Find the smallest λ that yields a strictly monotonic spline
    λ, coefs, status = tune_lambda_for_monotonic_spline(
        fit_spline_with_penalty, B, rA_norm_min, rA_norm_max,
        λ_min, λ_max, logλ_tolerance, max_lambda_iters, monotonicity_grid_size)

    # Abort if solver did not return an optimal solution
    status == MOI.OPTIMAL || throw(OptimizationError(λ, status))

    # Construct the final fitted spline
    spline = Spline(B, coefs)

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
        converted_extras
    )
end

function is_strictly_monotonic(spline::Spline, x::AbstractVector{<:Real}; tol=1e-8)
    y = spline.(x)
    all(diff(y) .> tol)
end

function try_fit_and_check_monotonicity(
    λ::Real, 
    fit_spline_with_penalty::F, 
    B::AbstractBSplineBasis, 
    check_grid::AbstractVector{<:Real}
    ) where {F}

    # Fit the spline coefficients using the provided fitting function with penalty λ
    coefs, status = fit_spline_with_penalty(λ)

    # Construct a spline object from the basis and the fitted coefficients
    spline = Spline(B, coefs)

    # Check if the optimization was successful (status == MOI.OPTIMAL)
    # and if the resulting spline is strictly monotonic over the check_grid points
    is_monotonic = status == MOI.OPTIMAL && is_strictly_monotonic(spline, check_grid)

    # Return whether the spline is monotonic, the coefficients, and the solver status
    is_monotonic, coefs, status
end

function tune_lambda_for_monotonic_spline(
    fit_spline_with_penalty::F, 
    B::AbstractBSplineBasis, 
    rA_norm_min::Real, 
    rA_norm_max::Real, 
    λ_min::Real, 
    λ_max::Real, 
    logλ_tolerance::Real, 
    max_lambda_iters::Integer,
    monotonicity_grid_size::Integer
    ) where {F}

    # Create a dense grid of points over the normalized retention time domain
    # where monotonicity will be checked.
    check_grid = LinRange(rA_norm_min, rA_norm_max, monotonicity_grid_size)

    # Step 1: Check if the spline fitted with the highest smoothing penalty (λ_max) produces
    # a monotonic spline. This confirms that monotonicity is achievable at the upper bound.
    λ = λ_max
    is_feasible, _, _ = try_fit_and_check_monotonicity(λ, fit_spline_with_penalty, B, 
        check_grid)
    if !is_feasible
        throw(ArgumentError(
            "λ_max = $(λ) does not produce a feasible monotonic spline. Increase λ_max."))
    end

    # Step 2: Check if the spline fitted with the lowest smoothing penalty (λ_min) is non-
    # monotonic. This confirms the lower bound is meaningful for the search.
    λ = λ_min
    is_feasible, _, _ = try_fit_and_check_monotonicity(λ, fit_spline_with_penalty, B, 
        check_grid)
    # if is_feasible
    #     throw(ArgumentError(
    #         "λ_min = $(λ) already yields a feasible monotonic spline. Decrease λ_min."))
    # end

    # Initialize variables for storing best spline coefficients and optimization status
    coefs, status = nothing, nothing

    # Convert λ_min and λ_max to logarithmic scale for numerical stability and efficient 
    # searching
    logλ_max = log10(λ_max)
    logλ_min = log10(λ_min)
    converged = false

    # Binary search loop to find the minimal λ that yields a monotonic spline
    for i in 1:max_lambda_iters
        # Midpoint of the current search interval in log space
        logλ_mid = (logλ_min + logλ_max) / 2
        # Convert midpoint back to original scale
        λ_mid = 10.0^logλ_mid

        # Fit spline with current λ_mid and check monotonicity
        is_feasible, coefs_tmp, status_tmp = try_fit_and_check_monotonicity(λ_mid, 
            fit_spline_with_penalty, B, check_grid)

        # If monotonic, narrow the upper bound to λ_mid to find a smaller λ
        if is_feasible
            logλ_max = logλ_mid
            coefs, status = coefs_tmp, status_tmp

        # If monotonic, narrow the upper bound to λ_mid to find a smaller λ
        else
            logλ_min = logλ_mid
        end

        if abs(logλ_max - logλ_min) < logλ_tolerance
            converged = true
            break
        end
    end

    # Issue a warning if binary search did not converge within max iterations
    if !converged
        @warn ("Reached maximum iterations ($max_lambda_iters) without satisfying " * "
            logλ_tolerance = $logλ_tolerance")
    end

    # Convert best found λ back to original scale
    best_λ = 10.0^logλ_max

    # Return the best found λ along with corresponding spline coefficients and optimization 
    # status
    best_λ, coefs, status
end
