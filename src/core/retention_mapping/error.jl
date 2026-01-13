struct OptimizationError <: Exception
    λ::Real
    status
end

Base.show(io::IO, e::OptimizationError) = 
    print(io, "OptimizationError: Spline optimization failed for λ = $(e.λ). " * 
        "Solver status: $(e.status)")

Base.showerror(io::IO, e::OptimizationError) = 
    print(io, "Spline optimization failed for λ = $(e.λ). Solver status: $(e.status). " *
              "Try increasing λ_max (more smoothing), relaxing monotonicity_grid_size, " *
              "or verifying that anchor points are strictly increasing.")
