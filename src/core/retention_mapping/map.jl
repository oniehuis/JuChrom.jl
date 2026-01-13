# ── applymap ──────────────────────────────────────────────────────────────────────────────

function _compute_retention_mapping(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T=1e-8,
    warn::Bool=false
) where {T<:Real}

    if retention isa AbstractQuantity
        rA = ustrip(rm.rA_unit, retention)
    else
        rA = retention
    end
    
    # Normalize rt to [0,1] scale based on original rt range
    rA_norm = minmax_scale(rA, rm.rA_min, rm.rA_max)
    
    # Handle extrapolation logic (your existing code)
    if rA_norm ≤ rm.rA_norm_min
        if warn && rA_norm ≤ rm.rA_norm_min - domain_boundary_threshold
            @info ("rA=$rA is below the spline domain boundary (rA_min=$(rm.rA_min)): " * 
                   "extrapolating left.")
        end
        spline′ = Derivative(1) * rm.spline
        slope = spline′(rm.rA_norm_min)
        delta = rA_norm - rm.rA_norm_min
        rB_norm = rm.spline(rm.rA_norm_min) + slope * delta
    elseif rA_norm ≥ rm.rA_norm_max
        if warn && rA_norm ≥ rm.rA_norm_max + domain_boundary_threshold
            @info ("rA=$rA is above the spline domain boundary (rA_max=$(rm.rA_max)): " * 
                   "extrapolating right.")
        end
        spline′ = Derivative(1) * rm.spline
        slope = spline′(rm.rA_norm_max)
        delta = rA_norm - rm.rA_norm_max
        rB_norm = rm.spline(rm.rA_norm_max) + slope * delta
    else
        rB_norm = rm.spline(rA_norm)
    end
    
    # Return raw result
    inverse_minmax_scale(rB_norm, rm.rB_min, rm.rB_max)
end

"""
    applymap(
        rm::RetentionMapper,
        retention::Union{<:Real, <:AbstractQuantity{<:Real}}; 
        domain_boundary_threshold::Real=1e-8,
        warn::Bool=false, 
        unit::Union{<:Nothing, Unitful.Units}=retentionunit_B(rm)
    )

Apply the fitted mapping to transform a retention value from domain A to domain B using the 
fitted spline within the `RetentionMapper`.

Inputs are a fitted `rm`, a retention value (unitful or unitless), and `unit` for output
conversion. `domain_boundary_threshold` sets how far beyond the normalized spline
domain a value must lie before logging an extrapolation warning when `warn=true`.

Out-of-domain values are extrapolated linearly using the boundary slope, and the result
is returned on the original scale with appropriate units.

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
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
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]u"Th"
       mapper = fitmap(retention_times, retention_indices);

julia> applymap(mapper, 5.0u"minute") ≈ 338.71144090340465u"Th"
true

julia> applymap(mapper, 0.5u"minute", warn=true) ≈ 44.68816899287087u"Th"
[ Info: rA=0.5 is below the spline domain boundary (rA_min=1.2): extrapolating left.
true
```
"""
function applymap(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T1=1e-8,
    warn::Bool=false,
    unit::T2=retentionunit_B(rm)
    ) where {T1<:Real, T2<:Union{Nothing, Unitful.Units}}

    # Validate input unit compatibility (ADD THIS)
    if rm.rA_unit isa Unitful.Units && !(retention isa AbstractQuantity)
        throw(ArgumentError("Input must have units when rm.rA_unit is a Unitful unit. " *
                           "Expected units compatible with $(rm.rA_unit)."))
    elseif rm.rA_unit isa Nothing && retention isa AbstractQuantity
        throw(ArgumentError("Input cannot have units when rm.rA_unit is Nothing. " *
                           "Expected unitless input."))
    end
    
    raw_result = _compute_retention_mapping(rm, retention; 
        domain_boundary_threshold=domain_boundary_threshold, warn=warn)
    rm.rB_unit isa Unitful.Units ? uconvert(unit, raw_result * rm.rB_unit) : raw_result
end


"""
    applymap(
        rmap::RetentionMapper,
        series::MassScanSeries;
        domain_boundary_threshold::Real=1e-8,
        unit::Union{<:Nothing, Unitful.Units}=retentionunit_B(rmap),
        warn::Bool=false
    ) -> MassScanSeries

Apply retention mapping `rmap` to `series`, transforming retention coordinates and scaling 
intensities by the Jacobian to preserve peak areas.

For parameter descriptions and basic behavior, see [`applymap`](@ref) for single retention 
values.

Returns a new `MassScanSeries` with mapped retentions and Jacobian-corrected intensities.

`applymap` treats intensities as densities with respect to retention coordinates. Monotonic 
reparameterization A→B rescales intensities by 1/(dB/dA) to preserve peak areas. If your 
data are integrated counts per scan with variable dwell times, convert to densities first 
using dwell time before calling `applymap`. Only numerical values are transformed — units 
remain unchanged. If intensities are densities with retention units in the denominator, 
handle unit transformation separately.

# Examples
"""
function applymap(
    rmap::RetentionMapper,
    series::MassScanSeries;
    domain_boundary_threshold::T1=1e-8,
    unit::T2=retentionunit_B(rmap),
    warn::Bool=false,
    ) where {T1<:Real, T2<:Union{Nothing, Unitful.Units}}

    @assert scancount(series) > 0 "Cannot transform empty series."
    
    # Use the first scan as a reference for type inference
    reference_scan = first(scans(series))
    
    # Infer types from the reference scan
    retention_type = Float64
    retention_unit_type = typeof(unit)
    mzvalues_type = typeof(rawmzvalues(reference_scan))
    mz_unit_type = typeof(mzunit(reference_scan))
    intensity_unit_type = typeof(intensityunit(reference_scan))
    level_type = typeof(level(reference_scan))
    metadata_type = typeof(attrs(reference_scan))
    
    # Define concrete MassScan type - always use Vector{Float64} for intensities
    MassScanT = MassScan{
        retention_type,
        retention_unit_type,
        mzvalues_type,
        mz_unit_type,
        Vector{Float64},
        intensity_unit_type,
        level_type,
        metadata_type
    }

    # Preallocate new scan array
    new_mscans = Vector{MassScanT}(undef, scancount(series))
    
    # Iterate over each scan in the series
    for (i, scan) in enumerate(series)
        new_retention = applymap(rmap, retention(scan);
                                 domain_boundary_threshold=domain_boundary_threshold,
                                 unit=unit,
                                 warn=warn)
        
        
        # Compute Jacobian
        jacobian = rawderivmap(rmap, retention(scan), 
                               domain_boundary_threshold=domain_boundary_threshold,
                               rB_unit=unit,
                               warn=false)

        # Validate Jacobian
        isfinite(jacobian) && jacobian > 0 || throw(
            ArgumentError("Jacobian must be finite and positive."))

        # Correct intensities
        new_ints = Float64.(intensities(scan) ./ jacobian)
        
        # Construct new MassScan with transformed retention and corrected intensities
        mscan = MassScan(new_retention,
            deepcopy(unit),
            deepcopy(rawmzvalues(scan)),
            mzunit(scan),
            new_ints,
            intensityunit(scan),
            level=deepcopy(level(scan)),
            attrs=deepcopy(attrs(scan))
        )
        
        # Store new MassScan
        new_mscans[i] = mscan
    end
    
    # Construct and return new MassScanSeries
    MassScanSeries(
        new_mscans,
        instrument = deepcopy(instrument(series)),
        acquisition = deepcopy(acquisition(series)),
        user = deepcopy(user(series)),
        sample = deepcopy(sample(series)),
        extras = deepcopy(extras(series))
    )
end

# """
#     applymap(rmap::RetentionMapper,
#              msm::MassScanMatrix,
#              variances::AbstractMatrix{V};
#              domain_boundary_threshold::Real=1e-8,
#              unit=rmap.rB_unit,
#              warn::Bool=false)
#         -> (MassScanMatrix)

# Apply retention mapping `rmap` to `msm`, transforming retention coordinates and scaling 
# intensities by the Jacobian to preserve peak areas.

# For parameter descriptions and basic behavior, see [`applymap`](@ref) for single retention 
# values.

# # Returns
# - New `MassScanMatrix` with mapped retentions and Jacobian-corrected intensities

# # Important Notes
# `applymap` treats intensities as densities with respect to retention coordinates. Monotonic 
# reparameterization A→B rescales intensities by 1/(dB/dA) to preserve peak areas. If your 
# data are integrated counts per scan with variable dwell times, convert to densities first 
# using dwell time before calling `applymap`. Only numerical values are transformed — units 
# remain unchanged. If intensities are densities with retention units in the denominator, 
# handle unit transformation separately.
# """
# function applymap(
#     rmap::RetentionMapper,
#     msm::MassScanMatrix;
#     domain_boundary_threshold::T1=1e-8,
#     unit::T2=rmap.rB_unit,
#     warn::Bool=false
# ) where {
#         T1<:Real,
#         T2<:Union{Nothing, Unitful.Units}
#     }
    
#     # Get original retention coordinates
#     ret = retentions(msm)

#     # Map retentions to new coordinates
#     new_rets = applymap.(rmap, ret;
#                          domain_boundary_threshold=domain_boundary_threshold,
#                          unit=unit, warn=warn)

#     # Determine new retention unit and strip units
#     new_retention_unit = 
#         first(new_rets) isa AbstractQuantity ? Unitful.unit(first(new_rets)) : nothing
#     new_retention_unitfree = ustrip.(new_rets)

#     # Compute Jacobians for each retention coordinate
#     jacobians = rawderivmap.(rmap, ret;
#                              domain_boundary_threshold = domain_boundary_threshold,
#                              rB_unit = unit,
#                              warn = false)

#     # Validate Jacobians
#     for j in jacobians
#         (isfinite(j) && j > 0) || throw(
#             ArgumentError("Jacobian values must be finite and positive."))
#     end

#     # Extract raw intensities and allocate output arrays
#     rawints = rawintensities(msm)
#     new_rawints = similar(rawints, Float64)
#     @inbounds @simd for c in axes(rawints, 2)
#         new_rawints[:, c] = rawints[:, c] ./ jacobians
#     end

#     # Construct and return MassScanMatrix
#     MassScanMatrix(
#         new_retention_unitfree,
#         new_retention_unit,
#         deepcopy(rawmzvalues(msm)),
#         mzunit(msm),
#         new_rawints,
#         intensityunit(msm),
#         level=deepcopy(level(msm)),
#         instrument=deepcopy(instrument(msm)),
#         acquisition=deepcopy(acquisition(msm)),
#         user=deepcopy(user(msm)),
#         sample=deepcopy(sample(msm)),
#         extras=deepcopy(extras(msm))
#     )
# end

"""
    applymap(
        rmap::RetentionMapper,
        msm::MassScanMatrix;
        domain_boundary_threshold::Real=1e-8,
        unit::Union{<:Nothing, Unitful.Units}=retentionunit_B(rmap),
        warn::Bool=false
    ) -> MassScanMatrix

Apply retention mapper `rmap` to `msm`, transforming retention coordinates and
scaling intensities by the Jacobian to preserve peak areas.

For additional parameter descriptions and basic behavior, see [`applymap`](@ref) for 
single retention values.

Returns a new `MassScanMatrix` with mapped retentions and Jacobian-corrected intensities.

`applymap` treats intensities as densities with respect to retention coordinates. Monotonic 
reparameterization A→B rescales intensities by 1/(dB/dA) to preserve peak areas. If your 
data are integrated counts per scan with variable dwell times, convert to densities first 
using dwell time before calling `applymap`. Only numerical values are transformed — units 
remain unchanged. If intensities are densities with retention units in the denominator, 
handle unit transformation separately.
"""
function applymap(
    rmap::RetentionMapper,
    msm::MassScanMatrix;
    domain_boundary_threshold::T1=1e-8,
    unit::T2=retentionunit_B(rmap),
    warn::Bool=false
    ) where {T1<:Real, T2<:Union{<:Nothing, Unitful.Units}}

    # Get original retention coordinates
    ret = retentions(msm)

    # Map retentions to new coordinates
    new_rets = applymap.(rmap, ret;
                         domain_boundary_threshold=domain_boundary_threshold,
                         unit=unit, warn=warn)

    # Determine new retention unit and strip units
    new_retention_unit = 
        first(new_rets) isa AbstractQuantity ? Unitful.unit(first(new_rets)) : nothing
    new_retention_unitfree = ustrip.(new_rets)

    # Compute Jacobians for each retention coordinate
    jacobians = rawderivmap.(rmap, ret;
                             domain_boundary_threshold = domain_boundary_threshold,
                             rB_unit = unit,
                             warn = false)

    # Validate Jacobians
    for j in jacobians
        (isfinite(j) && j > 0) || throw(
            ArgumentError("Jacobian values must be finite and positive."))
    end

    # Extract raw intensities and allocate output arrays
    rawints = rawintensities(msm)
    new_rawints = similar(rawints, Float64)

    # Correct intensities for each scan
    @inbounds @simd for c in axes(rawints, 2)
        new_rawints[:, c] = rawints[:, c] ./ jacobians
    end

    # Return new MassScanMatrix
    MassScanMatrix(
        new_retention_unitfree,
        new_retention_unit,
        deepcopy(rawmzvalues(msm)),
        mzunit(msm),
        new_rawints,
        intensityunit(msm),
        level=deepcopy(level(msm)),
        instrument=deepcopy(instrument(msm)),
        acquisition=deepcopy(acquisition(msm)),
        user=deepcopy(user(msm)),
        sample=deepcopy(sample(msm)),
        extras=deepcopy(extras(msm))
    )
end

# ── derivinvmap ───────────────────────────────────────────────────────────────────────────

function _compute_inverse_derivative_mapping(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T=1e-8,
    tol::Real=1e-8,
    warn::Bool=false
) where {T<:Real}

    if retention isa AbstractQuantity
        rB = ustrip(rm.rB_unit, retention)
    else
        rB = retention
    end

    # Normalize retention index to spline domain
    rB_norm = minmax_scale(rB, rm.rB_min, rm.rB_max)

    # Construct first derivative of the fitted spline
    spline′ = Derivative(1) * rm.spline

    # Handle extrapolation if input is below spline domain
    if rB_norm ≤ rm.rB_pred_norm_min
        if warn && rB_norm ≤ rm.rB_pred_norm_min - domain_boundary_threshold
            rB_min_str = Printf.@sprintf("%.6f", rm.rB_pred_min)
            @info ("rB=$rB is below the inverse spline domain boundary " *
                   "(rB_min=$rB_min_str): extrapolating left.")
        end
        d_rA_d_rB_norm = spline′(rm.rA_norm_min)

    # Handle extrapolation if input is above spline domain
    elseif rB_norm ≥ rm.rB_pred_norm_max
        if warn && rB_norm ≥ rm.rB_pred_norm_max + domain_boundary_threshold
            rB_max_str = Printf.@sprintf("%.6f", rm.rB_pred_max)
            @info ("rB=$rB is above the inverse spline domain boundary " *
                   "(rB_max=$rB_max_str): extrapolating right.")
        end
        d_rA_d_rB_norm = spline′(rm.rA_norm_max)
    else
        # If within domain, find normalized time using root finding
        g(x) = rm.spline(x) - rB_norm
        rA_norm = find_zero(g, (rm.rA_norm_min, rm.rA_norm_max), Bisection(), tol=tol)
        d_rA_d_rB_norm = spline′(rA_norm)
    end

    # Invert and rescale derivative to obtain d(rA)/d(rB)
    return (1 / d_rA_d_rB_norm) * ((rm.rA_max - rm.rA_min) / (rm.rB_max - rm.rB_min))
end

"""
    derivinvmap(
        rmap::RetentionMapper,
        retention::Union{<:Real, <:AbstractQuantity{<:Real}};
        domain_boundary_threshold::Real=1e-8,
        rA_unit::Union{Nothing, Unitful.Units}=retentionunit_A(rmap)),
        rB_unit::Union{Nothing, Unitful.Units}=retentionunit_B(rmap)),
        tol::Real=1e-8, 
        warn::Bool=false
    )

Compute the derivative of the inverse mapping, i.e. d(A)/d(B), at a given output value 
(domain B).

This function computes the local derivative `d(A)/d(B)` at a given output value using 
the inverse of the normalized B-spline approximation stored in the `RetentionMapper`. 
This derivative describes how a small change in the output (domain B) would translate 
back to a change in the input (domain A), under the smoothed, monotonic mapping. If the 
output value is outside the supported domain, extrapolation is performed using the boundary 
derivative. The derivative is adjusted to account for the normalization scales applied 
during spline fitting.

Inputs are a fitted `rm`, an output value in domain B, optional `rA_unit`/`rB_unit` for
unit conversion, and `tol` for root-finding (absolute bisection tolerance when solving the 
inverse; smaller values increase accuracy at higher cost). `domain_boundary_threshold` sets 
how far beyond the normalized inverse spline domain a value must lie before logging an
extrapolation warning when `warn=true`.

Returns the derivative d(A)/d(B) with units of `rA_unit/rB_unit`. This is the reciprocal
of [`derivmap`](@ref), is always positive under the monotonic constraint, and uses
root-finding to locate the corresponding input value. Out-of-domain values use the
boundary derivative.

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivmap`](@ref),
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
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]u"Th"
       mapper = fitmap(retention_times, retention_indices);

julia> derivinvmap(mapper, 350.0u"Th") ≈ 0.02857099580171666u"minute/Th"
true

julia> derivinvmap(mapper, 350.0)
ERROR: ArgumentError: Input must have units when rm.rB_unit is a Unitful unit. Expected units compatible with Th.

julia> test_indices = [250.0, 600.0]u"Th";

julia> derivinvmap.(mapper, test_indices) ≈ [0.015606366497378708, 0.03240806524912834]u"minute/Th"
true

julia> derivinvmap(mapper, 350.0u"Th", rA_unit=u"s", rB_unit=u"kTh") ≈ 1714.2597481029995u"s/kTh"
true

julia> derivinvmap(mapper, 50.0u"Th", warn=true) ≈ 0.012655503630152843u"minute/Th"
[ Info: rB=50.0 is below the inverse spline domain boundary (rB_min=100.000073): extrapolating left.
true
```
"""
function derivinvmap(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T1=1e-8,
    rA_unit::T2=retentionunit_A(rm),
    rB_unit::T3=retentionunit_B(rm),
    tol::T4=1e-8,
    warn::Bool=false
    ) where {
        T1<:Real,
        T2<:Union{Nothing, Unitful.Units},
        T3<:Union{Nothing, Unitful.Units},
        T4<:Real
    }
    # Validate input unit compatibility (FIXED - was missing this check)
    if rm.rB_unit isa Unitful.Units && !(retention isa AbstractQuantity)
        throw(ArgumentError("Input must have units when rm.rB_unit is a Unitful unit. " *
                           "Expected units compatible with $(rm.rB_unit)."))
    elseif rm.rB_unit isa Nothing && retention isa AbstractQuantity
        throw(ArgumentError("Input cannot have units when rm.rB_unit is Nothing. " *
                           "Expected unitless input."))
    end

    # Existing parameter validation (FIXED - was incomplete)
    if rm.rA_unit isa Unitful.Units && rA_unit isa Nothing
        throw(ArgumentError("rA_unit cannot be Nothing if rm.rA_unit is a Unitful unit."))
    elseif rm.rA_unit isa Nothing && rA_unit isa Unitful.Units
        throw(ArgumentError("rA_unit cannot be Unitful unit if rm.rA_unit is Nothing."))
    elseif rm.rB_unit isa Unitful.Units && rB_unit isa Nothing
        throw(ArgumentError("rB_unit cannot be Nothing if rm.rB_unit is a Unitful unit."))
    elseif rm.rB_unit isa Nothing && rB_unit isa Unitful.Units
        throw(ArgumentError("rB_unit cannot be Unitful unit if rm.rB_unit is a Nothing."))
    end
    
    result = _compute_inverse_derivative_mapping(rm, retention; 
        domain_boundary_threshold=domain_boundary_threshold, tol=tol, warn=warn)
    
    # Rescale to physical units and return converted derivative in requested units
    if rA_unit isa Unitful.Units && rB_unit isa Unitful.Units
        result *= rm.rA_unit / rm.rB_unit
        return uconvert((rA_unit / rB_unit), result)
    elseif rA_unit isa Unitful.Units && rB_unit isa Nothing
        result *= rm.rA_unit
        return uconvert(rA_unit, result)
    elseif rA_unit isa Nothing && rB_unit isa Unitful.Units
        result /= rm.rB_unit
        return uconvert(inverse(rB_unit), result)
    else
        return result
    end
end

# ── derivmap ──────────────────────────────────────────────────────────────────────────────

# Core computation logic - always returns unitless derivative
function _compute_derivative_mapping(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T=1e-8,
    warn::Bool=false
    ) where {T<:Real}

    if retention isa AbstractQuantity
        rA = ustrip(rm.rA_unit, retention)
    else
        rA = retention
    end
    
    # Normalize retention time to spline domain
    rA_norm = minmax_scale(rA, rm.rA_min, rm.rA_max)
    
    # Construct first derivative of the fitted spline
    spline′ = Derivative(1) * rm.spline
    
    # Handle extrapolation if input is below spline domain
    if rA_norm ≤ rm.rA_norm_min
        if warn && rA_norm ≤ rm.rA_norm_min - domain_boundary_threshold
            @info ("rA=$rA is below the spline domain boundary (rA_min=$(rm.rA_min)): " * 
                   "extrapolating left.")
        end
        d_rB_d_rA_norm = spline′(rm.rA_norm_min)
    # Handle extrapolation if input is above spline domain
    elseif rA_norm ≥ rm.rA_norm_max
        if warn && rA_norm ≥ rm.rA_norm_max + domain_boundary_threshold
            @info ("rA=$rA is above the spline domain boundary (rA_max=$(rm.rA_max)): " * 
                   "extrapolating right.")
        end
        d_rB_d_rA_norm = spline′(rm.rA_norm_max)
    else
        # If within domain, evaluate derivative in normalized space
        d_rB_d_rA_norm = spline′(rA_norm)
    end
    
    # Transform derivative from normalized to original units
    d_rB_d_rA_norm * (rm.rB_max - rm.rB_min) / (rm.rA_max - rm.rA_min)
end

"""
    derivmap(
        rm::RetentionMapper,
        retention::Union{<:Real, <:AbstractQuantity{<:Real}};
        rA_unit::Union{Nothing, Unitful.Units}=retentionunit_A(rm),
        rB_unit::Union{Nothing, Unitful.Units}=retentionunit_B(rm),
        warn::Bool=false
    )

Compute the first derivative of the fitted mapping with respect to the input (domain A) 
at a given input value.

This function evaluates the derivative of the fitted monotonic B-spline, extrapolating 
using boundary derivatives if the input lies outside the spline domain. The derivative is 
adjusted to account for the normalization scales applied during spline fitting.

Inputs are a fitted `rm`, an input value in domain A, optional `rA_unit`/`rB_unit` for
unit conversion, and `warn` to log extrapolation messages. `domain_boundary_threshold`
sets how far beyond the normalized spline domain a value must lie before logging an
extrapolation warning when `warn=true`.

Returns the derivative d(B)/d(A) with units of `rB_unit/rA_unit`. For retention
time/retention index mappings, this gives the local “retention index per unit time” rate.
The derivative is always positive due to the monotonicity constraint, and extrapolation
uses constant boundary derivatives outside the fitted domain.

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
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
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]
       mapper = fitmap(retention_times, retention_indices);

julia> derivmap(mapper, 5.0u"minute") ≈ 37.42885149584893u"minute^-1"
true

julia> test_times = [3.0, 11.0]u"minute";

julia> derivmap.(mapper, test_times) ≈ [67.10235298491122, 34.639458067769134]u"minute^-1"
true

julia> derivmap(mapper, 5.0u"minute", rA_unit=u"s") ≈ 0.6238141915974822u"s^-1"
true

julia> derivmap(mapper, 0.8u"minute", warn=true) ≈ 79.01700550402535u"minute^-1"
[ Info: rA=0.8 is below the spline domain boundary (rA_min=1.2): extrapolating left.
true
```
"""
function derivmap(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T1=1e-8,
    rA_unit::T2=retentionunit_A(rm),
    rB_unit::T3=retentionunit_B(rm),
    warn::Bool=false
) where {
    T1<:Real,
    T2<:Union{Nothing, Unitful.Units},
    T3<:Union{Nothing, Unitful.Units},
}
    # Validate input unit compatibility (ADD THIS)
    if rm.rA_unit isa Unitful.Units && !(retention isa AbstractQuantity)
        throw(ArgumentError("Input must have units when rm.rA_unit is a Unitful unit. " *
                           "Expected units compatible with $(rm.rA_unit)."))
    elseif rm.rA_unit isa Nothing && retention isa AbstractQuantity
        throw(ArgumentError("Input cannot have units when rm.rA_unit is Nothing. " *
                           "Expected unitless input."))
    end
    
    # Existing parameter validation logic
    if rm.rA_unit isa Unitful.Units && rA_unit isa Nothing
        throw(ArgumentError("rA_unit cannot be Nothing if rm.rA_unit is a Unitful unit."))
    elseif rm.rA_unit isa Nothing && rA_unit isa Unitful.Units
        throw(ArgumentError("rA_unit cannot be Unitful unit if rm.rA_unit is Nothing."))
    elseif rm.rB_unit isa Unitful.Units && rB_unit isa Nothing
        throw(ArgumentError("rB_unit cannot be Nothing if rm.rB_unit is a Unitful unit."))
    elseif rm.rB_unit isa Nothing && rB_unit isa Unitful.Units
        throw(ArgumentError("rB_unit cannot be Unitful unit if rm.rB_unit is a Nothing."))
    end
    
    result = _compute_derivative_mapping(rm, retention; 
        domain_boundary_threshold=domain_boundary_threshold, warn=warn)
    
    # Rescale to physical units and return converted derivative in requested units
    if rA_unit isa Unitful.Units && rB_unit isa Unitful.Units
        result *= rm.rB_unit / rm.rA_unit
        return uconvert((rB_unit / rA_unit), result)
    elseif rA_unit isa Unitful.Units && rB_unit isa Nothing
        result /= rm.rA_unit
        return uconvert(inverse(rA_unit), result)
    elseif rA_unit isa Nothing && rB_unit isa Unitful.Units
        result *= rm.rB_unit
        return uconvert(rB_unit, result)
    else
        return result
    end
end

# ── invmap ────────────────────────────────────────────────────────────────────────────────

function _compute_inverse_mapping(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T=1e-8,
    tol::Real=1e-8,
    warn::Bool=false
) where {T<:Real}

    if retention isa AbstractQuantity
        rB = ustrip(rm.rB_unit, retention)
    else
        rB = retention
    end

    # Normalize retention index to [0,1] based on min and max RI of data
    rB_norm = minmax_scale(rB, rm.rB_min, rm.rB_max)

    # Handle extrapolation if input is below spline domain
    if rB_norm ≤ rm.rB_pred_norm_min
        if warn && rB_norm ≤ rm.rB_pred_norm_min - domain_boundary_threshold
            rB_min_str = Printf.@sprintf("%.6f", rm.rB_pred_min)
            @info ("rB=$rB is below the inverse spline domain boundary " *
                   "(rB_min=$rB_min_str): extrapolating left.")
        end
        # Use spline derivative at left boundary to extend linearly
        spline′ = Derivative(1) * rm.spline
        slope = spline′(rm.rA_norm_min)
        ΔrB_norm = rB_norm - rm.rB_pred_norm_min
        rA_norm = rm.rA_norm_min + ΔrB_norm / slope
    
    # Handle extrapolation if input is above spline domain
    elseif rB_norm ≥ rm.rB_pred_norm_max
        if warn && rB_norm ≥ rm.rB_pred_norm_max + domain_boundary_threshold
            rB_max_str = Printf.@sprintf("%.6f", rm.rB_pred_max)
            @info ("rB=$rB is above the inverse spline domain boundary " *
                   "(rB_max=$rB_max_str): extrapolating right.")
        end
        # Use spline derivative at right boundary to extend linearly
        spline′ = Derivative(1) * rm.spline
        slope = spline′(rm.rA_norm_max)
        ΔrB_norm = rB_norm - rm.rB_pred_norm_max
        rA_norm = rm.rA_norm_max + ΔrB_norm / slope
    else
        # If within domain, use root finding (bisection) to invert spline
        g(x) = rm.spline(x) - rB_norm
        rA_norm = find_zero(g, (rm.rA_norm_min, rm.rA_norm_max), Bisection(), tol=tol)
    end

    # Return raw result
    inverse_minmax_scale(rA_norm, rm.rA_min, rm.rA_max)
end

"""
    invmap(
        rm::RetentionMapper,
        retention::Union{<:Real, <:AbstractQuantity{<:Real}};
        domain_boundary_threshold::Real=1e-8,
        unit::Union{Nothing, Unitful.Units}=retentionunit_A(rm),
        tol::Real=1e-8,
        warn::Bool=false
    )

Invert the fitted mapping by computing the input value (domain A) corresponding to a given 
output value (domain B).

Inputs are a fitted `rm`, an output value in domain B, an optional `unit` for output 
conversion, `tol` for root-finding (absolute bisection tolerance when solving the inverse; 
smaller values increase accuracy at higher cost), and `warn` to log extrapolation messages.

The inversion normalizes the input to the spline domain, extrapolates linearly using the
boundary derivative, uses bisection within the domain, and returns the result on the
original domain A scale with requested units.

This function leverages the strict monotonicity of the fitted spline to provide reliable 
inverse mapping, enabling bidirectional conversion between the two domains.

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`fitmap`](@ref),
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
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]u"Th"
       mapper = fitmap(retention_times, retention_indices);

julia> invmap(mapper, 350.0u"Th") ≈ 5.312425320906996u"minute"
true

julia> invmap(mapper, 350.0)
ERROR: ArgumentError: Input must have units when rm.rB_unit is a Unitful unit. Expected units compatible with Th.
[...]

julia> test_indices = [250.0u"Th", 450.0u"Th"];

julia> invmap.(mapper, test_indices) ≈ [3.2283166885897936, 8.08686182716771]u"minute"
true

julia> invmap(mapper, 400.0u"Th", unit=u"s") ≈ 407.99982102982375u"s"
true

julia> invmap(mapper, 50.0u"Th", warn=true) ≈ 0.5672238965934812u"minute"
[ Info: rB=50.0 is below the inverse spline domain boundary (rB_min=100.000073): extrapolating left.
true
```
"""
function invmap(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T1=1e-8,
    unit::T2=rm.rA_unit,
    tol::T3=1e-8,
    warn::Bool=false
    ) where {
        T1<:Real,
        T2<:Union{Nothing, Unitful.Units},
        T3<:Real
    }
    
    # Validate input unit compatibility
    if rm.rB_unit isa Unitful.Units && !(retention isa AbstractQuantity)
        throw(ArgumentError("Input must have units when rm.rB_unit is a Unitful unit. " *
                           "Expected units compatible with $(rm.rB_unit)."))
    elseif rm.rB_unit isa Nothing && retention isa AbstractQuantity
        throw(ArgumentError("Input cannot have units when rm.rB_unit is Nothing. " *
                           "Expected unitless input."))
    end
    
    raw_result = _compute_inverse_mapping(rm, retention; 
        domain_boundary_threshold=domain_boundary_threshold, tol=tol, warn=warn)
    
    # Add units back if needed
    if rm.rA_unit isa Unitful.Units
        return uconvert(unit, raw_result * rm.rA_unit)
    else
        return raw_result
    end
end

# ── invmapmin ─────────────────────────────────────────────────────────────────────────────

"""
    invmapmin(rm::AbstractRetentionMapper; unit::Union{Nothing, Unitful.Units}=nothing)

Return the minimum predicted value of domain B.

For unitful mappers, the result carries units and can be converted with `unit`. For
unitless mappers, `unit` must be `nothing`; otherwise `ArgumentError` is thrown with the
message "Cannot convert unitless retention B minimum to a unit".

See also [`rawinvmapmin`](@ref), [`invmapmax`](@ref), [`rawinvmapmax`](@ref), 
[`fitmap`](@ref).

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`fitmap`](@ref),
[`invmap`](@ref),
[`invmapmax`](@ref),
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
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> invmapmin(mapper) ≈ 100.00001207241763
true

julia> invmapmin(mapper) ≈ applymap(mapper, mapmin(mapper))
true

julia> (invmapmin(mapper), invmapmax(mapper)) .≈ (100.00001207241763, 499.9999793028839)
(true, true)

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> invmapmin(mapper_unitless) ≈ 100.00001207241763
true

julia> invmapmin(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A minimum to a unit
```
"""
@inline function invmapmin(
    rm::AbstractRetentionMapper{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rB_pred_min, unit, "retention A minimum")
end

@inline function invmapmin(
    rm::AbstractRetentionMapper{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert_scalar(rm.rB_pred_min, rm.rB_unit, unit)
end

# ── invmapmax ─────────────────────────────────────────────────────────────────────────────

"""
    invmapmax(rm::AbstractRetentionMapper; unit::Union{Nothing, Unitful.Units}=nothing)

Return the maximum predicted value of domain B.

For unitful mappers, the result carries units and can be converted with `unit`. For
unitless mappers, `unit` must be `nothing`; otherwise `ArgumentError` is thrown with the
message "Cannot convert unitless retention B maximum to a unit".

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`fitmap`](@ref),
[`invmap`](@ref),
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
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> invmapmax(mapper) ≈ 499.9999793028839
true

julia> invmapmax(mapper) ≈ applymap(mapper, mapmax(mapper))
true

julia> invmap(mapper, invmapmax(mapper)) ≈ mapmax(mapper)
true

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> invmapmax(mapper_unitless) ≈ 499.9999793028839
true

julia> invmapmax(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A minimum to a unit
```
"""
@inline function invmapmax(
    rm::AbstractRetentionMapper{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rB_pred_max, unit, "retention A minimum")
end

@inline function invmapmax(
    rm::AbstractRetentionMapper{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert_scalar(rm.rB_pred_max, rm.rB_unit, unit)
end

# ── mapmax ────────────────────────────────────────────────────────────────────────────────

"""
    mapmax(rm::AbstractRetentionMapper; unit::Union{Nothing, Unitful.Units}=nothing)

Return the maximum input value from domain A used during calibration.

For unitful mappers, the result carries units and can be converted with `unit`. For
unitless mappers, `unit` must be `nothing`; otherwise `ArgumentError` is thrown with the
message "Cannot convert unitless retention A maximum to a unit".

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`fitmap`](@ref),
[`invmap`](@ref),
[`invmapmax`](@ref),
[`invmapmin`](@ref),
[`mapmin`](@ref),
[`rawapplymap`](@ref),
[`rawderivinvmap`](@ref),
[`rawderivmap`](@ref),
[`rawinvmap`](@ref),
[`rawinvmapmax`](@ref),
[`rawinvmapmin`](@ref),
[`rawmapmax`](@ref),
[`rawmapmin`](@ref),
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> mapmax(mapper) ≈ 9.3u"minute"
true

julia> mapmax(mapper, unit=u"s") ≈ 558.0u"s"
true

julia> (mapmin(mapper), mapmax(mapper)) .≈ (1.2u"minute", 9.3u"minute")
(true, true)

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> mapmax(mapper_unitless) ≈ 9.3
true

julia> mapmax(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A minimum to a unit
```
"""
@inline function mapmax(
    rm::AbstractRetentionMapper{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rA_max, unit, "retention A minimum")
end

@inline function mapmax(
    rm::AbstractRetentionMapper{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert_scalar(rm.rA_max, rm.rA_unit, unit)
end

# ── mapmin ────────────────────────────────────────────────────────────────────────────────

"""
    mapmin(rm::AbstractRetentionMapper; unit::Union{Nothing, Unitful.Units}=nothing)

Return the minimum input value from domain A used during calibration.

For unitful mappers, the result carries units and can be converted with `unit`. For
unitless mappers, `unit` must be `nothing`; otherwise `ArgumentError` is thrown with the
message "Cannot convert unitless retention A minimum to a unit".

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`fitmap`](@ref),
[`invmap`](@ref),
[`invmapmax`](@ref),
[`invmapmin`](@ref),
[`mapmax`](@ref),
[`rawapplymap`](@ref),
[`rawderivinvmap`](@ref),
[`rawderivmap`](@ref),
[`rawinvmap`](@ref),
[`rawinvmapmax`](@ref),
[`rawinvmapmin`](@ref),
[`rawmapmax`](@ref),
[`rawmapmin`](@ref),
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> mapmin(mapper) ≈ 1.2u"minute"
true

julia> mapmin(mapper, unit=u"s") ≈ 72.0u"s"
true

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> mapmin(mapper_unitless) ≈ 1.2
true

julia> mapmin(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A minimum to a unit
```
"""
@inline function mapmin(
    rm::AbstractRetentionMapper{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rA_min, unit, "retention A minimum")
end

@inline function mapmin(
    rm::AbstractRetentionMapper{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_convert_scalar(rm.rA_min, rm.rA_unit, unit)
end


# ── rawapplymap ───────────────────────────────────────────────────────────────────────────

"""
    rawapplymap(
        rm::RetentionMapper,
        retention::Union{<:Real, <:AbstractQuantity{<:Real}};
        domain_boundary_threshold::Real=1e-8,
        unit::Union{<:Nothing, Unitful.Units}=retentionunit_B(rm),
        warn::Bool=false
    )

Apply the fitted mapping to transform an input value from domain A to domain B using the
fitted spline within the `RetentionMapper`, returning a unitless numeric result.

Inputs are a fitted `rm`, a unitful retention value, `warn` for extrapolation logs, and
`unit` for output conversion before stripping.

The input is normalized before spline evaluation. Out-of-domain values are extrapolated
linearly using the boundary slope, and the result is returned as a unitless number in the
requested unit system. This method requires unitful input (use [`applymap`](@ref) for
unitless inputs).

This function provides the same transformation as [`applymap`](@ref) but strips units from 
the result, which is useful for downstream calculations that require raw numeric values.

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`fitmap`](@ref),
[`invmap`](@ref),
[`invmapmax`](@ref),
[`invmapmin`](@ref),
[`mapmax`](@ref),
[`mapmin`](@ref),
[`rawderivinvmap`](@ref),
[`rawderivmap`](@ref),
[`rawinvmap`](@ref),
[`rawinvmapmax`](@ref),
[`rawinvmapmin`](@ref),
[`rawmapmax`](@ref),
[`rawmapmin`](@ref),
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]u"Th"
       mapper = fitmap(retention_times, retention_indices);

julia> rawapplymap(mapper, 5.0u"minute") ≈ 338.71144090340465
true

julia> test_times = [3.0, 11.0]u"minute";

julia> rawapplymap.(mapper, test_times) ≈ [235.0195638002354, 564.0302630521933]
true

julia> rawapplymap(mapper, 0.5u"minute", warn=true) ≈ 44.68816899287087
[ Info: rA=0.5 is below the spline domain boundary (rA_min=1.2): extrapolating left.
true
```
"""
function rawapplymap(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T=1e-8,
    unit::Union{<:Nothing, Unitful.Units}=retentionunit_B(rm),
    warn::Bool=false
    ) where {T<:Real}

    # Validate input unit compatibility
    if rm.rA_unit isa Unitful.Units && !(retention isa AbstractQuantity)
        throw(ArgumentError("Input must have units when rm.rA_unit is a Unitful unit. " *
                           "Expected units compatible with $(rm.rA_unit)."))
    elseif rm.rA_unit isa Nothing && retention isa AbstractQuantity
        throw(ArgumentError("Input cannot have units when rm.rA_unit is Nothing. " *
                           "Expected unitless input."))
    end
    
    raw_result = _compute_retention_mapping(rm, retention; 
        domain_boundary_threshold=domain_boundary_threshold, warn=warn)
    return rm.rB_unit isa Unitful.Units ? ustrip(unit, raw_result * rm.rB_unit) : raw_result
end

# ── rawderivinvmap ────────────────────────────────────────────────────────────────────────

"""
    rawderivinvmap(
        rm::RetentionMapper,
        retention::Union{<:Real, <:AbstractQuantity{<:Real}};
        domain_boundary_threshold::Real=1e-8,
        rA_unit::Union{Nothing, Unitful.Units}=retentionunit_A(rm),
        rB_unit::Union{Nothing, Unitful.Units}=retentionunit_B(rm),
        tol::Real=1e-8, 
        warn::Bool=false
    )

Compute the derivative of the inverse mapping, i.e. d(A)/d(B), at a given output value 
(domain B), returning a unitless numeric result.

This function computes the local derivative `d(A)/d(B)` at a given output value using 
the inverse of the normalized B-spline approximation stored in the `RetentionMapper`. 
This derivative describes how a small change in the output (domain B) would translate 
back to a change in the input (domain A), under the smoothed, monotonic mapping. If the 
output value is outside the supported domain, extrapolation is performed using the boundary 
derivative. The derivative is adjusted to account for the normalization scales applied 
during spline fitting.

Inputs are a fitted `rm`, an output value in domain B, optional `rA_unit`/`rB_unit` for 
unit conversion, and `tol` for root-finding (absolute bisection tolerance when solving the 
inverse; smaller values increase accuracy at higher cost). `domain_boundary_threshold` sets 
how far beyond the normalized inverse spline domain a value must lie before logging an
extrapolation warning when `warn=true`.

Returns the derivative d(A)/d(B) as a unitless number in the `rA_unit/rB_unit` system
before stripping. This is the reciprocal of [`rawderivmap`](@ref), is always positive
under the monotonic constraint, and uses root-finding to locate the corresponding input
value. Out-of-domain values use the boundary derivative.

This function provides the same transformation as [`derivinvmap`](@ref) but strips units 
from the result, which is useful for downstream calculations that require raw numeric 
values.

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`fitmap`](@ref),
[`invmap`](@ref),
[`invmapmax`](@ref),
[`invmapmin`](@ref),
[`mapmax`](@ref),
[`mapmin`](@ref),
[`rawapplymap`](@ref),
[`rawderivmap`](@ref),
[`rawinvmap`](@ref),
[`rawinvmapmax`](@ref),
[`rawinvmapmin`](@ref),
[`rawmapmax`](@ref),
[`rawmapmin`](@ref),
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]u"Th"
       mapper = fitmap(retention_times, retention_indices);

julia> rawderivinvmap(mapper, 350.0u"Th") ≈ 0.02857099580171666
true

julia> rawderivinvmap(mapper, 350.0)
ERROR: ArgumentError: Input must have units when rm.rB_unit is a Unitful unit. Expected units compatible with Th.
[...]

julia> test_indices = [250.0, 600.0]u"Th";

julia> rawderivinvmap.(mapper, test_indices) ≈ [0.015606366497378708, 0.03240806524912834]
true

julia> rawderivinvmap(mapper, 350.0u"Th", rA_unit=u"s") ≈ 1.7142597481029995
true

julia> rawderivinvmap(mapper, 50.0u"Th", warn=true) ≈ 0.012655503630152843
[ Info: rB=50.0 is below the inverse spline domain boundary (rB_min=100.000073): extrapolating left.
true

julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]
       mapper = fitmap(retention_times, retention_indices);

julia> rawderivinvmap(mapper, 350.0) ≈ 0.02857099580171666
true

julia> rawderivinvmap(mapper, 350.0u"Th")
ERROR: ArgumentError: Input cannot have units when rm.rB_unit is Nothing. Expected unitless input.
[...]
```
"""
function rawderivinvmap(
   rm::RetentionMapper,
   retention::Union{<:Real, <:AbstractQuantity{<:Real}};
   domain_boundary_threshold::T1=1e-8,
   rA_unit::T2=retentionunit_A(rm),
   rB_unit::T3=retentionunit_B(rm),
   tol::T4=1e-8,
   warn::Bool=false
) where {
    T1<:Real,
    T2<:Union{Nothing, Unitful.Units},
    T3<:Union{Nothing, Unitful.Units},
    T4<:Real
}
   # Validate input unit compatibility
   if rm.rB_unit isa Unitful.Units && !(retention isa AbstractQuantity)
       throw(ArgumentError("Input must have units when rm.rB_unit is a Unitful unit. " *
                          "Expected units compatible with $(rm.rB_unit)."))
   elseif rm.rB_unit isa Nothing && retention isa AbstractQuantity
       throw(ArgumentError("Input cannot have units when rm.rB_unit is Nothing. " *
                          "Expected unitless input."))
   end

   # Same parameter validation as derivinvmap
   if rm.rA_unit isa Unitful.Units && rA_unit isa Nothing
       throw(ArgumentError("rA_unit cannot be Nothing if rm.rA_unit is a Unitful unit."))
   elseif rm.rA_unit isa Nothing && rA_unit isa Unitful.Units
       throw(ArgumentError("rA_unit cannot be Unitful unit if rm.rA_unit is Nothing."))
   elseif rm.rB_unit isa Unitful.Units && rB_unit isa Nothing
       throw(ArgumentError("rB_unit cannot be Nothing if rm.rB_unit is a Unitful unit."))
   elseif rm.rB_unit isa Nothing && rB_unit isa Unitful.Units
       throw(ArgumentError("rB_unit cannot be Unitful unit if rm.rB_unit is a Nothing."))
   end
   
   result = _compute_inverse_derivative_mapping(rm, retention; 
        domain_boundary_threshold=domain_boundary_threshold, tol=tol, warn=warn)
   
   # Apply unit conversion then strip units
   if rA_unit isa Unitful.Units && rB_unit isa Unitful.Units
       result *= rm.rA_unit / rm.rB_unit
       return ustrip((rA_unit / rB_unit), result)
   elseif rA_unit isa Unitful.Units && rB_unit isa Nothing
       result *= rm.rA_unit
       return ustrip(rA_unit, result)
   elseif rA_unit isa Nothing && rB_unit isa Unitful.Units
       result /= rm.rB_unit
       return ustrip(inverse(rB_unit), result)
   else
       return result
   end
end

# ── rawderivmap ───────────────────────────────────────────────────────────────────────────

"""
    rawderivmap(
        rm::RetentionMapper,
        retention::AbstractQuantity{<:Real};
        domain_boundary_threshold::Real=1e-8,
        rA_unit::Union{Nothing, Unitful.Units}=retentionunit_A(rm), 
        rB_unit::Union{Nothing, Unitful.Units}=retentionunit_B(rm), 
        warn::Bool=false
    )

Compute the first derivative of the fitted mapping with respect to the input (domain A)
at a given input value, returning a unitless numeric result.

This function evaluates the derivative of the fitted monotonic B-spline, extrapolating
using boundary derivatives if the input lies outside the spline domain. The derivative is
adjusted to account for the normalization scales applied during spline fitting.

Inputs are a fitted `rm`, a unitful input value in domain A, optional `rA_unit`/`rB_unit`
for unit conversion, and `warn` to log extrapolation messages. `domain_boundary_threshold`
sets how far beyond the normalized spline domain a value must lie before logging an
extrapolation warning when `warn=true`.

Returns the derivative d(B)/d(A) as a unitless number in the `rB_unit/rA_unit` system
before stripping. This method requires unitful input (use [`derivmap`](@ref) for unitless
inputs). The derivative is always positive under the monotonic constraint, and
extrapolation uses constant boundary derivatives outside the fitted domain.

This function provides the same transformation as [`derivmap`](@ref) but strips units from 
the result, which is useful for downstream calculations that require raw numeric values.

See also 
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`fitmap`](@ref),
[`invmap`](@ref),
[`invmapmax`](@ref),
[`invmapmin`](@ref),
[`mapmax`](@ref),
[`mapmin`](@ref),
[`rawapplymap`](@ref),
[`rawderivinvmap`](@ref),
[`rawinvmap`](@ref),
[`rawinvmapmax`](@ref),
[`rawinvmapmin`](@ref),
[`rawmapmax`](@ref),
[`rawmapmin`](@ref),
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]
       mapper = fitmap(retention_times, retention_indices);

julia> rawderivmap(mapper, 5.0u"minute") ≈ 37.42885149584893
true

julia> test_times = [3.0, 11.0]u"minute";

julia> rawderivmap.(mapper, test_times) ≈ [67.10235298491122, 34.639458067769134]
true

julia> rawderivmap(mapper, 5.0u"minute", rA_unit=u"s") ≈ 0.6238141915974822
true

julia> rawderivmap(mapper, 0.8u"minute", warn=true) ≈ 79.01700550402535
[ Info: rA=0.8 is below the spline domain boundary (rA_min=1.2): extrapolating left.
true
```
"""
function rawderivmap(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T1=1e-8,
    rA_unit::T2=retentionunit_A(rm),
    rB_unit::T3=retentionunit_B(rm),
    warn::Bool=false
    ) where {T1<:Real,
        T2<:Union{Nothing, Unitful.Units},
        T3<:Union{Nothing, Unitful.Units},
    }
    # Validate input unit compatibility
    if rm.rA_unit isa Unitful.Units && !(retention isa AbstractQuantity)
        throw(ArgumentError("Input must have units when rm.rA_unit is a Unitful unit. " *
                           "Expected units compatible with $(rm.rA_unit)."))
    elseif rm.rA_unit isa Nothing && retention isa AbstractQuantity
        throw(ArgumentError("Input cannot have units when rm.rA_unit is Nothing. " *
                           "Expected unitless input."))
    end
    
    # Same validation as derivmap for unit parameters
    if rm.rA_unit isa Unitful.Units && rA_unit isa Nothing
        throw(ArgumentError("rA_unit cannot be Nothing if rm.rA_unit is a Unitful unit."))
    elseif rm.rA_unit isa Nothing && rA_unit isa Unitful.Units
        throw(ArgumentError("rA_unit cannot be Unitful unit if rm.rA_unit is Nothing."))
    elseif rm.rB_unit isa Unitful.Units && rB_unit isa Nothing
        throw(ArgumentError("rB_unit cannot be Nothing if rm.rB_unit is a Unitful unit."))
    elseif rm.rB_unit isa Nothing && rB_unit isa Unitful.Units
        throw(ArgumentError("rB_unit cannot be Unitful unit if rm.rB_unit is a Nothing."))
    end
    
    result = _compute_derivative_mapping(rm, 
        domain_boundary_threshold=domain_boundary_threshold, retention; warn=warn)
    
    # Apply unit conversion then strip units
    if rA_unit isa Unitful.Units && rB_unit isa Unitful.Units
        result *= rm.rB_unit / rm.rA_unit
        return ustrip((rB_unit / rA_unit), result)
    elseif rA_unit isa Unitful.Units && rB_unit isa Nothing
        result /= rm.rA_unit
        return ustrip(inverse(rA_unit), result)
    elseif rA_unit isa Nothing && rB_unit isa Unitful.Units
        result *= rm.rB_unit
        return ustrip(rB_unit, result)
    else
        return result
    end
end


# ── rawinvmap ─────────────────────────────────────────────────────────────────────────────

"""
    rawinvmap(
        rm::RetentionMapper, 
        retention::Union{<:Real, <:AbstractQuantity{<:Real}};
        domain_boundary_threshold::Real=1e-8, 
        unit::Union{Nothing, Unitful.Units}=retentionunit_A(rm), 
        tol::Real=1e-8,
        warn::Bool=false
    )

Invert the fitted mapping by computing the input value (domain A) corresponding to a given 
output value (domain B), returning a unitless numeric result.

Inputs are a fitted `rm`, a unitful output value in domain B, an optional `unit` for
conversion before stripping, `tol` for root-finding (absolute bisection tolerance when 
solving the inverse; smaller values increase accuracy at higher cost), and `warn` to log 
extrapolation messages.

The inversion normalizes the input to the spline domain, extrapolates linearly using the
boundary derivative, uses bisection within the domain, and returns a unitless result in
the requested unit system. This method requires unitful input (use [`invmap`](@ref) for
unitless inputs).

This function leverages the strict monotonicity of the fitted spline to provide reliable 
inverse mapping, enabling bidirectional conversion between the two domains while returning
raw numeric values useful for downstream calculations.

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
[`fitmap`](@ref),
[`invmap`](@ref),
[`invmapmax`](@ref),
[`invmapmin`](@ref),
[`mapmax`](@ref),
[`mapmin`](@ref),
[`rawapplymap`](@ref),
[`rawderivinvmap`](@ref),
[`rawderivmap`](@ref),
[`rawinvmapmax`](@ref),
[`rawinvmapmin`](@ref),
[`rawmapmax`](@ref),
[`rawmapmin`](@ref),
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]u"Th"
       mapper = fitmap(retention_times, retention_indices);

julia> rawinvmap(mapper, 350.0u"Th") ≈ 5.312425320906996
true

julia> rawinvmap(mapper, 350.0) ≈ 5.312425320906996
ERROR: ArgumentError: Input must have units when rm.rB_unit is a Unitful unit. Expected units compatible with Th.
[...]

julia> test_indices = [250.0, 450.0]u"Th";

julia> rawinvmap.(mapper, test_indices) ≈ [3.2283166885897936, 8.08686182716771]
true

julia> rawinvmap(mapper, 400.0u"Th", unit=u"s") ≈ 407.99982102982375
true

julia> rawinvmap(mapper, 50.0u"Th", warn=true) ≈ 0.5672238965934812
[ Info: rB=50.0 is below the inverse spline domain boundary (rB_min=100.000073): extrapolating left.
true

julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3, 12.1, 15.7]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]
       mapper = fitmap(retention_times, retention_indices);

julia> rawinvmap(mapper, 350.0) ≈ 5.312425320906996
true

julia> rawinvmap(mapper, 350.0u"Th") ≈ 5.312425320906996
ERROR: ArgumentError: Input cannot have units when rm.rB_unit is Nothing. Expected unitless input.
[...]
```
"""
function rawinvmap(
    rm::RetentionMapper,
    retention::Union{<:Real, <:AbstractQuantity{<:Real}};
    domain_boundary_threshold::T1=1e-8,
    unit::T2=retentionunit_A(rm),
    tol::T3=1e-8,
    warn::Bool=false
) where {
    T1<:Real,
    T2<:Union{Nothing, Unitful.Units},
    T3<:Real
}
    # Validate input unit compatibility (same as invmap)
    if rm.rB_unit isa Unitful.Units && !(retention isa AbstractQuantity)
        throw(ArgumentError("Input must have units when rm.rB_unit is a Unitful unit. " *
                           "Expected units compatible with $(rm.rB_unit)."))
    elseif rm.rB_unit isa Nothing && retention isa AbstractQuantity
        throw(ArgumentError("Input cannot have units when rm.rB_unit is Nothing. " *
                           "Expected unitless input."))
    end
    
    raw_result = _compute_inverse_mapping(rm, retention; 
        domain_boundary_threshold=domain_boundary_threshold, tol=tol, warn=warn)
    
    # Apply unit conversion then strip units
    if rm.rA_unit isa Unitful.Units
        return ustrip(unit, raw_result * rm.rA_unit)
    else
        return raw_result
    end
end

# ── rawinvmapmax ──────────────────────────────────────────────────────────────────────────

"""
    rawinvmapmax(rm::AbstractRetentionMapper; unit::Union{Nothing, Unitful.Units}=nothing)

Return the maximum predicted value of domain B with units stripped.

For unitful mappers, `unit` converts values before stripping. For unitless mappers,
`unit` must be `nothing`; otherwise `ArgumentError` is thrown with the message
"Cannot convert unitless retention B maximum to a unit".

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
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
[`rawinvmapmin`](@ref),
[`rawmapmax`](@ref),
[`rawmapmin`](@ref),
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> rawinvmapmax(mapper) ≈ 499.9999793028839
true

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> rawinvmapmax(mapper_unitless) ≈ 499.9999793028839
true

julia> rawinvmapmax(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A minimum to a unit
```
"""
@inline function rawinvmapmax(
    rm::AbstractRetentionMapper{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rB_pred_max, unit, "retention A minimum")
end

@inline function rawinvmapmax(
    rm::AbstractRetentionMapper{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip_scalar(rm.rB_pred_max, rm.rB_unit, unit)
end

# ── rawinvmapmin ──────────────────────────────────────────────────────────────────────────

"""
    rawinvmapmin(rm::AbstractRetentionMapper; unit::Union{Nothing, Unitful.Units}=nothing)

Return the minimum predicted value of domain B with units stripped.

For unitful mappers, `unit` converts values before stripping. For unitless mappers,
`unit` must be `nothing`; otherwise `ArgumentError` is thrown with the message
"Cannot convert unitless retention B minimum to a unit".

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
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
[`rawmapmax`](@ref),
[`rawmapmin`](@ref),
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> rawinvmapmin(mapper) ≈ 100.00001207241763
true

julia> (rawinvmapmin(mapper), rawinvmapmax(mapper)) .≈ (100.00001207241763, 499.9999793028839)
(true, true)

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> rawinvmapmin(mapper_unitless) ≈ 100.00001207241763
true

julia> rawinvmapmin(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A minimum to a unit
```
"""
@inline function rawinvmapmin(
    rm::AbstractRetentionMapper{<:Any, Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rB_pred_min, unit, "retention A minimum")
end

@inline function rawinvmapmin(
    rm::AbstractRetentionMapper{<:Any, <:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip_scalar(rm.rB_pred_min, rm.rB_unit, unit)
end

# ── rawmapmax ─────────────────────────────────────────────────────────────────────────────

"""
    rawmapmax(rm::AbstractRetentionMapper; unit::Union{Nothing, Unitful.Units}=nothing)

Return the maximum input value of domain A with units stripped.

For unitful mappers, `unit` converts values before stripping. For unitless mappers,
`unit` must be `nothing`; otherwise `ArgumentError` is thrown with the message
"Cannot convert unitless retention A maximum to a unit".

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
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
[`rawmapmin`](@ref),
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> rawmapmax(mapper) ≈ 9.3
true

julia> rawmapmax(mapper, unit=u"s") ≈ 558.0
true

julia> (rawmapmin(mapper), rawmapmax(mapper)) .≈ (1.2, 9.3)
(true, true)

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> rawmapmax(mapper_unitless) ≈ 9.3
true

julia> rawmapmax(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A minimum to a unit
```
"""
@inline function rawmapmax(
    rm::AbstractRetentionMapper{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rA_max, unit, "retention A minimum")
end

@inline function rawmapmax(
    rm::AbstractRetentionMapper{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip_scalar(rm.rA_max, rm.rA_unit, unit)
end

# ── rawmapmin ─────────────────────────────────────────────────────────────────────────────

"""
    rawmapmin(rm::AbstractRetentionMapper; unit::Union{Nothing, Unitful.Units}=nothing)

Return the minimum input value of domain A with units stripped.

For unitful mappers, `unit` converts values before stripping. For unitless mappers,
`unit` must be `nothing`; otherwise `ArgumentError` is thrown with the message
"Cannot convert unitless retention A minimum to a unit".

See also
[`AbstractRetentionMapper`](@ref JuChrom.AbstractRetentionMapper), 
[`RetentionMapper`](@ref JuChrom.RetentionMapper),
[`applymap`](@ref),
[`derivinvmap`](@ref),
[`derivmap`](@ref),
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
[`retentionunit_A`](@ref), 
[`retentionunit_B`](@ref).

# Examples
```jldoctest
julia> retention_times = [1.2, 2.5, 4.1, 6.8, 9.3]u"minute"
       retention_indices = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper = fitmap(retention_times, retention_indices);

julia> rawmapmin(mapper) ≈ 1.2
true

julia> rawmapmin(mapper, unit=u"s") ≈ 72.0
true

julia> inputs = [1.2, 2.5, 4.1, 6.8, 9.3]
       outputs = [100.0, 200.0, 300.0, 400.0, 500.0]
       mapper_unitless = fitmap(inputs, outputs);

julia> rawmapmin(mapper_unitless) ≈ 1.2
true

julia> rawmapmin(mapper_unitless, unit=u"s")
ERROR: ArgumentError: Cannot convert unitless retention A minimum to a unit
```
"""
@inline function rawmapmin(
    rm::AbstractRetentionMapper{Nothing};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitless(rm.rA_min, unit, "retention A minimum")
end

@inline function rawmapmin(
    rm::AbstractRetentionMapper{<:Unitful.Units};
    unit::Union{Nothing, Unitful.Units}=nothing)
    _handle_unitful_strip_scalar(rm.rA_min, rm.rA_unit, unit)
end
