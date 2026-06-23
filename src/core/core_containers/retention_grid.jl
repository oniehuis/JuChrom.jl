"""
    RetentionGrid{T<:Real,U}

Regular retention-bin grid returned by [`densestgrid`](@ref JuChrom.densestgrid).

`binedges`, `binwidth`, `tolerance`, `overlap_min`, and `overlap_max` are stored as raw
numeric values in the scale indicated by `retentionunit`. `retentionunit` is either a
`Unitful.Units` value or `nothing` for an explicitly unitless retention axis. `extras`
stores optional unstructured metadata.

Use [`binedges`](@ref JuChrom.binedges), [`binwidth`](@ref JuChrom.binwidth),
[`overlapmin`](@ref JuChrom.overlapmin), and [`overlapmax`](@ref JuChrom.overlapmax)
to retrieve unit-attached values when the grid has a retention unit.
"""
struct RetentionGrid{T<:Real,U}
    binedges::Vector{T}
    binwidth::T
    retentionunit::U
    tolerance::T
    overlap_min::T
    overlap_max::T
    extras::Dict{String,Any}

    function RetentionGrid{T,U}(
        binedges::Vector{T},
        binwidth::T,
        retentionunit::U,
        tolerance::T,
        overlap_min::T,
        overlap_max::T,
        extras::Dict{String,Any},
    ) where {T<:Real,U<:Union{Nothing,Unitful.Units}}

        length(binedges) >= 2 || throw(ArgumentError("RetentionGrid needs at least two bin edges."))
        all(isfinite, binedges) || throw(ArgumentError("RetentionGrid bin edges must be finite."))
        isfinite(binwidth) || throw(ArgumentError("RetentionGrid binwidth must be finite."))
        isfinite(tolerance) || throw(ArgumentError("RetentionGrid tolerance must be finite."))
        isfinite(overlap_min) || throw(ArgumentError("RetentionGrid overlap_min must be finite."))
        isfinite(overlap_max) || throw(ArgumentError("RetentionGrid overlap_max must be finite."))
        binwidth > zero(T) || throw(ArgumentError("RetentionGrid binwidth must be positive."))
        tolerance >= zero(T) || throw(ArgumentError("RetentionGrid tolerance must be non-negative."))
        overlap_max > overlap_min ||
            throw(ArgumentError("RetentionGrid overlap_max must be greater than overlap_min."))
        all(diff(binedges) .> zero(T)) ||
            throw(ArgumentError("RetentionGrid bin edges must be strictly increasing."))
        first(binedges) <= overlap_min + tolerance ||
            throw(ArgumentError("RetentionGrid bin edges do not cover overlap_min."))
        last(binedges) + tolerance >= overlap_max ||
            throw(ArgumentError("RetentionGrid bin edges do not cover overlap_max."))

        new{T,U}(binedges, binwidth, retentionunit, tolerance, overlap_min, overlap_max, extras)
    end
end

function RetentionGrid(
    binedges::AbstractVector{<:Real},
    binwidth::Real,
    retentionunit::Union{Nothing,Unitful.Units},
    tolerance::Real,
    overlap_min::Real,
    overlap_max::Real;
    extras=Dict{String,Any}(),
)
    T = promote_type(
        eltype(binedges),
        typeof(binwidth),
        typeof(tolerance),
        typeof(overlap_min),
        typeof(overlap_max),
    )
    converted_extras = Dict{String,Any}(string(k) => v for (k, v) in extras)
    RetentionGrid{T,typeof(retentionunit)}(
        Vector{T}(binedges),
        T(binwidth),
        retentionunit,
        T(tolerance),
        T(overlap_min),
        T(overlap_max),
        converted_extras,
    )
end

Base.broadcastable(grid::RetentionGrid) = Base.RefValue(grid)

_retention_grid_unit_label(unit) = unit === nothing ? "unitless" : string(unit)
_retention_grid_value_label(value, unit) = unit === nothing ? string(value) :
    string(value * unit)
_retention_grid_extras_label(grid::RetentionGrid) =
    isempty(grid.extras) ? "none" :
    "$(length(grid.extras)) $(length(grid.extras) == 1 ? "entry" : "entries")"

function Base.show(io::IO, grid::RetentionGrid)
    unit = grid.retentionunit
    print(io, "RetentionGrid(")
    print(io, length(grid.binedges) - 1, " bins")
    print(io, ", edges=", _retention_grid_value_label(first(grid.binedges), unit))
    print(io, " to ", _retention_grid_value_label(last(grid.binedges), unit))
    print(io, ", width=", _retention_grid_value_label(grid.binwidth, unit))
    print(io, ", unit=", _retention_grid_unit_label(unit))
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", grid::RetentionGrid)
    unit = grid.retentionunit
    println(io, "RetentionGrid with $(length(grid.binedges) - 1) bins")
    println(io, "  Retention unit: $(_retention_grid_unit_label(unit))")
    println(io, "  Bin edges: $(_retention_grid_value_label(first(grid.binedges), unit)) to " *
        "$(_retention_grid_value_label(last(grid.binedges), unit)) ($(length(grid.binedges)) edges)")
    println(io, "  Bin width: $(_retention_grid_value_label(grid.binwidth, unit))")
    println(io, "  Tolerance: $(_retention_grid_value_label(grid.tolerance, unit))")
    println(io, "  Overlap: $(_retention_grid_value_label(grid.overlap_min, unit)) to " *
        "$(_retention_grid_value_label(grid.overlap_max, unit))")
    println(io, "  Raw edge type: $(eltype(grid.binedges))")
    println(io, "  Extras: $(_retention_grid_extras_label(grid))")
end
