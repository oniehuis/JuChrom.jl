"""
    extras(grid::RetentionGrid) -> Dict{String, Any}

Return the unstructured metadata dictionary attached to a retention grid.
"""
extras(grid::RetentionGrid) = grid.extras

"""
    retentionunit(grid::RetentionGrid) -> Union{Unitful.Units, Nothing}

Return the unit associated with the retention-grid coordinates, or `nothing` for an
explicitly unitless grid.
"""
retentionunit(grid::RetentionGrid) = grid.retentionunit

"""
    rawbinedges(grid::RetentionGrid) -> Vector{<:Real}

Return the raw numeric bin edges in the grid's stored retention scale.
"""
rawbinedges(grid::RetentionGrid) = grid.binedges

"""
    binedges(grid::RetentionGrid; unit=retentionunit(grid))

Return bin edges, attaching or converting retention units when the grid has a unit.
"""
@inline function binedges(
    grid::RetentionGrid;
    unit::Union{Nothing,Unitful.Units}=retentionunit(grid))

    runit = retentionunit(grid)
    isnothing(runit) && return _handle_unitless(grid.binedges, unit, "bin edges")
    _handle_unitful_convert(grid.binedges, runit, unit)
end

"""
    rawbinwidth(grid::RetentionGrid) -> Real

Return the raw numeric bin width in the grid's stored retention scale.
"""
rawbinwidth(grid::RetentionGrid) = grid.binwidth

"""
    binwidth(grid::RetentionGrid; unit=retentionunit(grid))

Return the bin width, attaching or converting retention units when the grid has a unit.
"""
@inline function binwidth(
    grid::RetentionGrid;
    unit::Union{Nothing,Unitful.Units}=retentionunit(grid))

    runit = retentionunit(grid)
    isnothing(runit) && return _handle_unitless(grid.binwidth, unit, "bin width")
    _handle_unitful_convert_scalar(grid.binwidth, runit, unit)
end

"""
    rawtolerance(grid::RetentionGrid) -> Real

Return the raw numeric right-edge tolerance in the grid's stored retention scale.
"""
rawtolerance(grid::RetentionGrid) = grid.tolerance

"""
    tolerance(grid::RetentionGrid; unit=retentionunit(grid))

Return the right-edge tolerance, attaching or converting retention units when the grid has
a unit.
"""
@inline function tolerance(
    grid::RetentionGrid;
    unit::Union{Nothing,Unitful.Units}=retentionunit(grid))

    runit = retentionunit(grid)
    isnothing(runit) && return _handle_unitless(grid.tolerance, unit, "tolerance")
    _handle_unitful_convert_scalar(grid.tolerance, runit, unit)
end

"""
    rawoverlapmin(grid::RetentionGrid) -> Real

Return the raw numeric lower bound of the shared retention overlap used to build the grid.
"""
rawoverlapmin(grid::RetentionGrid) = grid.overlap_min

"""
    rawoverlapmax(grid::RetentionGrid) -> Real

Return the raw numeric upper bound of the shared retention overlap used to build the grid.
"""
rawoverlapmax(grid::RetentionGrid) = grid.overlap_max

"""
    overlapmin(grid::RetentionGrid; unit=retentionunit(grid))

Return the lower bound of the shared retention overlap, attaching or converting retention
units when the grid has a unit.
"""
@inline function overlapmin(
    grid::RetentionGrid;
    unit::Union{Nothing,Unitful.Units}=retentionunit(grid))

    runit = retentionunit(grid)
    isnothing(runit) && return _handle_unitless(grid.overlap_min, unit, "overlap minimum")
    _handle_unitful_convert_scalar(grid.overlap_min, runit, unit)
end

"""
    overlapmax(grid::RetentionGrid; unit=retentionunit(grid))

Return the upper bound of the shared retention overlap, attaching or converting retention
units when the grid has a unit.
"""
@inline function overlapmax(
    grid::RetentionGrid;
    unit::Union{Nothing,Unitful.Units}=retentionunit(grid))

    runit = retentionunit(grid)
    isnothing(runit) && return _handle_unitless(grid.overlap_max, unit, "overlap maximum")
    _handle_unitful_convert_scalar(grid.overlap_max, runit, unit)
end

function Base.length(::RetentionGrid)
    Int(2)
end

Base.iterate(grid::RetentionGrid) = (binedges(grid), Val(:binwidth))
Base.iterate(grid::RetentionGrid, ::Val{:binwidth}) = (binwidth(grid), Val(:done))
Base.iterate(::RetentionGrid, ::Val{:done}) = nothing
