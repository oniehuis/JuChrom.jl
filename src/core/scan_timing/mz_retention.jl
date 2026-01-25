"""
    mzretention(
        retention::Union{Real, Unitful.Quantity};
        mzindex::Union{Int, Nothing}=nothing,
        mzvalue::Union{Real, Unitful.Quantity, Nothing}=nothing,
        mzvalues::Union{AbstractVector{<:Union{Real, Unitful.Quantity}}, Nothing}=nothing,
        retention_ref::Symbol=:start,
        dwell_ref::Symbol=:middle,
        dwell::Symbol=:homogeneous,
        dwell_retention::Union{Real, Unitful.Quantity, Nothing}=nothing,
        scan_interval::Union{Real, Unitful.Quantity, Nothing}=nothing,
        dwell_retentions::Union{AbstractVector{<:Union{Real, Unitful.Quantity}}, Nothing}=nothing,
        mzcount::Union{Int, Nothing}=nothing,
        order::Symbol=:ascending,
        validate_span::Bool=true
    )

Return the retention coordinate at which a specific m/z (or m/z index) is sampled within a 
scan.

A scan spans a finite retention interval of width `scan_span`. The reported scan-level 
`retention` typically refers to a particular reference point of that interval 
(start/middle/end). Within that interval, m/z values are sampled sequentially, each over a 
dwell interval. This function maps the scan-level `retention` to the retention coordinate 
associated with one specific m/z dwell interval, using how `retention` is referenced within 
the scan interval (`retention_ref`), the dwell allocation model (`dwell`, 
`dwell_retention`/`scan_interval` or `dwell_retentions`), the m/z acquisition order 
(`order`), and which point within the dwell interval to return (`dwell_ref`).

`retention` and all dwell/interval inputs may be unitless (`Real`) or unitful 
(`Unitful.Quantity`), as long as they are mutually compatible for addition and subtraction. 
The target ion is specified either by index (`mzindex`) or by value (`mzvalue` together with 
`mzvalues`), in which case the index is found via `findfirst(==(mzvalue), mzvalues)`. The 
resolved index refers to the order in which `mzvalues` or `dwell_retentions` are provided; 
when the instrument acquires in descending m/z order but the vectors are stored in ascending 
order, set `order = :descending` to obtain the correct within-scan position.

Keyword arguments are interpreted as follows. `mzindex` is a 1-based index of the target m/z 
in the provided list and is inferred from `mzvalue` and `mzvalues` when omitted. `mzcount` 
is the total number of m/z values sampled in the scan and is inferred from 
`dwell_retentions` for heterogeneous dwell or from `mzvalues` otherwise; if neither is 
available, it must be provided explicitly. For homogeneous dwell, you must supply 
`dwell_retention` or `scan_interval` (from which `dwell_retention` is derived), and 
consistency is checked when both are provided. For heterogeneous dwell, `dwell_retentions` 
must be provided, must have length `mzcount`, and must be strictly positive; if 
`scan_interval` is provided and `validate_span` is `true`, the function checks that 
`sum(dwell_retentions) ≈ scan_interval` within a relative tolerance of `1e-6`. 
The scan-level reference is controlled by `retention_ref` (`:start`, `:middle`, or `:end`), 
and the dwell reference is controlled by `dwell_ref` (`:start`, `:middle`, or `:end`), 
corresponding to the start, midpoint, or end of the target dwell interval.

The scan start is computed from the scan-level reference as `scan_start=retention` for
`retention_ref==:start`, `scan_start=retention - scan_span / 2` for 
`retention_ref==:middle`, and `scan_start=retention - scan_span` for `retention_ref==:end`.

The returned value is
```math
scan\\_start + \\sum_{j<k} \\delta_j + \\phi(\\delta_k)
```
where `k` is the acquisition-order index, `δ_j` are dwell widths in acquisition order, 
and `φ` is the offset chosen by `dwell_ref` (`0`, `δ/2`, or `δ`).

# Examples
```jldoctest
julia> # Homogeneous dwell times with mz target selected by index:

julia> r₁ = mzretention(12.345u"minute";       # scan-level retention
                       mzindex=3,              # return acquisition time of 3rd ion
                       mzcount=10,             # 10 ions
                       dwell=:homogeneous,     # homogeneously dwelled
                       order=:descending,      # in descending m/z order
                       scan_interval=50u"ms",  # across 50 ms
                       retention_ref=:start,   # scantime represents scan start
                       dwell_ref=:middle);     # return acquisition time at dwell midpoint

julia> r₁ ≈ 12.345625u"minute"
true

julia> # Heterogeneous dwell times with mz target selected by value:

julia> mzvals = [43, 57, 71, 73, 77]u"Th"; 

julia> dw = [0.002, 0.002, 0.003, 0.004, 0.002]u"s";

julia> r₂ = mzretention(100.0u"s";             # scan-level retention
                        mzvalue=73u"Th",       # return acquisition time of m/z 73 Thomson
                        mzvalues=mzvals,       # vector with 5 ions
                        dwell=:heterogeneous,  # heterogeneously dwelled
                        order=:descending,     # in descending m/z order
                        dwell_retentions=dw,   # dwell times per ion
                        retention_ref=:start,  # scantime represents scan start
                        dwell_ref=:middle);    # return acquisition time at dwell midpoint
        

julia> r₂ ≈ 100.004u"s"
true

julia> r₃ = mzretention.([100.0, 110.0]u"s";  # dot form broadcasts over the retentions
                         mzvalue=73u"Th",
                         mzvalues=mzvals,
                         dwell=:heterogeneous,
                         order=:descending,
                         dwell_retentions=dw,
                         retention_ref=:start,
                         dwell_ref=:middle);

julia> r₃ ≈ [100.004, 110.004]u"s"
true
```
"""
function mzretention(
    retention::Union{Real, Unitful.Quantity};
    mzindex::T1=nothing,
    mzvalue::T2=nothing,
    mzvalues::T3=nothing,
    retention_ref::Symbol=:start,
    dwell_ref::Symbol=:middle,
    dwell::Symbol=:homogeneous,
    dwell_retention::T4=nothing,
    scan_interval::T5=nothing,
    dwell_retentions::T6=nothing,
    mzcount::T7=nothing,
    order::Symbol=:ascending,
    validate_span::Bool=true
    ) where {
        T1<:Union{Integer, Nothing},
        T2<:Union{Real, Unitful.Quantity, Nothing},
        T3<:Union{AbstractVector{<:Union{Real, Unitful.Quantity}}, Nothing},
        T4<:Union{Real, Unitful.Quantity, Nothing},
        T5<:Union{Real, Unitful.Quantity, Nothing},
        T6<:Union{AbstractVector{<:Union{Real, Unitful.Quantity}}, Nothing},
        T7<:Union{Integer, Nothing}
    }

    retention_ref ∈ (:start, :middle, :end) || throw(ArgumentError(
        "retention_ref must be :start, :middle, or :end"))
    dwell_ref ∈ (:start, :middle, :end) || throw(ArgumentError(
        "dwell_ref must be :start, :middle, or :end"))
    dwell ∈ (:homogeneous, :heterogeneous) || throw(ArgumentError(
        "dwell must be :homogeneous or :heterogeneous"))
    order ∈ (:ascending, :descending) || throw(ArgumentError(
        "order must be :ascending or :descending"))

    # resolve imzvalue index
    if isnothing(mzindex)
        (!isnothing(mzvalue) && !isnothing(mzvalues)) || throw(ArgumentError(
            "Provide mzindex, or provide mzvalue and mzvalues"))
        idx = findfirst(==(mzvalue), mzvalues)
        isnothing(idx) && throw(ArgumentError("mzvalue not found in mzvalues"))
        mzindex = idx
    end
    mzindex ≥ 1 || throw(ArgumentError("mzindex must be ≥ 1"))

    # determine mzcount
    if isnothing(mzcount)
        if dwell == :heterogeneous && !isnothing(dwell_retentions)
            mzcount = length(dwell_retentions)
        elseif !isnothing(mzvalues)
            mzcount = length(mzvalues)
        else
            throw(ArgumentError("mzcount is required unless it can be inferred"))
        end
    end
    mzindex ≤ mzcount || throw(ArgumentError("mzindex exceeds mzcount"))

    # build dwell vector
    local dwells::Vector{Union{Real, Unitful.Quantity}}

    if dwell == :homogeneous
        if isnothing(dwell_retention)
            isnothing(scan_interval) && throw(ArgumentError(
                "Provide dwell_retention or scan_interval"))
            scan_interval > zero(scan_interval) || throw(ArgumentError(
                "scan_interval must be > 0"))
            dwell_retention = scan_interval / mzcount
        else
            dwell_retention > zero(dwell_retention) || throw(ArgumentError(
                "dwell_retention must be > 0"))
            if !isnothing(scan_interval)
                scan_interval > zero(scan_interval) || throw(ArgumentError(
                    "scan_interval must be > 0"))
                implied = dwell_retention * mzcount
                tol = abs(scan_interval) * 1e-6
                abs(implied - scan_interval) ≤ tol || throw(ArgumentError(
                    "dwell_retention*mzcount inconsistent with scan_interval"))
            end
        end
        dwells = fill(dwell_retention, mzcount)
    else
         isnothing(dwell_retentions) && throw(ArgumentError(
            "dwell_retentions is required"))
        length(dwell_retentions) == mzcount || throw(ArgumentError(
            "length(dwell_retentions) must equal mzcount"))
        all(d -> d > zero(d), dwell_retentions) || throw(ArgumentError(
            "all dwell_retentions must be > 0"))
        dwells = collect(dwell_retentions)
        if !isnothing(scan_interval) && validate_span
            scan_interval > zero(scan_interval) || throw(ArgumentError(
                "scan_interval must be > 0"))
            implied = sum(dwells)
            tol = abs(scan_interval) * 1e-6
            abs(implied - scan_interval) ≤ tol || throw(ArgumentError(
                "sum(dwell_retentions) inconsistent with scan_interval"))
        end
    end

    # acquisition order
    if order == :descending
        dwells = reverse(dwells)
        mzindex = mzcount - mzindex + 1
    end

    # scan start from retention reference
    scan_span = sum(dwells)
    scan_start =
        retention_ref == :start  ? retention :
        retention_ref == :middle ? (retention - scan_span / 2) :
                                   (retention - scan_span)

    # offset within scan
    z = zero(scan_span)
    offset_before = mzindex == 1 ? z : sum(@view dwells[1:mzindex-1])
    this_dwell = dwells[mzindex]

    dwell_offset =
        dwell_ref == :start  ? z :
        dwell_ref == :middle ? this_dwell/2 :
                               this_dwell

    result = scan_start + offset_before + dwell_offset
    retention isa Unitful.Quantity ? 
        Unitful.uconvert(Unitful.unit(retention), result) : result
end
