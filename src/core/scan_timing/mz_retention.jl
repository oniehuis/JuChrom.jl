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

Return the retention coordinate at which a specific m/z (or m/z index) is sampled within a scan.

A scan spans a finite retention interval of width `scan_span`. The reported scan-level `retention`
typically refers to a particular reference point of that interval (start/middle/end). Within that interval,
m/z values are sampled sequentially, each over a dwell interval. This function maps the scan-level
`retention` to the retention coordinate associated with one specific m/z dwell interval, using:

1) how `retention` is referenced within the scan interval (`retention_ref`),
2) the dwell allocation model (`dwell`, `dwell_retention`/`scan_interval` or `dwell_retentions`),
3) the m/z acquisition order (`order`), and
4) which point within the dwell interval to return (`dwell_ref`).

`retention` and all dwell/interval inputs may be unitless (`Real`) or unitful (`Unitful.Quantity`), as long
as they are mutually compatible for addition/subtraction.

# Selecting the target m/z
Exactly one of the following identification modes must be used:

- **By index**: provide `mzindex`.
- **By value**: provide both `mzvalue` and `mzvalues`. The index is found via `findfirst(==(mzvalue), mzvalues)`.

The resolved index is interpreted as an index into `mzvalues` / `dwell_retentions` *as provided* (see `order`).

## Keywords
- `mzindex::Union{Int,Nothing}=nothing`:
  1-based index of the target m/z in the provided m/z list / dwell list.
  If `nothing`, it is resolved from `mzvalue` and `mzvalues`.

- `mzvalue::Union{Real,Unitful.Quantity,Nothing}=nothing` and
  `mzvalues::Union{AbstractVector{<:Union{Real,Unitful.Quantity}},Nothing}=nothing`:
  Use these together to locate the target by equality match. This is useful when you have an explicit m/z
  list and want the function to look up the index.

- `mzcount::Union{Int,Nothing}=nothing`:
  Total number of m/z values sampled in the scan (`N`). Required unless it can be inferred.
  Inference rules:
  - if `dwell == :heterogeneous` and `dwell_retentions` is given: `mzcount = length(dwell_retentions)`
  - else if `mzvalues` is given: `mzcount = length(mzvalues)`
  If neither is available, `mzcount` must be provided explicitly.

- `validate_span::Bool=true`:
  When `scan_interval` is provided for `dwell = :heterogeneous`, check that
  `sum(dwell_retentions) ≈ scan_interval` (relative tolerance `1e-6`).
  Set to `false` to skip this consistency check.

# Scan-level reference: `retention_ref`
`retention_ref` specifies what the input `retention` means with respect to the scan interval.

- `:start`: `retention` is the start of the scan interval.
- `:middle`: `retention` is the midpoint of the scan interval.
- `:end`: `retention` is the end of the scan interval.

Let `scan_span` be the total scan interval width (derived from dwell information). The scan start
coordinate used internally is then:

- `scan_start = retention`                                 if `retention_ref == :start`
- `scan_start = retention - scan_span/2`                   if `retention_ref == :middle`
- `scan_start = retention - scan_span`                     if `retention_ref == :end`

# Dwell allocation: `dwell`, `dwell_retention`, `scan_interval`, `dwell_retentions`
The scan span is defined as the sum of all dwell widths. Two allocation modes are supported:

## `dwell = :homogeneous`
All m/z values have the same dwell width.

You must provide either:
- `dwell_retention`: per-m/z dwell width, or
- `scan_interval`: total scan interval width. In this case `dwell_retention` is derived as
  `dwell_retention = scan_interval / mzcount`.

If both `dwell_retention` and `scan_interval` are provided, the function checks consistency:
`dwell_retention * mzcount ≈ scan_interval` within a relative tolerance of `1e-6`.

## `dwell = :heterogeneous`
Each m/z has its own dwell width.

- `dwell_retentions` must be provided and must have length `mzcount`.
- All elements of `dwell_retentions` must be strictly positive.
- If `scan_interval` is provided and `validate_span == true`, the function checks
  `sum(dwell_retentions) ≈ scan_interval` within a relative tolerance of `1e-6`.

# Acquisition order: `order`
`order` specifies the *acquisition order within the scan interval*:

- `:ascending`: the dwell sequence is used as provided (m/z sampled in increasing order).
- `:descending`: the dwell sequence is reversed and the provided index is mapped as
  `mzindex_acq = mzcount - mzindex + 1`.

Important interaction:
- `mzindex` (or the index resolved from `mzvalue`/`mzvalues`) is assumed to refer to the order
  in which `mzvalues` / `dwell_retentions` are provided. If the instrument acquired in descending order
  but your vectors are stored in ascending order (common), set `order = :descending` to obtain the correct
  within-scan position.

# Dwell reference: `dwell_ref`
`dwell_ref` selects which point of the target dwell interval to return:

- `:start`: return the start of the dwell interval.
- `:middle`: return the midpoint of the dwell interval.
- `:end`: return the end of the dwell interval.

If the target dwell width is `δ` and its start offset within the scan is `D`, then:
- `:start`  → offset `D`
- `:middle` → offset `D + δ/2`
- `:end`    → offset `D + δ`

# Return value
Returns:
\\[
scan\\_start + \\sum_{j<k} \\delta_j + \\phi(\\delta_k)
\\]
where `k` is the acquisition-order index, `δ_j` are dwell widths in acquisition order, and `φ` is the
offset chosen by `dwell_ref` (`0`, `δ/2`, or `δ`).

# Examples
```jldoctest
julia> # Homogeneous dwell times with mz target selected by index:

julia> r = mzretention(12.345u"minute";          # scan-level retention
                       mzindex=3,                # return acquisition time of 3rd ion
                       mzcount=10,               # 10 ions
                       dwell=:homogeneous,       # homogeneously dwelled
                       order=:descending,        # in descending m/z order
                       scan_interval=0.75u"ms",  # across 0.75 seconds
                       retention_ref=:start,     # scantime represents scan start
                       dwell_ref=:middle);       # return acquisition time at dwell mid point

julia> r ≈ 12.345009375u"minute"
true

julia> # Heterogeneous dwell times with mz target selected by value:

julia> mzvals = [43, 57, 71, 73, 77]u"Th"; 

julia> dw = [0.002, 0.002, 0.003, 0.004, 0.002]u"s";

julia> r₇₃ = mzretention(
        100.0u"s";
        mzvalue=73u"Th",
        mzvalues=mzvals,
        dwell=:heterogeneous,
        dwell_retentions=dw,
        retention_ref=:start,
        dwell_ref=:middle, 
        order=:descending);

julia> r₇₃ ≈ 100.004u"s"
true

julia> r₇₃ = mzretention.(
        [100.0, 110.0]u"s";
        mzvalue=73u"Th",
        mzvalues=mzvals,
        dwell=:heterogeneous,
        dwell_retentions=dw,
        retention_ref=:start,
        dwell_ref=:middle, 
        order=:descending);

julia> r₇₃ ≈ [100.004, 110.004]u"s"
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
    retention isa Unitful.Quantity ? Unitful.uconvert(Unitful.unit(retention), result) : result
end
