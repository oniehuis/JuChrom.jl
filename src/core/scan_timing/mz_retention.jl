"""
    mzretention(
        retention::Union{Real, Unitful.Quantity};
        mzindex::Union{Int, Nothing}=nothing,
        mzvalue::Union{Real, Unitful.Quantity, Nothing}=nothing,
        mzvalues::Union{AbstractVector{<:Union{Real, Unitful.Quantity}}, Nothing}=nothing,
        retentionref::Symbol=:start,
        dwellref::Symbol=:middle,
        dwell::Symbol=:homogeneous,
        dwellretention::Union{Real, Unitful.Quantity, Nothing}=nothing,
        scaninterval::Union{Real, Unitful.Quantity, Nothing}=nothing,
        dwellretentions::Union{AbstractVector{<:Union{Real, Unitful.Quantity}}, Nothing}=nothing,
        mzcount::Union{Int, Nothing}=nothing,
        order::Symbol=:ascending,
        validate_span::Bool=true,
        rtol::Real=1e-6,
        atol::Union{Nothing, Real, Unitful.Quantity}=nothing
    )

Return the retention coordinate at which a specific m/z (or m/z index) is sampled within a 
scan.

A scan spans a finite retention interval of width `scan_span`. The reported scan-level 
`retention` typically refers to a particular reference point of that interval 
(start/middle/end). Within that interval, m/z values are sampled sequentially, each over a 
dwell interval. This function maps the scan-level `retention` to the retention coordinate 
associated with one specific m/z dwell interval, using how `retention` is referenced within 
the scan interval (`retentionref`), the dwell allocation model (`dwell`, 
`dwellretention`/`scaninterval` or `dwellretentions`), the m/z acquisition order 
(`order`), and which point within the dwell interval to return (`dwellref`).

`retention` and all dwell/interval inputs may be unitless (`Real`) or unitful 
(`Unitful.Quantity`), as long as they are mutually compatible for addition and subtraction. 
The target ion is specified either by index (`mzindex`) or by value (`mzvalue` together with 
`mzvalues`), in which case the index is found via `findfirst(==(mzvalue), mzvalues)`. The 
resolved index refers to the order in which `mzvalues` or `dwellretentions` are provided; 
when the instrument acquires in descending m/z order but the vectors are stored in ascending 
order, set `order = :descending` to obtain the correct within-scan position.

Keyword arguments are interpreted as follows. `mzindex` is a 1-based index of the target m/z 
in the provided list and is inferred from `mzvalue` and `mzvalues` when omitted. `mzcount` 
is the total number of m/z values sampled in the scan and is inferred from 
`dwellretentions` for heterogeneous dwell or from `mzvalues` otherwise; if neither is 
available, it must be provided explicitly. For homogeneous dwell, you must supply 
`dwellretention` or `scaninterval` (from which `dwellretention` is derived), and 
consistency is checked when both are provided. For heterogeneous dwell, `dwellretentions` 
must be provided, must have length `mzcount`, and must be strictly positive; if 
`scaninterval` is provided and `validate_span` is `true`, the function checks that 
`sum(dwellretentions) ≈ scaninterval` using `rtol` and `atol`.
The scan-level reference is controlled by `retentionref` (`:start`, `:middle`, or `:end`), 
and the dwell reference is controlled by `dwellref` (`:start`, `:middle`, or `:end`), 
corresponding to the start, midpoint, or end of the target dwell interval.

The scan start is computed from the scan-level reference as `scan_start=retention` for
`retentionref==:start`, `scan_start=retention - scan_span / 2` for 
`retentionref==:middle`, and `scan_start=retention - scan_span` for `retentionref==:end`.

The returned value is
```math
scan\\_start + \\sum_{j<k} \\delta_j + \\phi(\\delta_k)
```
where `k` is the acquisition-order index, `δ_j` are dwell widths in acquisition order, 
and `φ` is the offset chosen by `dwellref` (`0`, `δ/2`, or `δ`).

# Examples
```jldoctest
julia> # Homogeneous dwell times with mz target selected by index:

julia> r₁ = mzretention(12.345u"minute";      # scan-level retention
                       retentionref=:start,   # scantime represents scan start
                       scaninterval=50u"ms",  # whole scan took 50 ms
                       mzindex=3,             # return acquisition time of 3rd ion
                       mzcount=10,            # total of 10 ions were scanned
                       order=:descending,     # scan was in descending m/z order
                       dwellref=:middle,      # return acquisition time at dwell midpoint
                       dwell=:homogeneous);   # ions were homogeneously dwelled;     

julia> r₁ ≈ 12.345625u"minute"
true

julia> # Heterogeneous dwell times with mz target selected by value:

julia> mzvals = [43, 57, 71, 73, 77]u"Th"; 

julia> dw = [0.002, 0.002, 0.003, 0.004, 0.002]u"s";

julia> r₂ = mzretention(100.0u"s";              # scan-level retention
                        retentionref=:start,    # scantime represents scan start
                        mzvalue=73u"Th",        # return acquisition time of m/z 73 Thomson
                        mzvalues=mzvals,        # vector with 5 ions
                        order=:descending,      # ions were scanned in descending m/z order
                        dwellretentions=dw,     # dwell times per ion
                        dwellref=:middle,       # return acquisition time at dwell midpoint
                        dwell=:heterogeneous);  # ions were heterogeneously dwelled
        

julia> r₂ ≈ 100.004u"s"
true

julia> r₃ = mzretention.([100.0, 110.0]u"s";  # dot form broadcasts over the retentions
                         retentionref=:start,
                         mzvalue=73u"Th",
                         mzvalues=mzvals,
                         order=:descending,
                         dwellretentions=dw,
                         dwellref=:middle,
                         dwell=:heterogeneous);

julia> r₃ ≈ [100.004, 110.004]u"s"
true
```
"""
function mzretention(
    retention::Union{Real, Unitful.Quantity};
    retentionref::Symbol=:start,
    scaninterval::T1=nothing,
    mzindex::T2=nothing,
    mzvalue::T3=nothing,
    mzvalues::T4=nothing,
    mzcount::T5=nothing,
    order::Symbol=:ascending,
    dwellretention::T6=nothing,
    dwellretentions::T7=nothing,
    dwellref::Symbol=:middle,
    dwell::Symbol=:homogeneous,
    validate_span::Bool=true,
    atol::T8=nothing,
    rtol::T9=1e-6
    ) where {
        T1<:Union{Real, Unitful.Quantity, Nothing},
        T2<:Union{Integer, Nothing},
        T3<:Union{Real, Unitful.Quantity, Nothing},
        T4<:Union{AbstractVector{<:Union{Real, Unitful.Quantity}}, Nothing},
        T5<:Union{Integer, Nothing},
        T6<:Union{Real, Unitful.Quantity, Nothing},
        T7<:Union{AbstractVector{<:Union{Real, Unitful.Quantity}}, Nothing},
        T8<:Union{Nothing, Real, Unitful.Quantity},
        T9<:Real
    }

    retentionref ∈ (:start, :middle, :end) || throw(ArgumentError(
        "retentionref must be :start, :middle, or :end"))
    dwellref ∈ (:start, :middle, :end) || throw(ArgumentError(
        "dwellref must be :start, :middle, or :end"))
    dwell ∈ (:homogeneous, :heterogeneous) || throw(ArgumentError(
        "dwell must be :homogeneous or :heterogeneous"))
    order ∈ (:ascending, :descending) || throw(ArgumentError(
        "order must be :ascending or :descending"))

    # Resolve imzvalue index
    if isnothing(mzindex)
        (!isnothing(mzvalue) && !isnothing(mzvalues)) || throw(ArgumentError(
            "Provide mzindex, or provide mzvalue and mzvalues"))
        idx = findfirst(==(mzvalue), mzvalues)
        isnothing(idx) && throw(ArgumentError("mzvalue not found in mzvalues"))
        mzindex = idx
    end
    mzindex ≥ 1 || throw(ArgumentError("mzindex must be ≥ 1"))

    # Determine mzcount
    if isnothing(mzcount)
        if dwell == :heterogeneous && !isnothing(dwellretentions)
            mzcount = length(dwellretentions)
        elseif !isnothing(mzvalues)
            mzcount = length(mzvalues)
        else
            throw(ArgumentError("mzcount is required unless it can be inferred"))
        end
    end
    mzindex ≤ mzcount || throw(ArgumentError("mzindex exceeds mzcount"))

    if dwell == :homogeneous
        if isnothing(dwellretention)
            isnothing(scaninterval) && throw(ArgumentError(
                "Provide dwellretention or scaninterval"))
            scaninterval > zero(scaninterval) || throw(ArgumentError(
                "scaninterval must be > 0"))
            dwellretention = scaninterval / mzcount
        else
            dwellretention > zero(dwellretention) || throw(ArgumentError(
                "dwellretention must be > 0"))
            if !isnothing(scaninterval)
                scaninterval > zero(scaninterval) || throw(ArgumentError(
                    "scaninterval must be > 0"))
                implied = dwellretention * mzcount
                atol_eff = isnothing(atol) ? zero(scaninterval) : atol
                isapprox(implied, scaninterval; rtol=rtol, atol=atol_eff) || throw(
                    ArgumentError("dwellretention*mzcount inconsistent with scaninterval"))
            end
        end
        dwells = fill(dwellretention, mzcount)
    else
         isnothing(dwellretentions) && throw(ArgumentError(
            "dwellretentions is required"))
        length(dwellretentions) == mzcount || throw(ArgumentError(
            "length(dwellretentions) must equal mzcount"))
        all(d -> d > zero(d), dwellretentions) || throw(ArgumentError(
            "all dwellretentions must be > 0"))
        dwells = collect(dwellretentions)
        if !isnothing(scaninterval) && validate_span
            scaninterval > zero(scaninterval) || throw(ArgumentError(
                "scaninterval must be > 0"))
            implied = sum(dwells)
            atol_eff = isnothing(atol) ? zero(scaninterval) : atol
            isapprox(implied, scaninterval; rtol=rtol, atol=atol_eff) || throw(ArgumentError(
                "sum(dwellretentions) inconsistent with scaninterval"))
        end
    end

    # Acquisition order
    if order == :descending
        dwells = reverse(dwells)
        mzindex = mzcount - mzindex + 1
    end

    # Scan start from retention reference
    scan_span = sum(dwells)
    scan_start =
        retentionref == :start  ? retention :
        retentionref == :middle ? (retention - scan_span / 2) :
                                  (retention - scan_span)

    # Offset within scan
    z = zero(scan_span)
    offset_before = mzindex == 1 ? z : sum(@view dwells[1:mzindex-1])
    this_dwell = dwells[mzindex]

    dwell_offset =
        dwellref == :start  ? z :
        dwellref == :middle ? this_dwell/2 :
                              this_dwell

    result = scan_start + offset_before + dwell_offset
    retention isa Unitful.Quantity ? 
        Unitful.uconvert(Unitful.unit(retention), result) : result
end
