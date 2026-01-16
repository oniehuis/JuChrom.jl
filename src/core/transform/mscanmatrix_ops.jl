
# ── subtract ──────────────────────────────────────────────────────────────────────────────

# DIE IMPLEMENTIERUNG AKZEPTIERT KEINE UNITS. IST DAS BEABSICHTIGT? ÜBERHAUPT IST DIE SUBTRACTION NUR
# FÜR FÄLLTE WIE BASELINE MÖGLICH. BEABSICHTIGT?

# function retentiontrim(
#     msmatrix::MassScanMatrix{<:Any, <:Unitful.Units}; 
#     start::T1=first(retentions(msmatrix)), 
#     stop::T2=last(retentions(msmatrix))
#     ) where {T1<:Unitful.AbstractQuantity, T2<:Unitful.AbstractQuantity}

#     _retentiontrim(msmatrix::MassScanMatrix, start, stop)
# end

# function retentiontrim(
#     msmatrix::MassScanMatrix{<:Any, <:Nothing}; 
#     start::T1=first(retentions(msmatrix)), 
#     stop::T2=last(retentions(msmatrix))
#     ) where {T1<:Real, T2<:Real}

#     _retentiontrim(msmatrix::MassScanMatrix, start, stop)
# end

# function _retentiontrim(msmatrix::MassScanMatrix, start, stop)

#         # Validate that the time range is in correct order
#     start ≤ stop || throw(ArgumentError("`start` must be less than or equal to `stop`."))

#     # Filter scans to include only those within the [start, stop] range
#     i_start = searchsortedfirst(retentions(msmatrix), start)
#     i_stop = searchsortedlast(retentions(msmatrix), stop)
#     indices = i_start:i_stop

#     # Ensure that trimming did not remove all scans
#     isempty(indices) && throw(ArgumentError("No scans left after trimming."))

#     # Create and return a new MassScanMatrix with the filtered scans 
#     MassScanMatrix(
#         rawretentions(msmatrix)[indices],
#         retentionunit(msmatrix),
#         mzvalues(msmatrix),
#         mzunit(msmatrix),
#         rawintensities(msmatrix)[indices, :],
#         intensityunit(msmatrix),
#         level=level(msmatrix),
#         instrument=deepcopy(instrument(msmatrix)),
#         acquisition=deepcopy(acquisition(msmatrix)),
#         user=deepcopy(user(msmatrix)),
#         sample=deepcopy(sample(msmatrix)),
#         extras=deepcopy(extras(msmatrix))
#     )
# end

function subtract(msmatrix_A::MassScanMatrix{<:Any, <:Nothing}, msmatrix_B::MassScanMatrix{<:Any, <:Nothing}, floor::Real)
    retentions(msmatrix_A) == retentions(msmatrix_B) ||
        throw(ArgumentError("MassScanMatrices must have identical retentions."))
    mzvalues(msmatrix_A) == mzvalues(msmatrix_B) ||
        throw(ArgumentError("MassScanMatrices must have identical m/z values."))
    size(rawintensities(msmatrix_A)) == size(rawintensities(msmatrix_B)) ||
        throw(DimensionMismatch("MassScanMatrices must have matching intensities dimensions."))
    level(msmatrix_A) == level(msmatrix_B) ||
        throw(ArgumentError("MassScanMatrices must have matching MS levels."))
    instrument(msmatrix_A) == instrument(msmatrix_B) ||
        throw(ArgumentError("MassScanMatrices must have matching instrument data."))
    acquisition(msmatrix_A) == acquisition(msmatrix_B) ||
        throw(ArgumentError("MassScanMatrices must have matching acquisition data."))
    user(msmatrix_A) == user(msmatrix_B) ||
        throw(ArgumentError("MassScanMatrices must have matching user data."))
    sample(msmatrix_A) == sample(msmatrix_B) ||
        throw(ArgumentError("MassScanMatrices must have matching sample data."))
    extras(msmatrix_A) == extras(msmatrix_B) ||
        throw(ArgumentError("MassScanMatrices must have matching extras."))

    resulting_intensities = intensities(msmatrix_A) .- intensities(msmatrix_B)
    resulting_intensities[resulting_intensities .< floor] .= floor

    MassScanMatrix(
            copy(rawretentions(msmatrix_A)),
            retentionunit(msmatrix_A),
            copy(rawmzvalues(msmatrix_A)),
            mzunit(msmatrix_A),
            resulting_intensities,
            intensityunit(msmatrix_A),
            level=level(msmatrix_A),
            instrument=deepcopy(instrument(msmatrix_A)),
            acquisition=deepcopy(acquisition(msmatrix_A)),
            user=deepcopy(user(msmatrix_A)),
            sample=deepcopy(sample(msmatrix_A)),
            extras=deepcopy(extras(msmatrix_A))
        )
end
