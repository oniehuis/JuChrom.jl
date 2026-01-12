
# ── subtract ──────────────────────────────────────────────────────────────────────────────

# DIE IMPLEMENTIERUNG AKZEPTIERT KEINE UNITS. IST DAS BEABSICHTIGT? ÜBERHAUPT IST DIE SUBTRACTION NUR
# FÜR FÄLLTE WIE BASELINE MÖGLICH. BEABSICHTIGT?

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
