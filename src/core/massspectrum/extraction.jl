"""
    extractmassspectrum(series::AbstractMassScanSeries, scanindices;
                        subtract=Int[]) -> MassSpectrum
    extractmassspectrum(matrix::AbstractMassScanMatrix, scanindices;
                        subtract=Int[]) -> MassSpectrum

Average the supplied scans and subtract the average of the scans in `subtract`.
Match m/z values exactly, treating absent entries as zero, and divide each group's
sum by its total number of scans, including scans without that m/z. Floor the final
intensities at zero. An empty group contributes zero; at least one group must
contain a scan. Each index contributes at most once to each group. An index in
both groups contributes to both averages with their respective group weights.

Binning is an explicit preprocessing step: use [`binmzvalues`](@ref) when nearby
m/z values should represent the same channel. No rounding, interpolation, or
base-peak normalization is performed here. Compatible units are converted to those of the
first contributing scan. Accumulation uses at least Float64 precision.

The result contains the sorted union of contributing m/z values, including zero
channels, preserves their numeric type, and records the contributing scan indices
in its attributes. Matrices use their shared m/z axis and row indices as scan
indices. Variance-bearing matrices are accepted; this operation extracts their
intensities without propagating variances into the `MassSpectrum`. The input is unchanged.

# Examples
```julia
prepared = binmzvalues(series)
result = extractmassspectrum(prepared, [10, 12, 13]; subtract=[5, 6])
```
"""
function extractmassspectrum(series::AbstractMassScanSeries,
    scanindices::AbstractVector{<:Integer}; subtract::AbstractVector{<:Integer}=Int[])
    positive, negative, allindices = extraction_scan_groups(series, scanindices, subtract)

    reference = scan(series, first(allindices))
    mz_unit = mzunit(reference)
    intensity_unit = intensityunit(reference)
    mz_type = eltype(rawmzvalues(reference; unit=mz_unit))
    intensity_type = promote_type(Float64, eltype(rawintensities(reference; unit=intensity_unit)))
    totals = Dict{mz_type, Tuple{intensity_type, intensity_type}}()
    zero_intensity = zero(intensity_type)
    positive_set, negative_set = Set(positive), Set(negative)
    for index in allindices
        current = scan(series, index)
        adding, subtracting = index in positive_set, index in negative_set
        mzs = rawmzvalues(current; unit=mz_unit)
        values = rawintensities(current; unit=intensity_unit)
        for (mz, value) in zip(mzs, values)
            positive_sum, negative_sum = get(totals, mz, (zero_intensity, zero_intensity))
            intensity = intensity_type(value)
            totals[mz] = (positive_sum + (adding ? intensity : zero_intensity),
                negative_sum + (subtracting ? intensity : zero_intensity))
        end
    end

    masses = sort!(collect(keys(totals)))
    # Empty groups have a zero sum. Every selected scan counts in the denominator,
    # including those where a particular m/z is absent.
    positive_count = max(length(positive), 1)
    negative_count = max(length(negative), 1)
    values = [max(zero_intensity,
        totals[mz][1] / positive_count - totals[mz][2] / negative_count) for mz in masses]
    return MassSpectrum(masses, mz_unit, values, intensity_unit;
        attrs=(; source=:scan_extraction, aggregation=:mean,
            added_scans=positive, subtracted_scans=negative))
end

function extractmassspectrum(matrix::AbstractMassScanMatrix,
    scanindices::AbstractVector{<:Integer}; subtract::AbstractVector{<:Integer}=Int[])
    positive, negative, _ = extraction_scan_groups(matrix, scanindices, subtract)
    intensities = rawintensities(matrix)
    intensity_type = promote_type(Float64, eltype(intensities))
    zero_intensity = zero(intensity_type)
    positive_count = max(length(positive), 1)
    negative_count = max(length(negative), 1)
    values = Vector{intensity_type}(undef, size(intensities, 2))
    for column in axes(intensities, 2)
        positive_sum = sum((intensity_type(intensities[i, column]) for i in positive);
            init=zero_intensity)
        negative_sum = sum((intensity_type(intensities[i, column]) for i in negative);
            init=zero_intensity)
        values[column] = max(zero_intensity,
            positive_sum / positive_count - negative_sum / negative_count)
    end
    return MassSpectrum(copy(rawmzvalues(matrix)), mzunit(matrix), values, intensityunit(matrix);
        attrs=(; source=:scan_extraction, aggregation=:mean,
            added_scans=positive, subtracted_scans=negative))
end

function extraction_scan_groups(data, scanindices, subtract)
    positive = sort!(unique(Int.(scanindices)))
    negative = sort!(unique(Int.(subtract)))
    allindices = sort!(union(positive, negative))
    isempty(allindices) && throw(ArgumentError("at least one contributing scan is required"))
    for index in allindices
        1 <= index <= scancount(data) || throw(BoundsError(data, index))
    end
    return positive, negative, allindices
end
