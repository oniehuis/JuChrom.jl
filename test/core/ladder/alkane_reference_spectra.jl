using Test
using JuChrom

@testset "alkane reference spectrum constructor" begin
    mzs = Int16[43, 57, 71]
    relativeintensities = Float32[0.25, 1.0, 0.5]
    label = SubString("pentadecane reference", 1, 11)

    spectrum = JuChrom.alkane_reference_spectrum(
        Int8(15),
        label,
        1500//1,
        mzs,
        relativeintensities,
    )

    @test spectrum isa MassSpectrum
    @test mzvalues(spectrum) == [43, 57, 71]
    @test eltype(mzvalues(spectrum)) == Int
    @test intensities(spectrum) == [0.25, 1.0, 0.5]
    @test eltype(intensities(spectrum)) == Float64

    spectrum_attrs = attrs(spectrum)
    @test spectrum_attrs.ri === 1500.0
    @test spectrum_attrs.order === 15
    @test spectrum_attrs.label == "pentadecane"
    @test spectrum_attrs.label isa String
    @test spectrum_attrs.normalization == :base_peak
    @test spectrum_attrs.minrelativeintensity == 0.0

    mzs[1] = 41
    relativeintensities[1] = 0.75
    @test mzvalues(spectrum) == [43, 57, 71]
    @test intensities(spectrum) == [0.25, 1.0, 0.5]

    @test_throws DimensionMismatch JuChrom.alkane_reference_spectrum(
        15,
        "invalid",
        1500.0,
        [43, 57],
        [1.0],
    )
end

@testset "built-in alkane reference spectra" begin
    spectra = alkanereferencespectra()

    @test length(spectra) == 33
    @test all(spectrum -> spectrum isa MassSpectrum, spectra)
    @test [attrs(spectrum).order for spectrum in spectra] == collect(8:40)
    @test [attrs(spectrum).ri for spectrum in spectra] == Float64.(800:100:4000)
    @test attrs(first(spectra)).label == "octane"
    @test attrs(last(spectra)).label == "tetracontane"

    for spectrum in spectra
        spectrum_attrs = attrs(spectrum)

        @test spectrum_attrs.normalization == :base_peak
        @test spectrum_attrs.minrelativeintensity == 0.0
        @test eltype(mzvalues(spectrum)) == Int
        @test eltype(intensities(spectrum)) == Float64
        @test length(mzvalues(spectrum)) == length(intensities(spectrum))
        @test !isempty(mzvalues(spectrum))
        @test all(isfinite, mzvalues(spectrum))
        @test all(isfinite, intensities(spectrum))
        @test all(>(0), mzvalues(spectrum))
        @test all(>=(0), intensities(spectrum))
        @test all(>(0), diff(mzvalues(spectrum)))
        @test maximum(intensities(spectrum)) ≈ 100.0
        @test all(x -> x == round(x; sigdigits=4), intensities(spectrum))
    end
end

@testset "alkane reference spectrum accessors" begin
    c8 = alkanereferencespectrum(8)
    c40 = alkanereferencespectrum(40)

    @test c8 isa MassSpectrum
    @test attrs(c8).order == 8
    @test attrs(c8).label == "octane"
    @test first(mzvalues(c8)) == 29

    @test c40 isa MassSpectrum
    @test attrs(c40).order == 40
    @test attrs(c40).label == "tetracontane"
    @test last(mzvalues(c40)) == 562

    @test_throws ArgumentError alkanereferencespectrum(7)
    @test_throws ArgumentError alkanereferencespectrum(41)
end

@testset "alkane reference spectra are returned as copies" begin
    spectra = alkanereferencespectra()
    original_c8_first_intensity = first(intensities(alkanereferencespectrum(8)))

    intensities(first(spectra))[1] = 999.0
    @test first(intensities(alkanereferencespectrum(8))) == original_c8_first_intensity

    c8 = alkanereferencespectrum(8)
    intensities(c8)[1] = 999.0
    @test first(intensities(alkanereferencespectrum(8))) == original_c8_first_intensity
end

@testset "default alkane standard" begin
    standard = defaultalkanestandard()

    @test standard isa AlkaneStandard
    @test standard.name == "n-alkane C8-C40 reference spectra"
    @test standard.metadata.carbon_range == 8:40
    @test standard.metadata.normalization == :base_peak
    @test standard.metadata.minrelativeintensity == 0.0
    @test length(standard.spectra) == 33
    @test [attrs(spectrum).order for spectrum in standard.spectra] == collect(8:40)

    original_c8_first_intensity = first(intensities(alkanereferencespectrum(8)))
    intensities(first(standard.spectra))[1] = 999.0
    @test first(intensities(defaultalkanestandard().spectra[1])) ==
        original_c8_first_intensity
end
