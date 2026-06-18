__precompile__()

module JuChrom

using Pkg.Artifacts
using Reexport: @reexport
using Unitful: @unit, AbstractQuantity, uconvert, ustrip
@reexport using Unitful

# Additional package dependencies
using BSplineKit: AbstractBSplineBasis, BSplineBasis, BSplineOrder, Derivative, Natural,
    RecombinedBSplineBasis, Spline, collocation_matrix
import HiGHS
using JuMP: Model, MOI, optimize!, set_silent, termination_status, value, @constraint,
    @objective, @variable
using LinearAlgebra: Diagonal, I, cholesky, dot, eigen, mul!, norm, svd, Symmetric,
    diagind
using OSQP: OSQP
using Roots: Bisection, find_zero
using SHA: SHA256_CTX, bytes2hex, digest!, update!
using SparseArrays: SparseMatrixCSC, sparse, spdiagm, spzeros, nnz
using Statistics: cor, mean, median, quantile, std, var
import Base: show, summary
import Random
using Printf
using Base.Threads: @threads

const agilent = artifact"Agilent"
const andi = artifact"ANDI"
const calibration = artifact"calibration"
const shimadzu = artifact"Shimadzu"

include("core/core_containers/matrices.jl")
include("core/core_containers/variance_matrices.jl")
include("core/core_containers/scans.jl")
include("core/core_containers/series.jl")
include("core/core_getters/utils.jl")
include("core/core_getters/matrices.jl")
include("core/core_getters/variance_matrices.jl")
include("core/core_getters/scans.jl")
include("core/core_getters/series.jl")
include("core/retention_mapping/container.jl")
include("core/retention_mapping/error.jl")
include("core/retention_mapping/fit.jl")
include("core/retention_mapping/map.jl")
include("core/alignment/alignment.jl")
include("core/baseline/arpls.jl")
include("core/convert/mscanmatrix.jl")
include("core/convert/mzchrom.jl")
include("core/variances/countvariances.jl")
include("core/variances/quadvarmodel/containers.jl")
include("core/variances/quadvarmodel/fit.jl")
include("core/scan_timing/mz_retention.jl")
include("core/transform/binning.jl")
include("core/transform/clr.jl")
include("core/transform/gridding.jl")
include("core/transform/intensity_units.jl")
include("core/transform/mscanmatrix_ops.jl")
include("core/transform/series_ops.jl")
include("core/transform/whitening.jl")
include("core/utils/collections.jl")
include("core/utils/math.jl")
include("core/utils/units.jl")
include("core/deconvolution/unimodalfit.jl")
include("core/deconvolution/parafac2.jl")
include("core/massspectrum/massspectrum.jl")
include("core/ladder/alkane_reference_spectra.jl")
include("core/ladder/alkaneseries.jl")
include("core/ladder/abundances.jl")
include("core/ladder/molecular_ions.jl")
include("core/ladder/path.jl")
include("core/ladder/apexes.jl")
include("core/ladder/additions.jl")
include("core/ladder/mass_spectra.jl")

export AbstractChromScan
export AbstractChromScanSeries
export AbstractMassScan
export AbstractMassScanSeries
export AbstractMassScanMatrix
export AbstractMassSpectrum
export AbstractScan
export AbstractScanSeries
export AbstractVarianceMassScanMatrix
export AlkaneAbundanceInfo
export AlkaneAbundanceWindow
export AlkaneChannelInfo
export AlkaneLadderCalibrationPoint
export AlkaneLadderStep
export AlkaneMolecularIonContrast
export AlkaneMolecularIonInfo
export AlkaneReferenceChannels
export AlkaneSeriesResult
export AlkaneStandard
export ChromScan
export ChromScanSeries
export CountVarianceEstimate
export MassScan
export MassScanMatrix
export MassScanSeries
export MassSpectrum
export OptimizationError
export Parafac2Fit
export QuadVarParams
export QuadVarFit
export RetentionMapper
export VarianceMassScanMatrix

export attrs
export acquisition
export arpls
export applymap
export binmzvalues
export binretentions
export clr
export cosdis
export cossim
export countvariances
export DENSE
export densestgrid
export derivinvmap
export derivmap
export dwellnormalize
export extras
export findalkanes
export findalkaneseries
export alkaneladdercalibrationpoints
export alkaneladdersteps
export alkaneladderpath
export alkaneladderapex
export alkaneladderapexes
export alkaneladderadditions
export alkaneladderscanorder
export alkaneladdermassspectrum
export alkaneladdermassspectra
export alkanereferencespectra
export alkanereferencespectrum
export defaultalkanestandard
export findclosest
export fitmap
export fitquadvarmodel
export gapalign
export indextrim
export indextrim!
export instrument
export integer
export intensities
export intensity
export intensityunit
export invmap
export invmapmax
export invmapmin
export ioncount
export ions
export level
export levels
export levelscans
export mapmin
export mapmax
export mscanmatrix
export mzchrom
export mzcount
export mzindex
export mzretention
export mzunit
export mzvalues
export residautocorr
export sample
export scan
export scancount
export scans
export SPARSE
export rawapplymap
export rawderivinvmap
export rawderivmap
export rawintensity
export rawintensities
export rawinvmap
export rawinvmapmax
export rawinvmapmin
export rawmapmax
export rawmapmin
export rawmzvalues
export rawretention
export rawretentions
export rawretentions_A
export rawretentions_B
export retention
export retentions
export retentions_A
export retentions_B
export retentiontrim
export retentiontrim!
export retentionunit
export retentionunit_A
export retentionunit_B
export typify
export unimodalfit
export unimodalfit_apexsearch
export uniquemzvalues
export user
export varpred
export varpredbias
export varianceunit
export variances
export vif
export whiten
export withintensityunit
export rawvariances
export parafac2
export parafac2abundances
export parafac2apexes
export parafac2fitpercent
export parafac2loss
export parafac2profilediagnostics
export parafac2profileminima
export parafac2reconstruct
export parafac2residuals
export parafac2scores
export parafac2spectra

include("IO/InputOutput.jl")
using .InputOutput

export FileFormatError
export load
export set_verbosity

include("IO/Loaders/AgilentFIDLoader.jl")
include("IO/Loaders/ChemStationMSLoader.jl")
include("IO/Loaders/MassHunterMSLoader.jl")

# Makie extension
export massspectrum
export massspectrum!
export tictrace
export tictrace!
function massspectrum end
function massspectrum! end
function tictrace end
function tictrace! end

# Unitful extension
function __init__()
    Unitful.register(@__MODULE__)
    Unitful.register(JuChromUnits)
end

end
