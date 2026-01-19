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
using LinearAlgebra: Diagonal, I, cholesky, dot, mul!, norm, Symmetric
using OSQP: OSQP
using Roots: Bisection, find_zero
using SparseArrays: SparseMatrixCSC, sparse, spdiagm, spzeros
using Statistics: cor, mean, median, quantile, std, var
import Base: show, summary
using Printf
using Base.Threads: @threads

const agilent = artifact"Agilent"
const andi = artifact"ANDI"
const calibration = artifact"calibration"
const shimadzu = artifact"Shimadzu"

include("core/core_containers/matrices.jl")
include("core/core_containers/scans.jl")
include("core/core_containers/series.jl")
include("core/core_getters/utils.jl")
include("core/core_getters/matrices.jl")
include("core/core_getters/scans.jl")
include("core/core_getters/series.jl")
include("core/retention_mapping/container.jl")
include("core/retention_mapping/error.jl")
include("core/retention_mapping/fit.jl")
include("core/retention_mapping/map.jl")
include("core/alignment/alignment.jl")
include("core/baseline/baseline.jl")
include("core/convert/mscanmatrix.jl")
include("core/convert/mzchrom.jl")
include("core/quadvar_model/containers.jl")
include("core/quadvar_model/fit.jl")
include("core/transform/binning.jl")
include("core/transform/clr.jl")
include("core/transform/gridding.jl")
include("core/transform/mscanmatrix_ops.jl")
include("core/transform/series_ops.jl")
include("core/transform/whitening.jl")
include("core/utils/collections.jl")
include("core/utils/math.jl")
include("core/utils/units.jl")
include("core/deconvolution/unimodalfit.jl")

export AbstractChromScan
export AbstractChromScanSeries
export AbstractMassScan
export AbstractMassScanSeries
export AbstractMassScanMatrix
export AbstractScan
export AbstractScanSeries
export ChromScan
export ChromScanSeries
export MassScan
export MassScanMatrix
export MassScanSeries
export OptimizationError
export QuadVarParams
export QuadVarFit
export RetentionMapper

export attrs
export acquisition
export airpls
export applymap
export binmzvalues
export binretentions
export clr
export cosdis
export cossim
export DENSE
export densestgrid
export derivinvmap
export derivmap
export extras
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
export unimodalfit_t0
export uniquemzvalues
export user
export varpred
export varpredbias
export vif
export whiten

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
function massspectrum end
function massspectrum! end

# Unitful extension
function __init__()
    Unitful.register(@__MODULE__)
end

end
