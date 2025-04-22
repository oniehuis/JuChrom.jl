module JuChrom

using Reexport

using Pkg.Artifacts
import BasicInterpolators
import Optim
import BSplineKit
import BSplineKit: diff
import Roots

@reexport using Unitful

const agilent = artifact"Agilent"
const andi = artifact"ANDI"
const calibration = artifact"calibration"
const shimadzu = artifact"Shimadzu"

include("base.jl")
include("retentionindices.jl")
include("massspectra.jl")

export AbstractChromatogram
export AbstractChrom
export AbstractChromMS
export Chrom
export ChromMS
export IonScanOrder
export LinearAscending
export LinearDescending
export AbstractRiMapper
export RiMapper
export PolationMethod
export NaturalCubicBSpline
export Linear
export PiecewiseLinear
export AbstractMassSpectrum
export MassSpectrum

export binions
export cosine
export extrapolationmethod
export integer
export intensities
export intensity
export intensitydifferences
export intensitysums
export interpolationmethod
export ion
export ioncount
export ionindex
export ions
export ionscantime
export ionscantimeindex
export ionscantimes
export ionscantimeshift
export massspectrum
export maxintensity
export maxion
export maxretentionindex
export maxretentiontime
export maxscantime
export meanintensities
export metadata
export minintensity
export minion
export minretentionindex
export minretentiontime
export minscantime
export retentionindex
export retentionindices
export retentionindexname
export retentiontimes
export retentiontime
export rimapper
export rimapper!
export runduration
export scancount
export scanduration
export scantime
export scantimeindex
export scantimes
export sharedions
export similarity
export totalionchromatogram

include("IO/InputOutput.jl")
using .InputOutput

export FileFormat
export AgilentFID
export ANDI
export ChemStationMS
export DelimitedText
export Excel
export MassHunterMS
export ShimadzuMS
export exportdata
export importdata

end