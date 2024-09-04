module JuChrom

using Reexport

using LinearAlgebra
using Pkg.Artifacts
using Polynomials
using PrecompileTools: @setup_workload, @compile_workload
@reexport using Unitful

agilent = artifact"Agilent"
andi = artifact"ANDI"

include("base.jl")
include("retentionindices.jl")

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

export binions
export cosine
export extrapolationmethod
export integer
export intensities
export intensity
export interpolationmethod
export ion
export ioncount
export ionindex
export ions
export ionscantime
export ionscantimeindex
export ionscantimes
export ionscantimeshift
export maxintensity
export maxion
export maxretentionindex
export maxretentiontime
export maxscantime
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
export rimapper
export runduration
export scancount
export scanduration
export scantime
export scantimeindex
export scantimes
export totalionchromatogram

include("IO/InputOutput.jl")
using .InputOutput

export FileFormat
export AgilentFID
export ANDI
export ChemStationMS
export DelimitedText
#export Excel
export MassHunterMS
export exportdata
export importdata

end