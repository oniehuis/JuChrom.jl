module JuChrom

import BSplineKit: BSplineOrder, Linear, extrapolate, interpolate
using Reexport
using Pkg.Artifacts
using PrecompileTools: @setup_workload, @compile_workload

@reexport using Unitful

agilent = artifact"Agilent"
andi = artifact"ANDI"

include("base.jl")

export AbstractChromatogram
export AbstractFID
export AbstractGC
export AbstractGCMS
export AbstractTIC
export FID
export RiFID
export GCMS
export RiGCMS
export TIC
export RiTIC
export IonScanOrder
export LinearAscending
export LinearDescending

export binions
export bsplineinterpolation
export cosine
export integer
export intensities
export intensity
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
export maxscantime
export metadata
export minintensity
export minion
export minretentionindex
export minscantime
export retentionindex
export retentionindices
export retentionindexname
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
#export DelimitedText
#export Excel
export MassHunterMS
#export exportdata
export importdata

end