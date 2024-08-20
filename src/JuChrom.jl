module JuChrom

using PrecompileTools: @setup_workload, @compile_workload

include("base.jl")

export AbstractChromatogram
export AbstractFID
export AbstractGC
export AbstractGCMS
export AbstractTIC
export FID
export RiFID
export GCMS
export TIC
export IonScanOrder
export LinearAscending
export LinearDescending

export binions
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
export retentionindices
export retentionindexname
export runduration
export scancount
export scanduration
export scantime
export scantimeindex
export scantimes
export totalionchromatogram

end