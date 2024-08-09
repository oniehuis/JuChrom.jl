var documenterSearchIndex = {"docs":
[{"location":"man/register/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"man/register/","page":"Index","title":"Index","text":"","category":"page"},{"location":"man/basics/#Basics","page":"Basics","title":"Basics","text":"","category":"section"},{"location":"man/export/#Data-export","page":"Data export","title":"Data export","text":"","category":"section"},{"location":"man/deconvolution/#Deconvolution","page":"Deconvolution","title":"Deconvolution","text":"","category":"section"},{"location":"man/explorer/#GUI-explorer","page":"GUI explorer","title":"GUI explorer","text":"","category":"section"},{"location":"man/import/#Data-import","page":"Data import","title":"Data import","text":"","category":"section"},{"location":"#JuChrom.jl","page":"Home","title":"JuChrom.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for JuChrom.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [JuChrom]\nOrder   = [:type, :function]","category":"page"},{"location":"#JuChrom.AbstractChromatogram","page":"Home","title":"JuChrom.AbstractChromatogram","text":"AbstractChromatogram\n\nThe AbstractChromatogram type is the supertype of all chromatogram implementations in  JuChrom. All subtypes of AbstractChromatogram (e.g., FID, GCMS, TIC) contain  scantimes, intensities, and source information.\n\nSee also AbstractGC, AbstractGCMS, AbstractFID,  AbstractTIC, FID, GCMS, TIC, intensities,  scantimes, source.\n\n\n\n\n\n","category":"type"},{"location":"#JuChrom.AbstractFID","page":"Home","title":"JuChrom.AbstractFID","text":"AbstractFID <: AbstractGC\n\nThe AbstractFID type is the supertype of all flame ionization detector chromatogram  implementations in JuMS (e.g., FID).\n\nSee also AbstractChromatogram, AbstractGC, FID,  intensities, scantimes, source.\n\n\n\n\n\n","category":"type"},{"location":"#JuChrom.AbstractGC","page":"Home","title":"JuChrom.AbstractGC","text":"AbstractGC <: AbstractChromatogram\n\nThe AbstractGC type is the supertype of all chromatogram implementations that have no  mass-charge ratio (m/z) data (= ions) and therefore have a single intensity value  associated with a given scantime in JuChrom (e.g., FID, TIC). The intensities are  stored in a vector whose index corresponds to the index of the associated scantime.\n\nSee also AbstractChromatogram, AbstractFID, AbstractTIC,  FID, TIC, intensities, scantimes, source.\n\n\n\n\n\n","category":"type"},{"location":"#JuChrom.AbstractGCMS","page":"Home","title":"JuChrom.AbstractGCMS","text":"AbstractGCMS <: AbstractChromatogram\n\nThe AbstractGCMS type is the supertype of all chromatogram implementations that contain  mass-charge ratio (m/z) data (= ions) and associated abundance values (= intensities)  and thus can have one or more ion intensity values associated with a given scan time in  JuChrom (e.g., GCMS). The intensities are stored in a matrix where the row index  corresponds to that of the associated scantime and where the column index corresponds  to that of the associated ion.\n\nSee also AbstractChromatogram, GCMS, intensities,  ions, scantimes, source.\n\n\n\n\n\n","category":"type"},{"location":"#JuChrom.AbstractTIC","page":"Home","title":"JuChrom.AbstractTIC","text":"AbstractTIC <: AbstractGC\n\nThe AbstractTIC type is the supertype of all total ion chromatogram implementations in  JuChrom (e.g., TIC).\n\nSee also AbstractChromatogram, AbstractGC, TIC,  intensities, scantimes, source.\n\n\n\n\n\n","category":"type"},{"location":"#JuChrom.FID-Union{Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{T1, T2}, Tuple{T1, T2, T3}} where {T1<:(AbstractVector{var\"#s4\"} where var\"#s4\"<:(Union{Quantity{T, 𝐓, U}, Level{L, S, Quantity{T, 𝐓, U}} where {L, S}} where {T, U})), T2<:(AbstractVector{var\"#s10\"} where var\"#s10\"<:Real), T3<:Union{Nothing, AbstractString}}","page":"Home","title":"JuChrom.FID","text":"FID(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, \nsource::Union{Nothing, AbstractString}=nothing)\n\nReturn a FID object consisting of scan times, intensities, and optionally data source  information.\n\nSee also AbstractChromatogram, AbstractGC, AbstractFID,  scantimes, intensities, source.\n\nIn the following examples, the number types of the arrays passed to the object constructor  are explicitly annotated to illustrate that the FID object preserves the types.\n\nExample\n\njulia> FID(Int64[1, 2, 3]u\"s\", Int32[12, 956, 1])\nFID {scantimes: Int64, intensities: Int32}\nSource: nothing\n3 scans; time range: 1 s - 3 s\nintensity range: 1 - 956\n\njulia> FID(Int32[1, 2, 3]u\"s\", Float64[12.0, 956.0, 1.0], \"example data\")\nFID {scantimes: Int32, intensities: Float64}\nSource: example data\n3 scans; time range: 1 s - 3 s\nintensity range: 1.0 - 956.0\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.GCMS-Union{Tuple{T4}, Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{T1, T2, T3}, Tuple{T1, T2, T3, T4}} where {T1<:(AbstractVector{var\"#s4\"} where var\"#s4\"<:(Union{Quantity{T, 𝐓, U}, Level{L, S, Quantity{T, 𝐓, U}} where {L, S}} where {T, U})), T2<:(AbstractVector{var\"#s14\"} where var\"#s14\"<:Real), T3<:(AbstractMatrix{var\"#s15\"} where var\"#s15\"<:Real), T4<:Union{Nothing, AbstractString}}","page":"Home","title":"JuChrom.GCMS","text":"GCMS(scantimes::AbstractVector{<:Unitful.Time}, ions::AbstractVector{<:Real}, \nintensities::AbstractMatrix{<:Real}; source::Union{Nothing, AbstractString})\n\nReturn a GCMS object consisting of scantimes, ions, intensities, and optionally data  source information.\n\nSee also AbstractChromatogram, AbstractGCMS, intensities,  ions, scantimes, source.\n\nIn the following examples, the number types of the arrays passed to the object constructor  are explicitly annotated to illustrate that the GCMS object preserves the types.\n\nExample\n\njulia> GCMS(Int32[1, 2, 3]u\"s\", Int64[85, 100], Int32[0 12; 34 956; 23 1])\nGCMS {scantimes: Int32, ions: Int64, intensities: Int32}\nSource: nothing\n3 scans; time range: 1 s - 3 s\n2 ions; range: m/z 85 - 100\nintensity range: 0 - 956\n\njulia> GCMS([1.1f0, 2.1f0]u\"s\", Float32[35.1, 76.2], Int64[0 12; 34 956], \"example data\")\nGCMS {scantimes: Float32, ions: Float32, intensities: Int64}\nSource: example data\n2 scans; time range: 1.1f0 s - 2.1f0 s\n2 ions; range: m/z 35.1 - 76.2\nintensity range: 0 - 956\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.TIC-Union{Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{T1, T2}, Tuple{T1, T2, T3}} where {T1<:(AbstractVector{var\"#s8\"} where var\"#s8\"<:(Union{Quantity{T, 𝐓, U}, Level{L, S, Quantity{T, 𝐓, U}} where {L, S}} where {T, U})), T2<:(AbstractVector{var\"#s7\"} where var\"#s7\"<:Real), T3<:Union{Nothing, AbstractString}}","page":"Home","title":"JuChrom.TIC","text":"TIC(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, \nsource::Union{Nothing, AbstractString}=nothing)\n\nReturn a TIC object consisting of scan times, intensities, and optionally data source  information.\n\nSee also AbstractChromatogram, AbstractGC, AbstractTIC,  scantimes, intensities, source.\n\nIn the following examples, the number types of the arrays passed to the object constructor  are explicitly annotated to illustrate that the TIC object preserves the types.\n\nExample\n\njulia> TIC(Int64[1, 2, 3]u\"s\", Int32[12, 956, 1])\nTIC {scantimes: Int64, intensities: Int32}\nSource: nothing\n3 scans; time range: 1 s - 3 s\nintensity range: 1 - 956\n\njulia> TIC(Int32[1, 2, 3]u\"s\", Float64[12.0, 956.0, 1.0], \"example data\")\nTIC {scantimes: Int32, intensities: Float64}\nSource: example data\n3 scans; time range: 1 s - 3 s\nintensity range: 1.0 - 956.0\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.intensities-Tuple{AbstractChromatogram}","page":"Home","title":"JuChrom.intensities","text":"intensities(chrom::AbstractChromatogram)\n\nReturn the intensities of the AbstractChromatogram subtype object.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  minintensity, maxintensity.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", [85, 100], Int64[0 12; 34 956; 23 1]);\n\njulia> intensities(gcms)\n3×2 Matrix{Int64}:\n  0   12\n 34  956\n 23    1\n\njulia> fid = FID([1, 2, 3]u\"s\", Float32[12.0, 956.0, 23.0]);\n\njulia> intensities(fid)\n3-element Vector{Float32}:\n  12.0\n 956.0\n  23.0\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.ioncount-Tuple{AbstractGCMS}","page":"Home","title":"JuChrom.ioncount","text":"ioncount(gcms::AbstractGCMS) -> Int\n\nReturn the number of ions in the AbstractGCMS subtype object.\n\nSee also AbstractGCMS, GCMS, ions.\n\nExample\n\njulia> gcms = GCMS([1.0, 2.0, 3.0]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> ioncount(gcms)\n2\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.ions-Tuple{AbstractGCMS}","page":"Home","title":"JuChrom.ions","text":"ions(gcms::AbstractGCMS)\n\nReturn the ions of the AbstractGCMS subtype object.\n\nSee also AbstractGCMS, GCMS, minion, maxion,  ioncount.\n\nIn the following example, the number type of ions passed to the object constructor  is explicitly annotated to illustrate that the GCMS object preserves the type.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", Int64[85, 100], [0 12; 34 956; 23 1]);\n\njulia> ions(gcms)\n2-element Vector{Int64}:\n  85\n 100\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.maxintensity-Tuple{AbstractChromatogram}","page":"Home","title":"JuChrom.maxintensity","text":"maxintensity(chrom::AbstractChromatogram)\n\nReturn the maximum intensity in the AbstractChromatogram subtype object.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  minintensity, intensities.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", [85, 100], Int64[0 12; 34 956; 23 1]);\n\njulia> maxintensity(gcms)\n956\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.maxion-Tuple{AbstractGCMS}","page":"Home","title":"JuChrom.maxion","text":"maxion(gcms::AbstractGCMS)\n\nReturn the largest ion in the AbstractGCMS subtype object.\n\nSee also AbstractGCMS, GCMS, minion, ions,  ioncount.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", Int64[85, 100], [0 12; 34 956; 23 1]);\n\njulia> maxion(gcms)\n100\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.maxscantime-Tuple{AbstractChromatogram}","page":"Home","title":"JuChrom.maxscantime","text":"maxscantime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, \nustripped::Bool)\n\nReturn the time of the last scan from the AbstractChromatogram subtype object. The optional  keyword argument timeunit allows you to change the unit of the returned scan time. All  time units defined in the package Unitful.jl  (e.g., u\"s\", u\"minute\") are supported. The optional keyword argument ustripped allows  you to specify whether the unit is stripped from the returned value.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  minscantime, scantimes, scancount.\n\nExample\n\njulia> gcms = GCMS(Float64[1.0, 2.0, 3.0]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> maxscantime(gcms)\n3.0 s\n\njulia> maxscantime(gcms, timeunit=u\"minute\")\n0.05 minute\n\njulia> maxscantime(gcms, timeunit=u\"minute\", ustripped=true)\n0.05\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.minintensity-Tuple{AbstractChromatogram}","page":"Home","title":"JuChrom.minintensity","text":"minintensity(chrom::AbstractChromatogram)\n\nReturn the minimum intensity in the AbstractChromatogram subtype object.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  maxintensity, intensities.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", [85, 100], Int64[0 12; 34 956; 23 1]);\n\njulia> minintensity(gcms)\n0\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.minion-Tuple{AbstractGCMS}","page":"Home","title":"JuChrom.minion","text":"minion(gcms::AbstractGCMS)\n\nReturn the smallest ion in the AbstractGCMS subtype object.\n\nSee also AbstractGCMS, GCMS, ions, ioncount.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", Int64[85, 100], [0 12; 34 956; 23 1]);\n\njulia> minion(gcms)\n85\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.minscantime-Tuple{AbstractChromatogram}","page":"Home","title":"JuChrom.minscantime","text":"minscantime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, \nustripped::Bool=false)\n\nReturn the time of the first scan from the AbstractChromatogram subtype object. The  optional keyword argument timeunit allows you to change the unit of the returned scan  time. All time units defined in the package  Unitful.jl (e.g., u\"s\", u\"minute\") are  supported. The optional keyword argument ustripped allows you to specify whether the unit  is stripped from the returned value.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  maxscantime, scantimes, scancount.\n\nExample\n\njulia> gcms = GCMS(Float64[1.0, 2.0, 3.0]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> minscantime(gcms)\n1.0 s\n\njulia> minscantime(gcms, timeunit=u\"minute\")\n0.016666666666666666 minute\n\njulia> minscantime(gcms, timeunit=u\"minute\", ustripped=true)\n0.016666666666666666\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.scancount-Tuple{AbstractChromatogram}","page":"Home","title":"JuChrom.scancount","text":"scancount(chrom::AbstractChromatogram) -> Int\n\nReturn the number of scans from the AbstractChromatogram subtype object.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  scantimes.\n\nExample\n\njulia> fid = FID([1, 2, 3]u\"s\", [12, 956, 23]);\n\njulia> scancount(fid)\n3\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.scantimes-Tuple{AbstractChromatogram}","page":"Home","title":"JuChrom.scantimes","text":"scantimes(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, \nustripped::Bool=false)\n\nReturn the scantimes of the AbstractChromatogram subtype object. The optional keyword  argument timeunit allows you to change the unit of the returned scantimes. All time  units defined in the package Unitful.jl  (e.g., u\"s\", u\"minute\") are supported. The optional keyword argument ustripped  allows you to specify whether the unit is stripped from the returned values.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  minscantime, maxscantime.\n\nIn the following example, the element type of the scantime vector passed to the object  constructor is explicitly annotated to illustrate that the GCMS object preserves the type.\n\nExample\n\njulia> gcms = GCMS(Float32[1.0, 2.0, 3.0]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> scantimes(gcms)\n3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:\n 1.0f0 s\n 2.0f0 s\n 3.0f0 s\n\njulia> scantimes(gcms, timeunit=u\"minute\")\n3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:\n 0.016666668f0 minute\n 0.033333335f0 minute\n 0.050000004f0 minute\n\njulia> scantimes(gcms, timeunit=u\"minute\", ustripped=true)\n3-element Vector{Float32}:\n 0.016666668\n 0.033333335\n 0.050000004\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.source-Tuple{AbstractChromatogram}","page":"Home","title":"JuChrom.source","text":"source(chrom::AbstractChromatogram) -> Type{<:Union{AbstractString, Nothing}}\n\nReturn the source information of the AbstractChromatogram subtype object.\n\nSee also AbstractChromatogram, FID, GCMS, TIC.\n\nExample\n\njulia> gcms₁ = GCMS([1, 2, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1], \"example data\");\n\njulia> source(gcms₁)\n\"example data\"\n\njulia> gcms₂ = GCMS([1, 2, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> source(gcms₂)\n\njulia> isnothing(source(gcms₂))\ntrue\n\n\n\n\n\n","category":"method"},{"location":"man/installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"}]
}
