var documenterSearchIndex = {"docs":
[{"location":"man/register/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"man/register/","page":"Index","title":"Index","text":"","category":"page"},{"location":"man/basics/#Basics","page":"Basics","title":"Basics","text":"","category":"section"},{"location":"man/basics/","page":"Basics","title":"Basics","text":"AbstractChromatogram\nAbstractFID\nAbstractGC\nAbstractGCMS\nAbstractTIC\nFID\nGCMS\nTIC\nbinions\ninteger\nintensities\nion\nioncount\nionindex\nions\nmaxintensity\nmaxion\nmaxscantime\nmetadata\nminintensity\nminion\nminscantime\nrunduration\nscancount\nscantime\nscantimeindex\nscantimes\ntotalionchromatogram","category":"page"},{"location":"man/basics/#JuChrom.AbstractChromatogram","page":"Basics","title":"JuChrom.AbstractChromatogram","text":"AbstractChromatogram\n\nSupertype of all chromatogram implementations. All subtypes (e.g., FID, GCMS, TIC)  contain scantimes, intensities, and metadata.\n\nSee also AbstractGC, AbstractFID, AbstractGCMS,  AbstractTIC, FID, GCMS, TIC, scantimes,  intensities, metadata.\n\n\n\n\n\n","category":"type"},{"location":"man/basics/#JuChrom.AbstractFID","page":"Basics","title":"JuChrom.AbstractFID","text":"AbstractFID <: AbstractGC\n\nSupertype of all flame ionization detector chromatogram implementations in JuMS (e.g.,  FID).\n\nSee also AbstractChromatogram, AbstractFID, AbstractGC,  FID, scantimes, intensities, metadata.\n\n\n\n\n\n","category":"type"},{"location":"man/basics/#JuChrom.AbstractGC","page":"Basics","title":"JuChrom.AbstractGC","text":"AbstractGC <: AbstractChromatogram\n\nSupertype of all chromatogram implementations that have no mass-charge ratio (m/z) data  (= ions) and therefore have a single intensity value associated with a given scantime  (e.g., FID, TIC). The intensities are stored in a vector whose index corresponds to  the index of the associated scantime.\n\nSee also AbstractChromatogram, AbstractFID, AbstractTIC,  FID, TIC, scantimes, intensities,  metadata.\n\n\n\n\n\n","category":"type"},{"location":"man/basics/#JuChrom.AbstractGCMS","page":"Basics","title":"JuChrom.AbstractGCMS","text":"AbstractGCMS <: AbstractChromatogram\n\nSupertype of all chromatogram implementations that contain mass-charge ratio (m/z) data  (= ions) and associated abundance values (= intensities) and thus can have one or more  ion intensity values associated with a given scan time (e.g., GCMS). The intensities  are stored in a matrix where the row index corresponds to that of the associated scantime  and where the column index corresponds to that of the associated ion.\n\nSee also AbstractChromatogram, AbstractGCMS, GCMS,  scantimes, ions, intensities, metadata.\n\n\n\n\n\n","category":"type"},{"location":"man/basics/#JuChrom.AbstractTIC","page":"Basics","title":"JuChrom.AbstractTIC","text":"AbstractTIC <: AbstractGC\n\nSupertype of all total ion chromatogram implementations (e.g., TIC).\n\nSee also AbstractChromatogram, AbstractGC, AbstractTIC,  TIC, scantimes, intensities, metadata.\n\n\n\n\n\n","category":"type"},{"location":"man/basics/#JuChrom.FID","page":"Basics","title":"JuChrom.FID","text":"FID(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, \nmetadata::Dict=Dict{Any, Any})\n\nConstruct a FID object consisting of scantimes, intensities, and metadata. Note that  the scantimes must be in ascending order and the intensities must not contain values  less than zero.\n\nSee also AbstractChromatogram, AbstractGC, AbstractFID,  scantimes, intensities, metadata.\n\nIn the following examples, the number types of the arrays passed to the object constructor  are explicitly annotated to illustrate that the FID object preserves the types.\n\nExample\n\njulia> FID(Int64[1, 2, 3]u\"s\", Int32[12, 956, 1])\nFID {scantimes: Int64, intensities: Int32}\n3 scans; time range: 1 s - 3 s\nintensity range: 1 - 956\nmetadata: 0\n\njulia> FID(Int32[1, 2, 3]u\"s\", Float64[12.0, 956.0, 1.0], Dict(\"name\" => \"sample\"))\nFID {scantimes: Int32, intensities: Float64}\n3 scans; time range: 1 s - 3 s\nintensity range: 1.0 - 956.0\nmetadata: 1\n\njulia> FID([2, 1, 3]u\"s\", [12.0, 956.0, 1.0])\nERROR: ArgumentError: scantimes not in ascending order\n[...]\n\njulia> FID([1, 2, 3]u\"s\", [-12.0, 956.0, 1.0])\nERROR: ArgumentError: intensity values contain values less than zero\n[...]\n\n\n\n\n\n","category":"type"},{"location":"man/basics/#JuChrom.GCMS","page":"Basics","title":"JuChrom.GCMS","text":"GCMS(scantimes::AbstractVector{<:Unitful.Time}, ions::AbstractVector{<:Real}, \nintensities::AbstractMatrix{<:Real}; metadata::Dict=Dict{Any, Any}())\n\nConstruct a GCMS object consisting of scantimes, ions, intensities, and metadata.  Note that the scantimes and the ions must be in ascending order and the intensities  must not contain values less than zero.\n\nSee also AbstractChromatogram, AbstractGCMS, scantimes,  ions, intensities, metadata, totalionchromatogram.\n\nIn the following examples, the number types of the arrays passed to the object constructor  are explicitly annotated to illustrate that the GCMS object preserves the types.\n\nExample\n\njulia> GCMS(Int32[1, 2, 3]u\"s\", Int64[85, 100], Int32[0 12; 34 956; 23 1])\nGCMS {scantimes: Int32, ions: Int64, intensities: Int32}\n3 scans; time range: 1 s - 3 s\n2 ions; range: m/z 85 - 100\nintensity range: 0 - 956\nmetadata: 0\n\njulia> GCMS([1.1f0, 2.1f0]u\"s\", [35.1f0, 76.2f0], Int64[0 12; 34 956], Dict(:id => 1))\nGCMS {scantimes: Float32, ions: Float32, intensities: Int64}\n2 scans; time range: 1.1f0 s - 2.1f0 s\n2 ions; range: m/z 35.1 - 76.2\nintensity range: 0 - 956\nmetadata: 1\n\njulia> GCMS([2, 1, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1])\nERROR: ArgumentError: scantimes not in ascending order\n[...]\n\njulia> GCMS([1, 2, 3]u\"s\", [100, 85], [0 12; 34 956; 23 1])\nERROR: ArgumentError: ions not in ascending order\n[...]\n\njulia> GCMS([1, 2, 3]u\"s\", [85, 100], [0 -12; 34 956; 23 1])\nERROR: ArgumentError: intensity values contain at least one value less than zero\n[...]\n\n\n\n\n\n","category":"type"},{"location":"man/basics/#JuChrom.TIC","page":"Basics","title":"JuChrom.TIC","text":"TIC(scantimes::AbstractVector{<:Unitful.Time}, intensities::AbstractVector{<:Real}, \nmetadata::Dict=Dict{Any, Any}())\n\nConstruct a TIC object consisting of scantimes, intensities, and metadata. Note that  the scantimes must be in ascending order and the intensities must not contain values  less than zero.\n\nSee also AbstractChromatogram, AbstractGC, AbstractTIC,  scantimes, intensities, metadata,  totalionchromatogram.\n\nIn the following examples, the number types of the arrays passed to the object constructor  are explicitly annotated to illustrate that the TIC object preserves the types.\n\nExample\n\njulia> TIC(Int64[1, 2, 3]u\"s\", Int32[12, 956, 1])\nTIC {scantimes: Int64, intensities: Int32}\n3 scans; time range: 1 s - 3 s\nintensity range: 1 - 956\nmetadata: 0\n\njulia> TIC(Int32[1, 2, 3]u\"s\", Float64[12.0, 956.0, 1.0], Dict(\"name\" => \"sample\"))\nTIC {scantimes: Int32, intensities: Float64}\n3 scans; time range: 1 s - 3 s\nintensity range: 1.0 - 956.0\nmetadata: 1\n\njulia> TIC([2, 1, 3]u\"s\", [12.0, 956.0, 1.0])\nERROR: ArgumentError: scantimes not in ascending order\n[...]\n\njulia> TIC([1, 2, 3]u\"s\", [-12.0, 956.0, 1.0])\nERROR: ArgumentError: intensity values contain at least one value less than zero\n[...]\n\n\n\n\n\n","category":"type"},{"location":"man/basics/#JuChrom.binions","page":"Basics","title":"JuChrom.binions","text":"binions(gcms::AbstractGCMS; ionbin::Function=integerion)\n\nReturn a GCMS object in which the ions are binned according to the ionbin function (default  function is integer) and the intensity values of the binned ions are summed.\n\nSee also AbstractGCMS, GCMS, integer, intensities,  ions, ioncount.\n\nExample\n\njulia> gcms = GCMS((1:3)u\"s\", [84.8, 85.2, 100.9], [0 24 12; 0 0 956; 23 0 1])\nGCMS {scantimes: Int64, ions: Float64, intensities: Int64}\n3 scans; time range: 1 s - 3 s\n3 ions; range: m/z 84.8 - 100.9\nintensity range: 0 - 956\nmetadata: 0\n\njulia> intensities(gcms)\n3×3 Matrix{Int64}:\n  0  24   12\n  0   0  956\n 23   0    1\n\njulia> gcmsᵢ = binions(gcms)\nGCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n3 scans; time range: 1 s - 3 s\n2 ions; range: m/z 85 - 101\nintensity range: 0 - 956\nmetadata: 0\n\njulia> ions(gcmsᵢ)\n2-element Vector{Int64}:\n  85\n 101\n\njulia> intensities(gcmsᵢ)\n3×2 Matrix{Int64}:\n 24   12\n  0  956\n 23    1\n\n\njulia> custom_ionbin(ion) = integer(ion, start=0.9);\n\njulia> gcmsᵢ = binions(gcms, ionbin=custom_ionbin)\nGCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n3 scans; time range: 1 s - 3 s\n3 ions; range: m/z 84 - 101\nintensity range: 0 - 956\nmetadata: 0\n\njulia> ions(gcmsᵢ)\n3-element Vector{Int64}:\n  84\n  85\n 101\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.integer","page":"Basics","title":"JuChrom.integer","text":"integer(value:::Real; start::Real=0.7) -> Int\n\nReturn the integer for the given value that satisfies the following condition:  integer - 1 + start ≤ value < integer + start, where 0 ≤ start < 1.\n\nSee also AbstractGCMS, GCMS, binions, ions.\n\nExample\n\njulia> integer(29.7)\n30\n\njulia> integer(30.0)\n30\n\njulia> integer(30.69)\n30\n\njulia> integer(29.7, start=0.8)\n29\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.intensities","page":"Basics","title":"JuChrom.intensities","text":"intensities(chrom::AbstractChromatogram)\n\nReturn the intensities.\n\nSee also AbstractChromatogram, AbstractGC, AbstractGCMS,  FID, GCMS, TIC, minintensity,  maxintensity, scantimes, scancount, ions,  ioncount.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", [85, 100], Int64[0 12; 34 956; 23 1]);\n\njulia> intensities(gcms)\n3×2 Matrix{Int64}:\n  0   12\n 34  956\n 23    1\n\n\njulia> intensities(gcms)[3, :]  # all intensites of 3rd scan\n2-element Vector{Int64}:\n 23\n  1\n\njulia> intensities(gcms)[:, 1]  # all intensites of 1st ion\n3-element Vector{Int64}:\n  0\n 34\n 23\n\njulia> intensities(gcms)[2, 1]  # intensity of 2nd scan, 1st ion\n34\n\njulia> fid = FID([1, 2, 3]u\"s\", Float32[12.0, 956.0, 23.0]);\n\njulia> intensities(fid)\n3-element Vector{Float32}:\n  12.0\n 956.0\n  23.0\n\njulia> intensities(fid)[2]  # intensity of 2nd scan\n956.0f0\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.ion","page":"Basics","title":"JuChrom.ion","text":"ion(gcms::AbstractGCMS, index::Integer)\n\nReturn the ion at the specified index.\n\nSee also AbstractGCMS, ions, ionindex, minion,  maxion, ioncount.\n\nExample\n\njulia> gcms = GCMS((1:3)u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> ion(gcms, 1)\n85\n\njulia> ion(gcms, 2)\n100\n\njulia> ion(gcms, 3)\nERROR: BoundsError: attempt to access 2-element Vector{Int64} at index [3]\n[...]\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.ioncount","page":"Basics","title":"JuChrom.ioncount","text":"ioncount(gcms::AbstractGCMS) -> Int\n\nReturn the number of ions.\n\nSee also AbstractGCMS, GCMS, ions, ion,  minion, maxion.\n\nExample\n\njulia> gcms = GCMS([1.0, 2.0, 3.0]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> ioncount(gcms)\n2\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.ionindex","page":"Basics","title":"JuChrom.ionindex","text":"ionindex(gcms::AbstractGCMS, ion::Real) -> Int\n\nReturn the index of the ion. If the ion does not exist, an error is thrown.\n\nSee also AbstractGCMS, ions, ion, minion,  maxion, ioncount.\n\nExample\n\njulia> gcms = GCMS((1:3)u\"s\", [85.2f0, 100.1f0], [0 12; 34 956; 23 1]);\n\njulia> ionindex(gcms, 100.1)\n2\n\njulia> ionindex(gcms, 201.1)\nERROR: ArgumentError: ion 201.1 does not exist\n[...]\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.ions","page":"Basics","title":"JuChrom.ions","text":"ions(gcms::AbstractGCMS)\n\nReturn the ions.\n\nSee also AbstractGCMS, GCMS, ion, minion,  maxion, ionindex, ioncount.\n\nIn the following example, the number type of ions passed to the object constructor  is explicitly annotated to illustrate that the GCMS object preserves the type.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", Int64[85, 100], [0 12; 34 956; 23 1]);\n\njulia> ions(gcms)\n2-element Vector{Int64}:\n  85\n 100\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.maxintensity","page":"Basics","title":"JuChrom.maxintensity","text":"maxintensity(chrom::AbstractChromatogram)\n\nReturn the maximum intensity.\n\nSee also AbstractChromatogram, AbstractGC, AbstractGCMS,  FID, GCMS, TIC, minintensity,  intensities, scantimes, scancount, ions,  ioncount.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> maxintensity(gcms)\n956\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.maxion","page":"Basics","title":"JuChrom.maxion","text":"maxion(gcms::AbstractGCMS)\n\nReturn the largest ion.\n\nSee also AbstractGCMS, GCMS, minion, ions,  ion, ioncount.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> maxion(gcms)\n100\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.maxscantime","page":"Basics","title":"JuChrom.maxscantime","text":"maxscantime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, \nustripped::Bool)\n\nReturn the time of the last scan. The optional keyword argument timeunit allows you to  change the unit of the returned scan time. All time units defined in the package  Unitful.jl (e.g., u\"s\", u\"minute\") are  supported. The optional keyword argument ustripped allows you to specify whether the unit  is stripped from the returned value.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  minscantime, scantimes, scancount.\n\nExample\n\njulia> gcms = GCMS([1.0, 2.0, 3.0]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> maxscantime(gcms)\n3.0 s\n\njulia> maxscantime(gcms, timeunit=u\"minute\")\n0.05 minute\n\njulia> maxscantime(gcms, timeunit=u\"minute\", ustripped=true)\n0.05\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.metadata","page":"Basics","title":"JuChrom.metadata","text":"metadata(chrom::AbstractChromatogram) -> Dict{Any, Any}\n\nReturn the metadata.\n\nSee also AbstractChromatogram, FID, GCMS, TIC.\n\nExample\n\njulia> gcms₁ = GCMS(Int64[1, 2]u\"s\", Int64[85, 100], Int64[0 12; 34 956], Dict(:id => 1))\nGCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n2 scans; time range: 1 s - 2 s\n2 ions; range: m/z 85 - 100\nintensity range: 0 - 956\nmetadata: 1\n\njulia> metadata(gcms₁)\nDict{Any, Any} with 1 entry:\n  :id => 1\n\njulia> gcms₂ = GCMS(Int64[1, 2]u\"s\", Int64[85, 100], Int64[0 12; 34 956])\nGCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n2 scans; time range: 1 s - 2 s\n2 ions; range: m/z 85 - 100\nintensity range: 0 - 956\nmetadata: 0\n\njulia> metadata(gcms₂)\nDict{Any, Any}()\n\njulia> metadata(gcms₂)[\"name\"] = \"sample\"\n\"sample\"\n\njulia> metadata(gcms₂)[:id] = 123\n123\n\njulia> metadata(gcms₂)\nDict{Any, Any} with 2 entries:\n  \"name\" => \"sample\"\n  :id    => 123\n\njulia> gcms₂\nGCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n2 scans; time range: 1 s - 2 s\n2 ions; range: m/z 85 - 100\nintensity range: 0 - 956\nmetadata: 2\n\njulia> delete!(metadata(gcms₂), \"name\")\nDict{Any, Any} with 1 entry:\n  :id => 123\n\njulia> gcms₂\nGCMS {scantimes: Int64, ions: Int64, intensities: Int64}\n2 scans; time range: 1 s - 2 s\n2 ions; range: m/z 85 - 100\nintensity range: 0 - 956\nmetadata: 1\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.minintensity","page":"Basics","title":"JuChrom.minintensity","text":"minintensity(chrom::AbstractChromatogram)\n\nReturn the minimum intensity.\n\nSee also AbstractChromatogram, AbstractGC, AbstractGCMS,  FID, GCMS, TIC, maxintensity,  intensities, scantimes, scancount, ions,  ioncount.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> minintensity(gcms)\n0\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.minion","page":"Basics","title":"JuChrom.minion","text":"minion(gcms::AbstractGCMS)\n\nReturn the smallest ion.\n\nSee also AbstractGCMS, GCMS, maxion, ions,  ion, ioncount.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> minion(gcms)\n85\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.minscantime","page":"Basics","title":"JuChrom.minscantime","text":"minscantime(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, \nustripped::Bool=false)\n\nReturn the time of the first scan. The optional keyword argument timeunit allows you to  change the unit of the returned scan time. All time units defined in the package  Unitful.jl (e.g., u\"s\", u\"minute\") are  supported. The optional keyword argument ustripped allows you to specify whether the unit  is stripped from the returned value.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  maxscantime, scantimes, scancount.\n\nExample\n\njulia> gcms = GCMS([1.0, 2.0, 3.0]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> minscantime(gcms)\n1.0 s\n\njulia> minscantime(gcms, timeunit=u\"minute\")\n0.016666666666666666 minute\n\njulia> minscantime(gcms, timeunit=u\"minute\", ustripped=true)\n0.016666666666666666\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.runduration","page":"Basics","title":"JuChrom.runduration","text":"runduration(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, \nustripped::Bool=false)\n\nReturn the runduration. The optional keyword argument timeunit allows you to change the  unit of the returned time interval. All time units defined in the package  Unitful.jl(e.g., u\"s\", u\"minute\") are  supported. The optional keyword argument ustripped allows you to specify whether the unit  is stripped from the returned value. \n\nSee also AbstractChromatogram, scantimes, minscantime,  maxscantime, scancount.\n\nExample\n\njulia> fid = FID([30.1u\"minute\", 40.8u\"minute\", 51.5u\"minute\"], [12, 956, 23])\nFID {scantimes: Float64, intensities: Int64}\n3 scans; time range: 30.1 minute - 51.5 minute\nintensity range: 12 - 956\nmetadata: 0\n\njulia> runduration(fid)\n21.4 minute\n\njulia> runduration(fid, timeunit=u\"s\")\n1284.0 s\n\njulia> runduration(fid, timeunit=u\"s\", ustripped=true)\n1284.0\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.scancount","page":"Basics","title":"JuChrom.scancount","text":"scancount(chrom::AbstractChromatogram) -> Int\n\nReturn the number of scans.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  scantimes, minscantime, maxscantime.\n\nExample\n\njulia> fid = FID([1, 2, 3]u\"s\", [12, 956, 23]);\n\njulia> scancount(fid)\n3\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.scantime","page":"Basics","title":"JuChrom.scantime","text":"scantime(chrom::AbstractChromatogram, index::Integer; timeunit::Unitful.TimeUnits, \nustripped::Bool=false)\n\nReturn the scantime by specifying the scan index. The optional parameter  timeunit allows you to change the unit of the returned time value. All time units defined  in the package Unitful.jl (e.g., u\"s\",  u\"minute\") are supported. The optional keyword argument ustripped allows you to specify  whether the unit is stripped from the returned value.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  scantimes, minscantime, maxscantime, scancount.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> scantime(gcms, 2)\n2 s\n\njulia> scantime(gcms, 2, timeunit=u\"minute\")\n1//30 minute\n\njulia> scantime(gcms, 2, timeunit=u\"minute\", ustripped=true)\n1//30\n\njulia> scantime(gcms, 5, timeunit=u\"minute\", ustripped=true)\nERROR: BoundsError: attempt to access 3-element Vector{Quantity{Int64, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}} at index [5]\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.scantimeindex","page":"Basics","title":"JuChrom.scantimeindex","text":"scantimeindex(gcms::AbstractGCMS, time::Unitful.Time; precisetime::Bool=false) -> Int\n\nReturn the index of the scantime in the scantimes that is closest to the given time.  All time units defined in the package  Unitful.jl (e.g., u\"s\", u\"minute\") are supported. If there is a tie, the larger scan time is used.  If the optional parameterprecisetimeis set totrue, the specified time must exist inscantimes`, otherwise an error is thrown.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  scantimes, scantime, minscantime, maxscantime,  scancount, ionindex.\n\nExample\n\njulia> gcms = GCMS([1.1f0, 2.1f0, 3.1f0]u\"s\", [85, 100], [0 12; 34 956; 23 1])\nGCMS {scantimes: Float32, ions: Int64, intensities: Int64}\n3 scans; time range: 1.1f0 s - 3.1f0 s\n2 ions; range: m/z 85 - 100\nintensity range: 0 - 956\nmetadata: 0\n\njulia> scantimeindex(gcms, 1.1f0u\"s\", precisetime=true)\n1\n\njulia> scantimeindex(gcms, 2.1u\"s\", precisetime=true)\n2\n\njulia> scantimeindex(gcms, 2.2u\"s\", precisetime=true)\nERROR: ArgumentError: scantime 2.2 s does not exist\n[...]\n\njulia> scantimeindex(gcms, 2.2u\"s\")\n2\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.scantimes","page":"Basics","title":"JuChrom.scantimes","text":"scantimes(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, \nustripped::Bool=false)\n\nReturn the scantimes. The optional keyword argument timeunit allows you to change the  unit of the returned scantimes. All time units defined in the package  Unitful.jl (e.g., u\"s\", u\"minute\") are  supported. The optional keyword argument ustripped allows you to specify whether the unit  is stripped from the returned values.\n\nSee also AbstractChromatogram, FID, GCMS, TIC,  scantime, minscantime, maxscantime, scancount.\n\nIn the following example, the element type of the scantimes vector passed to the object  constructor is explicitly annotated to illustrate that the GCMS object preserves the type.\n\nExample\n\njulia> gcms = GCMS(Float32[1.0, 2.0, 3.0]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> scantimes(gcms)\n3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:\n 1.0f0 s\n 2.0f0 s\n 3.0f0 s\n\njulia> scantimes(gcms, timeunit=u\"minute\")\n3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:\n 0.016666668f0 minute\n 0.033333335f0 minute\n 0.050000004f0 minute\n\njulia> scantimes(gcms, timeunit=u\"minute\", ustripped=true)\n3-element Vector{Float32}:\n 0.016666668\n 0.033333335\n 0.050000004\n\n\n\n\n\n","category":"function"},{"location":"man/basics/#JuChrom.totalionchromatogram","page":"Basics","title":"JuChrom.totalionchromatogram","text":"totalionchromatogram(gcms::GCMS)\n\nCompute the total ion chromatrogram.\n\nSee also AbstractChromatogram, AbstractGC, GCMS,  AbstractTIC, TIC, scantimes, intensities,  metadata.\n\nIn the following example, the element type of intensities passed to the object  constructor  is explicitly annotated to show that the TIC object has the same type.\n\nExample\n\njulia> gcms = GCMS((1.1:3.1)u\"s\", [85, 100], Int64[0 12; 34 956; 23 1]);\n\njulia> tic = totalionchromatogram(gcms)\nTIC {scantimes: Float64, intensities: Int64}\n3 scans; time range: 1.1 s - 3.1 s\nintensity range: 12 - 990\nmetadata: 0\n\njulia> intensities(tic)\n3-element Vector{Int64}:\n  12\n 990\n  24\n\njulia> gcms = GCMS((1.1:3.1)u\"s\", [85, 100], Float64[0 12; 34 956; 23 1]);\n\njulia> tic = totalionchromatogram(gcms)\nTIC {scantimes: Float64, intensities: Float64}\n3 scans; time range: 1.1 s - 3.1 s\nintensity range: 12.0 - 990.0\nmetadata: 0\n\njulia> intensities(tic)\n3-element Vector{Float64}:\n  12.0\n 990.0\n  24.0\n\n\n\n\n\n","category":"function"},{"location":"man/export/#Data-export","page":"Data export","title":"Data export","text":"","category":"section"},{"location":"man/internals/#Internals","page":"Internals","title":"Internals","text":"","category":"section"},{"location":"man/internals/","page":"Internals","title":"Internals","text":"JuChrom.findclosest\nJuChrom.copy_with_eltype","category":"page"},{"location":"man/internals/#JuChrom.findclosest","page":"Internals","title":"JuChrom.findclosest","text":"JuChrom.findclosest(A::AbstractVector{<:Number}, x::Number) -> Int\n\nReturn the index of the number closest to number x in a list of numbers sorted in ascending  order. If there is a tie, the index of the larger number is returned.\n\nExample\n\njulia> JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 0)\n3\n\njulia> JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], 1.5)\n5\n\njulia> JuChrom.findclosest([-2, -1, 0, 1, 2, 3, 4, 5], -1.5)\n2\n\n\n\n\n\n","category":"function"},{"location":"man/internals/#JuChrom.copy_with_eltype","page":"Internals","title":"JuChrom.copy_with_eltype","text":"JuChrom.copy_with_eltype(array::AbstractArray, elementtype::Type)\n\nCreate a mutable copy of the array with the type of its elements converted to  elementtype.\n\nExample\n\njulia> JuChrom.copy_with_eltype(Int[1, 2, 3, 4, 5, 6], Float64)\n6-element Vector{Float64}:\n 1.0\n 2.0\n 3.0\n 4.0\n 5.0\n 6.0\n\njulia> JuChrom.copy_with_eltype(Float64[1, 2, 3, 4, 5, 6], Int32)\n6-element Vector{Int32}:\n 1\n 2\n 3\n 4\n 5\n 6\n\njulia> JuChrom.copy_with_eltype(Float64[1.1, 2, 3, 4, 5, 6], Int32)\nERROR: InexactError: Int32(1.1)\n[...]\n\n\n\n\n\n","category":"function"},{"location":"man/deconvolution/#Deconvolution","page":"Deconvolution","title":"Deconvolution","text":"","category":"section"},{"location":"man/explorer/#GUI-explorer","page":"GUI explorer","title":"GUI explorer","text":"","category":"section"},{"location":"man/import/#Data-import","page":"Data import","title":"Data import","text":"","category":"section"},{"location":"#JuChrom.jl","page":"Home","title":"JuChrom.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for JuChrom.jl","category":"page"},{"location":"man/installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"}]
}
