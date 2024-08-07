var documenterSearchIndex = {"docs":
[{"location":"#JuChrom.jl","page":"JuChrom.jl","title":"JuChrom.jl","text":"","category":"section"},{"location":"","page":"JuChrom.jl","title":"JuChrom.jl","text":"Documentation for JuChrom.jl","category":"page"},{"location":"","page":"JuChrom.jl","title":"JuChrom.jl","text":"Modules = [JuChrom]\nOrder   = [:type, :function]","category":"page"},{"location":"#JuChrom.AbstractChromatogram","page":"JuChrom.jl","title":"JuChrom.AbstractChromatogram","text":"AbstractChromatogram\n\nThe AbstractChromatogram type is the supertype of all chromatogram implementations in  JuChrom. All subtypes of AbstractChromatogram (e.g., FID, GCMS, TIC) contain  scantimes, intensities, and source information.\n\nSee also AbstractGCMS.\n\n\n\n\n\n","category":"type"},{"location":"#JuChrom.AbstractGCMS","page":"JuChrom.jl","title":"JuChrom.AbstractGCMS","text":"AbstractGCMS <: AbstractChromatogram\n\nThe AbstractGCMS type is the supertype of all chromatogram implementations that contain  mass-charge ratio (m/z) data (= ions) and associated abundance values (= intensities)  and thus can have one or more ion intensity values associated with a given scan time  in JuChrom (e.g., GCMS).\n\nSee also AbstractChromatogram, GCMS.\n\n\n\n\n\n","category":"type"},{"location":"#JuChrom.GCMS-Union{Tuple{T4}, Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{T1, T2, T3}, Tuple{T1, T2, T3, T4}} where {T1<:(AbstractVector{var\"#s3\"} where var\"#s3\"<:(Union{Quantity{T, 𝐓, U}, Level{L, S, Quantity{T, 𝐓, U}} where {L, S}} where {T, U})), T2<:(AbstractVector{var\"#s11\"} where var\"#s11\"<:Real), T3<:(AbstractMatrix{var\"#s12\"} where var\"#s12\"<:Real), T4<:Union{Nothing, AbstractString}}","page":"JuChrom.jl","title":"JuChrom.GCMS","text":"GCMS(scantimes::AbstractVector{<:Unitful.Time}, ions::AbstractVector{<:Real}, \nintensities::AbstractMatrix{<:Real}; source::Union{Nothing, AbstractString})\n\nReturn a GCMS object consisting of scantimes, ions, intensities, and optionally data  source information.\n\nSee also AbstractGCMS, intensities, ions, scantimes.\n\nIn the following examples, the number types of the arrays passed to the object constructor are explicitly annotated to illustrate that the GCMS object preserves the types.\n\nExample\n\njulia> GCMS(Int32[1, 2, 3]u\"s\", Int64[85, 100], Int32[0 12; 34 956; 23 1])\nGCMS {scantimes: Int32, ions: Int64, intensities: Int32}\nSource: nothing\n3 scans; time range: 1 s - 3 s\n2 ions; range: m/z 85 - 100\nintensity range: 0 - 956\n\njulia> GCMS([1.1f0, 2.1f0]u\"s\", Float32[35.1, 76.2], Int64[0 12; 34 956], \"example data\")\nGCMS {scantimes: Float32, ions: Float32, intensities: Int64}\nSource: example data\n2 scans; time range: 1.1f0 s - 2.1f0 s\n2 ions; range: m/z 35.1 - 76.2\nintensity range: 0 - 956\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.intensities-Tuple{AbstractChromatogram}","page":"JuChrom.jl","title":"JuChrom.intensities","text":"intensities(chrom::AbstractChromatogram)\n\nReturn the intensities of an AbstractChromatogram subtype object.\n\nSee also AbstractChromatogram.\n\nExample\n\njulia> gcms = GCMS([1u\"s\", 2u\"s\", 3u\"s\"], [85, 100], Int64[0 12; 34 956; 23 1]);\n\njulia> intensities(gcms)\n3×2 Matrix{Int64}:\n  0   12\n 34  956\n 23    1\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.ions-Tuple{AbstractGCMS}","page":"JuChrom.jl","title":"JuChrom.ions","text":"ions(gcms::AbstractGCMS)\n\nReturn the ions of an AbstractGCMS subtype object (e.g., GCMS`).\n\nSee also AbstractGCMS.\n\nIn the following example, the number type of ions passed to the object  constructor  is explicitly annotated to illustrate that the GCMS object preserves the type.\n\nExample\n\njulia> gcms = GCMS([1, 2, 3]u\"s\", Int64[85, 100], [0 12; 34 956; 23 1]);\n\njulia> ions(gcms)\n2-element Vector{Int64}:\n  85\n 100\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.scantimes-Tuple{AbstractChromatogram}","page":"JuChrom.jl","title":"JuChrom.scantimes","text":"scantimes(chrom::AbstractChromatogram; timeunit::Unitful.TimeUnits, \nustripped::Bool=false)\n\nReturn the scantimes of an AbstractChromatogram subtype object (e.g., FID, GCMS,  TIC). The  optional keyword  argument timeunit allows you to change the unit of the  returned  scantimes. All time  units defined in the package Unitful.jl  (e.g., u\"s\", u\"minute\") are supported. The optional keyword argument ustripped  allows  you to specify whether the unit is stripped from the returned values. \n\nSee also AbstractChromatogram.\n\nIn the following example, the element type of the scantime vector passed to the object  constructor is explicitly annotated to illustrate that the GCMS object preserves the type.\n\nExample\n\njulia> gcms = GCMS(Float32[1.0, 2.0, 3.0]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> scantimes(gcms)\n3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(s,), 𝐓, nothing}}}:\n 1.0f0 s\n 2.0f0 s\n 3.0f0 s\n\njulia> scantimes(gcms, timeunit=u\"minute\")\n3-element Vector{Quantity{Float32, 𝐓, Unitful.FreeUnits{(minute,), 𝐓, nothing}}}:\n 0.016666668f0 minute\n 0.033333335f0 minute\n 0.050000004f0 minute\n\njulia> scantimes(gcms, timeunit=u\"minute\", ustripped=true)\n3-element Vector{Float32}:\n 0.016666668\n 0.033333335\n 0.050000004\n\n\n\n\n\n","category":"method"},{"location":"#JuChrom.source-Tuple{AbstractChromatogram}","page":"JuChrom.jl","title":"JuChrom.source","text":"source(chrom::AbstractChromatogram) -> Type{<:Union{AbstractString, Nothing}}\n\nReturn the source information of the AbstractChromatogram subtype object.\n\nSee also AbstractChromatogram.\n\nExample\n\njulia> gcms₁ = GCMS([1, 2, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1], \"example data\");\n\njulia> source(gcms₁)\n\"example data\"\n\njulia> gcms₂ = GCMS([1, 2, 3]u\"s\", [85, 100], [0 12; 34 956; 23 1]);\n\njulia> source(gcms₂)\n\njulia> isnothing(source(gcms₂))\ntrue\n\n\n\n\n\n","category":"method"}]
}
