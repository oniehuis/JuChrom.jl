using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
if get(ENV, "CI", "false") == "true"
    Pkg.instantiate()
end

using Documenter
using PyCall
using JuChrom
Base.retry_load_extensions()
let ext = Base.get_extension(JuChrom, :PyCallExtension)
    if ext === nothing
        @info "PyCallExtension not loaded; forcing load for docs"
        include(joinpath(@__DIR__, "..", "ext", "PyCallExtension.jl"))
        ext = PyCallExtension
    end
    if !isdefined(JuChrom, :ShimadzuMSLoader)
        @eval JuChrom const ShimadzuMSLoader = $(ext.ShimadzuMSLoader)
    end
end
if isdefined(JuChrom, :ShimadzuMSLoader)
    using JuChrom.ShimadzuMSLoader
end

DocMeta.setdocmeta!(JuChrom, :DocTestSetup, quote
    using JuChrom
    using SparseArrays
end; recursive = true)

withenv("UNITFUL_FANCY_EXPONENTS" => "false") do
    makedocs(
	    sitename = "JuChrom.jl",
        format = Documenter.HTML(),
        modules = isdefined(JuChrom, :ShimadzuMSLoader) ? [JuChrom, JuChrom.ShimadzuMSLoader] : [JuChrom],
	    authors = "Oliver Niehuis",
        root = @__DIR__,
        pages = [
            "Home" => "index.md",
            "Manual" => Any[
                "Core containers" => Any[
                    "Design overview" => "man/containers.md",
                    "man/scans.md",
                    "man/scanseries.md",
                    "man/scanmatrices.md",
                ],
                "Units and Quantities" => "man/units_quantities.md",
                "Loaders" => Any[
                    "Loader API" => "man/loaderapi.md",
                    "man/agilentfid.md",
                    "man/chemstationms.md",
                    "man/masshunterms.md",
                    "man/shimadzums.md",
                ],
                "Data processing" => Any[
                    "Alignment" => "man/alignment.md",
                    "Baseline" => "man/baseline.md",
                    "Binning and Gridding" => Any[
                        "m/z binning" => "man/mz_binning.md",
                        "Retention binning and gridding" => "man/binning_gridding.md",
                    ],
                    "Deconvolution" => "man/deconvolution.md",
                    "Format conversion" => "man/convert.md",
                    "Retention mapping" => Any[
                        "Overview" => "man/mapping_overview.md",
                        "Mapping tools" => "man/mapping_tools.md",
                    ],
                    "Scan timing" => "man/scan_timing.md",
                    "Transformation" => "man/transformation.md",
                    "Trimming and Filtering" => "man/transform.md",
                    "Variance modeling" => "man/quadvar_model.md",
                ]
			    ],
            "Internals" => Any[
                "internals/deconvolution.md",
                "internals/transform.md",
                "internals/utils.md"
		        ],
            "Index" => "man/register.md"
        ]
    )
end

deploydocs(
	repo = "github.com/oniehuis/JuChrom.jl"
)
