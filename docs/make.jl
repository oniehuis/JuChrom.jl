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
            "Home" => "Index.md",
            "Manual" => Any[
                "Core containers" => Any[
                    "Design overview" => "man/Containers.md",
                    "man/Scans.md",
                    "man/ScanSeries.md",
                    "man/ScanMatrices.md",
                ],
                "Loaders" => Any[
                    "Loader API" => "man/LoaderAPI.md",
                    "man/AgilentFID.md",
                    "man/ChemStationMS.md",
                    "man/MassHunterMS.md",
                    "man/ShimadzuMS.md",
                ],
                "Data processing" => Any[
                    "Retention mapping" => Any[
                        "Mapping overview" => "man/Mapping_overview.md",
                        "Mapping tools" => "man/Mapping_tools.md",
                    ],
                    "Alignment" => "man/Alignment.md",
                    "Baseline" => "man/Baseline.md",
                    "Binning and Gridding" => "man/Binning_gridding.md",
                    "Format conversion" => "man/Convert.md",
                    "Transformation" => "man/Transformation.md",
                    "Trimming and Filtering" => "man/Transform.md",
                    "Variance modeling" => "man/Quadvar_model.md",
                ]
			    ],
            "Internals" => Any[
                "internals/Transform.md",
                "internals/Utils.md"
		        ],
            "Index" => "man/Register.md"
        ]
    )
end

deploydocs(
	repo = "github.com/oniehuis/JuChrom.jl"
)
