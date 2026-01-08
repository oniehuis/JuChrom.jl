using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
if get(ENV, "CI", "false") == "true"
    Pkg.instantiate()
end

using Documenter, JuChrom
using PyCall
Base.retry_load_extensions()

DocMeta.setdocmeta!(JuChrom, :DocTestSetup, quote
    using JuChrom
    using SparseArrays
end; recursive = true)

withenv("UNITFUL_FANCY_EXPONENTS" => "false") do
    makedocs(
	    sitename = "JuChrom",
        format = Documenter.HTML(),
        modules = [JuChrom],
	    authors = "Oliver Niehuis",
        root = @__DIR__,
        pages = [
            "Home" => "index.md",
            "Manual" => Any[
                "man/Scans.md",
                "man/ScanSeries.md",
                "man/ScanMatrices.md",
                "man/loaders.md",
                "man/convert.md",
                "man/transform.md",
                "man/retention_mapper.md",
                "man/quadvar_model.md",
                "man/baseline.md",
                "man/alignment.md"
			    ],
            "Internals" => Any[
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
