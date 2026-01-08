using Documenter, JuChrom
using PyCall
Base.retry_load_extensions()

DocMeta.setdocmeta!(JuChrom, :DocTestSetup, quote
    using JuChrom
    using SparseArrays
end; recursive = true)

makedocs(
	sitename = "JuChrom",
    format = Documenter.HTML(),
    modules = [JuChrom],
	authors = "Oliver Niehuis",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "man/scans.md",
            "man/scanseries.md",
            "man/scanmatrices.md",
            "man/convert.md",
            "man/loader.md",
            "man/jld2.md",
            "man/retention_mapper.md",
            "man/baseline.md",
            "man/quadvar_model.md",
            "man/alignment.md",
            "man/transform.md"
			],
        "Internals" => Any[
            "internals/transform.md",
            "internals/utils.md"
		    ],
        # "Index" => "man/register.md"
    ]
)


deploydocs(
	repo = "github.com/oniehuis/JuChrom.jl"
)
