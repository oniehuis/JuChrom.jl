using Documenter, JuChrom

DocMeta.setdocmeta!(JuChrom, :DocTestSetup, :(using JuChrom); recursive=true)

makedocs(
	sitename = "JuChrom",
    format = Documenter.HTML(),
    modules = [JuChrom],
	authors = "Oliver Niehuis",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "man/basics.md",
            "man/retentionindices.md",
            "man/massspectra.md",
            "man/import.md",
            "man/export.md"
			],
        "Internals" => Any[
            "internals/base.md",
            # "internals/deconvolution.md",
            "internals/explorer.md",
            "internals/inputoutput.md"
			],
        "Tutorials" => Any[
            "tutorials/retentionindices.md"],
        "Index" => "man/register.md"
    ]
)


deploydocs(
	repo = "github.com/oniehuis/JuChrom.jl"
)