using Documenter, JuChrom

DocMeta.setdocmeta!(JuChrom, :DocTestSetup, :(using JuChrom); recursive=true)

makedocs(
    authors = "Oliver Niehuis",
    sitename = "JuChrom.jl",
    format = Documenter.HTML(),
    modules = [JuChrom],
    pages = [
        "Home" => "index.md",
        "Manual" => Any["man/installation.md",
            "man/basics.md",
            "man/import.md",
            "man/export.md",
            "man/explorer.md",
            "man/deconvolution.md"],
        "Index" => "man/register.md"
    ]
)

deploydocs(
    repo = "github.com/oniehuis/JuChrom.jl"
)
