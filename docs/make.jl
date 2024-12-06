using Documenter, JuChrom

DocMeta.setdocmeta!(JuChrom, :DocTestSetup, :(using JuChrom); recursive=true)

makedocs(
    authors = "Oliver Niehuis",
    sitename = "JuChrom.jl",
    format = Documenter.HTML(),
    modules = [JuChrom],
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "man/basics.md",
            "man/retentionindices.md",
            "man/massspectra.md",
            "man/import.md",
            "man/export.md"],
        "Internals" => Any[
            "internals/base.md",
            "internals/deconvolution.md",
            "internals/explorer.md",
            "internals/inputoutput.md"],
        "Tutorials" => Any[
            "tutorials/retentionindices.md"],
        "Index" => "man/register.md"
    ]
)

deploydocs(
    repo = "github.com/oniehuis/JuChrom.jl"
)
