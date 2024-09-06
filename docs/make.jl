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
            "man/import.md",
            "man/export.md"],
        "Tutorials" => Any[
            "tutorials/retentionindices.md"],
        "Developers" => Any["man/internals.md"],
        "Index" => "man/register.md"
    ]
)

deploydocs(
    repo = "github.com/oniehuis/JuChrom.jl"
)
