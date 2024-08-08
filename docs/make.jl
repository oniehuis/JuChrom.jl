using Documenter, JuChrom

DocMeta.setdocmeta!(JuChrom, :DocTestSetup, :(using JuChrom); recursive=true)

makedocs(
    sitename = "JuChrom.jl",
    format = Documenter.HTML(),
    modules = [JuChrom]
)

deploydocs(
    repo = "github.com/oniehuis/JuChrom.jl"
)
