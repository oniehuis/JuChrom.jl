using Documenter
using JuChrom

makedocs(
    sitename = "JuChrom",
    format = Documenter.HTML(),
    modules = [JuChrom]
)

deploydocs(
    repo = "github.com/oniehuis/JuChrom.jl"
)
