using Documenter, JuChrom

DocMeta.setdocmeta!(JuChrom, :DocTestSetup, :(using JuChrom); recursive=true)

makedocs(
	sitename = "JuChrom",
	format = Documenter.HTML(),
	modules = [JuChrom]
)

deploydocs(
	repo = "github.com/oniehuis/JuChrom.jl"
)