using Documenter, Juscan

if haskey(ENV, "DOCSARGS")
  for arg in split(ENV["DOCSARGS"])
    (arg in ARGS) || push!(ARGS, arg)
  end
end

makedocs(
  sitename="Juscan.jl",
  authors="lihuax(LiZehua)",
  modules=[Juscan],
  format=if "pdf" in ARGS
    Documenter.LaTeX()
  else
    Documenter.HTML(
      assets=["assets/favicon.ico"],
      footer="Copyright © lihuax (LiZehua) Juscan.jl 2025",
    )
  end,
  pages=[
    "Home" => "index.md",
    "Tutorial" => Any["Quick Start" => "tutorial/quickstart.md"],
    "API" => Any[
      "Anndata Utils" => "api/anndata.md",
      "Preprocessing" => "api/preprocessing.md",
      "Tools" => "api/tools.md",
    ],
  ],
)

if get(ENV, "CI", "false") == "true"
  deploydocs(repo="https://github.com/zehua0417/Juscan.jl.git")
end
