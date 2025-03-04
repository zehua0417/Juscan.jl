using Documenter, Juscan

makedocs(
  sitename="Juscan.jl",
  modules=[Juscan],
  pages=[
    "Home" => "index.md",
    "Tutorial" => Any["Quick Start" => "tutorial/quickstart.md"],
    "API" => Any["Anndata Utils" => "api/anndata.md", "Preprocessing" => "api/preprocessing.md"],
  ],
)

if get(ENV, "CI", "false") == "true"
  deploydocs(repo="https://github.com/zehua0417/Juscan.jl.git")
end
