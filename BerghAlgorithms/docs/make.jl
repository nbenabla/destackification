using Documenter

push!(LOAD_PATH,"../src/")

using BerghAlgorithms

DocMeta.setdocmeta!(BerghAlgorithms, :DocTestSetup, :(using BerghAlgorithms); recursive = true)


makedocs(
    modules=[BerghAlgorithms],
    sitename="Bergh Documentation",
    pages = [
        "Home" => "index.md",
        "Functions" => "functions.md",
        "Contributors" => "contributors.md"
    ]
)

deploydocs(;
    repo="github.com/nbenabla/destackification",
)