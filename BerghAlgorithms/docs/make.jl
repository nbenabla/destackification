using Documenter

push!(LOAD_PATH,"../src/")

using BerghAlgorithms

makedocs(
    modules=[BerghAlgorithms],
    sitename="Bergh Documentation",
    pages = [
        "Home" => "index.md",
        "Functions" => "functions.md"
    ]
)

deploydocs(;
    repo="github.com/nbenabla/destackification",
)