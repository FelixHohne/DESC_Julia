using DESC
using Documenter

DocMeta.setdocmeta!(DESC, :DocTestSetup, :(using DESC); recursive=true)

makedocs(;
    modules=[DESC],
    authors="Benjamin Faber <bfaber@wisc.edu> and contributors",
    repo="https://gitlab.com/wistell/DESC.jl/blob/{commit}{path}#{line}",
    sitename="DESC.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wistell.gitlab.io/DESC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
