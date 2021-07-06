using IndirectTrajOpt
using Documenter

DocMeta.setdocmeta!(IndirectTrajOpt, :DocTestSetup, :(using IndirectTrajOpt); recursive=true)

makedocs(;
    modules=[IndirectTrajOpt],
    authors="Grant Hecht",
    repo="https://github.com/GrantHecht/IndirectTrajOpt.jl/blob/{commit}{path}#{line}",
    sitename="IndirectTrajOpt.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GrantHecht.github.io/IndirectTrajOpt.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GrantHecht/IndirectTrajOpt.jl",
    devbranch="main",
)
