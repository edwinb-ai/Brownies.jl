using Documenter, Brownies

makedocs(;
    modules=[Brownies],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/edwinb-ai/Brownies.jl/blob/{commit}{path}#L{line}",
    sitename="Brownies.jl",
    authors="Edwin Bedolla, University of Guanajuato",
    assets=String[],
)

deploydocs(;
    repo="github.com/edwinb-ai/Brownies.jl",
)
