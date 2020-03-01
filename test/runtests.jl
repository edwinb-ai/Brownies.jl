using Brownies
using Test

@testset "Brownies.jl" begin
    include("typestest.jl")
    # TODO: potentials.jl
    # TODO: energy.jl
    # TODO: utils.jl
    # TODO: move.jl
    # TODO: algorithms/ermak.jl
    # TODO: observables/rdf.jl
    include("structuretests.jl")
    include("thermodynamictests.jl")
    include("dynamicstests.jl")
end
