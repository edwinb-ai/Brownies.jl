using Brownies
using Test

@testset "Brownies.jl" begin
    include("typestest.jl")
    include("structuretests.jl")
    include("thermodynamictests.jl")
    include("dynamicstests.jl")
    include("utilstests.jl")
end
