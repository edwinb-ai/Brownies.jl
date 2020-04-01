@testset "Utils" begin
    # Use a random seed
    p = Parameters(0.4, 1.4737, 7^3, 1e-5, Float64)
    @test p.seed == 0

    # And change the diameter of the particles
    s = SimulationSystem(p, 3, Float64)
    potential = PseudoHS(50.0, 1.5, p.kT)
    @test potential.σ == 1.5
    # And just check that everything works as expected
    initialize!(s)
    # As well as particles can be moved
    move!(1000, s, potential; interval = 500)
end

@testset "DataFrames" begin
    p = Parameters(0.4, 1.4737, 7^3, 1e-5, Float64)
    s = SimulationSystem(p, 3, Float64)
    potential = PseudoHS(50.0, p.kT)
    initialize!(s)

    move!(1000, s, potential; interval = 500)
    # Save positions and forces to DataFrames
    posdf, forcesdf = todataframe(s.system.positions, s.system.forces)

    # Check consistency for positions
    @test posdf.x == s.system.positions[:, 1]
    @test posdf.y == s.system.positions[:, 2]
    @test posdf.z == s.system.positions[:, 3]

    # Check consistency for forces
    @test forcesdf.x == s.system.forces[:, 1]
    @test forcesdf.y == s.system.forces[:, 2]
    @test forcesdf.z == s.system.forces[:, 3]
end
