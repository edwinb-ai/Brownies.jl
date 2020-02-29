using UnPack

@testset "Parameters" begin
    # First constructor
    p = Parameters(0.4, 1.4737, 1e-5, Float32)
    @unpack ϕ, kT, N, τ, seed = p
    @test ϕ == 0.4f0
    @test kT == 1.4737f0
    @test N == 512
    @test τ == 0.00001f0
    @test seed == 393216

    # Second constructor, define number of particles
    num_particles = 256
    p1 = Parameters(0.4, 1.4737, num_particles, 1e-5, Float32)
    @unpack ϕ, kT, N, τ, seed = p1
    @test ϕ == 0.4f0
    @test kT == 1.4737f0
    @test N == num_particles
    @test τ == 0.00001f0
    @test seed == 393216

    # Third constructor, define number of particles and random seed
    num_particles = 256
    new_seed = 804533
    p2 = Parameters(0.4, 1.4737, num_particles, 1e-5, new_seed, Float32)
    @unpack ϕ, kT, N, τ, seed = p2
    @test ϕ == 0.4f0
    @test kT == 1.4737f0
    @test N == num_particles
    @test τ == 0.00001f0
    @test seed == new_seed
end
@testset "SimulationSystem" begin
    # Create the parameters and simulation system from it
    p = Parameters(0.4, 1.4737, 1e-5, Float64)
    # 3 is the dimension of the system, a 3D simulation system
    dims = 3
    simul_sys = SimulationSystem(p, dims, Float64)
    @unpack params, ρ, boxl, rc, system, energy = simul_sys
    # Check the parameters
    @unpack ϕ, kT, N, τ, seed = p
    @test ϕ == 0.4
    @test kT == 1.4737
    @test N == 512
    @test τ == 0.00001
    @test seed == 393216
    # Now check the system parameters
    density = 6.0 * params.ϕ / π
    @test density ≈ ρ
    @test boxl == cbrt(params.N / density)
    @test rc == boxl * 0.5
    @test energy == 0.0
    # Check that the positions and forces are correctly initialized
    total_size = (512, dims)
    @test size(system.positions) == total_size
    @test size(system.forces) == total_size
    @test all(system.positions .== 0.0)
    @test all(system.forces .== 0.0)
end
