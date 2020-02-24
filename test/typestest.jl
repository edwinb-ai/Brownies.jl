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
