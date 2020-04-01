using Statistics: mean

@testset "MeanSquareDisplacement" begin
    p = Parameters(0.4, 1.4737, 7^3, 1e-5, 393216, Float64)
    s = SimulationSystem(p, 3, Float64)
    potential = PseudoHS(50.0, p.kT)
    initialize!(s)
    move!(1000, s, potential; interval = 100)
    msd = MeanSquaredDisplacement(s, 500, 10)
    move!(5000, s, potential, msd; interval = 500)

    # Check that the value is "almost" the self-diffusion coefficient
    # which is 6 for a 3D system
    true_diffusion = 6.0
    slope = (msd.wt[2, 2] - msd.wt[1, 2]) / (msd.wt[2, 1] - msd.wt[1, 1])
    relative_error = abs(true_diffusion - slope) / true_diffusion
    # Only accept 15% deviation
    @test relative_error <= 0.15
end

@testset "IntermediateSelfScatteringFunction" begin
    # Create a simulation system
    p = Parameters(0.4, 1.4737, 7^3, 1e-5, 393216, Float64)
    s = SimulationSystem(p, 3, Float64)
    potential = PseudoHS(50.0, p.kT)
    # Initialize it and equilibrate it
    initialize!(s)
    move!(1000, s, potential; interval = 100)

    # Define the MeanSquareDisplacement and SelfScatteringFunction objects
    msd = MeanSquaredDisplacement(s, 500, 10)
    sft = SelfScatteringFunction(s, 500, 10, 6.6)
    # Move the particles around to obtain useful statistics
    move!(5000, s, potential, msd, sft; interval = 500)

    # Check that the starting values are close to one, as expected from
    # the scattering function
    mean_value = mean(sft.dft[1:20, 2]) # These are closer to 0.9 actually, numerical issues
    @show mean_value
    true_value = 0.99
    relative_error = (true_value - mean_value) / true_value
    # Only accept 15% deviation
    @test relative_error <= 0.15
end
