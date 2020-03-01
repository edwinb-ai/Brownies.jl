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
    slope = (msd.wt[2, 2] - msd.wt[1, 2]) / (msd.wt[2, 1] - msd.wt[1, 1])
    relative_error = abs(6.0 - slope) / 6
    # Only accept 15% deviation
    @test relative_error <= 0.15
end
