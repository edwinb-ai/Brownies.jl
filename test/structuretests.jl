using Statistics

@testset "StructureFactor" begin
    p = Parameters(0.4, 1.4737, 7^3, 1e-5, Float64)
    s = SimulationSystem(p, 3, Float64)
    potential = PseudoHS(50.0, p.kT)
    initialize!(s)
    move!(10, s, potential; interval = 10)
    gfunc = PairDistributionFunction(s, 2^9)
    move!(3, s, potential, gfunc; interval = 1)
    sqfactor = StructureFactor(s, gfunc)
    compute_structure!(sqfactor, gfunc, s)

    # Should be "almost" 1
    @test isapprox(mean(sqfactor.sq[:, 2]), 1.0; atol=1e-1)
    # Check spacing
    dq = π / s.rc
    @test sqfactor.sq[2, 1] == dq
    # Check the normalizing constant
    const_val = 4.0 * π * s.ρ
    @test sqfactor.norm_const == const_val
    # And just check consistency
    @test sqfactor.nm == gfunc.nm
end
