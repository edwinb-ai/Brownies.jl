@testset "CompressibilityFactor" begin
    p = Parameters(0.4, 1.4737, 7^3, 1e-5, 393216, Float64)
    s = SimulationSystem(p, 3, Float64)
    potential = PseudoHS(50.0, p.kT)
    initialize!(s)
    move!(1000, s, potential; interval = 500)
    zfac = ZFactor(0.0, 0.0)
    move!(1000, s, potential, zfac; interval = 500)

    # Check consistency
    @test zfac.naverage == 1000.0
    # Compute the Carnahan-Starling value
    true_zcomp = 8.24498784197415
    # Check value, larger means that it needs more steps (which is the case)
    @test zfac.zval ≈ true_zcomp
end

@testset "CompressibilityFactorRDF" begin
    p = Parameters(0.4, 1.4737, 7^3, 1e-5, 393216, Float64)
    s = SimulationSystem(p, 3, Float64)
    potential = PseudoHS(50.0, p.kT)
    initialize!(s)
    move!(1000, s, potential; interval = 500)

    gfunc = PairDistributionFunction(s, 2^9)
    zfac = ZFactor(0.0, 0.0)
    move!(1000, s, potential, gfunc, zfac; interval = 500)

    # Check that the RDF is physically consistent
    rdf_values = @view gfunc.gofr[:, 2]
    half_idx = fld1(length(rdf_values), 2)
    @test isapprox(mean(rdf_values[half_idx:end]), 1.0; atol=1e-1)

    # Check consistency
    @test zfac.naverage == 1000.0
    # Compute the Carnahan-Starling value
    true_zcomp = 8.24498784197415
    # Check value, larger means that it needs more steps (which is the case)
    @test zfac.zval ≈ true_zcomp
end
