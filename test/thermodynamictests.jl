@testset "CompressibilityFactor" begin
    p = Parameters(0.4, 1.4737, 7^3, 1e-5, 393216, Float64)
    s = SimulationSystem(p, 3, Float64)
    potential = PseudoHS(50.0, p.kT)
    initialize!(s)
    move!(1000, s, potential; interval = 500)
    zfac = ZFactor(0.0, 0.0)
    move!(1000, s, potential, zfac; interval = 500)
    @show zfac.zval

    # Check consistency
    @test zfac.naverage == 1000.0
    # Check value, NOT physically representative
    @test zfac.zval ≈ 8.24498784197415
end
