using UnPack

@testset "Parameters" begin
    p = Parameters(0.4, 1.4737, 1e-5, Float32)
    @unpack ϕ, kT, N, τ, seed = p
    @test ϕ == 0.4f0
    @test kT == 1.4737f0
    @test N == 512
    @test τ == 0.00001f0
    @test seed == 393216
end
