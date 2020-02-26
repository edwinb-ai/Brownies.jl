function _create_rngs(num_rngs::Integer; seed::Integer = nothing)
    # Create the RNG to create seeds
    if isnothing(seed)
        # Without a given seed
        rng_master = PCG.PCGStateOneseq()
    else
        # With a given seed
        rng_master = PCG.PCGStateOneseq(seed)
    end
    # From this RNG, create `num_rngs` seeds
    seed_list = zeros(UInt64, num_rngs)
    rand!(rng_master, seed_list)
    # With these seeds, seed the new RNG's
    rng_list = map(Xorshifts.Xoroshiro128Plus, seed_list)

    return (seed_list, rng_list)
end

function rng_matrix!(rnd_matrix::AbstractArray, rng_list::AbstractArray)
    rnd_type = eltype(rnd_matrix)
    for i in axes(rnd_matrix, 2)
        rnd_matrix[:, i] = randn(rng_list[i], rnd_type, size(rnd_matrix, 1))
    end
end

function savetofile(s::SimulationSystem, energies)
    # Save positions and forces
    @unpack positions, forces = s.system
    @save "positions-$(s.params.ϕ)-$(s.params.N).jld2" positions
    @save "forces-$(s.params.ϕ)-$(s.params.N).jld2" forces
    # Save the computed energies as well
    @save "energy-$(s.params.ϕ)-$(s.params.N).jld2" energies
end
