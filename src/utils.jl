function _create_rngs(num_rngs::Integer; seed::Integer = 0)
    # Create the RNG to create seeds
    if seed == 0
        # Without a given seed
        rng_master = PCG.PCGStateOneseq()
    else
        # With a given seed
        rng_master = PCG.PCGStateOneseq(seed)
    end
    # From this RNG, create `num_rngs` seeds
    seed_list = rand(rng_master, UInt64, num_rngs)
    # With these seeds, seed the new RNG's
    rng_list = map(Xorshifts.Xoroshiro128Plus, seed_list)

    # Convert to static vector for faster performance
    seed_list = SVector{num_rngs}(seed_list)
    rng_list = SVector{num_rngs}(rng_list)

    return (seed_list, rng_list)
end

function rng_matrix!(rnd_matrix::AbstractArray, rng_list::AbstractArray)
    rnd_type = eltype(rnd_matrix)
    for i in axes(rnd_matrix, 2)
        rnd_matrix[:, i] = randn(rng_list[i], rnd_type, size(rnd_matrix, 1))
    end
end
