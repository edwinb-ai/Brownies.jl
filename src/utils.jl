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
    N, M = size(rnd_matrix)
    rnd_type = eltype(rnd_matrix)
    for i = 1:M
        @inbounds for j = 1:N
            rnd_matrix[j, i] = randn(rng_list[i], rnd_type)
        end
    end
end