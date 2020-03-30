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

function _save_positions(positions, forces, ϕ, N, seed; move = false)
    move_string = "-average.csv"
    filename = ".csv"

    # Save positions
    pos_df = DataFrame()
    pos_df.x = positions[:, 1]
    pos_df.y = positions[:, 2]
    pos_df.z = positions[:, 3]
    positions_name = "positions-$ϕ-$N-$seed"

    # Save forces
    forces_df = DataFrame()
    forces_df.x = forces[:, 1]
    forces_df.y = forces[:, 2]
    forces_df.z = forces[:, 3]
    forces_name = "forces-$ϕ-$N-$seed"

    if move
        positions_name *= move_string
        forces_name *= move_string
    else
        positions_name *= filename
        forces_name *= filename
    end

    # Save to files
    CSV.write(positions_name, pos_df)
    CSV.write(forces_name, forces_df)
end

function savetofile(s::SimulationSystem; move = false)
    @unpack positions, forces = s.system
    _save_positions(
        positions,
        forces,
        s.params.ϕ,
        s.params.N,
        s.params.seed;
        move = move,
    )
end

function savetofile(
    s::SimulationSystem,
    grobject::PairDistributionFunction;
    move = false,
)
    @unpack positions, forces = s.system
    _save_positions(
        positions,
        forces,
        s.params.ϕ,
        s.params.N,
        s.params.seed;
        move = move,
    )

    gofr_df = DataFrame()
    gofr_df.r = grobject.gofr[:, 1]
    gofr_df.gr = grobject.gofr[:, 2]
    gr_filename = "gr-$(s.params.ϕ)-$(s.params.N)-$(s.params.seed).csv"
    CSV.write(gr_filename, gofr_df)
end

function savetofile(
    s::SimulationSystem,
    msd::MeanSquaredDisplacement;
    move = false,
)
    @unpack positions, forces = s.system
    _save_positions(
        positions,
        forces,
        s.params.ϕ,
        s.params.N,
        s.params.seed;
        move = false,
    )

    msd_df = DataFrame()
    msd_df.t = msd.wt[:, 1]
    msd_df.xt = msd.wt[:, 2]
    msd_filename = "msd-$(s.params.ϕ)-$(s.params.N)-$(s.params.seed).csv"
    CSV.write(msd_filename, msd_df)
end
