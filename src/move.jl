function _move_loop!(
    N::Integer,
    positions::AbstractArray,
    forces::AbstractArray,
    pot::PairwisePotential,
    energies::AbstractArray,
    interval::Integer,
    params;
    rdf = nothing,
    zfactor = nothing,
    msd = nothing,
)
    # * Main loop
    for i = 1:N
        # Always re-compute random numbers
        rng_matrix!(params.random_matrix, params.rng_list)

        if !isnothing(msd)
            ermak!(
                positions,
                forces,
                params.τ,
                params.boxl,
                params.random_matrix;
                pbc = false,
            )
            total_energy = energy_force!(positions, forces, params, pot)
            if i % msd.interval == 0
                record_positions!(msd, positions, params.τ)
            end
        else

            ermak!(
                positions,
                forces,
                params.τ,
                params.boxl,
                params.random_matrix,
            )
            if isnothing(rdf) & isnothing(zfactor)
                total_energy = energy_force!(positions, forces, params, pot)
            elseif !isnothing(rdf)
                total_energy =
                    energy_force!(positions, forces, params, pot; gofr = rdf)
                rdf.naverage += 1
            elseif !isnothing(zfactor)
                total_energy = energy_force!(
                    positions,
                    forces,
                    params,
                    pot;
                    zfactor = zfactor,
                )
                zfactor.naverage += 1
            end
        end

        if i % interval == 0
            @show total_energy
            idx = Int(i / interval)
            energies[idx] = total_energy
        end
    end
    total_energy = energy_force!(positions, forces, params, pot)

    return total_energy
end

function _prepare(s::SimulationSystem, N::Integer, interval::Integer)
    # Allocate information for random values
    (_, rng_list) = _create_rngs(s.dims; seed = s.params.seed)
    random_matrix = zeros(typeof(s.boxl), size(s.system.positions))
    # Create an energy array to store values

    energy_size = Int(N / interval)
    # 2 because we only need the energy value and the movement number
    energies = zeros(energy_size, 2)

    # Some useful collections
    params = (
        N = s.params.N,
        rc = s.rc,
        boxl = s.boxl,
        rc2 = s.rc * s.rc,
        random_matrix = random_matrix,
        rng_list = rng_list,
        τ = s.params.τ,
    )

    return energies, params
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential;
    interval = 1000,
    tofiles = false,
)
    # Retrieve system information
    @unpack positions, forces = s.system
    (energies, params) = _prepare(s, N, interval)
    total_energy =
        _move_loop!(N, positions, forces, pot, energies, interval, params)
    s.energy = total_energy
    if tofiles
        savetofile(s)
    end
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential,
    grobject::PairDistributionFunction;
    interval = 1000,
    tofiles = false,
)
    # Retrieve system information
    @unpack positions, forces = s.system
    (energies, params) = _prepare(s, N, interval)
    total_energy = _move_loop!(
        N,
        positions,
        forces,
        pot,
        energies,
        interval,
        params;
        rdf = grobject,
    )
    s.energy = total_energy
    grobject.naverage = N
    compute_rdf!(grobject, s)
    if tofiles
        savetofile(s, grobject; move = true)
    end
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential,
    grobject::PairDistributionFunction,
    zfactor::ZFactor;
    interval = 1000,
    tofiles = false,
)
    # Retrieve system information
    @unpack positions, forces = s.system
    (energies, params) = _prepare(s, N, interval)
    total_energy = _move_loop!(
        N,
        positions,
        forces,
        pot,
        energies,
        interval,
        params;
        rdf = grobject,
        zfactor = zfactor,
    )
    grobject.naverage = N
    compute_rdf!(grobject, s)
    s.energy = total_energy
    if tofiles
        savetofile(s, grobject; move = true)
    end
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential,
    zfactor::ZFactor;
    interval = 1000,
    tofiles = false,
)
    # Retrieve system information
    @unpack positions, forces = s.system
    (energies, params) = _prepare(s, N, interval)
    total_energy = _move_loop!(
        N,
        positions,
        forces,
        pot,
        energies,
        interval,
        params;
        zfactor = zfactor,
    )
    zfactor.zval = 1.0 - (zfactor.zval / (3.0 * zfactor.naverage * s.params.N))
    s.energy = total_energy
    if tofiles
        savetofile(s; move = true)
    end
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential,
    msd::MeanSquaredDisplacement;
    interval = 1000,
    tofiles = false,
)
    # Retrieve system information
    @unpack positions, forces = s.system
    (energies, params) = _prepare(s, N, interval)
    _move_loop!(
        N,
        positions,
        forces,
        pot,
        energies,
        interval,
        params;
        msd = msd,
    )
    difusion!(msd)
    if tofiles
        savetofile(s, msd; move = true)
    end
end
