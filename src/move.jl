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

        # We check for the MeanSquaredDisplacement first as it needs special treatment
        # with the periodic boundary conditions
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
            # We save positions in time here as needed by the MeanSquaredDisplacement
            if i % msd.interval == 0
                record_positions!(msd, positions, params.τ)
            end
        else
            # If it happens that we don't need the MeanSquaredDisplacement, then we proceed
            # moving the particles around
            ermak!(positions, forces, params.τ, params.boxl, params.random_matrix)

            # We check to see if we need to compute physical observables, if not
            # we just move particles around and check their energy
            if isnothing(rdf) && isnothing(zfactor)
                total_energy = energy_force!(positions, forces, params, pot)

                # We now check to see if the PairDistributionFunction is needed, if so
                # we step in and compute it
            elseif !isnothing(rdf)
                total_energy = energy_force!(positions, forces, params, pot; gofr = rdf)
                rdf.naverage += 1
                # If it so happens that the CompressibilityFactor needs to be computed
                # as well, we check to see if the object is passed and compute it
                if !isnothing(zfactor)
                    total_energy =
                        energy_force!(positions, forces, params, pot; zfactor = zfactor)
                    zfactor.naverage += 1
                end

                # If nothing but the CompressibilityFactor is needed, we compute it here
            elseif !isnothing(zfactor)
                total_energy =
                    energy_force!(positions, forces, params, pot; zfactor = zfactor)
                zfactor.naverage += 1
            end
        end

        # Finally, we check the interval and print the total energy
        # as well as save it
        if i % interval == 0
            @show total_energy
            idx = Int(i / interval)
            energies[idx] = total_energy
        end
    end

    # Lastly, we update the energy for the full system
    total_energy = energy_force!(positions, forces, params, pot)

    return total_energy
end

function _prepare(s::SimulationSystem, N::Integer, interval::Integer)
    # Allocate information for random values
    (_, rng_list) = _create_rngs(s.dims; seed = s.params.seed)
    random_matrix = zeros(typeof(s.boxl), size(s.system.positions))

    # Create an energy array to store values
    energy_size = Int(N / interval)

    # 2 dimensions because we only need the energy value and the movement number
    energies = zeros(energy_size, 2)

    # A collection of multiple parameters, it makes it easier to move it around
    # functions
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

function move!(N::Integer, s::SimulationSystem, pot::PairwisePotential; interval = 1000)
    # Retrieve system information
    @unpack positions, forces = s.system

    # Prepare all the arrays
    (energies, params) = _prepare(s, N, interval)

    # And then we move all the particles
    total_energy = _move_loop!(N, positions, forces, pot, energies, interval, params)
    # Lastly, we assign the system's energy
    s.energy = total_energy
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential,
    grobject::PairDistributionFunction;
    interval = 1000,
)
    # Retrieve system information
    @unpack positions, forces = s.system
    (energies, params) = _prepare(s, N, interval)

    # Pass in the PairDistributionFunction object and move particles
    total_energy =
        _move_loop!(N, positions, forces, pot, energies, interval, params; rdf = grobject)
    s.energy = total_energy

    # After all the statistics have been collected, we normalize the
    # PairDistributionFunction object
    grobject.naverage = N
    compute_rdf!(grobject, s)
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential,
    grobject::PairDistributionFunction,
    zfactor::ZFactor;
    interval = 1000,
)
    # Retrieve system information
    @unpack positions, forces = s.system
    (energies, params) = _prepare(s, N, interval)

    # Pass in the PairDistributionFunction and CompressibilityFactor objects
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

    # We compute the full CompressibilityFactor here
    zfactor.zval = compute_z(zfactor, s.params.N)

    # After all the statistics have been collected, we normalize the
    # PairDistributionFunction object
    grobject.naverage = N
    compute_rdf!(grobject, s)

    # Assign the full system's energy
    s.energy = total_energy
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential,
    zfactor::ZFactor;
    interval = 1000,
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
    # We compute the full CompressibilityFactor here
    zfactor.zval = compute_z(zfactor, s.params.N)

    # Assign the full system's energy
    s.energy = total_energy
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential,
    msd::MeanSquaredDisplacement;
    interval = 1000,
)
    # Retrieve system information
    @unpack positions, forces = s.system
    (energies, params) = _prepare(s, N, interval)
    _move_loop!(N, positions, forces, pot, energies, interval, params; msd = msd)

    # Compute and normalize the MeanSquaredDisplacement with all the recorded positions
    difusion!(msd)
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential,
    msd::MeanSquaredDisplacement,
    sft::SelfScatteringFunction;
    interval = 1000,
)
    # Retrieve system information
    @unpack positions, forces = s.system
    (energies, params) = _prepare(s, N, interval)
    _move_loop!(N, positions, forces, pot, energies, interval, params; msd = msd)

    # Copy the displacement vector from the MeanSquaredDisplacement to the
    # SelfScatteringFunction
    sft.dft[:, 1] = copy(msd.wt[:, 1])

    # Compute and normalize the MeanSquaredDisplacement and the SelfScatteringFunction
    # with all the recorded positions
    difusion!(msd; sft = sft)
end
