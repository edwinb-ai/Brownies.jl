function _move_loop!(
    N::Integer,
    positions::AbstractArray,
    forces::AbstractArray,
    pot::PairwisePotential,
    energies::AbstractArray,
    interval::Integer,
    params;
    average = false,
)
    naverage = 0
    # * Main loop
    for i = 1:N
        # Always re-compute random numbers
        rng_matrix!(params.random_matrix, params.rng_list)
        # Move particles and update positions
        ermak!(positions, forces, params.τ, params.boxl, params.random_matrix)
        # Re-compute the energy and forces from the system
        total_energy = energy_force!(positions, forces, params, pot)
        if i % interval == 0
            @show total_energy
            idx = Int(i / interval)
            energies[idx] = total_energy
        end
        if average
            if i % 100 == 0
                naverage += 1
            end
        end
    end
    total_energy = energy_force!(positions, forces, params, pot)
    if average
        return total_energy, naverage
    else
        return total_energy
    end
end

function _move_rdf!(
    N::Integer,
    positions::AbstractArray,
    forces::AbstractArray,
    pot::PairwisePotential,
    energies::AbstractArray,
    interval::Integer,
    params,
    gofr;
    average = false,
)
    naverage = 0
    # * Main loop
    for i = 1:N
        # Always re-compute random numbers
        rng_matrix!(params.random_matrix, params.rng_list)
        # Move particles and update positions
        ermak!(positions, forces, params.τ, params.boxl, params.random_matrix)
        # Re-compute the energy and forces from the system
        total_energy = energy_force!(
            positions,
            forces,
            params,
            pot;
            rdf = true,
            gofr = gofr,
        )
        if i % interval == 0
            @show total_energy
            idx = Int(i / interval)
            energies[idx] = total_energy
        end
        if average
            if i % 100 == 0
                naverage += 1
            end
        end
    end
    total_energy = energy_force!(positions, forces, params, pot)
    if average
        return total_energy, naverage
    else
        return total_energy
    end
end

function thermalize!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential;
    interval = 1000,
    tofiles = false,
)
    # Retrieve system information
    @unpack positions, forces = s.system

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

    # Check first energy
    total_energy = zero(s.energy)
    total_energy = energy_force!(positions, forces, params, pot)
    println("First energy ", total_energy)

    # Move the particles !
    total_energy =
        _move_loop!(N, positions, forces, pot, energies, interval, params)

    # Save all relevant values
    s.energy = total_energy
    @pack! s.system = positions, forces

    # Store in results files
    if tofiles
        savetofile(s, energies)
    end
end

function move!(
    N::Integer,
    s::SimulationSystem,
    pot::PairwisePotential;
    interval = 1000,
    gdr = false,
    nm = 2^9,
    tofiles = false,
)
    # Retrieve system information
    @unpack positions, forces = s.system

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
    #
    if gdr
        gfunc = PairDistributionFunction(s, nm)
        total_energy, naverage = _move_rdf!(
            N,
            positions,
            forces,
            pot,
            energies,
            interval,
            params,
            gfunc;
            average = true,
        )
        gfunc.naverage = naverage
        compute_rdf!(gfunc, s)
        if tofiles
            @save "gr-$(s.params.ϕ)-$(s.params.N).jld2" gfunc
        end
        return gfunc
    else
        total_energy =
            _move_loop!(N, positions, forces, pot, energies, interval, params)
    end
    # Store in results files
    if tofiles
        savetofile(s, energies; move = true)
    end
end
