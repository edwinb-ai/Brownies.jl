function _move_loop!(
    N::Integer,
    positions::AbstractArray,
    forces::AbstractArray,
    pot::PairwisePotential,
    energies::AbstractArray,
    interval::Integer,
)
    # * Main loop
    for i = 1:N
        # Always re-compute random numbers
        rng_matrix!(random_matrix, rng_list)
        # Move particles and update positions
        ermakv!(positions, forces, s.params.τ, s.boxl, random_matrix)
        # Re-compute the energy and forces from the system
        total_energy = energy_force!(positions, forces, params, pot)
        if i % interval == 0
            @show total_energy
            energies[i] = total_energy
        end
    end
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

    # Some useful collections
    params = (N = s.params.N, rc = s.rc, boxl = s.boxl, rc2 = s.rc * s.rc)

    # Check first energy
    total_energy = zero(s.energy)
    total_energy = energy_force!(positions, forces, params, pot)
    println("First energy ", total_energy)

    # Allocate information for random values
    (_, rng_list) = _create_rngs(s.dims; seed = s.params.seed)
    random_matrix = zeros(typeof(s.boxl), size(s.system.positions))
    # Create an energy array to store values
    energy_size = Int(N / interval)
    # 2 because we only need the energy value and the movement number
    energies = zeros(energy_size, 2)

    # Move the particles !
    _move_loop(N, positions, forces, pot, energies, interval)

    # Save all relevant values
    s.energy = total_energy
    @pack! s.system = positions, forces

    # Store in results files
    if tofiles
        savetofile(s, energies)
    end
end
