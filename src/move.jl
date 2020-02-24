function _move_loop!(N::Integer, s::SimulationSystem, pot::PairwisePotential)
    # Retrieve system information
    @unpack positions, forces = s.system

    # Some useful collections
    params = (N = s.params.N, rc = s.rc, boxl = s.boxl, rc2 = s.rc * s.rc)

    # Check first energy
    total_energy = zero(s.energy)
    total_energy = energy_force!(positions, forces, params, pot)
    println("First energy ", total_energy, typeof(total_energy))

    # Allocate information for random values
    (_, rng_list) = _create_rngs(3; seed = s.params.seed)
    random_matrix = zeros(typeof(s.boxl), size(s.system.positions))

    # * Main loop
    for i = 1:N
        # Always re-compute random numbers
        rng_matrix!(random_matrix, rng_list)
        ermakv!(positions, forces, s.params.τ, s.boxl, random_matrix)
        total_energy = energy_force!(positions, forces, params, pot)
        if i % 1000 == 0
            @show total_energy
        end
    end
    s.energy = total_energy
end

function move!(N::Integer, s::SimulationSystem, pot::PairwisePotential)
    _move_loop(N, s, pot)
end
