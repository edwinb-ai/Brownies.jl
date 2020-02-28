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
)
    # * Main loop
    for i = 1:N
        # Always re-compute random numbers
        rng_matrix!(params.random_matrix, params.rng_list)
        # Move particles and update positions
        ermak!(positions, forces, params.τ, params.boxl, params.random_matrix)
        # Re-compute the energy and forces from the system

        if isnothing(rdf) & isnothing(zfactor)
            total_energy = energy_force!(positions, forces, params, pot)
        elseif !isnothing(rdf)
            total_energy = energy_force!(positions, forces, params, pot; gofr = rdf)
        elseif !isnothing(zfactor)
            total_energy = energy_force!(positions, forces, params, pot; zfactor = zfactor)
        end

        if i % interval == 0
            @show total_energy
            idx = Int(i / interval)
            energies[idx] = total_energy
        end
    end
    total_energy = energy_force!(positions, forces, params, pot)
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
    _move_loop!(N, positions, forces, pot, energies, interval, params)
    if tofiles
        savetofile(s, energies)
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
    _move_loop!(N, positions, forces, pot, energies, interval, params; rdf = grobject)
    grobject.naverage = N
    compute_rdf!(grobject, s)
    if tofiles
        savetofile(s, energies; move = true)
        @save "gr-$(s.ρ)-$(s.params.N).jld2" grobject.gofr
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
    _move_loop!(
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
    if tofiles
        savetofile(s, energies; move = true)
        @save "gr-$(s.ϕ)-$(s.params.N).jld2" grobject.gofr
    end
end
