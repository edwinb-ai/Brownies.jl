function _lattice_configuration!(x, y, z, dist, dist_half, rc, N)
    # Define the first positions
    x[1] = -rc + dist_half
    y[1] = -rc + dist_half
    z[1] = -rc + dist_half
    # Create a complete lattice
    @inbounds for i = 2:N
        x[i] = x[i - 1] + dist
        y[i] = y[i - 1]
        z[i] = z[i - 1]

        if x[i] > rc
            x[i] = x[1]
            y[i] = y[i - 1] + dist

            if y[i] > rc
                x[i] = x[1]
                y[i] = y[1]
                z[i] = z[i - 1] + dist
            end
        end
    end
end

@doc raw"""
    initialize!(sys::SimulationSystem)

Creates a meshgrid-like initial configuration for the given `sys` simulation
system. It uses the internal parameters to make use of the density, cut-off
radius and total number of particles. It modifies `sys` positions by using
views of the array.

# Arguments
- `sys::SimulationSystem`: The system to initialize. It modifies the system's
positions.
"""
function initialize!(sys::SimulationSystem)
    # Get the array views from the system positions
    x = view(sys.system.positions, :, 1)
    y = view(sys.system.positions, :, 2)
    z = view(sys.system.positions, :, 3)
    # Define the inter-particle distance
    distance = oftype(sys.ρ, ∛(1.0 / sys.ρ))
    # And then half of that
    distance_half = oftype(sys.ρ, distance * 0.5)
    _lattice_configuration!(x, y, z, distance, distance_half, sys.rc, sys.params.N)
end
