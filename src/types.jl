"""
    ParticleSystem{T <: AbstractVector}

A handy type to manipulate the positions of a set of given particles.

# Fields
- `x::AbstractVector{T}`: stores the positions of the particles in the x
    directions, in rectangular coordinates
- `y::AbstractVector{T}`: stores the positions of the particles in the x
    directions, in rectangular coordinates
- `z::AbstractVector{T}`: stores the positions of the particles in the x
    directions, in rectangular coordinates
"""

mutable struct ParticleSystem{T<:AbstractArray}
    positions::T
    forces::T
end

function ParticleSystem(N::V, ::Type{T}) where {T<:AbstractFloat,V<:Integer}
    num_coord = 3
    pos = zeros(T, N, num_coord)
    forces = zeros(T, N, num_coord)
    return ParticleSystem{AbstractArray{T}}(pos, forces)
end

function ParticleSystem(
    N::V,
    M::V,
    ::Type{T},
) where {T<:AbstractFloat,V<:Integer}
    pos = zeros(T, N, M)
    forces = zeros(T, N, M)
    return ParticleSystem{AbstractArray{T}}(pos, forces)
end

"""
"""
mutable struct Parameters{U<:AbstractFloat,V<:Integer}
    ϕ::U
    kT::U
    N::V
    τ::U
    seed::V
end

# Just define packing fraction, temperature and time step
function Parameters(ϕ, kT, τ, ::Type{U}) where {U<:AbstractFloat}
    convert(U, ϕ)
    convert(U, kT)
    convert(U, τ)
    return Parameters{U,Integer}(ϕ, kT, 512, τ, 393216)
end

# Keep default seed, define total number of particles
function Parameters(
    ϕ,
    kT,
    N::Integer,
    τ,
    ::Type{U},
) where {U<:AbstractFloat,V<:Integer}
    return Parameters{U,Integer}(U(ϕ), U(kT), N, U(τ), 393216)
end

# Define a different seed
function Parameters(
    ϕ,
    kT,
    N::Integer,
    τ,
    s::Integer,
    ::Type{U},
) where {U<:AbstractFloat,V<:Integer}
    return Parameters{U,Integer}(U(ϕ), U(kT), N, U(τ), s)
end

"""
    SimulationSystem{U<:AbstractFloat} <: PhysicalSystem

A type that stores the most important parameters to perform a simulation.
These values are computed through the `Parameters` type.

# Fields
- `ρ::U`: the packing fraction of the system
"""
mutable struct SimulationSystem{U<:AbstractFloat}
    params::Parameters
    ρ::U
    boxl::U
    rc::U
    system::ParticleSystem
    energy::U
end

function SimulationSystem(
    params::Parameters,
    ::Type{T},
) where {T<:AbstractFloat}
    # Assign density for a 3D system
    ρ = params.ϕ * 6.0 / π
    boxl = ∛(params.N / ρ)
    rc = boxl * 0.5
    # Build a system of particles with the total number of particles
    psys = ParticleSystem(params.N, T)
    return SimulationSystem{T}(params, T(ρ), T(boxl), T(rc), psys, T(0.0))
end

"""
    Overload the scatter method from Plots to enable inspection
of the system.
"""
function scatter(s::SimulationSystem)
    x = view(s.system.positions, :, 1)
    y = view(s.system.positions, :, 2)
    z = view(s.system.positions, :, 3)
    scatter(x, y, z)
end
