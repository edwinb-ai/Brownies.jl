# Some useful constants
const sixoverpi = 6.0 / π
const fourpi = 4.0 * π

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

function ParticleSystem(N::V, M::V, ::Type{T}) where {T<:AbstractFloat,V<:Integer}
    pos = zeros(T, N, M)
    forces = zeros(T, N, M)
    return ParticleSystem{AbstractArray{T}}(pos, forces)
end

"""
"""
mutable struct Parameters{U<:AbstractFloat}
    ϕ::U
    kT::U
    N::Integer
    τ::U
    seed::Integer
end

# Just define packing fraction, temperature and time step
function Parameters(
    ϕ::AbstractFloat,
    kT::AbstractFloat,
    τ::AbstractFloat,
    ::Type{U},
) where {U}
    return Parameters{U}(U(ϕ), U(kT), 512, U(τ), 0)
end

# Keep default seed, define total number of particles
function Parameters(
    ϕ::AbstractFloat,
    kT::AbstractFloat,
    N::Integer,
    τ::AbstractFloat,
    ::Type{U},
) where {U}
    return Parameters{U}(U(ϕ), U(kT), N, U(τ), 0)
end

# Define a different seed
function Parameters(
    ϕ::AbstractFloat,
    kT::AbstractFloat,
    N::Integer,
    τ::AbstractFloat,
    s::Integer,
    ::Type{U},
) where {U}
    return Parameters{U}(U(ϕ), U(kT), N, U(τ), s)
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
    dims::Integer
end

function SimulationSystem(params::Parameters, dims::Integer, ::Type{T}) where {T}
    # Assign density for a 3D system
    ρ = params.ϕ * sixoverpi
    boxl = ∛(params.N / ρ)
    rc = boxl * 0.5
    # Build a system of particles with the total number of particles
    psys = ParticleSystem(params.N, dims, T)
    return SimulationSystem{T}(params, T(ρ), T(boxl), T(rc), psys, T(0.0), dims)
end

abstract type Structure end
mutable struct PairDistributionFunction <: Structure
    gofr::AbstractArray
    nm::Integer
    dr::AbstractFloat
    naverage::AbstractFloat
    norm_const::AbstractFloat
end

function PairDistributionFunction(s::SimulationSystem, nm::Integer)
    normalizing_constant = fourpi * s.ρ
    gofr = zeros(typeof(s.ρ), (nm, 2))
    dr = s.rc / nm
    return PairDistributionFunction(gofr, nm, dr, 0.0, normalizing_constant)
end

mutable struct StructureFactor{T<:AbstractFloat} <: Structure
    sq::AbstractArray
    dq::T
    nm::Integer
    naverage::T
    norm_const::T
end

function StructureFactor(s::SimulationSystem, rdf::PairDistributionFunction)
    compute_type = typeof(s.ρ)
    normalizing_constant = fourpi * s.ρ
    sq = zeros(compute_type, (rdf.nm, 2))
    dq = 0.1 * π / s.rc
    return StructureFactor{compute_type}(sq, dq, rdf.nm, rdf.naverage, normalizing_constant)
end

abstract type Thermodynamics end
mutable struct ZFactor{T<:AbstractFloat} <: Thermodynamics
    zval::T
    naverage::T
end

function compute_z(zfac::ZFactor, N::Real)
    return 1.0 - (zfac.zval / (3.0 * zfac.naverage * N))
end

abstract type Dynamics end
mutable struct MeanSquaredDisplacement{V<:Integer} <: Dynamics
    displacement::AbstractArray
    mt::V
    wt::AbstractArray
    naverage::V
    interval::V
end
function MeanSquaredDisplacement(s::SimulationSystem, mt::Integer, interval::Integer)
    displacement = zeros(mt, s.params.N, s.dims)
    wt = zeros(mt, 2)
    return MeanSquaredDisplacement(displacement, mt, wt, 0, interval)
end
