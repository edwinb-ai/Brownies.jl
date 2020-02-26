abstract type PairwisePotential end

struct PseudoHS{T<:AbstractFloat} <: PairwisePotential
    λ::T
    σ::T
    a::T
    b::T
    temp::T
end

function _pseudoHS(λ::T, σ::T, temp::T) where {T<:AbstractFloat}
    b = λ / (λ - 1.0)
    b = convert(T, b)
    a = λ * b^(λ - 1.0)
    a = convert(T, a)
    return PseudoHS{T}(λ, σ, a, b, temp)
end

PseudoHS(λ::T, σ::T, temp::T) where {T<:AbstractFloat} = _pseudoHS(λ, σ, temp)
function PseudoHS(λ::T, temp::T) where {T<:AbstractFloat}
    return _pseudoHS(λ, convert(T, 1.0), temp)
end

@inline function _potential_func(
    x,
    λ::T,
    a::T,
    temp::T,
) where {T<:AbstractFloat}
    energy = (a / temp) * (x^-λ - x^-(λ - 1.0))
    return energy + (1.0 / temp)
end

@inline function _force_func(x, λ::T, a::T, temp::T) where {T<:AbstractFloat}
    force = (λ * x^-(λ + 1.0)) - ((λ - 1.0) * x^-λ)
    return force * (a / temp)
end  # function _force_func

function _apply!(x, λ, a, b, temp)
    if x < b
        energy = _potential_func(x, λ, a, temp)
        energy = oftype(b, energy)
        force = _force_func(x, λ, a, temp)
        force = oftype(b, force)
    else
        energy = zero(b)
        force = zero(b)
    end
    return energy, force
end

function apply!(p::PseudoHS, x::AbstractFloat)
    _apply!(x, p.λ, p.a, p.b, p.temp)
end
# TODO: Work on a cleaner interface for custom potentials
