function _add_forces!(
    f::Tuple,
    p,
    idx::Integer,
    jdx::Integer,
    fp::AbstractFloat,
    rp::AbstractFloat,
)
    @inbounds for (k, c) in zip(p, f)
        @fastmath c[idx] += (fp * k) / rp
        @fastmath c[jdx] -= (fp * k) / rp
    end
end

function _compute_distance(
    x::SubArray,
    y::SubArray,
    z::SubArray,
    i::Integer,
    j::Integer,
    boxl::Real,
)
    xij = x[i] - x[j]
    yij = y[i] - y[j]
    zij = z[i] - z[j]

    xij -= boxl * round(xij / boxl)
    yij -= boxl * round(yij / boxl)
    zij -= boxl * round(zij / boxl)

    rij2 = xij * xij + yij * yij + zij * zij

    # st_positions = @SVector [xij, yij, zij]
    st_positions = (xij, yij, zij)

    return st_positions, rij2
end

function _compressibilityz(pos, fp::AbstractFloat, rij::AbstractFloat, total_sum)
    for p in pos
        total_sum += (p^2 * -fp) / rij
    end
end

function _compute_energy!(
    positions::AbstractArray,
    forces::AbstractArray,
    params::NamedTuple,
    pot::PairwisePotential;
    rdfobj = nothing,
    zfactor = nothing,
)
    # Initialize necessary variables
    energy = zero(params.rc2)

    x = view(positions, :, 1)
    y = view(positions, :, 2)
    z = view(positions, :, 3)
    fx = view(forces, :, 1)
    fy = view(forces, :, 2)
    fz = view(forces, :, 3)

    for i = 1:params.N-1
        @inbounds @fastmath for j = (i+1):params.N
            (static_positions, rij2) = _compute_distance(x, y, z, i, j, params.boxl)

            if rij2 < params.rc2
                rij2 = √rij2
                u_pair, f_pair = apply!(pot, rij2)
                _add_forces!((fx, fy, fz), static_positions, i, j, f_pair, rij2)
                energy += u_pair
                if !isnothing(rdfobj)
                    simple_rdf!(rdfobj, rij2)
                end
                if !isnothing(zfactor)
                    # tmp_sum = @. (static_positions^2.0 * -f_pair) / rij2
                    # zfactor.zval += sum(tmp_sum)
                    for p in static_positions
                        zfactor.zval += (p * -f_pair) / rij2
                    end
                end
            end
        end
    end
    return energy / params.N
end

function energy_force!(
    positions::AbstractArray,
    forces::AbstractArray,
    params::NamedTuple,
    pot::PairwisePotential;
    gofr = nothing,
    zfactor = nothing,
)
    # Retrieve data from the system
    fill!(forces, zero(params.boxl))
    energy =
        _compute_energy!(positions, forces, params, pot; rdfobj = gofr, zfactor = zfactor)
end
