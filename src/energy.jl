function _add_forces!(
    f::Tuple,
    p::Tuple,
    idx::Integer,
    jdx::Integer,
    fp::AbstractFloat,
    rp::AbstractFloat,
)
    @inbounds for (k, c) in zip(p, f)
        c[idx] += (fp * k) / rp
        c[jdx] -= (fp * k) / rp
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

    return (xij, yij, zij, rij2)
end

function _compute_energy!(
    positions::AbstractArray,
    forces::AbstractArray,
    params::NamedTuple,
    pot::PairwisePotential;
    rdf = false,
    rdfobj = nothing,
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
        @inbounds for j = (i+1):params.N
            (xij, yij, zij, rij2) =
                _compute_distance(x, y, z, i, j, params.boxl)

            if rij2 < params.rc2
                rij2 = sqrt(rij2)
                u_pair, f_pair = apply!(pot, rij2)
                _add_forces!((fx, fy, fz), (xij, yij, zij), i, j, f_pair, rij2)
                energy += u_pair
                if rdf
                    nbin = round(rij2 / rdfobj.dr) + 1
                    nbin = Int(nbin)
                    if nbin <= rdfobj.nm
                        rdfobj.gofr[nbin, 2] += 2.0
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
    rdf = false,
    gofr = nothing,
)
    # Retrieve data from the system
    fill!(forces, zero(params.boxl))
    energy = _compute_energy!(
        positions,
        forces,
        params,
        pot;
        rdf = rdf,
        rdfobj = gofr,
    )
end
