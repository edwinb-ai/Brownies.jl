function rdf(
    positions::AbstractArray,
    boxl::AbstractFloat,
    gofr::AbstractArray,
    rc::AbstractFloat,
    dr,
    nm,
)
    N = size(positions, 1)

    for i = 1:N-1
        @inbounds @fastmath for j = (i+1):N
            xij = positions[i, 1] - positions[j, 1]
            yij = positions[i, 2] - positions[j, 2]
            zij = positions[i, 3] - positions[j, 3]

            xij -= boxl * round(xij / boxl)
            yij -= boxl * round(yij / boxl)
            zij -= boxl * round(zij / boxl)

            rij2 = xij * xij + yij * yij + zij * zij

            if rij2 < rc
                rij2 = √rij2
                nbin = round(rij2 / dr) + 1
                nbin = Int(nbin)
                if nbin <= nm
                    gofr[nbin, 2] += 2.0
                end
            end
        end
    end
end

function normalize(
    gofr::AbstractArray,
    nm::Integer,
    dr::AbstractFloat,
    cnst
)
    for i = 1:nm
        rpos[i, 1] = dr * (i - 0.5)
        r_lo = (i - 1) * dr
        r_hi = r_lo + dr
        h_id = cnst * (r_hi^3 - r_lo^3)
        gofr[i, 2] /= h_id
    end
end
