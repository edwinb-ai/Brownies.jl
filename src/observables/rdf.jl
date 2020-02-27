const pi4 = 4.0 * π

function _rdf!(
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

            rij2 = sqrt(xij * xij + yij * yij + zij * zij)

            if rij2 < rc
                nbin = round(rij2 / dr) + 1
                nbin = Int(nbin)
                if nbin <= nm
                    gofr[nbin, 2] += 2.0
                end
            end
        end
    end
end

function _normalize!(
    gofr::AbstractArray,
    nm::Integer,
    dr::AbstractFloat,
    cnst,
    N,
    nstep,
)
    for i = 1:nm
        gofr[i, 1] = dr * (i - 1.0)
        dv = pi4 .* gofr[i, 1].^2 .* dr
        gofr[i, 2] ./= dv .* nstep .* N
    end
end

function compute_rdf!(grobject::PairDistributionFunction, s::SimulationSystem)
    # _rdf!(
    #     s.system.positions,
    #     s.boxl,
    #     grobject.gofr,
    #     s.rc,
    #     grobject.dr,
    #     grobject.nm,
    # )
    _normalize!(
        grobject.gofr,
        grobject.nm,
        grobject.dr,
        grobject.norm_const,
        s.params.N,
        grobject.naverage,
    )
end
