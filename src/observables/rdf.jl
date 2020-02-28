function simple_rdf!(rdfobj::PairDistributionFunction, rij::AbstractFloat)
    nbin = round(rij / rdfobj.dr) + 1
    nbin = Int(nbin)
    if nbin <= rdfobj.nm
        rdfobj.gofr[nbin, 2] += 2.0
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
    for i = 2:nm
        gofr[i, 1] = dr * (i - 1.0)
        dv = cnst * gofr[i, 1]^2 * dr
        gofr[i, 2] /= dv * nstep * N
    end
end

function compute_rdf!(grobject::PairDistributionFunction, s::SimulationSystem)
    _normalize!(
        grobject.gofr,
        grobject.nm,
        grobject.dr,
        grobject.norm_const,
        s.params.N,
        grobject.naverage,
    )
end
