function _fourier_structure!(
    gofr::AbstractArray,
    rdfpos::AbstractArray,
    sqobj::StructureFactor,
    dr::AbstractFloat,
)
    sum_total = 0.0

    for i = 1:sqobj.nm
        sqobj.sq[i, 1] = (i - 1.0) * sqobj.dq
    end

    for i = 2:sqobj.nm
        integrand =
            @. rdfpos^2 * gofr * sin(sqobj.sq[i, 1] * rdfpos) / (sqobj.sq[i, 1] * rdfpos)
        sum_total = sum(integrand)
        sum_total *= sqobj.norm_const * dr
        sqobj.sq[i, 2] = 1.0 + sum_total
    end
end

function compute_structure!(
    sqobj::StructureFactor,
    gofr::PairDistributionFunction,
    s::SimulationSystem,
)
    rdfpos = @view gofr.gofr[2:end, 1]
    hfunc = gofr.gofr[2:end, 2] .- 1.0
    _fourier_structure!(hfunc, rdfpos, sqobj, gofr.dr)
end
