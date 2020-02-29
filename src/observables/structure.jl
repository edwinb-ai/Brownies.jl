function _fourier_structure!(gofr, rdfpos, sqobj, dr, dq, nm)
    sum_total = 0.0

    for i = 1:sqobj.nm
        sqobj.sq[i, 1] = (i - 1.0) * dq
    end

    for i = 2:sqobj.nm
        sum_total = sum(
            rdfpos .^ 2 .* gofr .* sin.(sqobj.sq[i, 1] * rdfpos) ./
                (sqobj.sq[i, 1] * rdfpos)
        )
        sum_total *= sqobj.norm_const * dr
        sqobj.sq[i, 2] = 1.0 + sum_total
    end
end

function compute_structure!(sqobj, gofr, s)
    rdfpos = @view gofr.gofr[2:end, 1]
    hfunc = gofr.gofr[2:end, 2] .- 1.0
    _fourier_structure!(hfunc, rdfpos, sqobj, gofr.dr, sqobj.dq, sqobj.nm)
end
