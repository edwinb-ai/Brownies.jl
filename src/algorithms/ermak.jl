function _ermak_loop!(
    pos::AbstractArray,
    forces::AbstractArray,
    rnd_matrix::AbstractArray,
    ttime::AbstractFloat,
    boxl::AbstractFloat;
    pbc::Bool = true,
)
    for j  = axes(pos, 2)
        @inbounds for i = axes(pos, 1)
            pos[i, j] += (forces[i, j] * ttime) + rnd_matrix[i, j]
            if pbc
                pos[i, j] -= boxl * round(pos[i, j] / boxl)
            end
        end
    end
end

function ermak!(
    positions::AbstractArray,
    forces::AbstractArray,
    ttime::AbstractFloat,
    boxl::AbstractFloat,
    rnd_matrix::AbstractArray;
    pbc::Bool = true,
)
    𝝈 = √(2.0 * ttime)
    𝝈 = oftype(ttime, 𝝈)
    rnd_matrix .*= 𝝈
    _ermak_loop!(positions, forces, rnd_matrix, ttime, boxl; pbc = pbc)
    # @. positions += (forces * ttime) + rnd_matrix
    # if pbc
    #     @. positions -= boxl * round(positions - boxl)
    # end
end
