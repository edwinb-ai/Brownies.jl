function _ermak_loop!(
    pos::AbstractArray,
    forces::AbstractArray,
    rnd_matrix::AbstractArray,
    ttime::AbstractFloat,
    boxl::AbstractFloat;
    pbc::Bool = true,
)
    N = size(pos, 1)
    for j  = axes(pos, 2)
        @inbounds @fastmath for i = 1:N
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
end
