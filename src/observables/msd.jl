function record_positions!(
    msd::MeanSquaredDisplacement,
    positions::AbstractArray,
    ttime::AbstractFloat,
)
    msd.naverage += 1
    msd.wt[msd.naverage, 1] = ttime * msd.interval * (msd.naverage - 1.0)

    for i in axes(msd.displacement, 2)
        msd.displacement[msd.naverage, i, :] = view(positions, i, :)
    end
end

function difusion!(msd::MeanSquaredDisplacement)
    total_sum = 0.0
    num_particle = size(msd.displacement, 2)
    norm_val = 0.0

    for i = 1:(msd.naverage-1)
        total_sum = 0.0
        @inbounds @fastmath for j = 1:(msd.naverage-i)
            ddiff = msd.displacement[j+i, :, :] - msd.displacement[j, :, :]
            total_sum += sum(ddiff .^ 2)
        end
        norm_val = num_particle * (msd.naverage - i)
        total_sum /= norm_val
        msd.wt[i, 2] = total_sum
    end
end
