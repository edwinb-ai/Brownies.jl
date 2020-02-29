function record_positions!(msd::MeanSquaredDisplacement, positions::AbstractArray, ttime::AbstractFloat)
    msd.naverage += 1
    msd.timearray[msd.naverage] = ttime * msd.interval * (msd.naverage - 1.0)

    for i = axes(msd.displacement, 2)
        msd.displacement[msd.naverage, i, :] = positions[i, :]
    end
end

function difusion(msd::MeanSquaredDisplacement)
    total_sum = 0.0
    num_particle = size(msd.displacement, 2)
    norm_val = 0.0

    for i = 1:(msd.naverage-1)
        for j = 1:(msd.naverage - i)
            for k = axes(msd.displacement, 2)
                ddiff = msd.displacement[j+i, k, :] - msd.displacement[j, k, :]
                total_sum += sum(ddiff .^ 2)
            end
        end
        norm_val = num_particle * (msd.naverage - i)
        ddiff /= norm_val
        msd.wt[i] += ddiff
    end
end
