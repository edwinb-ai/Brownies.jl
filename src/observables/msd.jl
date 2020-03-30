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

function difusion!(msd::MeanSquaredDisplacement; sft::SelfScatteringFunction = nothing)
    total_sum = 0.0
    num_particle = size(msd.displacement, 2)
    norm_val = 0.0
    total_scattering = 0.0

    for i = 1:msd.naverage
        total_sum = 0.0
        @inbounds @fastmath for j = 1:(msd.naverage-i)
            ddiff = msd.displacement[j+i, :, :] - msd.displacement[j, :, :]
            ddiff = @. ddiff^2
            total_sum += sum(ddiff)

            # Is needed, compute the scattering values for an isotropic system
            if !isnothing(sft)
                total_scattering += _scattering(ddiff, sft.κ)
            end
        end

        # Normalize the MeanSquaredDisplacement
        norm_val = num_particle * (msd.naverage - i)
        total_sum /= norm_val
        msd.wt[i, 2] = total_sum

        # Normalize the SelfScatteringFunction
        if !isnothing(sft)
            total_scattering /= norm_val
            sft.dft[i, 2] = total_scattering
        end
    end
end

function _scattering(ddiff::AbstractArray, κ::AbstractFloat)
    distance = @. √(ddiff)
    scattering = @. sin(distance * κ) / (distance * κ)

    return scattering
end
