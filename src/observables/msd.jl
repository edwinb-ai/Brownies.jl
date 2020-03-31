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
    num_particle = size(msd.displacement, 2)

    for i = 1:msd.naverage
        total_msd = 0.0

        @inbounds for j = 1:(msd.naverage-i)
            ddiff = msd.displacement[j+i, :, :] - msd.displacement[j, :, :]
            ddiff = sum(ddiff.^2)
            total_msd += ddiff
        end

        # Normalize the MeanSquaredDisplacement
        norm_val = num_particle * (msd.naverage - i)
        total_msd /= norm_val
        msd.wt[i, 2] = total_msd
    end
end

function difusion!(msd::MeanSquaredDisplacement, sft::SelfScatteringFunction)
    num_particle = size(msd.displacement, 2)

    for i = 1:msd.naverage
        total_msd = 0.0
        total_scattering = 0.0

        @inbounds for j = 1:(msd.naverage-i)
            ddiff = msd.displacement[j+i, :, :] - msd.displacement[j, :, :]
            ddiff = norm(ddiff)
            total_msd += ddiff^2

            # If needed, compute the scattering values for an isotropic system
            total_scattering += sin(ddiff * sft.κ) / (ddiff * sft.κ)
        end

        # Normalize the MeanSquaredDisplacement
        norm_val = num_particle * (msd.naverage - i)
        total_msd /= norm_val
        msd.wt[i, 2] = total_msd

        # Normalize the SelfScatteringFunction
        total_scattering /= norm_val
        sft.dft[i, 2] = total_scattering
    end
end
