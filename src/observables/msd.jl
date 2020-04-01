function record_positions!(
    msd::MeanSquaredDisplacement,
    positions::AbstractArray,
    ttime::AbstractFloat,
)
    # If we record positions, we update the counter here
    msd.naverage += 1
    # We keep record of the time the positions are collected
    msd.wt[msd.naverage, 1] = ttime * msd.interval * (msd.naverage - 1.0)

    # For each average value, we record the positions of all the particles
    # in that specific point in time
    for i in axes(msd.displacement, 2)
        msd.displacement[msd.naverage, i, :] = @view positions[i, :]
    end
end

function difusion!(msd::MeanSquaredDisplacement, sft::SelfScatteringFunction)
    num_particle = size(msd.displacement, 2)

    for i = 1:msd.naverage
        total_msd = 0.0
        total_scattering = 0.0

        for j = 1:msd.naverage-i
            # Obtain the distance between particles along all spatial dimensions
            ddiff = msd.displacement[(j+i), :, :] - msd.displacement[j, :, :]
            # We sum only the dimensions, not the particles
            ddiff = sum(ddiff .^ 2; dims = 2)

            # For the MeanSquareDisplacement, we actually need the sum of all the particles
            total_msd += sum(ddiff)

            # For the SelfScatteringFunction, we need the actual Euclidean distance
            # per particle, so we need to loop and extract the distance for each particle
            @inbounds @fastmath for d in ddiff
                # We already have the sum of square, just need the sqrt to complete
                # the Euclidean distance computation
                distance = √d
                # Given that the system is isotropic, we compute it as such
                total_scattering += sin(distance * sft.κ) / (distance * sft.κ)
            end
        end

        # Normalize the MeanSquaredDisplacement
        norm_val = num_particle * (msd.naverage - i)
        msd.wt[i, 2] = total_msd / norm_val

        # Normalize the SelfScatteringFunction
        sft.dft[i, 2] = total_scattering / norm_val
    end
end

function difusion!(msd::MeanSquaredDisplacement)
    num_particle = size(msd.displacement, 2)

    for i = 1:msd.naverage
        total_msd = 0.0

        @inbounds @fastmath for j = 1:msd.naverage-i
            # Obtain the distance between particles along all spatial dimensions
            ddiff = msd.displacement[(j+i), :, :] - msd.displacement[j, :, :]
            # And we can sum both the dimensions and the particles
            total_msd += sum(ddiff .^ 2)
        end

        # Normalize the MeanSquaredDisplacement
        norm_val = num_particle * (msd.naverage - i)
        msd.wt[i, 2] = total_msd / norm_val
    end
end
