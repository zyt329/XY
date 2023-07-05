


"""
    Wolff Alogrithm 
"""
function MC_Wolff(;
    T::Float64,
    measurement_runs::Int64,
    warmup_runs::Int64,
    init_micstate::Microstate
)

    # pass in initial microstate and parameters
    micstate = init_micstate

    # interaction strengths
    J_xy = (micstate.J1, micstate.J2)
    Jz = micstate.Jz

    # lattice information
    lattice_info = init_micstate.lattice_info

    L1 = lattice_info.L1
    L2 = lattice_info.L2
    xy_neibs = lattice_info.xy_neibs
    z_neibs = lattice_info.z_neibs
    all_neibs = lattice_info.all_neibs

    # all possible coordinates
    coords = [(z, i, j) for z in 1:2, i in 1:L1, j in 1:L2]

    # initialize sample holders
    samples = Sample()

    # perform sweeps
    for run_ind in 1:(measurement_runs+warmup_runs)

        # Choose a random reflection direction r
        θ_r = rand(Uniform(0.0, 2 * π))
        r = Spin(θ_r)

        # Choose a random site (x_init) to start building cluster
        x_init = rand(coords)
        sp_x_init = micstate.spins[x_init...]

        # flip x_init
        θ_new = mod1(π - sp_x_init.θ + 2θ_r, 2π)
        sp_x_new = Spin(θ_new)
        update(micstate, x_init, sp_x_new)

        # mark (hold) flipped(marked) site x_init in c
        c = [x_init]

        # loop over all sites in C
        # Note c is being constructed concurrently
        for x in c

            sp_x = micstate.spins[x...]

            # sp_x before flip
            θ_sp_x_original = mod1(π - sp_x.θ + 2θ_r, 2π)
            sp_x_original = Spin(θ_sp_x_original)

            # check all neighbors (y)
            for y in all_neibs[x]

                # spin value of y
                sp_y = micstate.spins[y...]

                if y[1] == x[1] # if x and y in the same layer
                    layer = x[1]
                    probability = 1 - exp(min(0, 2 * J_xy[layer] / T * (r * sp_x_original) * (r * sp_y)))
                else # if not in the same layer
                    probability = 1 - exp(min(0, 2 * Jz / T * (r * sp_x_original) * (r * sp_y)))
                end

                if rand() < probability
                    # flip spin Y
                    θ_new = mod1(π - sp_y.θ + 2θ_r, 2π)
                    sp_y_new = Spin(θ_new)
                    update(micstate, y, sp_y_new)

                    # mark spin Y
                    push!(c, y)
                end

            end
        end

        # Do not take sample while thermalizing
        # ======================================

        run_ind < warmup_runs + 1 && continue

        # =========== Take samples =============

        push!(samples.E, micstate.E)
        push!(samples.M, copy(micstate.M))

    end


    return samples, micstate
end
