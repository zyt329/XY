using XY
using BinningAnalysis
using Printf
using HDF5
using Statistics


function simulate()

    # model parameters
    J1 = -1.0
    J2 = -1.0
    Jz = 0.0

    L1 = 250
    L2 = 250

    Ts = range(0.4, 0.4, 1)

    # Monte Carlo parameters
    measurement_sweeps = 5*10^3
    warmup_sweeps = 1000

    # initialize spins
    single_flip_init_micstate = Microstate(J1, J2, Jz, L1, L2)
    Wolff_init_micstate = Microstate(J1, J2, Jz, L1, L2)

    # ===================================================#
    # =====   make folder and file to store data   ======#
    # ===================================================#

    # make folder to hold thermal avg data
    folder_prefix = "benchmark_single_layer"

    folder_full_path = make_indexed_folder(folder_prefix=folder_prefix, folder_path="/nfs/home/zyt329/metal_Hubbard_layer/XY_bilayer_runs/")

    single_flip_data = joinpath(folder_full_path, "single_flip.csv")


    # file header
    open(single_flip_data, "a+") do file
        write(file, "T    E    STD  Binned_STD    tau \n")
    end

    Wolff_data = joinpath(folder_full_path, "Wolff.csv")


    # file header
    open(Wolff_data, "a+") do file
        write(file, "T    E    E_STD    Binned_STD    tau \n")
    end

    for T in Ts

        print("simulating T = $T \r")
        # ================ MC ==================

        # single flip
        # use the configuration form last T for the new T
        single_flip_samples, single_flip_init_micstate = MC_single_flip(
            T=T,
            measurement_sweeps=measurement_sweeps,
            warmup_sweeps=warmup_sweeps,
            init_micstate=single_flip_init_micstate,
        )

        # Wolff Alogrithm
        # use the configuration form last T for the new T
        Wolff_samples, Wolff_init_micstate = MC_Wolff(
            T=T,
            measurement_runs=measurement_sweeps,
            warmup_runs=warmup_sweeps,
            init_micstate=Wolff_init_micstate,
        )

        # ============ Take Averages ============

        samples = (single_flip_samples, Wolff_samples)
        averages = (binning_anal(single_flip_samples), binning_anal(Wolff_samples))

       """ # Average quantities
        single_flip_samples_E = mean(single_flip_samples.E)
        Wolff_samples_E = mean(Wolff_samples.E)

        # Error estimation
        single_flip_std = std_error(single_flip_samples.E, method=:jackknife)
        Wolff_std = std_error(Wolff_samples.E, method=:jackknife)
        # Error Estimation by hand
        num_meas = length(single_flip_samples.E)
        num_bin = 10
        num_per_bin = Int(floor(num_meas / num_bin))
        single_flip_std_byhand = std([mean(single_flip_samples.E[((i-1)*num_per_bin+1):i*num_per_bin]) for i in 1:num_bin])
        Wolff_std_byhand = std([mean(Wolff_samples.E[((i-1)*num_per_bin+1):i*num_per_bin]) for i in 1:num_bin])


        single_flip_LB_E = LogBinner(single_flip_samples.E)
        Wolff_LB_E = LogBinner(Wolff_samples.E)

        single_flip_tau_E = tau(single_flip_LB_E)
        Wolff_tau_E = tau(Wolff_LB_E)"""

        # ============== save data ==============

        # prefixes
        prefixes = ("single_flip_", "Wolff_")

        for ind, prefix in enumerate(prefixes)
            # save averages
            save_name = prefix*"thermal_averages.h5"
            thermal_avg_path = joinpath(folder_full_path, save_name)
            thermal_avg_file = h5open(thermal_avg_path, "cw")

            # create data group for hdf5 file for the T
            group_name =@sprintf "T = %.4f" T
            g = create_group(thermal_avg_file, group_name)

            # =========  write averages to the group =========

            # write quantities
            for (key,val) in averages[ind]
                write(g, key, val)
            end
            close(thermal_avg_file)

            # =========  write spins configurations to file =========
            # save spins configurations
            save_name = prefix*"screen_shots.h5"
            screen_shot_path = joinpath(folder_full_path, save_name)
            screen_shot_file = h5open(screen_shot_path, "cw")

            write(screen_shot_file, group_name, samples[ind])
            close(screen_shot_file)
        end

        """# save averages and STD
        open(single_flip_data, "a+") do file
            write(file, "$T    $single_flip_samples_E    $single_flip_std   $single_flip_std_byhand    $single_flip_tau_E \n")
        end

        open(Wolff_data, "a+") do file
            write(file, "$T    $Wolff_samples_E    $Wolff_std   $Wolff_std_byhand    $Wolff_tau_E \n")
        end"""

        """# save screen shots
        save_name_T = @sprintf "spins-T=%.2f" T
        single_flip_screen_shot = joinpath(folder_full_path, "single_flip_screen_shot.h5")
        single_flip_screen_shot_file = h5open(single_flip_screen_shot, "cw")

        Wolff_screen_shot = joinpath(folder_full_path, "Wolff_screen_shot.h5")
        Wolff_screen_shot_file = h5open(Wolff_screen_shot, "cw")

        write(single_flip_screen_shot_file, save_name_T, spins2array(single_flip_init_micstate))
        write(Wolff_screen_shot_file, save_name_T, spins2array(Wolff_init_micstate))
        # close HDF5 files
        close(single_flip_screen_shot_file)
        close(Wolff_screen_shot_file)"""

    end

end


simulate()