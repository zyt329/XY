"""
Deciding whether to accept a proposed change or not.

Parameter:
micstate::Microstate : The current state.
numbering::Int64 : A tuple of spin index to flip.

Return:
True of false. True means to flip false means not to.
"""
function accept(micstate::Microstate, coord::Tuple{Int64,Int64,Int64}, Sp_new, T::Float64)
    ΔE = Energy_Diff(micstate, coord, Sp_new)

    p = exp(-ΔE / T) / (1 + exp(-ΔE / T))
    r = rand()
    if r < p
        return true
    else
        return false
    end
end



"""
Function to generate and store samples for a single temperature.

Parameter:
T::Float64 : Temperature
steps::Int64 : steps of the Monte Carlo Simulation
cutoff::Int64 : After the cutoff do we start to take down the samples.
σ::Float64 : Variation of the proposed new value for the Potts value each step in the simulation.

Output:
Samples::Array{Microstate,1} : An Array of samples, each entry is of type Microstate, representing one Microstate.
"""
function MC_single_flip(;
    T::Float64,
    measurement_sweeps::Int64,
    warmup_sweeps::Int64,
    init_micstate::Microstate
)

    # pass in initial microstate and parameters
    micstate = init_micstate

    L1 = init_micstate.lattice_info.L1
    L2 = init_micstate.lattice_info.L2

    # initialize sample holders
    samples = Sample()

    # perform sweeps
    for sweep_ind in 1:(measurement_sweeps+warmup_sweeps)

        # sweep over all sites
        coords = [(z, i, j) for z in 1:2, i in 1:L1, j in 1:L2]
        for coord in coords

            # propose a new spin
            θ_new = rand(Uniform(0.0, 2 * π))
            sp_new = Spin(θ_new)

            #val = mod1(conf.conf[m,n] + Int(ceil(rand(Normal(0, σ)))), conf.Q)
            if accept(micstate, coord, sp_new, T)
                update(micstate, coord, sp_new)
            end
        end

        # Do not take sample while thermalizing
        # =====================================================
        sweep_ind < warmup_sweeps + 1 && continue
        # =====================================================

        # =====================================================
        # ===========       Taking Samples!      ==============
        # =====================================================

        push!(samples.E, micstate.E)
        push!(samples.M, copy(micstate.M))

    end

    return samples, micstate
end



