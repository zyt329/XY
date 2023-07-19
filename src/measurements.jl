"""
Container to store value of measurement for each microstate in the Markov Chain.
"""
struct Sample
    E::Vector{Float64}
    E²::Vector{Float64}
    M::Vector{Vector{Float64}}
    M²::Vector{Float64}
    γ₁::Vector{Float64}
    γ₂::Vector{Float64}
    γ₁₂::Vector{Float64}
    spins::Vector{Array{Float64,3}}
    function Sample()
        new(Float64[], Float64[], Vector{Float64}[], Float64[], Float64[], Float64[], Float64[], Array{Float64,3}[])
    end
end

"""
    calculate helicity modulus γ of 1,2 and 1+2 layers of a single microstate.
        γ=β/L²∑ₓcos(θₓ-θₓ₊₁) - β²/L²(∑ₓsin(θₓ-θₓ₊₁))² = β/L²term1 - β²/L²term2²
"""
function γs(micstate, T)
    # load information
    β = 1 / T
    lattice_info = micstate.lattice_info

    L1 = lattice_info.L1
    L2 = lattice_info.L2
    xy_neibs = lattice_info.xy_neibs

    # all possible coordinates in a layer
    coords = [(z, i, j) for z in 1:2, i in 1:L1, j in 1:L2]

    # initialize γ
    γ₁ = 0
    γ₂ = 0
    γ₁₂ = 0

    # layer 1
    z = 1
    term1 = 0
    term2 = 0
    for j in 1:L2
        for i in 1:L1
            term1 += cos(micstate.spins[1, i, j] * micstate.spins[1, mod1(i + 1, L1), j])
            term2 += cos(micstate.spins[1, i, j] ⊗ micstate.spins[1, mod1(i + 1, L1), j])
        end
    end
    γ₁ = β / (L1 * L2) * term1 - β^2 / (L1 * L2) * term2^2

    # layer 2
    z = 2
    term1 = 0
    term2 = 0
    for j in 1:L2
        for i in 1:L1
            term1 += cos(micstate.spins[1, i, j] * micstate.spins[1, mod1(i + 1, L1), j])
            term2 += cos(micstate.spins[1, i, j] ⊗ micstate.spins[1, mod1(i + 1, L1), j])
        end
    end
    γ₂ = β / (L1 * L2) * term1 - β^2 / (L1 * L2) * term2^2

    return γ₁, γ₂, (γ₁ + γ₂)
end

function Eavg_site(E_sample::Array; L1, L2)
    samplelength = length(E_sample)
    Eavg = 0
    for i = 1:samplelength
        Eavg += E_sample[i]
    end
    Eavg = Eavg / samplelength / (L1 * L2)
    return Eavg
end

function C_site(E_sample::Array; T, L1, L2)
    samplelength = length(E_sample)
    Eavg = Eavg_site(E_sample, L1=L1, L2=L2)
    E2avg = 0
    for i = 1:samplelength
        E2avg += (E_sample[i] / (L1 * L2))^2
    end
    E2avg = E2avg / samplelength
    C = 1 / T * (E2avg - Eavg^2) * (L1 * L2)
    return C
end


"""
    Average of S^2 of the whole lattice, without normalizing it with respect to system size N. <S^2_α> =  ∑_i,j <S_α,i S_α,j> (α=x,y,z being three directions)

    return [<S^2_x>, <S^2_y>, <S^2_z>]
"""
function S2xyz(S_sample::Array; L1, L2)
    samplelength = length(S_sample)
    S2xyz_avg = [0, 0, 0]
    for i = 1:samplelength
        S2xyz_avg += S_sample[i] .* S_sample[i]
    end
    S2xyz_avg = S2xyz_avg / samplelength
    return S2xyz_avg
end

function χ_xy(S_sample::Array; L1, L2, Temp)
    samplelength = length(S_sample)
    S2xyz_avg = S2xyz(S_sample; L1=L1, L2=L2)
    S_parallel = 0 #parallel part of S on lattice
    for i = 1:samplelength
        S_parallel += sqrt(S_sample[i][1]^2 + S_sample[i][2]^2)
    end
    S_parallel = S_parallel / samplelength
    χ_xy = (S2xyz_avg[1] + S2xyz_avg[2] - S_parallel^2) / (L1 * L2)# * 1/Temp
    return χ_xy
end

"""
    Average of staggered S_stag^2 of the whole lattice, without normalizing it with respect to system size N. <S^2_α> =  ∑_i,j (-1)^(i+j) <S_α,i S_α,j> (α=x,y,z being three directions)

    return [<S_stag^2_x>, <S_stag^2_y>, <S_stag^2_z>]
"""
function S2xyz_stag(S_stag_sample::Array; L1, L2)
    samplelength = length(S_stag_sample)
    S2xyz_stag_avg = [0, 0, 0]
    for i = 1:samplelength
        S2xyz_stag_avg += S_stag_sample[i] .* S_stag_sample[i]
    end
    S2xyz_stag_avg = S2xyz_stag_avg / samplelength
    return S2xyz_stag_avg
end

function χ_xy_stag(S_stag_sample::Array; L1, L2, Temp)
    samplelength = length(S_stag_sample)
    S2xyz_stag_avg = S2xyz_stag(S_stag_sample; L1=L1, L2=L2)
    S_stag_parallel = 0 #parallel part of S on lattice
    for i = 1:samplelength
        S_stag_parallel += sqrt(S_stag_sample[i][1]^2 + S_stag_sample[i][2]^2)
    end
    S_stag_parallel = S_stag_parallel / samplelength
    χ_xy_stag = (S2xyz_stag_avg[1] + S2xyz_stag_avg[2] - S_stag_parallel^2) / (L1 * L2)# * 1/Temp
    return χ_xy_stag
end

function B(S_sample::Array; S2xyz_stag, L1, L2)
    samplelength = length(S_sample)
    S4_avg = 0
    for i = 1:samplelength
        S4_avg += (S_sample[i][1] * S_sample[i][1] + S_sample[i][2] * S_sample[i][2] + S_sample[i][3] * S_sample[i][3])^2
    end
    S4_avg = S4_avg / samplelength
    S2_avg = S2xyz_stag[1] + S2xyz_stag[2] + S2xyz_stag[3]
    B = 1 - S4_avg / ((1 + 2 / (3 - 1)) * S2_avg^2)
    return S4_avg
end

function E_hsty_sampling(E_sample::Array; L1, L2, sweep_btw_sample=10^3)
    samplelength = length(E_sample)
    @assert(samplelength >= sweep_btw_sample, "sample length is shorter than sweep_btw_sample, need longer sample length or shorter sweep_btw_sample")
    E_samples = []
    for i in 1:Int(floor(samplelength / 10^3))
        push!(E_samples, E_sample[i*sweep_btw_sample] / (L1 * L2))
    end
    return E_samples
end
