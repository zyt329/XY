
"""
    Binning analysis of data, without using external software
"""
function binning_anal(sample::Sample; num_bin=10)
    # E,E²,M,M² parameters
    num_meas = length(sample.E)
    num_per_bin = Int(floor(num_meas / num_bin))

    # E
    E_mean = mean(sample.E)
    E_std = std([mean(sample.E[((i-1)*num_per_bin+1):i*num_per_bin]) for i in 1:num_bin])

    # E²
    E²_mean = mean(sample.E²)
    E²_std = std([mean(sample.E²[((i-1)*num_per_bin+1):i*num_per_bin]) for i in 1:num_bin])

    # M
    M_mean = mean(sample.M)
    M_std = std([mean(sample.M[((i-1)*num_per_bin+1):i*num_per_bin]) for i in 1:num_bin])

    # M²
    M²_mean = mean(sample.M²)
    M²_std = std([mean(sample.M²[((i-1)*num_per_bin+1):i*num_per_bin]) for i in 1:num_bin])

    # γ parameters
    γ_num_meas = length(sample.γ₁)
    γ_num_per_bin = Int(floor(γ_num_meas / num_bin))

    # γ₁
    γ₁_mean = mean(sample.γ₁)
    γ₁_std = std([mean(sample.γ₁[((i-1)*γ_num_per_bin+1):i*γ_num_per_bin]) for i in 1:num_bin])

    # γ₂
    γ₂_mean = mean(sample.γ₂)
    γ₂_std = std([mean(sample.γ₂[((i-1)*γ_num_per_bin+1):i*γ_num_per_bin]) for i in 1:num_bin])

    # γ₁₂
    γ₁₂_mean = mean(sample.γ₁₂)
    γ₁₂_std = std([mean(sample.γ₁₂[((i-1)*γ_num_per_bin+1):i*γ_num_per_bin]) for i in 1:num_bin])

    return Dict("E_mean" => E_mean, "E_std" => E_std, "E2_mean" => E²_mean, "E2_std" => E²_std, "M_mean" => M_mean, "M_std" => M_std, "M2_mean" => M²_mean, "M2_std" => M²_std, "gamma1_mean" => γ₁_mean, "gamma1_std" => γ₁_std, "gamma2_mean" => γ₂_mean, "gamma2_std" => γ₂_std, "gamma12_mean" => γ₁₂_mean, "gamma12_std" => γ₁₂_std)
end