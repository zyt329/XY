
"""
Calculate energy of a certain configuration. The convention for the cross product is that spin at odd sites cross spin at even site.
"""
function Energy(J1, J2, Jz, lattice_info, spins::Array{Spin,3})

    # interaction strength
    J_xy = (J1, J2)

    # lattice information
    xy_neibs = lattice_info.xy_neibs
    z_neibs = lattice_info.z_neibs

    L1 = lattice_info.L1
    L2 = lattice_info.L2
    N = L1 * L2

    # ========  start calculating energy =========
    E = 0

    coords = [(z, i, j) for z in 1:2, i in 1:L1, j in 1:L2]
    # loop over sites
    for coord in coords
        sp = spins[coord...]

        # layer of the spin
        layer = coord[1]

        for neib in xy_neibs[coord...]
            E += J_xy[layer] * (sp * spins[neib...])
        end

        for neib in z_neibs[coord...]
            E += Jz * (sp * spins[neib...])
        end

    end

    # avoid double counting
    E = E / 2

    return E
end



"""
    Magnetization
"""
function Magnetization(spins::Array{Spin,3})

    Mag = [0.0, 0.0]
    for spin in spins
        Mag[1] += cos(spin.θ)
        Mag[2] += sin(spin.θ)
    end

    return Mag

end



"""
    Initialize spins.
"""
function init_spins_config_gen(L1, L2)
    spins = Array{Spin,3}(undef, 2, L1, L2)
    for z in 1:2, i in 1:L1, j in 1:L2
        ϕ_init = 0#rand(Uniform(0, 2π))
        spins[z, i, j] = Spin(ϕ_init)
    end
    return spins
end



"""
Represent a Microstate.

fields:  T, L, steps, dimension, spin
"""
mutable struct Microstate

    # interaction strength
    J1::Real
    J2::Real
    Jz::Real

    # Lattice info
    lattice_info::Lattice

    # spin configuration
    spins::Array{Spin,3}#(undef,3,3,3)

    # Basic measurements
    E::Real
    M::Vector{Float64}

    function Microstate(;
        J1::Real,
        J2::Real,
        Jz::Real,
        L1::Int64,
        L2::Int64
    )
        lattice_info = Lattice(L1, L2)

        spins = init_spins_config_gen(L1, L2)

        E = Energy(J1, J2, Jz, lattice_info, spins)

        M = Magnetization(spins)

        new(J1, J2, Jz, lattice_info, spins, E, M)
    end
end


"""
Calculate the energy difference after rotating one spin.

Parameter:
micstate::Microstate : Input the
numbering::Int64 : the numbering of the rotated spin
"""
function Energy_Diff(micstate::Microstate, coord::Tuple{Int64,Int64,Int64}, sp_new::Spin)

    # interaction strength
    J_xy = [micstate.J1, micstate.J2]
    Jz = micstate.Jz

    # lattice information
    lattice_info = micstate.lattice_info
    xy_neibs = lattice_info.xy_neibs
    z_neibs = lattice_info.z_neibs

    L1 = lattice_info.L1
    L2 = lattice_info.L2
    N = L1 * L2

    # spin configuration
    spins = micstate.spins
    sp_old = spins[coord...]

    ######Start calculating energy difference before and after
    ΔE = 0

    layer = coord[1]
    # loop over its xy neighbors
    for neib in xy_neibs[coord...]
        ΔE += J_xy[layer] * (spins[neib...] * sp_new - spins[neib...] * sp_old)
    end

    # loop over its z neighbors
    for neib in z_neibs[coord...]
        ΔE += Jz * (spins[neib...] * sp_new - spins[neib...] * sp_old)
    end

    return ΔE

end
@debug ΔE = Energy_Diff(micstate, 1, Spin(ϕ=π / 2))

function M_Diff(micstate::Microstate, coord::Tuple{Int64,Int64,Int64}, sp_new::Spin)

    # spin configuration
    spins = micstate.spins
    sp_old = spins[coord...]

    # return ΔM
    return [(cos(sp_new.θ) - cos(sp_old.θ)), (sin(sp_new.θ) - sin(sp_old.θ))]

end


"""
Updating the current Microstate to a new Microstate. Also change the energy and S.

Parameter:
micstate::Microstate : current configuration.
numbering::Int64 : numbering of the changed spin.
sp_new::Spin : updated spin value for the selected spin.
"""
function update(micstate::Microstate, coord::Tuple{Int64,Int64,Int64}, sp_new::Spin)

    micstate.E += Energy_Diff(micstate, coord, sp_new)
    micstate.M += M_Diff(micstate, coord, sp_new)

    micstate.spins[coord...] = sp_new

    nothing
end
# micstate = Microstate(J1=1,J2=1,Jz=1,L1=4,L2=4)
# update(micstate, 1, Spin(ϕ=π / 2))
