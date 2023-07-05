import Base.*

"""
Represent a XY spin on a site. Angles are in radian([0,2π))

fields:
θ::Float64 : θ angle the spin makes with z axis, takes value [0,π)
"""
mutable struct Spin
    θ::Float64
    function Spin(θ=0)
        new(θ)
    end
end

"""
Define spin dot product. Returns the dot product of two spins.
"""
function *(sp1::Spin, sp2::Spin)

    return cos(sp1.θ - sp2.θ)
end

"""
    Computing cross product of 2 spins
"""
function ⊗(sp1::Spin, sp2::Spin)
    return sin(sp2.θ - sp1.θ)
end


"""
    Rotate the spin by Δθ and Δϕ.
    Subtlety:
    1. when θ is rotated to be larger than π or smaller than 0, it should also add to ϕ an additional π.
"""
function rotate!(sp::Spin, Δθ::T) where {T<:Real}
    # update θ
    sp.θ += Δθ

    # mod 2π
    sp.θ = mod(sp.θ, 2π)

    nothing
end


