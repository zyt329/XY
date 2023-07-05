function coordinate(n::Int64; L1::Int64=2, L2::Int64=2)
    @assert((n ≤ L1 * L2) && (1 ≤ n), "The numbering (bit position) of a site shouldn't exceed the total number of sites $(L1 * L2), and should be bigger than 0.")
    i::Int64 = Int(ceil(n / L1))
    j::Int64 = mod1(n, L1)  #site i is at i-th row, j-th column
    return (i, j)
end

function numbering(coordinate::Tuple{Int64,Int64}; L1::Int64=2, L2::Int64=2)
    @assert((coordinate[1] ≤ L2) && (coordinate[2] ≤ L1), "The cooridnate should be within the range of the lattice size $L1 by $L2")
    n = (coordinate[1] - 1) * L1 + coordinate[2]
    return n
end

function neibs_xy(coord::Tuple{Int64,Int64,Int64}; L1::Int64=2, L2::Int64=2)
    neibs = Set{Tuple{Int64,Int64,Int64}}()

    # X, Y direction neighbors
    push!(neibs, (coord[1], mod1(coord[2] + 1, L2), coord[3]))
    push!(neibs, (coord[1], mod1(coord[2] - 1, L2), coord[3]))
    push!(neibs, (coord[1], coord[2], mod1(coord[3] + 1, L1)))
    push!(neibs, (coord[1], coord[2], mod1(coord[3] - 1, L1)))

    return neibs
end

function neibs_z(coord::Tuple{Int64,Int64,Int64}; L1::Int64=2, L2::Int64=2)
    neibs = Set{Tuple{Int64,Int64,Int64}}()

    # Z direction neighbors
    push!(neibs, (mod1(coord[1] + 1, 2), coord[2], coord[3]))

    return neibs
end


function neib_xy_list_gen(; L1::Int64=2, L2::Int64=2)
    neib_list = Dict{Tuple{Int64,Int64,Int64},Set{Tuple{Int64,Int64,Int64}}}()
    for z in 1:2, i in 1:L1, j in 1:L2
        coord = (z, i, j)
        push!(neib_list, coord => neibs_xy(coord, L1=L1, L2=L2))
    end
    return neib_list
end


function neib_z_list_gen(; L1::Int64=2, L2::Int64=2)
    neib_list = Dict{Tuple{Int64,Int64,Int64},Set{Tuple{Int64,Int64,Int64}}}()
    for z in 1:2, i in 1:L1, j in 1:L2
        coord = (z, i, j)
        push!(neib_list, coord => neibs_z(coord, L1=L1, L2=L2))
    end
    return neib_list
end

function neib_all_list_gen(; L1::Int64=2, L2::Int64=2)
    neib_list = Dict{Tuple{Int64,Int64,Int64},Set{Tuple{Int64,Int64,Int64}}}()
    for z in 1:2, i in 1:L1, j in 1:L2
        coord = (z, i, j)
        push!(neib_list, coord => union(neibs_xy(coord, L1=L1, L2=L2), neibs_z(coord, L1=L1, L2=L2)))
    end
    return neib_list
end




mutable struct Lattice
    L1::Int64
    L2::Int64
    xy_neibs::Dict{Tuple{Int64,Int64,Int64},Set{Tuple{Int64,Int64,Int64}}}
    z_neibs::Dict{Tuple{Int64,Int64,Int64},Set{Tuple{Int64,Int64,Int64}}}
    all_neibs::Dict{Tuple{Int64,Int64,Int64},Set{Tuple{Int64,Int64,Int64}}}
    function Lattice(
        L1, L2;
        xy_neibs=neib_xy_list_gen(L1=L1, L2=L2),
        z_neibs=neib_z_list_gen(L1=L1, L2=L2),
        all_neibs=neib_all_list_gen(L1=L1, L2=L2)
    )

        new(L1, L2,
            xy_neibs,
            z_neibs,
            all_neibs,
        )
    end
end

