module XY

# import external dependencies
using LinearAlgebra
using JLD2
using DelimitedFiles
using Dates
using Statistics
using Random
using Distributions
using OffsetArrays: Origin
using OffsetArrays
using Printf
using JSON
using BinningAnalysis

# ========================================================
# ===========     load utility functions     =============
# ========================================================

include("utilities.jl")
export make_indexed_folder

# ========================================================
# ===========     load spin definitions      =============
# ========================================================

include("spin.jl")
export *, Spin, ⊗, rotate!

# ========================================================
# ===========     load spin definitions      =============
# ========================================================

include("bilayer_square_lattice.jl")
export coordinate, numbering, neib, neibs_xy, neibs_z, neib_xy_list_gen, neib_z_list_gen, neib_all_list_gen, init_spins_config_gen, Lattice

# ==============================================================
# ===========     load  Microstate construction    =============
# ==============================================================

include("microstate.jl")
export Energy, Microstate, Energy_Diff, Magnetization, M_Diff, update, spins2array

# ==============================================================
# ===========      load measurement functions     ==============
# ==============================================================

include("measurements.jl")
export Sample

# ==============================================================
# ===========      load MC utility functions       =============
# ==============================================================

include("MC_single_flip.jl")
export accept, MC_single_flip

include("MC_Wolff.jl")
export MC_Wolff

# ===================================================================
# ===========      load binning analysis functions     ==============
# ===================================================================

include("data_analysis.jl")
export binning_anal







end
