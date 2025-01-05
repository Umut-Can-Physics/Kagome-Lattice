using LinearAlgebra
using OffsetArrays
using ProgressMeter
using QuantumOptics
using Revise
using SparseArrays
includet("../Scripts/FirstBandApproximation.jl")
includet("../Scripts/ManyBody.jl")
includet("Hofstadter_SP.jl")

Nx = 6
Ny = 6
p = 1
q = 6
pn = 2
U = 1

N = Nx*Ny
NPhi0 = Int(Nx*Ny*(p/q))

matrix = Hofstadter_SP(Nx, Ny, p / q, 0);

ϵ_sp, λ_sp = eigen(Matrix(matrix))

HardCore = true
U = 0

basis_sp = NLevelBasis(N)
H = Sp_Op(basis_sp, matrix)
basis_mb = get_Bosonic_MB_Basis(basis_sp, pn, true)

@btime H_full_hard_core, basis_mb = H_Hubbard(N, U, pn, matrix, HardCore);
H_full_hard_core, basis_mb = H_Hubbard(N, U, pn, matrix, HardCore);

#E_full_hard_core, psi_full_hard_core = eigenstates(dense(H_full_hard_core));

@btime Hhop_mb = get_mb_hopping(basis_mb, H)
Hhop_mb = get_mb_hopping(basis_mb, H)
