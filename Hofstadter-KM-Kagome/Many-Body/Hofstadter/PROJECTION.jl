using LinearAlgebra
using OffsetArrays
using Plots
using ProgressMeter
using QuantumOptics
using Revise
using SparseArrays
using LaTeXStrings
using Base.Threads
includet("../Scripts/FirstBandApproximation.jl")
includet("../Scripts/ManyBody.jl")
includet("../Scripts/Impurity.jl")
includet("Hofstadter_SP.jl")
includet("Hofstadter_SP(TwistedAngle).jl")

# PARAMETERS
Nx = 5
Ny = 5
p = 1
q = 5
ϕ = p/q
pn = 2
N = Nx*Ny
Vrand = 0
V = 2

N = Nx*Ny
b = NLevelBasis(N)
s = bosonstates(b, [pn])

# SINGLE PARTICLE WITH TWISTED ANGLE SECTION
Tx = Ty = 0
matrix_theta = HSP_T(Nx, Ny, ϕ, Tx, Ty)
matrix_thetaa = Operator(b, matrix_theta)
E, psi_sp_theta = eigenstates(matrix_thetaa)

ψ_θ = psi_sp_theta[1:5]
P_θ = ψ_θ'ψ_θ
Two_Body_Int = sparse(zeros(5, 5))


sub_basis = SubspaceBasis(b, ψ_θ)
P = projector(sub_basis, b)
Pt = dagger(P)

C_tilde_θ(i) = create()
for i in N
    for j in N
        
    end
end
