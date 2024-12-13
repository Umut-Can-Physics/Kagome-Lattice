##############################################
# PHYSICAL REVIEW LETTERS 122, 146601 (2019) #
##############################################

using Pkg
Pkg.status()
VERSION

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
Nx = 4
Ny = 4
p = 1
q = 4
ϕ = p/q
pn = 2
N = Nx*Ny
Vrand = 0
V = 2

# SINGLE PARTICLE SECTION 
matrix = Hofstadter_SP(Nx, Ny, p / q, 0)
ϵ_sp, psi_sp = eigen(Matrix(matrix))
scatter(real(ϵ_sp))

# SINGLE PARTICLE WITH TWISTED ANGLE SECTION
Tx = Ty = 0
matrix_theta = HSP_T(Nx, Ny, ϕ, Tx, Ty)
ϵ_sp_theta, psi_sp_theta = eigen(Matrix(matrix_theta))
scatter(real(ϵ_sp_theta))

# MANY-BODY SECTION
# HardCore = true
# H_hubbard, basis_mb = H_Hubbard(N, U, pn, matrix_theta, HardCore)

U = 1
HardCore = false
NPhi0 = Int(Nx*Ny*(p/q)); Cut_Off = Nx*Ny
H_hubbard, P, Pt, basis_cut_mb, _ = H_Hubbard_Projection(N, pn, U, matrix, Cut_Off, HardCore);
E_hubbard, Psi_hubbard = eigenstates(dense(H_hubbard))
scatter(E_hubbard, title="Hard-Core= $HardCore, U= $U")

serialize("Chern/H_hubbard.data", H_hubbard)

# IMPURITY SECTION
V0 = [V]; Imp_Site = [15]
NPin = 1
Impurity_Data = Impurity(V0, Imp_Site)
#Number_MB_Operator_List = [number(basis_mb, site) for site in 1:N]
# Hard-Core Limit (U>>Sp Gap)
num_sub_list = get_num_sub_list(N, P, Pt)
basis_cut_sp = NLevelBasis(Cut_Off)
Sub_Number_MB_Operator_List = get_num_mb_op(N, basis_cut_sp, num_sub_list, basis_cut_mb);

Impurity_H = Imp_H(Sub_Number_MB_Operator_List, Impurity_Data, Vrand)
Total_H = H_hubbard + Impurity_H
Degeneracy, nu_eff = ground_degeneracy(Nx, Ny, p, q, NPin, pn)
ϵ_imp, psi_imp = eigenstates(Total_H, Degeneracy+5)
scatter(real(ϵ_imp))

# MB CHERN NUMBER FOR GS MANIFOLD SECTION
Nθ_x = 10
Nθ_y = 10

function get_hard_core_GS_Chern_Num(Nx, Ny, N, pn, ϕ, U, Impurity_H, Degeneracy, Nθ_x, Nθ_y)

    dx = 2*pi/Nθ_x; dy = dx
    Tx = collect(0:dx:2*pi-dx); Ty = collect(0:dy:2*pi-dy)
    HardCore = false

    θ_Sum = 0

    @showprogress for tx in Tx
        for ty in Ty
            # (θx, θy)
            matrix_theta = HSP_T(Nx, Ny, ϕ, tx, ty)
            # H_hubbard, _ = H_Hubbard(N, U, pn, matrix_theta, HardCore)
            H_hubbard, _, _, _, _ = H_Hubbard_Projection(N, pn, U, matrix, Cut_Off, HardCore);
            Total_H = H_hubbard + Impurity_H
            _, v1 = eigenstates(Total_H, Degeneracy)
            v1 = hcat([v1[i].data for i in 1:Degeneracy] ...)
            
            # (θx+dx, θy)
            matrix_theta = HSP_T(Nx, Ny, ϕ, tx+dx, ty)
            #H_hubbard, _ = H_Hubbard(N, U, pn, matrix_theta, HardCore)
            H_hubbard, _, _, _, _ = H_Hubbard_Projection(N, pn, U, matrix, Cut_Off, HardCore);
            Total_H = H_hubbard + Impurity_H
            _, v2 = eigenstates(Total_H, Degeneracy)
            v2 = hcat([v2[i].data for i in 1:Degeneracy] ...)

            # (θx, θy+dy)
            matrix_theta = HSP_T(Nx, Ny, ϕ, tx, ty+dy)
            #H_hubbard, _ = H_Hubbard(N, U, pn, matrix_theta, HardCore)
            H_hubbard, _, _, _, _ = H_Hubbard_Projection(N, pn, U, matrix, Cut_Off, HardCore);
            Total_H = H_hubbard + Impurity_H
            _, v3 = eigenstates(Total_H, Degeneracy)
            v3 = hcat([v3[i].data for i in 1:Degeneracy] ...)

            # (θx+dx, θy+dy)
            matrix_theta = HSP_T(Nx, Ny, ϕ, tx+dx, ty+dy)
            #H_hubbard, _ = H_Hubbard(N, U, pn, matrix_theta, HardCore)
            H_hubbard, _, _, _, _ = H_Hubbard_Projection(N, pn, U, matrix, Cut_Off, HardCore);
            Total_H = H_hubbard + Impurity_H
            _, v4 = eigenstates(Total_H, Degeneracy)
            v4 = hcat([v4[i].data for i in 1:Degeneracy] ...)

            U1 = det(v1'*v2) / abs(det(v1'*v2))
            U2 = det(v2'*v4) / abs(det(v2'*v4))
            U3 = det(v3'*v4) / abs(det(v3'*v4))
            U4 = det(v1'*v3) / abs(det(v1'*v3))

            F=log(U1 * U2 * (1/U3) * (1/U4))
            θ_Sum += F
        end
    end
    return θ_Sum*1/(2*pi*1im)
end
ParameterInfo(NPin, pn, Nx, Ny, p, q)
println("\n Chern per state (ν_eff) = ", nu_eff)
CHERN = get_hard_core_GS_Chern_Num(Nx, Ny, N, pn, ϕ, U, Impurity_H, Degeneracy, Nθ_x, Nθ_y)
println("\n => MB Real Space Chern Number of GS Manifold is: ", CHERN)