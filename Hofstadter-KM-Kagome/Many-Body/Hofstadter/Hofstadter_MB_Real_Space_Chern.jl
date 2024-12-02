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
Nx = 8
Ny = 8
p = 1
q = 8
ϕ = p/q
pn = 3
N = Nx*Ny
Vrand = 0
V = 2
NPhi0 = Int(Nx*Ny*(p/q))

# SINGLE PARTICLE SECTION 
matrix = Hofstadter_SP(Nx, Ny, p / q, 0)
ϵ_sp, psi_sp = eigen(Matrix(matrix))
scatter(real(ϵ_sp))

# SINGLE PARTICLE WITH TWISTED ANGLE SECTION
Tx = Ty = 0
matrix_theta = HSP_T(Nx, Ny, ϕ, Tx, Ty)
ϵ_sp_theta, psi_sp_theta = eigen(Matrix(matrix_theta))
scatter(real(ϵ_sp_theta))

# MANY-BODY SECTOR
HardCore = true
H_hubbard, basis_mb = H_Hubbard(N, pn, matrix_theta, HardCore)
E_hubbard, Psi_hubbard = eigenstates(dense(H_hubbard))
scatter(E_hubbard)

# IMPURITY SECTOR
V0 = [V]; Imp_Site = [15]
NPin = 1
Impurity_Data = Impurity(V0, Imp_Site)
num_list = get_num_list(N)
basis_sp = NLevelBasis(N)
Number_MB_Operator_List = [number(basis_mb, site) for site in 1:N]
Impurity_H = Imp_H(Number_MB_Operator_List, Impurity_Data, Vrand)
Total_H = H_hubbard + Impurity_H
Degeneracy, nu_eff = ground_degeneracy(Nx, Ny, p, q, NPin, pn)
ϵ_imp, psi_imp = eigenstates(Total_H, Degeneracy+2, info=false)
scatter(real(ϵ_imp))

# MB CHERN NUMBER FOR GS MANIFOLD SECTION
Nθ_x = 10
Nθ_y = 10

function get_hard_core_GS_Chern_Num(Nx, Ny, N, pn, ϕ, Impurity_H, Degeneracy, Nθ_x, Nθ_y)
    dx = 2*pi/Nθ_x; dy = dx
    Tx = collect(0:dx:2*pi-dx); Ty = collect(0:dy:2*pi-dy)
    HardCore = true

    θ_Sum = 0

    @threads for tx in Tx
        for ty in Ty
            # (θx, θy)
            matrix_theta = HSP_T(Nx, Ny, ϕ, tx, ty)
            H_hubbard, _ = H_Hubbard(N, pn, matrix_theta, HardCore)
            Total_H = H_hubbard + Impurity_H
            _, v1 = eigenstates(Total_H, Degeneracy, info=false)
            v1 = hcat([v1[i].data for i in 1:Degeneracy] ...)
            
            # (θx+dx, θy)
            matrix_theta = HSP_T(Nx, Ny, ϕ, tx+dx, ty)
            H_hubbard, _ = H_Hubbard(N, pn, matrix_theta, HardCore)
            Total_H = H_hubbard + Impurity_H
            _, v2 = eigenstates(Total_H, Degeneracy, info=false)
            v2 = hcat([v2[i].data for i in 1:Degeneracy] ...)

            # (θx, θy+dy)
            matrix_theta = HSP_T(Nx, Ny, ϕ, tx, ty+dy)
            H_hubbard, _ = H_Hubbard(N, pn, matrix_theta, HardCore)
            Total_H = H_hubbard + Impurity_H
            _, v3 = eigenstates(Total_H, Degeneracy, info=false)
            v3 = hcat([v3[i].data for i in 1:Degeneracy] ...)

            # (θx+dx, θy+dy)
            matrix_theta = HSP_T(Nx, Ny, ϕ, tx+dx, ty+dy)
            H_hubbard, _ = H_Hubbard(N, pn, matrix_theta, HardCore)
            Total_H = H_hubbard + Impurity_H
            _, v4 = eigenstates(Total_H, Degeneracy, info=false)
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
CHERN = get_hard_core_GS_Chern_Num(Nx, Ny, N, pn, ϕ, Impurity_H, Degeneracy, Nθ_x, Nθ_y)
println("\n => MB Real Space Chern Number of GS Manifold is: ", CHERN)