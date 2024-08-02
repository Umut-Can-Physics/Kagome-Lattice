#   INTERACTING CASE
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

# executeme

# using NBInclude
# @nbinclude("Hofstadter Single Particle in Julia.ipynb")

include("Hofstadter_Single_Particle_in_Julia.jl")

using QuantumOptics

# Nx = 3; Ny = 3; N=Nx*Ny; q = 3; alpha=1/q
# PN = 2
# U = 2;

# Nx = 4; Ny = 4; N=Nx*Ny; q = Ny; alpha=1/q
# PN = 2
# U = 2;

# Nx = 4; Ny = 2; N=Nx*Ny; q=Ny; alpha=1/q
# PN = 2
# U = 2
# TSize = 15;

# Proj. yapmadan bunu hesaplamanın yöntemi yok! 88k state geliyor :D
# Nx = 8; Ny = 10; N=Nx*Ny; q = Ny; alpha=1/q
# PN = 3
# U = 2;

# executeme

f#= unction get_Bosonic_MB_Basis(N,PN)
   
    NBasis = NLevelBasis(N)
    NStates = bosonstates(NBasis, PN)
    
    NBasis_MB = ManyBodyBasis(NBasis, NStates)
    
    return NBasis_MB, NBasis
end =#

# basis_mb, basis = get_Bosonic_MB_Basis(N,PN)

# executeme

function get_Kinetic_Part(N, MB_Basis, Sp_Op)
    KT = SparseOperator(MB_Basis)
    for m in 1:N
        for n in 1:N
            KT = KT + Sp_Op[m,n] * transition(MB_Basis, m, n)
        end
    end
    
    return KT
end

# basis2 = basis ⊗ basis

# interaction : at_i at_i a_i a_i = at_i a_i at_i a_i - at_i a_i = n_i n_i - n_i
    
# Vint2 = SparseOperator(basis2)

# for n in 1:N
#     Vint2 += U/2*transition(basis,n,n)⊗transition(basis,n,n)
# end

# Vint_mb = manybodyoperator(basis_mb, Vint2)

# executeme

function get_Int_Part(N, MB_Basis, U)
    IT = SparseOperator(MB_Basis)
    for m in 1:N
        IT = IT + U/2 * number(MB_Basis, m) * ( number(MB_Basis, m) - identityoperator(MB_Basis) ) 
    end
    
    return IT
end

# Int_mb = get_Int_Part(N, basis_mb, U)

# executeme

function Hofstadter_Finite_U(Nx, Ny, alpha, PN, U)
    
    N = Nx*Ny
    
    MB_Basis, Basis = get_Bosonic_MB_Basis(N,PN)
    
    Sp_Op = Hofstadter_SP(Nx, Ny, alpha, 0)
    Kin = get_Kinetic_Part(N, MB_Basis, Sp_Op)
    
    Int = get_Int_Part(N, MB_Basis, U)
    
    H = Kin + Int
    
    return H
end

#     1. If particle number=1 and U=0, Hofstadter Finite U energies has to be
#        equal to Hofstadter Single Particle energies.

#U=0
#H_mb = Hofstadter_Finite_U(3, 3, 1/3, 2, U);

# using Plots
# scatter(eigenenergies(dense(H_mb)))

# eigenenergies(dense(Hofstadter_Finite_U(Nx,Ny,1/q,PN,U)))

# using LinearAlgebra
# eigen(Hofstadter_SP(4, 4, 1/4, 0))

function get_Fermionic_MB_Basis(N,PN)
    b_hard = NLevelBasis(N)
    states_hard = fermionstates(b_hard, [PN])
    b_mb_hard = ManyBodyBasis(b_hard, states_hard)
    
    return b_mb_hard
end

function Hofstadter_Hard_Core(Nx, Ny, alpha, PN, U)
    
    N = Nx*Ny
    
    MB_Basis = get_Fermionic_MB_Basis(N,PN)
    
    Sp_Op = Hofstadter_SP(Nx, Ny, alpha, 0)
    Kin = get_Kinetic_Part(N, MB_Basis, Sp_Op)
    
    Int = get_Int_Part(N, MB_Basis, U)
    
    H = Kin + Int
    
    return H
end

#   If U>>1, Hofstadter Finite U energies converges at Hofstadter Hard Core energies.

using Plots

# U = 100

# E1 = eigenenergies(dense(Hofstadter_Hard_Core(Nx, Ny, alpha, PN, U)))
# E2 = eigenenergies(dense(Hofstadter_Finite_U(Nx, Ny, alpha, PN, U)))

# plot(1:length(E1), E1, seriestype=:scatter, markershape=:star5, markersize=6, label="Hard-Core")
# plot!(1:length(E2), E2, seriestype=:scatter, label="Finite-U")

# xlabel!("n");ylabel!("E")

# E2

#   CHERN (INTERACTING CASE)
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

# @nbinclude("Hofstadter Single Particle in Theta Space.ipynb")

# n1 = 1
# n2 = 2

# function Hofstadter_Finite_U_C(Nx, Ny, alpha, Tx, Ty, PN, U)
    
#     N = Nx*Ny
    
#     MB_Basis, Basis = get_Bosonic_MB_Basis(N,PN)
    
#     Sp_Op = HSP_T(Nx, Ny, alpha, Tx, Ty, 0)
#     Kin = get_Kinetic_Part(N, MB_Basis, Sp_Op)
    
#     Int = get_Int_Part(N, MB_Basis, U)
    
#     H = Kin + Int
    
#     return H
# end

# # Twist Angle Parameter Space #
# dx=2*pi/TSize
# dy=dx
# Tx=collect(range(start=0, stop=2*pi-dx, step=dx))
# Ty=collect(range(start=0, stop=2*pi-dy, step=dy))
# # ---- #

# # Link Variable and Berry Curvature #
# Sum=0
# for tx in range(start=1, stop=length(Tx))
#     for ty in range(start=1, stop=length(Ty))
        
#         H_mb = Hofstadter_Finite_U_C(Nx, Ny, alpha, Tx[tx], Ty[ty], PN, U)
#         w1, v1 = eigen(dense(H_mb).data)
#         # i = sortperm(w1, by=real);w1 = w1[i];v1 = v1[:,i]
#         v1 = v1[:,n1:n2]  
#         #------------------------------------
#         H_mb = Hofstadter_Finite_U_C(Nx, Ny, alpha, Tx[tx]+dx, Ty[ty], PN, U)
#         w2, v2 = eigen(dense(H_mb).data)
#         #i = sortperm(w2, by=real);w2 = w2[i];v2 = v2[:,i]
#         v2 = v2[:,n1:n2]
#         #------------------------------------
#         H_mb = Hofstadter_Finite_U_C(Nx, Ny, alpha, Tx[tx], Ty[ty]+dy, PN, U)
#         w3, v3 = eigen(dense(H_mb).data)
#         #i = sortperm(w3, by=real);w3 = w3[i];v3 = v3[:,i]
#         v3 = v3[:,n1:n2]
#         #------------------------------------
#         H_mb = Hofstadter_Finite_U_C(Nx, Ny, alpha, Tx[tx]+dx, Ty[ty]+dy, PN, U)
#         w4, v4 = eigen(dense(H_mb).data)
#         #i = sortperm(w4, by=real);w4 = w4[i];v4 = v4[:,i]
#         v4 = v4[:,n1:n2]
#         #----------LINK VARIABLES------------
#         U1=det(adjoint(v1)*v2)
#         U1=U1/abs(U1)
#         U2=det(adjoint(v2)*v4)
#         U2=U2/abs(U2)
#         U3=det(adjoint(v3)*v4)
#         U3=U3/abs(U3)
#         U4=det(adjoint(v1)*v3)
#         U4=U4/abs(U4)
#         #----------BERRY CURVATURE-----------
#         F=log(U1*U2*1/U3*1/U4)
#         Sum=Sum+F 
#     end
# end
# # ---- #

# print("Chern Number is: ", 1/(2*pi*1im)*Sum)

# U = 2
# PN = 2
# H_Finite = Hofstadter_Finite_U_C(Nx, Ny, alpha, Tx[2], Ty[2], PN, U)
# Ee, Uu = eigenenergies(dense(H_Finite ))

# 4
# -3
# 6.0046910114196645
# -6.9939547502816435
# -3
# 4