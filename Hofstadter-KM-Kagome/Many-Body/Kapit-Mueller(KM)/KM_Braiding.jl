using QuantumOptics; using SparseArrays; using Plots; using LinearAlgebra; using ProgressMeter; using Revise; using DataFrames; using Interpolations; using LaTeXStrings

includet("../Scripts/FirstBandApproximation.jl"); includet("../Scripts/ManyBody.jl"); includet("KM_Model.jl"); includet("../Scripts/Impurity.jl"); includet("../Scripts/Braiding.jl")

Nx = 4
Ny = 4
p = 1
q = 2
ϕ=p/q
pn = 3
U = 4 # on-site interaction potential
V = 10 # Impurity potential
NPin = 2
N_Site = Nx*Ny
N = N_Site
t = 1
NPhi0 = Int(Nx*Ny*(p/q))
PN = vcat( ( [i] for i in 0:pn) ... )
Cut_Off = NPhi0

co, _ = plot_square_lattice(N_Site, Nx, Ny)

matrix = KM(Nx, Ny, t, p, q)

ϵ_sp, ψ_sp = eigen(Matrix((matrix+matrix')/2))

scatter(ϵ_sp)
heatmap(abs.(ψ_sp).^2)
heatmap(reshape(sum(abs.(ψ_sp[:,1:10]).^2,dims=2),(Nx,Ny)))
reshape(sum(abs.(ψ_sp[:,1:10]).^2,dims=2),(Nx,Ny))

V0 = [V, V]; Imp_Site = [13,7]
Impurity_Data = Impurity(V0, Imp_Site)
ParameterInfo(NPin, pn, Nx, Ny, p, q)

HardCore = true 
HTotal, basis_mb = H_Hubbard(N, pn, matrix, HardCore)
num_list = get_num_list(N)
basis_sp = NLevelBasis(N)
Number_MB_Operator_List = get_num_mb_op(N, basis_sp, num_list, basis_mb)
Impurity_H_hard_core = Imp_H(HTotal, Number_MB_Operator_List, Impurity_Data); 

ϵ, λ = eigenstates(dense(HTotal))
scatter(ϵ[1:100], title="Nx=$(Nx), Ny=$(Ny), ϕ=$(p//q), pn=$(pn)", label="Hard-Core=$(HardCore)")

ϵ, λ = eigenstates(Impurity_H_hard_core)
scatter(ϵ[1:100], title="Nx=$(Nx), Ny=$(Ny), p=$(p), q=$(q), pn=$(pn), NPin=$(NPin)", label="Hard-Core=$(HardCore)")

factor = 1
Degeneracy, _, _, _ = ground_degeneracy(Nx, Ny, p, q, NPin, pn)
plot_density(Nx, Ny, Degeneracy, N_Site, Number_MB_Operator_List, basis_mb, λ, factor)

moving_point = Imp_Site[2]
braid_path = braiding_path(moving_point, Nx, Ny, co)
fixed_point = 13

start_point_1 = Impurity_Data.Imp_Site[1]
lens_1 = [2, 2, 2, 2]
dirs_1 = [1,-Ny, -1, Ny]
#= lens_1 = [2, 2]
dirs_1 = [1,-Ny] =#
rec_path_1 = rectangular_path(start_point_1,lens_1,dirs_1)
#rec_path_1 = reverse(rec_path_1)
rec_path_1 = repeat([13],17)
#rec_path_1 = [13,14,15,11,7]