using Pkg
println("**Packages and VERSION info** \n",Pkg.status(), VERSION)

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

N_Site = Nx*Ny
N = N_Site
t = 1
NPhi0 = Int(Nx*Ny*(p/q))
PN = vcat( ( [i] for i in 0:pn) ... )
Cut_Off = NPhi0

co, _ = plot_square_lattice(N_Site, Nx, Ny)

matrix = KM(Nx, Ny, t, p, q)

# SP SPECTRUM
ϵ_sp, ψ_sp = eigen(Matrix((matrix+matrix')/2))
title_sp = latexstring("KM Single Particle Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
sp_spectrum = scatter(ϵ_sp, xlabel=L"n", ylabel=L"E", label=false, title=title_sp)
savefig(sp_spectrum,"Hofstadter-KM-Kagome/Many-Body/Kapit-Mueller(KM)/Braiding_Data/sp_spectrum.png")

V0 = [V, V]; Imp_Site = [10, 4]
Impurity_Data = Impurity(V0, Imp_Site)
NPin = 2 # It just defines # of degeneracy 

Degeneracy, _, _, _ = ground_degeneracy(Nx, Ny, p, q, NPin, pn)

println("**Parameters**")
print(ParameterInfo(NPin, pn, Nx, Ny, p, q))

# HARD-CORE IMPURTIY HAMILTONIAN
HardCore = true 
HTotal, basis_mb = H_Hubbard(N, pn, matrix, HardCore)
num_list = get_num_list(N)
basis_sp = NLevelBasis(N)
Number_MB_Operator_List = get_num_mb_op(N, basis_sp, num_list, basis_mb)
Impurity_H_hard_core = Imp_H(HTotal, Number_MB_Operator_List, Impurity_Data)

#= ϵϵ, λλ = eigenstates(dense(HTotal))
title_mb = latexstring("KM MB Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
mb_spectrum = scatter(ϵϵ, xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb)
savefig(mb_spectrum,"Hofstadter-KM-Kagome/Many-Body/Kapit-Mueller(KM)/Braiding_Data/mb_spectrum.png") =#

ϵ, λ = eigenstates(Impurity_H_hard_core)
title_mb_imp = latexstring("KM MB Impurity Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
mb_imp_spectrum = scatter(ϵ, xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
mb_imp_spectrum_zoom = scatter(ϵ[1:Degeneracy+5], xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
savefig(mb_imp_spectrum_zoom,"Hofstadter-KM-Kagome/Many-Body/Kapit-Mueller(KM)/Braiding_Data/mb_imp_spectrum_zoom.png")

Initial_Configuration = heatmap(Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Number_MB_Operator_List, basis_mb, λ)', title="Initial Configuration of QH",c=:roma)
savefig(Initial_Configuration, "Hofstadter-KM-Kagome/Many-Body/Kapit-Mueller(KM)/Braiding_Data/Initial_Configuration.png")

# EXCHANGE & DOUBLE EXCHANGE PATH
start_point_1 = 10
start_point_2 = 4
lens = [2, 2, 2, 2]
rec_path_1, rec_path_2 = exchange_path(lens, start_point_1, start_point_2)


# BRAIDING PATH
moving_point = 7
fixed_point = 13
path_1, path_2 = braiding_path(moving_point, fixed_point, Nx, Ny, co)

plot_paths(co, rec_path_1, rec_path_2)

V_initial = 0
step_break = 0.01
V_final = 1
#-step_break
STEP = V_initial:step_break:V_final
println("Step size between two sites:",length(STEP))

# BRAIDING
PATH1 = rec_path_1
PATH2 = rec_path_2
ψ_first, ψ_mat, Converge, ψ_op, imp_data, ψ_mat_list, ϵ_list = get_phases(Impurity_Data, PATH1, PATH2, basis_mb, STEP, HTotal, Number_MB_Operator_List, Degeneracy)

conv_plot = plot(Converge, title="Converge Values \n Step Size=$(length(STEP))", xlabel=L"n", ylabel=L"C = abs|<ψ|\tilde{\psi}>|", legend=false)
savefig(conv_plot, "Hofstadter-KM-Kagome/Many-Body/Kapit-Mueller(KM)/Braiding_Data/conv_plot.png")

BerryMatrix = ψ_mat'*ψ_first
BerryE, BerryU = eigen(BerryMatrix)
println(BerryE,angle.(BerryE)./pi)

# AB PHASE
#= N_mov = 2
number_of_plaq = 1
θ_AB, ex, charge = th_AB_phase(pn, p, q, NPin, N_mov, number_of_plaq)
θ_exchange = charge*pi
θ_AB =#

#= print(" :: EXCHANGE PATH ::
Eigenvalues: $(BerryE) \n 
SUM and AVARAGE of eigenvalues:
sum(BerryE)= $(sum(BerryE))
sum(BerryE)/$(Degeneracy) = $(sum(BerryE)/Degeneracy) \n
Phases of SUM of eigenvalues:
log(sum(BerryE))/(im*pi) = $(log(sum(BerryE))/(im*pi)) 
log(sum(BerryE))/(2*im*pi) = $(log(sum(BerryE))/(2*im*pi)) \n
Phases of AVARAGE of eigenvalues:
log(sum(BerryE)/$(Degeneracy))/(im*pi) = $(log(sum(BerryE)/Degeneracy)/(im*pi)) 
log(sum(BerryE)/$(Degeneracy))/(2*im*pi) = $(log(sum(BerryE)/Degeneracy)/(2*im*pi)) \n
Phases of eigenvalues:
log.(BerryE)./(im*pi) = $(log.(BerryE)./(im*pi)) 
log.(BerryE)./(2*im*pi) = $(log.(BerryE)./(2*im*pi)) \n
SUM and AVARAGE of Phases of eigenvalues:
sum(log.(BerryE)./(im*pi)) = $(sum(log.(BerryE)./(im*pi)))
sum(log.(BerryE)./(im*pi))/$(Degeneracy) = $(sum(log.(BerryE)./(im*pi))/Degeneracy) \n
----
") =#
#θ_AB= $(θ_AB) | exp(θ_AB)=$(ex) | Q=$(charge) | θ_exchange(Q*pi) = $(θ_exchange)