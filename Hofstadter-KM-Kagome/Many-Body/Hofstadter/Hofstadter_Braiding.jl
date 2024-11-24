using Pkg
println("\n",Pkg.status(), VERSION)

using QuantumOptics, SparseArrays, Plots, LinearAlgebra, ProgressMeter, Revise, OffsetArrays, Interpolations, LaTeXStrings, Base.Threads
includet("../Scripts/FirstBandApproximation.jl")
includet("../Scripts/ManyBody.jl")
includet("Hofstadter_SP.jl")
includet("../Scripts/Impurity.jl")
includet("../Scripts/Braiding.jl")

Nx = 6
Ny = 6
p = 1
q = 6
ϕ = p/q
pn = 2
U = 1
V = 0.5
Vrand = 1e-5

N_Site = Nx*Ny
N = N_Site
NPhi0 = Int(Nx*Ny*(p/q))
Cut_Off = NPhi0

co, lat_plot = plot_square_lattice(N_Site, Nx, Ny)
lat_plot
#savefig(lat_plot,"./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/Lattice.png")

matrix = Hofstadter_SP(Nx, Ny, p / q, 0)

V0 = [V, V]
Imp_Site = [13,34]
Impurity_Data = Impurity(V0, Imp_Site)
NPin = 2 # It just defines # of degeneracy 

# SP SPECTRUM
e1, psi1 = eigen(Matrix((matrix+matrix')/2))
title_sp = latexstring("KM Single Particle Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
sp_spectrum = scatter(e1, xlabel=L"n", ylabel=L"E", label=false, title=title_sp)

Degeneracy, nu_eff = ground_degeneracy(Nx, Ny, p, q, NPin, pn)
ParameterInfo(NPin, pn, Nx, Ny, p, q)

#error("STOP!")

# HARD-CORE IMPURTIY HAMILTONIAN
HardCore = true 
println("\n Constructing MB Operator...")
# NO PROJ
HHubbard, basis_mb, basis_sp = H_Hubbard(N, pn, matrix, HardCore)
#PROJ
#= Cut_Off = NPhi0
HHubbard, P, Pt, basis_cut_mb = H_Hubbard_Projection(N, pn, matrix, Cut_Off, HardCore) =#
println("\n End of the MB Operator")
# NO PROJ
#num_list = get_num_list(N)
#Number_MB_Operator_List = get_num_mb_op(N, basis_sp, num_list, basis_mb)
Number_MB_Operator_List = [number(basis_mb, site) for site in 1:N]
# PROJ
#= num_sub_list = get_num_sub_list(N_Site, P, Pt)
basis_cut_sp = NLevelBasis(Cut_Off)
Sub_Number_MB_Operator_List = get_num_mb_op(N_Site, basis_cut_sp, num_sub_list, basis_cut_mb);
 =#

Himpurity = Imp_H(Sub_Number_MB_Operator_List, Impurity_Data, Vrand)

HTotal = HHubbard + Himpurity

ϵ, UPsi = eigenstates(HTotal, 100, info=false)
title_mb_imp = latexstring("Hofstadter MB Impurity Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
mb_imp_spectrum = scatter(real(ϵ), xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
mb_imp_spectrum_zoom = scatter(real(ϵ)[1:Degeneracy+5], xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
#savefig(mb_imp_spectrum_zoom, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/mb_imp_spectrum_zoom.png")

UPsi_first = hcat([UPsi[i].data for i in 1:Degeneracy] ...)

Initial_Configuration = heatmap(Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Number_MB_Operator_List, basis_mb, UPsi)', title="Initial Configuration of QH",c=:roma)
#savefig(Initial_Configuration, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/Initial_Configuration.png")

####################
# BRAIDING SECTION #
####################

# EXCHANGE AND DOUBLE EXCHANGE
# Simultaneously
start_site_1 = Imp_Site[1] # top-left corner
start_site_2 = Imp_Site[2] # bottom-right corner
lens = [4, 4]
# 1.
path_no_1 = collect(start_site_1:1:start_site_1+lens[1])
for i in 1:lens[1]-1
    push!(path_no_1, path_no_1[end])
end
# 2.
path_no_2 = collect(start_site_2:-1:start_site_2-lens[1])
for i in 1:lens[1]
    pushfirst!(path_no_2, start_site_2)
end
for i in 1:lens[1]-1
    push!(path_no_2, path_no_2[end])
end
# 3.
path_no_3 = collect(path_no_1[end]:-Nx:path_no_1[end]-lens[2]*Nx)
for i in 1:lens[1]
    push!(path_no_3, path_no_3[end])
end
# 4.
path_no_4 = collect(path_no_2[end]:Nx:path_no_2[end]+lens[2]*Nx)
exch_path_1 = cat(path_no_1, path_no_3, dims=1)
exch_path_2 = cat(path_no_2, path_no_4, dims=1)

double_exch_path_1 = reduce(vcat, [exch_path_1, exch_path_2[2:end]])
double_exch_path_2 = reduce(vcat, [exch_path_2, exch_path_1[2:end]])
# Simultaneously End

# BRAIDING PATH
fixed_point = 73
moving_point = 33
path_1, path_2 = braiding_path(moving_point, fixed_point, Nx, Ny, co)
#path_1 = reverse([12,13,14,15,23,31,39,38,37,36,28,20,12])
#path_1 = reverse(path_1)
#path_1 = [41,50,59,68,77,5,14,23,32,41]
#path_1 = [1,2,3,4,5,6,7,8,17,26,35,44,53,62,71,70,69,68,67,66,65,64,55,46,37,28,19,10,1]
#path_1 = [41,42,43,44,45,37,38,39,40,41,32,23,14,5,77,68,59,50,41,40,39,38,37,45,44,43,42,41,50,59,68,77,5,14,23,32,41]
#path_1 = [41,42,43,44,45,37,38,39,40,41].-9
#path_1 = [58,59,60,61,52,43,34,33,32,31,40,49,58]
path_1 = [33,42,51,60,69,78,6,15,24,33, 32,31,30,29,28,36,35,34, 25,16,7,79,70,61,52,43,34, 35,36,28,29,30,31,32,33]

###########################
# T2^(-1)T1^(-1)T2T1 PATH #
###########################
#T2_path = reduce(vcat,[13,14,15,16,17,18,13, 13,13,13,13,13,13, 18,17,16,15,14,13, 13,13,13,13,13,13])
#T1_path = reduce(vcat,[34,34,34,34,34,34,34, 28,22,16,10,4, 34, 34,34,34,34,34,34,  4,10,16,22,28,34])
T2_path = reduce(vcat,[37,38,39,40,41,42,43,44,45,37, repeat([37], Ny), 45,44,43,42,41,40,39,38,37, repeat([37], Ny)])
T1_path = reduce(vcat,[repeat([77], Nx+1),   68,59,50,41,32,23,14,5,77, repeat([77], Nx),  5,14,23,32,41,50,59,68,77])

qh_paths = plot_paths(co, T2_path, T1_path)
#savefig(qh_paths, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/qh_paths.png")

V_initial = 0
step_break = 0.05
V_final = 1
STEP = V_initial:step_break:V_final
println("Step size between two sites:",length(STEP))

# CHOOSE A PATH !
PATH1 = T2_path
PATH2 = T1_path
V1 = V0[2] # path1 (moved)
V2 = V0[1] # path2 (fixed)
ψ_first, ψ_mat, Converge, ψ_op, imp_data, ψ_mat_list, V0_step_list, band_width, gap = get_phases(UPsi_first, V1, V2, Vrand, PATH1, PATH2, STEP, HHubbard, Number_MB_Operator_List, Degeneracy)

conv_plot = plot(Converge, title="Converge Values \n Step Size=$(length(STEP))", xlabel=L"n", ylabel=L"C = abs|<ψ|\tilde{\psi}>|", legend=false)
#savefig(conv_plot, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/conv_plot.png")

width_gap = plot(real(band_width./gap),title=L"δ/Δ", xlabel="Step")
#savefig(width_gap, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/width_gap.png")

ParamInfo, Degeneracy, nu_eff = ParameterInfo(NPin, pn, Nx, Ny, p, q)

BerryMatrix = ψ_mat'*ψ_first
BerryE, BerryU = eigen(BerryMatrix)

println("\n BerryE=", BerryE)

println("\n |BerryE|=", abs.(BerryE))

θ_berry_1 = angle.(BerryE)./pi

println("\n Phases=", θ_berry_1)

mod(sum(θ_berry_1)/14,2)

println("\n Average Phase (angle) =", sum(mod.(θ_berry_1,2))/Degeneracy)

θ_berry_2 = real(log.(BerryE)./(im*pi))

println("\n Berry Phase = ", sum(mod.(θ_berry_2,2))/Degeneracy)

N_enc = lens[1]*lens[2]
N_enc = 9*8
θ_AB = mod(2 * ϕ * nu_eff * N_enc, 2)

println("\n nu_eff = ",float(nu_eff))

println("\n 2*nu_eff = ",2*float(nu_eff))

println("\n AB Phase = ", θ_AB)

println("\n End of Hofstadter Braiding \n")

# error("STOP BEFORE MOVIE")

#########
# MOVIE #
#########
@gif for i in 1:length(ψ_op)
    heatmap(Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Number_MB_Operator_List, basis_mb, ψ_op[i])', title="Initial Configuration of QH",c=:roma)
end 