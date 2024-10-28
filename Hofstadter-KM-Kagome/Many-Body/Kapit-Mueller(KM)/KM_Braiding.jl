#using Pkg
#println("**Packages and VERSION info** \n",Pkg.status(), VERSION)

using QuantumOptics, SparseArrays, Plots, LinearAlgebra, ProgressMeter, Revise, DataFrames, LaTeXStrings

includet("../Scripts/FirstBandApproximation.jl"); 
includet("../Scripts/ManyBody.jl"); 
includet("KM_Model.jl"); 
includet("../Scripts/Impurity.jl"); 
includet("../Scripts/Braiding.jl")

Nx = 6
Ny = 6
p = 1
q = 3
ϕ=p/q
pn = 3
U = 1 # on-site interaction potential
V = 1 # Impurity potential
Vrand = 1e-6

N_Site = Nx*Ny
N = N_Site
t = 1
NPhi0 = Int(Nx*Ny*(p/q))
Cut_Off = NPhi0

co = plot_square_lattice(N_Site, Nx, Ny, path="./Braiding_Data/square_lattice.png")

matrix = KM(Nx, Ny, t, p, q)

V0 = [V, V]; 
Imp_Site = [13, 7];
Impurity_Data = Impurity(V0, Imp_Site);
NPin = 2 # It just defines # of degeneracy 

#for imp in 1:length(Impurity_Data.V0)
#    matrix[Impurity_Data.Imp_Site[imp],Impurity_Data.Imp_Site[imp]] += #Impurity_Data.V0[imp]
#end
#matrix += diagm(repeat([Vrand],N))

# SP SPECTRUM
e1, psi1 = eigen(Matrix((matrix+matrix')/2))
title_sp = latexstring("KM Single Particle Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
sp_spectrum = scatter(e1, xlabel=L"n", ylabel=L"E", label=false, title=title_sp)
#savefig(sp_spectrum,"./Braiding_Data/sp_spectrum.png")

#Degeneracy, _, _, _ = ground_degeneracy(Nx, Ny, p, q, NPin, pn)
Degeneracy = 50

#error("stopping...")

println("**Parameters**")
ParameterInfo(NPin, pn, Nx, Ny, p, q)

# HARD-CORE IMPURTIY HAMILTONIAN
HardCore = true 
HHubbard, basis_mb, basis_sp = H_Hubbard(N, pn, matrix, HardCore)
num_list = get_num_list(N)
#Number_MB_Operator_List = get_num_mb_op(N, basis_sp, num_list, basis_mb)
Number_MB_Operator_List = [number(basis_mb, site) for site in 1:N]

Himpurity = Imp_H(Number_MB_Operator_List, Impurity_Data, Vrand)

HTotal = HHubbard + Himpurity

#= ϵϵ, λλ = eigenstates(dense(HTotal))
title_mb = latexstring("KM MB Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
mb_spectrum = scatter(ϵϵ, xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb)
savefig(mb_spectrum,"Hofstadter-KM-Kagome/Many-Body/Kapit-Mueller(KM)/Braiding_Data/mb_spectrum.png") =#

ϵ, UPsi = eigenstates(HTotal, Degeneracy+5)
title_mb_imp = latexstring("KM MB Impurity Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
mb_imp_spectrum = scatter(ϵ, xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
mb_imp_spectrum_zoom = scatter(ϵ[1:Degeneracy+5], xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
#savefig(mb_imp_spectrum_zoom,"./Braiding_Data/mb_imp_spectrum_zoom.png")

UPsi_first = hcat([UPsi[i].data for i in 1:Degeneracy] ...)

Initial_Configuration = heatmap(Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Number_MB_Operator_List, basis_mb, UPsi)', title="Initial Configuration of QH",c=:roma)
#savefig(Initial_Configuration, "./Braiding_Data/Initial_Configuration.png")

#error("stopping...")

####################
# BRAIDING SECTION #
####################

# EXCHANGE & DOUBLE EXCHANGE PATH (One by One & Simultaneously)
start_point_1 = 13
start_point_2 = 7
lens = [2, 2, 2, 2]
rec_path_1, rec_path_2 = exchange_path(lens, start_point_1, start_point_2)
rec_path_1 = [13,14,15, 15,15, 11,7, 7,7]
rec_path_2 = [7, 7,7, 6,5, 5,5, 9,13] 
#= rec_path_1 = [13,14,15,15,15,11,7,7,7, 7,7,6,5,5,5,9,13]
rec_path_2 = [7,7,7,6,5,5,5,9,13, 14,15,15,15,11,7,7,7] =# 
#= rec_path_1 = reverse(rec_path_1)
rec_path_2 = reverse(rec_path_2) =#

# Simultaneously
start_site_1 = 13
start_site_2 = 7
lens = [2, 2]
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


#= two_neighbor_sites = [[13,14],[14,15],[7,6],[6,5],[15,11],[11,7],[5,9],[9,13]]
const_qhs = [7,7,15,15,5,5,7,7] =#

# BRAIDING PATH
moving_point = 7
fixed_point = 13
path_1, path_2 = braiding_path(moving_point, fixed_point, Nx, Ny, co)
path_1 = reverse(path_1)
path_2 = reverse(path_2)
path_1 = [7,11,15,16,12,8,4,3,7]
#path_1 = [13,14,15,15,15,11,7,7,7,7,7,7,6,5,5,5,9,13]
#path_2 = [7,7,7,6,5,5,5,9,13,13,14,15,15,15,11,7,7,7]
#path_1 = [31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31]
#path_2 = [22,28,34,10,16,22,21,20,19,24,23,22,16,10,34,28,22,23,24,19,29,21,22]
# plot_paths(co, rec_path_1, rec_path_2)

V_initial = 0
step_break = 0.05
V_final = 1
#-step_break
STEP = V_initial:step_break:V_final
println("Step size between two sites:",length(STEP))

PATH1 = exch_path_1
PATH2 = exch_path_2
ψ_first, ψ_mat, Converge, ψ_op, imp_data, ψ_mat_list, ϵ_list, V0_step_list, band_width, gap = get_phases(UPsi_first, V, Vrand, PATH1, PATH2, STEP, HHubbard, Number_MB_Operator_List, Degeneracy)
#ψ_first, ψ_mat, Converge, ψ_op, imp_data, ψ_mat_list, ϵ_list, V0_step_list = get_phases_simultaneously(λ_firstt, two_neighbor_sites, const_qhs, V, V_rand, STEP, HTotal, Number_MB_Operator_List, Degeneracy)

conv_plot = plot(Converge, title="Converge Values \n Step Size=$(length(STEP))", xlabel=L"n", ylabel=L"C = abs|<ψ|\tilde{\psi}>|", legend=false)
#savefig(conv_plot, "./Braiding_Data/conv_plot.png")

plot(real(band_width./gap),title=L"δ/Δ", xlabel="Step")

BerryMatrix = ψ_mat'*ψ_first
BerryE, BerryU = eigen(BerryMatrix)

println("\n", BerryE)

println("\n |BerryE|=", abs.(BerryE))

println("\n Phases=", angle.(BerryE)./pi)
#sprint(show,angle.(BerryE)./pi, context=:compact => true)

println("\n Avarage Phase=", sum(angle.(BerryE)./pi)/Degeneracy)

println("\n End of KM_Braiding \n")

# MOVIE
#    @gif for i in 1:800
#        heatmap(Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Number_MB_Operator_List, basis_mb, ψ_op[i])', title="Initial Configuration of QH",c=:roma)
#    end

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