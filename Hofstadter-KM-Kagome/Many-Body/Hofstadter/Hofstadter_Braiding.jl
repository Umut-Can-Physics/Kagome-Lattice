using Pkg
println("\n",Pkg.status(), VERSION)

using QuantumOptics, SparseArrays, Plots, LinearAlgebra, ProgressMeter, Revise, OffsetArrays, Interpolations, LaTeXStrings
includet("../Scripts/FirstBandApproximation.jl")
includet("../Scripts/ManyBody.jl")
includet("Hofstadter_SP.jl")
includet("../Scripts/Impurity.jl")
includet("../Scripts/Braiding.jl")

Nx = 8
Ny = 8
p = 1
q = 8
ϕ = p/q
pn = 2
U = 1
V = 0.5
NPin = 1
Vrand = 1e-14

N_Site = Nx*Ny
N = N_Site
NPhi0 = Int(Nx*Ny*(p/q))
Cut_Off = NPhi0

co, lat_plot = plot_square_lattice(N_Site, Nx, Ny)
lat_plot
savefig(lat_plot,"./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/Lattice.png")

matrix = Hofstadter_SP(Nx, Ny, p / q, 0)

V0 = [V, V]
Imp_Site = [57, 29]
Impurity_Data = Impurity(V0, Imp_Site)
NPin = 2 # It just defines # of degeneracy 

# SP SPECTRUM
e1, psi1 = eigen(Matrix((matrix+matrix')/2))
title_sp = latexstring("KM Single Particle Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
sp_spectrum = scatter(e1, xlabel=L"n", ylabel=L"E", label=false, title=title_sp)

ParamInfo, Degeneracy, nu_eff = ParameterInfo(NPin, pn, Nx, Ny, p, q)
#Degeneracy = 14

#error("STOP!")

# HARD-CORE IMPURTIY HAMILTONIAN
HardCore = true 
println("\n Constructing MB Operator...")
HHubbard, basis_mb, basis_sp = H_Hubbard(N, pn, matrix, HardCore)
println("\n End of the MB Operator")
num_list = get_num_list(N)
#Number_MB_Operator_List = get_num_mb_op(N, basis_sp, num_list, basis_mb)
Number_MB_Operator_List = [number(basis_mb, site) for site in 1:N]

Himpurity = Imp_H(Number_MB_Operator_List, Impurity_Data, Vrand)

HTotal = HHubbard + Himpurity

ϵ, UPsi = eigenstates(HTotal, 100)
title_mb_imp = latexstring("Hofstadter MB Impurity Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
mb_imp_spectrum = scatter(real(ϵ), xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
mb_imp_spectrum_zoom = scatter(real(ϵ)[1:Degeneracy+5], xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
savefig(mb_imp_spectrum_zoom, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/mb_imp_spectrum_zoom.png")

UPsi_first = hcat([UPsi[i].data for i in 1:Degeneracy] ...)

Initial_Configuration = heatmap(Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Number_MB_Operator_List, basis_mb, UPsi)', title="Initial Configuration of QH",c=:roma)
savefig(Initial_Configuration, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/Initial_Configuration.png")

####################
# BRAIDING SECTION #
####################

# EXCHANGE AND DOUBLE EXCHANGE
# Simultaneously
start_site_1 = 57 # top left corner
start_site_2 = 29 # bottom right corner
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
fixed_point = 57
moving_point = 29
path_1, path_2 = braiding_path(moving_point, fixed_point, Nx, Ny, co)
#path_1 = reverse([12,13,14,15,23,31,39,38,37,36,28,20,12])
#path_1 = reverse(path_1)

qh_paths = plot_paths(co, path_1, path_2)
savefig(qh_paths, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/qh_paths.png")

V_initial = 0
step_break = 0.05
V_final = 1
#-step_break
STEP = V_initial:step_break:V_final
println("Step size between two sites:",length(STEP))

PATH1 = path_1
PATH2 = path_2
ψ_first, ψ_mat, Converge, ψ_op, imp_data, ψ_mat_list, V0_step_list, band_width, gap = get_phases(UPsi_first, V, Vrand, PATH1, PATH2, STEP, HHubbard, Number_MB_Operator_List, Degeneracy)

conv_plot = plot(Converge, title="Converge Values \n Step Size=$(length(STEP))", xlabel=L"n", ylabel=L"C = abs|<ψ|\tilde{\psi}>|", legend=false)
savefig(conv_plot, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/conv_plot.png")

width_gap = plot(real(band_width./gap),title=L"δ/Δ", xlabel="Step")
savefig(width_gap, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/width_gap.png")

ParamInfo, Degeneracy, nu_eff = ParameterInfo(NPin, pn, Nx, Ny, p, q)

BerryMatrix = ψ_mat'*ψ_first
BerryE, BerryU = eigen(BerryMatrix)

println("\n BerryE=", BerryE)

println("\n |BerryE|=", abs.(BerryE))

θ_berry_1 = angle.(BerryE)./pi

println("\n Phases=", θ_berry_1)

println("\n Avarage Phase (angle) =", sum(mod.(θ_berry_1,2))/Degeneracy)

θ_berry_2 = real(log.(BerryE)./(im*pi))

println("\n Avarage Phase=", sum(mod.(θ_berry_2,2))/Degeneracy)

N_enc = lens[1]*lens[2]
θ_AB = 2 * ϕ * nu_eff * N_enc

println("\n AB Phase=", θ_AB)

println("\n End of Hofstadter Braiding \n")

# error("STOP BEFORE MOVIE")

# MOVIE
#= @gif for i in 1:800
    heatmap(Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Number_MB_Operator_List, basis_mb, ψ_op[i])', title="Initial Configuration of QH",c=:roma)
end =#