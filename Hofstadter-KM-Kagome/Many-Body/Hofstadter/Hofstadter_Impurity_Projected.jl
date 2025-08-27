cd(@__DIR__)

using Pkg
println("\n",Pkg.status(), VERSION)

using QuantumOptics, SparseArrays, Plots, LinearAlgebra, ProgressMeter, Revise, OffsetArrays, LaTeXStrings, BenchmarkTools, DelimitedFiles, Serialization
includet("../Scripts/FirstBandApproximation.jl")
includet("../Scripts/ManyBody.jl")
includet("Hofstadter_SP.jl")
includet("../Scripts/Impurity.jl")
includet("../Scripts/Braiding.jl")

Nx = 6;
Ny = 6;
p = 1;
q = 6;
ϕ = p/q;
pn = 2;
U = 1;
V = 0.5;
Vrand = 1e-5;

N_Site = Nx*Ny;
N = N_Site;
NPhi0 = Int(Nx*Ny*(p/q));

dir_name = mkdir("Projected__Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0)")

# SINGLE PARTICLE
matrix = Hofstadter_SP(Nx, Ny, p / q, 0);
e1, psi1 = eigen(Matrix((matrix+matrix')/2));
title_sp = latexstring("Hofstadter Single Particle Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$");
sp_spectrum = scatter(e1, xlabel=L"n", ylabel=L"E", label=false, title=title_sp)
savefig(sp_spectrum, joinpath(dir_name, "SingleParticleSpectrum.png"))

HardCore = false
Cut_Off = NPhi0
HHubbard_proj, P, Pt, basis_cut_mb = H_Hubbard_Projection(N, pn, U, matrix, Cut_Off, HardCore);
num_sub_list = get_num_sub_list(N_Site, P, Pt)
basis_cut_sp = NLevelBasis(Cut_Off)
Sub_Number_MB_Operator_List = get_num_mb_op(N_Site, basis_cut_sp, num_sub_list, basis_cut_mb)

serialize(joinpath(dir_name, "DataNumber.data"), Sub_Number_MB_Operator_List) 

co, lat_plot = plot_square_lattice(N_Site, Nx, Ny);
V0 = [V, V];
Imp_Site = [find_middle_site_x(Nx, co), find_middle_site_y(Ny, co)]; # [T2_position, T1_position]
Impurity_Data = Impurity(V0, Imp_Site);
NPin = 2 # It just defines # of degeneracy

Degeneracy, nu_eff = ground_degeneracy(Nx, Ny, p, q, NPin, pn);

Impurity_H_proj = Imp_H(Sub_Number_MB_Operator_List, Impurity_Data, Vrand);
HTotal_proj = HHubbard_proj + Impurity_H_proj;
ϵ_proj, UPsi_proj = eigenstates(HTotal_proj);

serialize(joinpath(dir_name, "DataHubbard.data"), HTotal_proj) 

title_mb_imp = latexstring("Hofstadter MB Impurity Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
mb_imp_spectrum_zoom = scatter(real(ϵ_proj)[1:Degeneracy+5], xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
savefig(mb_imp_spectrum_zoom, joinpath(dir_name, "MbImpuritySpectrumZoomed.png"))

UPsi_first_proj = hcat([UPsi_proj[i].data for i in 1:Degeneracy] ...)
serialize(joinpath(dir_name, "DataInitialState.data"), UPsi_first_proj) 

Initial_Configuration_proj = heatmap(Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Sub_Number_MB_Operator_List, basis_cut_mb, UPsi_proj)', title="Initial Configuration of QH",c=:roma)
savefig(Initial_Configuration_proj, joinpath(dir_name, "DataInitialConfiguration.png"))