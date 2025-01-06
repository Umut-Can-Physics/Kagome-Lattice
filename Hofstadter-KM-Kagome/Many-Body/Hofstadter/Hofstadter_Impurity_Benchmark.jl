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
Cut_Off = NPhi0;

co, lat_plot = plot_square_lattice(N_Site, Nx, Ny);
lat_plot
savefig(lat_plot,"./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/Lattice.png") 

# SP SPECTRUM
matrix = Hofstadter_SP(Nx, Ny, p / q, 0);
e1, psi1 = eigen(Matrix((matrix+matrix')/2));
title_sp = latexstring("KM Single Particle Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$");
sp_spectrum = scatter(e1, xlabel=L"n", ylabel=L"E", label=false, title=title_sp);
savefig(sp_spectrum,"./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/sp_spectrum.png") 

# HARD-CORE IMPURTIY HAMILTONIAN
HardCore = true; 
println("\n Constructing MB Operator...")
basis_sp = NLevelBasis(N);
H1_op = Sp_Op(basis_sp, matrix);
basis_mb = get_Bosonic_MB_Basis(basis_sp, pn, HardCore);
HHubbard = get_mb_hopping(basis_mb, H1_op);
println("\n End of the MB Operator")
Number_MB_Operator_List = [number(basis_mb, site) for site in 1:N];

V0 = [V, V];
Imp_Site = [find_middle_site_x(Nx, co), find_middle_site_y(Ny, co)]; # [T2_position, T1_position]
Impurity_Data = Impurity(V0, Imp_Site);
NPin = 2 # It just defines # of degeneracy

Degeneracy, nu_eff = ground_degeneracy(Nx, Ny, p, q, NPin, pn);
ParameterInfo(NPin, pn, Nx, Ny, p, q);
println("\n Vrand=",Vrand);
println("\n V=", V);

Himpurity = Imp_H(Number_MB_Operator_List, Impurity_Data, Vrand);

HTotal = HHubbard + Himpurity;

serialize("Hofstadter-KM-Kagome/Many-Body/Hofstadter/DataHamiltonian_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data", HTotal) 
#HTotal = deserialize("Hofstadter-KM-Kagome/Many-Body/Hofstadter/Hamiltonian_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data")

ϵ, UPsi = eigenstates(HTotal, 100, info=false);

title_mb_imp = latexstring("Hofstadter MB Impurity Spectrum \n \$ N_x=$(Nx), N_y=$(Ny), \\phi=$(p//q), N=$(pn), U=$(U), V=$(V) \$")
mb_imp_spectrum = scatter(real(ϵ), xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
mb_imp_spectrum_zoom = scatter(real(ϵ)[1:Degeneracy+5], xlabel=L"n", ylabel=L"E", label="Hard-Core=$(HardCore)", title=title_mb_imp)
savefig(mb_imp_spectrum_zoom, "./Hofstadter-KM-Kagome/Many-Body/Hofstadter/Braiding_Data/mb_imp_spectrum_zoom.png")

UPsi_first = hcat([UPsi[i].data for i in 1:Degeneracy] ...)
serialize("Hofstadter-KM-Kagome/Many-Body/Hofstadter/DataInitialState_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data", UPsi_first) 
#UPsi_first = deserialize("Hofstadter-KM-Kagome/Many-Body/Hofstadter/InitialState_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data")

#= open("Hofstadter-KM-Kagome/Many-Body/Hofstadter/geek.txt", "w") do file
    writedlm(file, [ParameterInfo(NPin, pn, Nx, Ny, p, q)])
end =#