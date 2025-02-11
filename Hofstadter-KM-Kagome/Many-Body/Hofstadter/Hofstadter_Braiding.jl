cd(@__DIR__)

using Pkg
println("\n",Pkg.status(), VERSION)

using QuantumOptics, SparseArrays, Plots, LinearAlgebra, ProgressMeter, Revise, OffsetArrays, LaTeXStrings, BenchmarkTools, Serialization, DelimitedFiles
includet("../Scripts/FirstBandApproximation.jl")
includet("../Scripts/ManyBody.jl")
includet("Hofstadter_SP.jl")
includet("../Scripts/Impurity.jl")
includet("../Scripts/Braiding.jl")

Nx = 10;
Ny = 10;
p = 1;
q = 10;
ϕ = p/q;
pn = 4;
U = 1;
V = 0.5;
Vrand = 1e-5;

N_Site = Nx*Ny;
N = N_Site;
NPhi0 = Int(Nx*Ny*(p/q));

co, lat_plot = plot_square_lattice(N_Site, Nx, Ny);

V0 = [V, V];
Imp_Site = [find_middle_site_x(Nx, co), find_middle_site_y(Ny, co)]; # [T2_position, T1_position]
Impurity_Data = Impurity(V0, Imp_Site);
NPin = 2 # It just defines # of degeneracy

Degeneracy, nu_eff = ground_degeneracy(Nx, Ny, p, q, NPin, pn);

T2_path = reduce(vcat,[Imp_Site[1]:Imp_Site[1]+Nx-1,Imp_Site[1], repeat([Imp_Site[1]], Ny), Imp_Site[1]+Nx-1:-1:Imp_Site[1], repeat([Imp_Site[1]], Ny)])
T1_path = reduce(vcat,[repeat([Imp_Site[2]], Nx+1),   push!(collect(Imp_Site[2]-Ny:-Ny:Imp_Site[2]-(Ny-1)*Ny),Imp_Site[2]), repeat([Imp_Site[2]], Nx),  Imp_Site[2]-(Ny-1)*Ny:Ny:Imp_Site[2]])

qh_paths = plot_paths(co, T2_path, T1_path)
savefig(qh_paths, "./Braiding_Data/qh_paths.png")

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
#= Number_MB_Operator_List = deserialize("DataNumber_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data") 
HHubbard = deserialize("DataHubbard_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data") 
UPsi_first = deserialize("DataInitialState_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data")
 =#
Sub_Number_MB_Operator_List = deserialize("DataNumberProjected_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data") 
HHubbard_proj = deserialize("DataHubbardProjected_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data") 
UPsi_first_proj = deserialize("DataInitialStateProjected_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data")
ψ_first, ψ_mat, Converge, ψ_op, imp_data, ψ_mat_list, V0_step_list, band_width, gap = get_phases(UPsi_first_proj, V1, V2, Vrand, PATH1, PATH2, STEP, HHubbard_proj, Sub_Number_MB_Operator_List, Degeneracy)

conv_plot_proj = plot(Converge, title="Converge Values \n Step Size=$(length(STEP))", xlabel=L"n", ylabel=L"C = abs|<ψ|\tilde{\psi}>|", legend=false)
ConvergeData_proj = serialize("DataConvergeProjected_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0)_1.data", conv_plot_proj)
savefig(conv_plot_proj, "./Braiding_Data/convProjected.png")

width_gap_proj = plot(real(band_width./gap),title=L"δ/Δ", xlabel="Step")
width_Data_proj = serialize("DataBandWidth_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0)_1.data", real(band_width./gap))
savefig(width_gap_proj, "./Braiding_Data/width_gapProjected.png")

BerryMatrix = ψ_mat'*ψ_first
BerryE, BerryU = eigen(BerryMatrix)

serialize("DataBerryE_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data", BerryE) 
serialize("DataBerryU_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).data", BerryU) 

open("DataBRAIDINGProjected_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).txt", "w") do file
    writedlm(file, "")
end

open("DataBRAIDINGProjected_Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0).txt", "a") do file

    writedlm(file, ["BRAIDING RESULTS FOR PROJECTED SPACE \n Cut-off = $(Cut_Off)"])

    # Write parameter info
    writedlm(file, [string(ParameterInfo(NPin, pn, Nx, Ny, p, q))])

    writedlm(file, ["\n U = $(U)"])

    writedlm(file, ["\n V = $(Vrand)"])

    writedlm(file, ["\n V_Rand = $(Vrand)"])

    # Write BerryE value
    writedlm(file, ["\n BerryE=$(BerryE)"])

    # Add requested lines
    writedlm(file, ["\n |BerryE|=$(abs.(BerryE))"])
    
    θ_berry = angle.(BerryE) ./ pi
    writedlm(file, ["\n Phases=$(θ_berry)"])
    
    writedlm(file, ["\n Average Phase (angle)=$(sum(mod.(θ_berry, 2)) / Degeneracy)"])
    
    writedlm(file, ["\n nu_eff=$(float(nu_eff))"])
    
    writedlm(file, ["\n 2*nu_eff=$(2 * float(nu_eff))"])
end

println("\n End of Hofstadter Braiding \n")