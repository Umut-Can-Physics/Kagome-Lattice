cd(@__DIR__)

using Pkg
println("\n",Pkg.status(), VERSION)

using QuantumOptics, SparseArrays, Plots, LinearAlgebra, ProgressMeter, Revise, OffsetArrays, LaTeXStrings, BenchmarkTools, Serialization, DelimitedFiles
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

HC = "HardCore_true"
FI = "HardCore_false"
PR = "Projected"

Int_type = PR

dir_name = Int_type*"__Nx_$(Nx)__Ny_$(Ny)__PN_$(pn)__NPhi_$(NPhi0)"

co, lat_plot = plot_square_lattice(N_Site, Nx, Ny);
V0 = [V, V];
Imp_Site = [find_middle_site_x(Nx, co), find_middle_site_y(Ny, co)]; # [T2_position, T1_position]
Impurity_Data = Impurity(V0, Imp_Site);
NPin = 2 # It just defines # of degeneracy

Degeneracy, nu_eff = ground_degeneracy(Nx, Ny, p, q, NPin, pn);

T2_path = reduce(vcat,[Imp_Site[1]:Imp_Site[1]+Nx-1,Imp_Site[1],repeat([Imp_Site[1]], Ny),Imp_Site[1]+Nx-1:-1:Imp_Site[1], repeat([Imp_Site[1]], Ny)])
T1_path = reduce(vcat,[repeat([Imp_Site[2]], Nx+1),push!(collect(Imp_Site[2]-Nx:-Nx:Imp_Site[2]-(Ny-1)*Nx),Imp_Site[2]), repeat([Imp_Site[2]], Nx),  Imp_Site[2]-(Ny-1)*Nx:Nx:Imp_Site[2]])

qh_paths = plot_paths(co, T2_path, T1_path)
savefig(qh_paths, joinpath(dir_name, "Path.png"))

V_initial = 0
step_break = 0.05
V_final = 1
STEP = V_initial:step_break:V_final
println("Step size between two sites:",length(STEP))

PATH1 = T2_path
PATH2 = T1_path
V1 = V0[2] 
V2 = V0[1] 

Number_MB_Operator_List = deserialize(joinpath(dir_name,"DataNumber.data")) 
HHubbard = deserialize(joinpath(dir_name, "DataHubbard.data")) 
UPsi_first = deserialize(joinpath(dir_name,"DataInitialState.data"))
ψ_first, ψ_mat, Converge, ψ_op, imp_data, ψ_mat_list, V0_step_list, band_width, gap = get_phases(UPsi_first, V1, V2, Vrand, PATH1, PATH2, STEP, HHubbard, Number_MB_Operator_List, Degeneracy)

conv_plot = plot(Converge, title="Converge Values \n Step Size=$(length(STEP))", xlabel=L"n", ylabel=L"C = abs|<ψ|\tilde{\psi}>|", legend=false)
serialize(joinpath(dir_name, "DataConverge.data"), Converge)
savefig(conv_plot, joinpath(dir_name, "ConvergePlot.png"))

width_gap = plot(real(band_width./gap),title=L"δ/Δ", xlabel="Step")
serialize(joinpath(dir_name, "DataBandWidthGap.data"), real(band_width./gap))
savefig(width_gap, joinpath(dir_name, "BandWidthGapPlot.png"))

BerryMatrix = ψ_mat'*ψ_first
BerryE, BerryU = eigen(BerryMatrix)

serialize(joinpath(dir_name, "DataBerryEnergies.data"), BerryE) 
serialize(joinpath(dir_name, "DataBerryStates.data"), BerryU) 

open(joinpath(dir_name, "BradidingResults.txt"), "w") do file
    writedlm(file, "")
end

open(joinpath(dir_name, "BradidingResults.txt"), "a") do file

    write(file, "BRAIDING OF ABELIAN ANYONS WITH A VANISHING AHORONOV-BOHM PHASE \n")

    if Int_type == PR
        writedlm(file, ["Cut-off = $(Cut_Off) \n"]) 
    elseif Int_type == HC
        write(file, "Hard-Core=TRUE \n")
    elseif Int_type == FI
        write(file, "Hard-Core=FALSE \n")
    end

    writedlm(file, [string(ParameterInfo(NPin, pn, Nx, Ny, p, q))])

    writedlm(file, ["\n U = $(U)"])

    writedlm(file, ["\n V = $(Vrand)"])

    writedlm(file, ["\n V_Rand = $(Vrand)"])

    writedlm(file, ["\n BerryE=$(BerryE)"])

    writedlm(file, ["\n |BerryE|=$(abs.(BerryE))"])
    
    θ_berry = angle.(BerryE) ./ pi
    writedlm(file, ["\n Phases=$(θ_berry)"])
    
    writedlm(file, ["\n Average Phase (angle)=$(sum(mod.(θ_berry, 2)) / Degeneracy)"])
    
    writedlm(file, ["\n nu_eff=$(float(nu_eff))"])
    
    writedlm(file, ["\n 2*nu_eff=$(2 * float(nu_eff))"])
end

println("\n End of Hofstadter Braiding \n")