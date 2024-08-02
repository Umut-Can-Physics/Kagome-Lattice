using NBInclude
using Plots
using LaTeXStrings
using SparseArrays
using DataFrames
using Interpolations
using LinearAlgebra
using ProgressMeter
include("Hofstadter_Single_Particle_in_Julia.jl")
include("First_Band_Approximation_Functions.jl")
include("Hofstadter_Many_Body_Utils.jl")
include("Braiding_Utils.jl")
include("Torus_Distance.jl")
include("Braiding_Phase_Depletion_Script.jl")
include("Hofstadter_Many_Body_Utils.jl")

# 2 QH

Nx=23
Ny=23
p=1
q=23
par_num = 7
U = 2
V_Strength = 4

N=Nx*Ny
NPhi = Int( Nx * Ny * (p/q) )
cut_off = NPhi
PN = [i for i in 1:par_num]

plot_square_lattice(Nx, Ny)

x = 1
ref_site = 265
V0 = [V_Strength,x*V_Strength,x*V_Strength,x*V_Strength,x*V_Strength]
Imp_Site = [ref_site,ref_site-1,ref_site+1,ref_site+Nx,ref_site-Nx]
N_Pin = 2

ϵ, Degeneracy, avg_density = impurity_control(V_Strength, V0, Imp_Site, N_Pin)

savefig(scatter(ϵ[1:5]),"First_4_energies_Npin_$(N_Pin).png")

savefig(scatter(ϵ[1:Degeneracy+1]),"Zoom_In_Energies_Npin_$(N_Pin).png")

ref_par_density = get_ref_prtcl_density(par_num,p,q,NPhi,N_Pin)
println("Summation (2qh):",sum(ref_par_density .- avg_density)) 

filling_frac = par_num / (NPhi-N_Pin)
frac_charge = (filling_frac)*N_Pin
println("Charge(2qh):",frac_charge)

Coords = get_coords_square(Nx, Ny);

dens_2 = avg_density'[1:end-1,:]

# 1 QH

Ny=Ny-1
N=Nx*Ny
NPhi = Int( Nx * Ny * (p/q) )
cut_off = NPhi

x = 0
V0 = [V_Strength,x*V_Strength,x*V_Strength,x*V_Strength,x*V_Strength]
N_Pin = 1

ϵ, Degeneracy, avg_density = impurity_control(V_Strength, V0, Imp_Site, N_Pin)

savefig(scatter(ϵ[1:5]),"First_4_energies_Npin_$(N_Pin).png")

savefig(scatter(ϵ[1:Degeneracy+1]),"Zoom_In_Energies_Npin_$(N_Pin).png")

ref_par_density = get_ref_prtcl_density(par_num,p,q,NPhi,N_Pin)
println("Summation(1qh):",sum(ref_par_density .- avg_density)) 

filling_frac = par_num / (NPhi-N_Pin)
frac_charge = (filling_frac)*N_Pin
println("Charge(1qh):",frac_charge)

Coords = get_coords_square(Nx, Ny);

distance_array = []
for j in 1:N
    push!(distance_array,distance_func(Coords, Coords[ref_site], Coords[j]))
end
distance_array = reshape(distance_array,Nx,Ny)'

dens_1 = avg_density'

a = 1 # lattica constant
α = p/q # flux per plaquette
l_b = a/sqrt(2*pi*α)
braid_phase = []
braiding_density = ((1/2)/l_b^2)*(dens_2 .- 2*dens_1).*(distance_array.^2)
savefig(heatmap(dens_1),"dens_1_heatmap.png")
savefig(heatmap(dens_2),"dens_2_heatmap.png")
savefig(heatmap(braiding_density),"braiding_density_heatmap.png")
R_max = sort(unique(distance_array))
println(R_max)
for r in R_max
    push!(braid_phase, sum(braiding_density .* (distance_array.<=r)))
end
av_braiding = sum(mod.(braid_phase,1))/length(braid_phase)
plot(R_max,mod.(braid_phase,1),xlabel=L"R_{max}[l_b]", ylabel=L"\varphi_{br}[2\pi]" ,title=L"\frac{\varphi_{br}}{2\pi} = \frac{1}{2l_b^2} \sum_j \left[d_{2qh}-2d_{1qh} \right] |\rho_j|^2", label="$(av_braiding)",marker=(:circle,5))
x=0;g(x)=frac_charge
P = plot!(g, x, 7, label="$(filling_frac)")
savefig(P,"Braiding_Phase.png")






