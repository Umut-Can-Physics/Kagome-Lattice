struct Impurity
    V0::Vector{Float64}
    Imp_Site::Vector{Int64}
end

function Imp_H(Number_MB, Impurity_Data, Vrand)
    # Random Impurtiy: To unsure that finding all degeneracy by studying sparse matrices
    Vimp = Impurity_Data.V0[1] * Number_MB[Impurity_Data.Imp_Site[1]] 
    for imp in 2:length(Impurity_Data.V0)
        Vimp += Impurity_Data.V0[imp] * Number_MB[Impurity_Data.Imp_Site[imp]] 
    end
    for Nop in Number_MB
        Vimp += Vrand * rand() * Nop
    end
    return (Vimp'+Vimp)/2
end

function ground_degeneracy(Nx, Ny, p, q, N_Pin, pn)
    NPhi0 = Nx*Ny*(p/q)
    NPhi = NPhi0-N_Pin
    nu0 = 1/2
    N_d = Int(NPhi - pn/nu0)
    nu_eff = pn/NPhi
    #if length(pn) == 1
    #    Degeneracy = 1
    #else
    Degeneracy = Int((factorial(N_d + pn - 1) / (factorial(N_d) * factorial(pn - 1))) * (NPhi / pn))
    #end
    return Degeneracy, nu_eff
end

function ParameterInfo(NPin, pn, Nx, Ny, p, q)
    nu0 = 1/2 # The quasihole numbers depend on Laughlin fraction
    NPhi = NPhi0-NPin
    N_d = Int(NPhi - pn/nu0)
    nu = pn//NPhi0
    nu_eff = pn//NPhi
    Degeneracy = Int((factorial(N_d + pn - 1) / (factorial(N_d) * factorial(pn - 1))) * (NPhi / pn))
    return println("\n Lattice Size: ",Nx,"x",Ny,
        "\n The Number of Flux per Unit-Cell (ϕ) = ",p//q,
        "\n Filling Fraction (ν) = ",nu," (The Number of Quasiholes = ",(NPin+N_d),
        ")\n The Number of Flux Quanta = ", NPhi,
        " (Delocalised Number = ", N_d,
        ")\nThe Total Number of Particle = ",pn,
        "\n The Effective Filling = ",nu_eff,
        "\nThe Number of Ground State Degeneracy = ",Degeneracy)
end

function Get_Density_Profile(N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States, index)
    Expectation_List = []
    for site in 1:N_Site
        push!(Expectation_List, expect(Sub_Number_MB_Operator_List[site], Fil_States[index]))
    end
    return real(Expectation_List)
end

function Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States)
    Avg_Density = spzeros(Nx,Ny)
    for index in 1:Degeneracy
        Avg_Density += reshape(Get_Density_Profile(N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States, index), Nx, Ny)
    end    
    return Avg_Density / Degeneracy
end

function plot_density(Nx, Ny, Degeneracy, N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States, factor)
    avg_density = Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States)
    #return Plots.heatmap(Interp(avg_density', factor), aspect_ratio=:equal)
    return Plots.heatmap(avg_density', aspect_ratio=:equal)
end

function movie(Nx, Ny, Degeneracy, N, Sub_Number_MB_Operator_List, basis_cut_mb, Eigen_List, factor)
    @gif for i in 1:length(Eigen_List)
        data1 = Get_Avg_Density(Nx, Ny, Degeneracy, N, Sub_Number_MB_Operator_List, basis_cut_mb, Eigen_List[i])'
        heatmap(Interp(data1, factor), aspect_ratio=:equal)
    end
    return nothing
end