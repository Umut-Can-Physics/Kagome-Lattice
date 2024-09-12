struct Impurity
    V0::Vector{Float64}
    Imp_Site::Vector{Int64}
end

function Imp_H(H, Sub_Number_MB_Operator_List, Impurity_Data)
    for imp in 1:length(Impurity_Data.V0)
        H += Impurity_Data.V0[imp] * Sub_Number_MB_Operator_List[Impurity_Data.Imp_Site[imp]]
    end
    return dense((H'+H)/2)
end

function ground_degeneracy(Nx, Ny, p, q, N_Pin, pn)
    NPhi0 = Nx*Ny*(p/q)
    NPhi = NPhi0-N_Pin
    nu0 = 1/2
    N_d = Int(NPhi - pn/nu0)
    if length(PN) == 1
        Degeneracy = 1
    else
        Degeneracy = Int((factorial(N_d + pn - 1) / (factorial(N_d) * factorial(pn - 1))) * (NPhi / pn))
    end
    return Degeneracy, pn, NPhi0, N_d
end

function ParameterInfo(NPin, pn, Nx, Ny, p, q)
    nu0 = 1/2 # The quasihole numbers depend on Laughlin fraction
    NPhi = NPhi0-NPin
    N_d = Int(NPhi - pn/nu0)
    nu = pn//NPhi0
    Degeneracy = Int((factorial(N_d + pn - 1) / (factorial(N_d) * factorial(pn - 1))) * (NPhi / pn))
    println("Lattice: ",Nx,"x",Ny,
        "\nThe Number of Flux per Unit-Cell (ϕ) = ",p//q,
        "\nFilling Fraction (ν) = ",nu," (The Number of Quasiholes = ",(NPin+N_d),
        ")\nThe Number of Flux Quanta = ", NPhi,
        " (Delocalised Number = ", N_d,
        ")\nThe Total Number of Particle = ",pn,
        "\nThe Number of State of Ground Degeneracy = ",Degeneracy)
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