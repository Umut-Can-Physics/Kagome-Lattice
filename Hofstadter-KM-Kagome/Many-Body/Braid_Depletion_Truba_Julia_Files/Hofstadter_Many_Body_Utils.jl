#   Module for Many Body Calculations with FBA
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function H_sp(N, Nx, Ny, p, q)
    sp_basis = NLevelBasis(N)
    periodicity = 0 #periodic (select 1 for hard-wall conditions)
    sp_matrix = Hofstadter_SP(Nx, Ny, p/q, periodicity)
    H1 = get_sp_op(sp_basis, sp_matrix)
    return H1, sp_basis
end

function H_sub(N, Nx, Ny, p, q, H1, cut_off, sp_basis)
    sub_states = get_sub_states(H1, cut_off)
    basis_sub, P, Pt = get_projector_op(sub_states, sp_basis)
    H1_sub = get_subspace_op(H1, P, Pt)
    return H1_sub, basis_sub, P, Pt
end

function H_Kin_MB(basis_sub, PN, H1_sub)
    states_mb = bosonstates(basis_sub, PN) 
    basis_mb = ManyBodyBasis(basis_sub, states_mb)
    H1_MB = get_mb_op(basis_mb, H1_sub)
    return H1_MB
end

function H_Total_Sub(basis_cut_mb, basis_sub, PN, H1_sub,P, Pt, cut_off, U)
    H_Kin = SparseOperator(basis_cut_mb)
    H1_MB = H_Kin_MB(basis_sub, PN, H1_sub)
    H_Kin.data = H1_MB.data
    H_Int = Hubbard_Interaction_op(P, Pt, basis_cut_mb, cut_off, U)
    return H_Kin, H_Kin + H_Int
end

function Number_sub(N, sp_basis, P, Pt, basis_cut_sp, basis_cut_mb, basis_sub)
    num_sub_list = get_num_sub_list(N, sp_basis, P, Pt)
    Sub_Number_MB_Operator_List = get_num_mb_op(N, basis_cut_sp, num_sub_list, basis_cut_mb, basis_sub)
    return Sub_Number_MB_Operator_List
end

function get_H_Impurity(N, sp_basis, P, Pt, basis_cut_sp, basis_cut_mb, basis_sub, Total_H, Sub_Number_MB_Operator_List, Impurity_Data)
    H_Impurity = Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data)
    H_Impurity = dense((H_Impurity'+H_Impurity)/2);
    return H_Impurity
end

function ground_degeneracy(Nx, Ny, p, q, N_Pin, PN)
    NPhi0 = Nx*Ny*(p/q)
    NPhi = NPhi0-N_Pin
    pn = maximum(PN)
    nu0 = 1/2
    N_d = Int(NPhi - pn/nu0)
    if length(PN) == 1
        Degeneracy = 1
    else
        Degeneracy = Int((factorial(N_d + pn - 1) / (factorial(N_d) * factorial(pn - 1))) * (NPhi / pn))
    end
    return Degeneracy, pn, NPhi0, N_d
end

function energies_imp(H_Impurity, PN, basis_cut_mb)
    E, V = eigenstates(H_Impurity)
    pn = maximum(PN)
    filtered_energies = get_filtered_energies(pn, E, V, basis_cut_mb)
    return filtered_energies, V
end

function plot_density(Nx, Ny, Degeneracy, N, Sub_Number_MB_Operator_List, basis_cut_mb, r_hubbard_states, factor)
    avg_density = Get_Avg_Density(Nx, Ny, Degeneracy, N, Sub_Number_MB_Operator_List, basis_cut_mb, r_hubbard_states)
    return Plots.heatmap(Interp(avg_density', factor), aspect_ratio=:equal)
end

function get_phases(Impurity_Data,rec_path_1,rec_path_2, Imp_Site, Total_H, Sub_Number_MB_Operator_List, Degeneracy)

    # Initial Configuration
    N_Pin = 4
    V1 = Impurity_Data.V0[1]
    V2 = Impurity_Data.V0[2]
    Imp_Site = [rec_path_1[1], rec_path_1[2], rec_path_2[1], rec_path_2[2]]
    V0 = [V1, 0, V2, 0]

    Impurity_Data = Impurity(V0, Imp_Site)
    Impurity_H = Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data)
    Impurity_H = dense((Impurity_H+Impurity_H')/2)
    E_Imp_0, ψ = eigenstates(Impurity_H, Degeneracy)
    
    ψ = hcat([ψ[i].data for i in 1:Degeneracy] ...)
    ψ_first = copy(ψ)

    Imp_Site_List = [ [imp, rec_path_1[idx+1], rec_path_2[idx], rec_path_2[idx+1] ] for (idx,imp) in (enumerate(rec_path_1[1:end-1])) ]
    V0_List = [ [V1*(1-step), V1*step, V2*(1-step), V2*step] for step in STEP ]
    
    Impurity_Data_List = [ [Impurity(V00, Imp_Sitee)] for Imp_Sitee in Imp_Site_List for V00 in V0_List ]
    Impurity_H_List = [ Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Dataa[1]) for Impurity_Dataa in Impurity_Data_List]
    
    Ψ_list = []
    Ψ_list_2 = []
    
    @showprogress for Impurity_HH in Impurity_H_List
        
        H = dense((Impurity_HH + Impurity_HH')/2)
        ϵ, ψ_tilde = eigenstates(H, Degeneracy) 

        push!(Ψ_list, ψ_tilde)
        
        ψ_tilde = hcat([ψ_tilde[i].data for i in 1:Degeneracy] ...)

        # KM Algorithm #
        # A = ψ'*ψ_tilde
        # A_inv = inv(A)
        # ψ = ψ_tilde*A_inv 
        # !!!
        
        # !!! Vanderbilt
        A = ψ' * ψ_tilde
        V, Σ, W = svd(A)
        M = V * W'
        ψ = ψ_tilde * M'
        # !!!

        # Method 1
        # for i in 1:Degeneracy
        #     norm = sqrt(ψ[:,i]'*ψ[:,i])
        #     ψ[:,i] = ψ[:,i] ./ norm
        # end

        # Method 2
        # ψ = qr(ψ).Q * Matrix(I, size(ψ)...) #Vanderbilt'de olacak mı?

        # x EVOLUTION OF GENERATED FUNCTION x #
        push!(Ψ_list_2, ψ)
        
    end

    BerryEnergies, BerryStates = eigen(ψ' * ψ_first)

    ϕ_tot = -imag(log(det(ψ' * ψ_first)))
    
    return ψ, ψ_first, Ψ_list, ϕ_tot, BerryEnergies
end

# using NBInclude
# @nbinclude("Hofstadter/Hofstadter MB in Julia.ipynb"; regex=r"#.*executeme");
include("Hofstadter_MB_in_Julia.jl")

function Get_MB(Nx, Ny, p, q, cut_off, PN, U, Impurity_Data, factor, N_Pin)
    N = Nx*Ny
    H1, sp_basis = H_sp(N, Nx, Ny, p, q)
    H1_sub, basis_sub, P, Pt = H_sub(N, Nx, Ny, p, q, H1, cut_off, sp_basis)
    basis_cut_mb, basis_cut_sp = get_Bosonic_MB_Basis(cut_off, PN)
    H_Kin, Total_H = H_Total_Sub(basis_cut_mb, basis_sub, PN, H1_sub,P, Pt, cut_off, U)
    Sub_Number_MB_Operator_List = Number_sub(N, sp_basis, P, Pt, basis_cut_sp, basis_cut_mb, basis_sub)
    H_Impurity = get_H_Impurity(N, sp_basis, P, Pt, basis_cut_sp, basis_cut_mb, basis_sub, Total_H, Sub_Number_MB_Operator_List, Impurity_Data)
    E, V = energies_imp(H_Impurity, PN, basis_cut_mb)
    filtered_energies = E
    r_hubbard_states = Restricted_Hubbard_States(V, filtered_energies)
    Degeneracy, pn, NPhi0 = ground_degeneracy(Nx, Ny, p, q, N_Pin, PN)
    Plot_1 = scatter(E,legend=false,title="Degeneracy=$(Degeneracy)")
    savefig(Plot_1, "Energies_Npin_$(N_Pin).png")
    display(Plot_1)
    Plot_2 = plot_density(Nx, Ny, Degeneracy, N, Sub_Number_MB_Operator_List, basis_cut_mb, r_hubbard_states, factor)
    savefig(Plot_2, "Density_Npin_$(N_Pin).png")
    display(Plot_2)
    return E, Sub_Number_MB_Operator_List, basis_cut_mb, Degeneracy, Total_H, Sub_Number_MB_Operator_List, r_hubbard_states, pn, NPhi0, H1
end

function movie(Nx, Ny, Degeneracy, N, Sub_Number_MB_Operator_List, basis_cut_mb, Eigen_List, factor)
    @gif for i in 1:length(Eigen_List)
        data1 = Get_Avg_Density(Nx, Ny, Degeneracy, N, Sub_Number_MB_Operator_List, basis_cut_mb, Eigen_List[i])'
        heatmap!(Interp(data1, factor), aspect_ratio=:equal)
    end
    return nothing
end

function get_braiding_path(Imp_Site, Nx, Ny, co)
    start_point = Imp_Site[1] # Moving site of quasihole
    First_Path = Int64[]
    mod_list = [mod(start_point,Nx),mod(start_point,Ny)]
    bottom_site_in_the_Ny_direction = mod_list[mod_list .> 0][1]
    last_site_in_the_Ny_direction = (bottom_site_in_the_Ny_direction+(Nx*(Ny-1)))
    site_number_from_upward = (last_site_in_the_Ny_direction-start_point)/Nx
    site_number_from_downward = (start_point-bottom_site_in_the_Ny_direction)/Nx
    for i in 0:site_number_from_upward
        step_site = start_point+i*Nx
        push!(First_Path, step_site)
        if step_site == last_site_in_the_Ny_direction
            for j in 0:site_number_from_downward
                push!(First_Path, bottom_site_in_the_Ny_direction+j*Nx)
            end
        end
    end
    A = findall(x->x==co[:,2][start_point],co[:,2])
    B = findall(x->x==start_point, A)[1]
    Second_Path = push!(vcat(reverse(A[1:B]),reverse(A[B+1:end])),start_point)
    Third_Path = reverse(First_Path)
    Fourth_Path = reverse(Second_Path)
    rec_path_braiding = vcat(First_Path,Second_Path,Third_Path,Fourth_Path)
    return rec_path_braiding
end

function th_AB_phase(pn, p, q, N_Pin, N_mov, number_of_plaq)
    NPhi = Int( Nx * Ny * (p/q) )
    charge = pn/(NPhi-N_Pin)
    θ_AB = N_mov * (p/q) * charge * number_of_plaq
    ex = exp(2*im*pi*θ_AB)
    return θ_AB, ex, charge
end

function chose_path(Path_Type, rec_path_exch_1, rec_path_exch_2, Imp_Site, Nx, Ny, co)
    # N_p = 1 or Closed Path
    if Path_Type == 1
        rec_path_1 = rec_path_exch_1
        rec_path_2 = rec_path_1        
    # x2 Exchange 
    elseif Path_Type == 2
        rec_path_1 = rec_path_exch_1
        rec_path_2 = rec_path_exch_2
    # Braiding
    rec_path_braiding = get_braiding_path(Imp_Site, Nx, Ny, co)
    elseif Path_Type == 3
        rec_path_1 = rec_path_braiding
        rec_path_2 = repeat([Imp_Site[2]],length(rec_path_1))
    # Manuel
    elseif Path_Type == 4
        rec_path_1 = [29,30,31,42,53,52,51,40,29]
        rec_path_2 = repeat([Imp_Site[2]],length(rec_path_1))
    end
    return rec_path_1, rec_path_2
end

function plot_paths(co, First_Path, Second_Path)
    M = Matrix{Int64}(undef,length(First_Path),2)
    N = Matrix{Int64}(undef,length(Second_Path),2)
    for (idx,value) in enumerate(First_Path)
            M[idx,:] = co[value,:]
    end
    for (idx,value) in enumerate(Second_Path)
            N[idx,:] = co[value,:]
    end
    p1 = scatter(co[:,1],co[:,2], series_annotations = text.([i for i in 1:Nx*Ny], :bottom), legend=false, aspect_ratio = :equal)
    p2 = plot!(M[:, 1], M[:, 2], linewidth=2, c=:blue,legend=:false)
    p3 = plot!(N[:, 1], N[:, 2], linewidth=2, c=:red,legend=:false)
    return display(p3)
end