function Impurity_Control(N, Nx, Ny, pn, t1, L1, t2, L2, cut_off, U, Impurity_Data)
    
    # SP LEVEL
    sp_basis = NLevelBasis(N)
    H1 = get_SP_H(Nx, Ny, t1, L1, t2, L2)
    sp_spectrum = eigenenergies(dense(H1))
    
    # PROJECTOR LEVEL
    sub_states = get_sub_states(H1, cut_off)
    basis_sub, P, Pt = get_projector_op(sub_states, sp_basis)
    H1_sub = get_subspace_op(H1, P, Pt)
    num_sub_list = get_num_sub_list(N, sp_basis, P, Pt)

    # MB LEVEL
    states_mb = bosonstates(basis_sub, PN) 
    basis_mb = ManyBodyBasis(basis_sub, states_mb)
    H1_MB = get_mb_op(basis_mb, H1_sub)
    basis_cut_mb, basis_cut_sp = get_Bosonic_MB_Basis(cut_off,PN)
    H_Int = Hubbard_Interaction(P, Pt, basis_cut_mb, cut_off, U) 
    H1cut = SparseOperator(basis_cut_mb)
    H1cut.data = H1_MB.data
    Total_H = H1cut + H_Int
    Total_H = (Total_H'+Total_H)/2
    mb_spectrum = eigenenergies(dense(Total_H))

    # IMPURITY LEVEL
    number_mb_list_operators = get_num_mb_op(N, basis_cut_sp, num_sub_list, basis_cut_mb, basis_sub)
    Impurity_H = Imp_H(Total_H, number_mb_list_operators, Impurity_Data)
    Impurity_H = dense((Impurity_H+Impurity_H')/2)
    E_full, V_full = eigenstates(dense(Impurity_H))
    filtered_energies = get_filtered_energies(pn, E_full, V_full, basis_cut_mb) 
    r_hubbard_states = Restricted_Hubbard_States(V_full, filtered_energies)
    
    return sp_spectrum, mb_spectrum, filtered_energies, basis_cut_mb, number_mb_list_operators, filtered_energies, r_hubbard_states
end

function ground_state_deg(Nx, Ny, p, q, N_Pin, PN)
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

function Avg_Density(Nx, Ny, N, p, q, N_Pin, PN, number_mb_list_operators, basis_cut_mb, r_hubbard_states)
    Degeneracy, pn, NPhi0, N_d = ground_state_deg(Nx, Ny, p, q, N_Pin, PN)
    av_density = Get_Avg_Density_List(Nx, Ny, Degeneracy, N, number_mb_list_operators, basis_cut_mb, r_hubbard_states)
    return Degeneracy, av_density
end

function Plot_Density_Kagome(Nx, Ny, a1_vec, a2_vec, Basis)
    N = Nx*Ny*3
    x_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[4], 1:N)
    y_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[5], 1:N)
    z = exp_list1.(x_co, y_co)
    z = collect(Iterators.flatten(z))
    density = exp_list1.(x_co, y_co)
    density = collect(Iterators.flatten(z))
    surface(x_co,y_co,z,xlabel=L"$x$",ylabel=L"$y$", camera = (0,90), size=(800,800), aspect_ratio=:equal)
    # Plot kagome as project to the plot of density profile
    return scatter!(x_co, y_co, 0. *density, camera=(0,90), legend=false)
end

function annotate_densities(Nx, Ny, a1_vec, a2_vec, Basis)
    N = Nx*Ny*3
    x_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[4], 1:N)
    y_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[5], 1:N)
    z = exp_list1.(x_co, y_co)
    z = collect(Iterators.flatten(z))
    density = exp_list1.(x_co, y_co)
    density = round.(collect(Iterators.flatten(z)), digits=4)
    annotate_size = 10
    return scatter(x_co, y_co, series_annotations = text.(density, :top, annotate_size))
end

function calculate_frac_charge_kagome(radius_list, ref_par_density, avg_density)
    density_list = []
    density_list_2 = []
    summ_list = []
    
    for i in radius_list
        
        summ = 0
    
        for j in Inner_Sites_Kagome(Nx, Ny, a1_vec, a2_vec, Basis, i)
            
            # r içindeki her bir sitenin ortalama yoğunluğu 
            push!(density_list, collect(Iterators.flatten((ref_par_density .- avg_density|>transpose)))[j])
            # Örgüde sadece bir tane quasi-parçacık boştayken gelen neredeyse ortalama yoğunluk 3/70'tir!
            
            # aynı yoğunluklar gelirse sadece onların birisini tut
            density_list_2 = unique!(density_list)
            
            # her biri biricik ve her bir site için olan yoğunlukları topla
            summ = sum(density_list_2)
        end
        
        # Sitelerin toplam yoğunluklarını bir listede her bir yarı-çap için biriktir
        push!(summ_list, summ)
        # println("Rho:",round(i,digits=3),"\t","Q_rho:",summ)
    end
    return summ_list
end

function plot_density_charge_dep_prof(radius_list, Charges,frac_charge)
    Plots.plot(radius_list, Charges, xlabel=L"\rho", ylabel=L"Q_{\rho}", title="Density (Charge) Depletion Profile", guidefontsize=17,legend=false, linewidth=3, m = (5, :white, stroke(1, :blue)))
    g(x)=frac_charge*N_Pin;x=0
    return Plots.plot!(g, x, length(radius_list), line=(:dot,2), xlim=(0,7))
end
