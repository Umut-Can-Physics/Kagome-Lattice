function path_line(point,len,direction)
    path = Array{Int64}(undef, 0)
    for i in 1:len+1
        append!(path,point)
        point += direction
    end
    return path
end

function rectangular_path(start_point,lens,dirs)
    paths = Array{Int64}(undef, 0)
    for (len,dir) in zip(lens,dirs)
        path = path_line(start_point,len,dir)
        append!(paths, path[1:end-1])
        start_point = last(path)
    end
    append!(paths, start_point)
    return vcat(paths)
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

function get_phases(Impurity_Data, rec_path_1, rec_path_2, basis_cut_mb, STEP, Total_H, Sub_Number_MB_Operator_List, Degeneracy)

    V1 = Impurity_Data.V0[1]
    V2 = Impurity_Data.V0[2]
    Imp_Site = [rec_path_1[1], rec_path_1[2], rec_path_2[1], rec_path_2[2]]
    V0 = [V1, 0, V2, 0] 

    Impurity_Data = Impurity(V0, Imp_Site)
    Impurity_H = Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data)
    E0, ψ = eigenstates(Impurity_H, Degeneracy)
    E0, ψ = fixed_pn_sector(pn, E0, ψ, basis_cut_mb) 
    
    ψ = hcat([ψ[i].data for i in 1:Degeneracy] ...) # matrix form
    # ψ = hcat([ψ[i].data for i in 1:length(E0)] ...)
    ψ_first = copy(ψ)

    Imp_Site_List = [ [imp, rec_path_1[idx+1], rec_path_2[idx], rec_path_2[idx+1] ] for (idx,imp) in (enumerate(rec_path_1[1:end-1])) ]
    V0_List = [ [V1*(1-step), V1*step, V2*(1-step), V2*step] for step in STEP ]
    
    Impurity_Data_List = [ [Impurity(V00, Imp_Sitee)] for Imp_Sitee in Imp_Site_List for V00 in V0_List ]
    Impurity_H_List = [ Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Dataa[1]) for Impurity_Dataa in Impurity_Data_List]
    
    ψ_tilde_list = []
    ψ_list = []
    A_inv_list = []
    A_list = []
    E_list = []

    @showprogress for H in Impurity_H_List
        
        ϵ, ψ_tilde = eigenstates(H, Degeneracy) 
        ϵ, ψ_tilde = fixed_pn_sector(pn, ϵ, ψ_tilde, basis_cut_mb) 
        push!(E_list, ϵ)
        
        ψ_tilde = hcat([ψ_tilde[i].data for i in 1:Degeneracy] ...)
        push!(ψ_tilde_list, ψ_tilde)

        # KM Algorithm #
        #= A = ψ'*ψ_tilde
        push!(A_list, A)
        A_inv = inv(A) 
        push!(A_inv_list, A_inv)
        ψ = ψ_tilde*A_inv =# 
        # !!!
        
        # !!! Vanderbilt
        A = ψ' * ψ_tilde
        push!(A_list, A)
        V, Σ, W = svd(A)
        M = V * W'
        ψ = ψ_tilde * M'
        # !!!

        # DOESN'T WORK !!!
        # Orthogonalize Method 1
        # DOESN'T WORK !!!
        # for i in 1:Degeneracy  # DOESN'T WORK !!!
        #     norm = sqrt(ψ[:,i]'*ψ[:,i])  # DOESN'T WORK !!!
        #     ψ[:,i] = ψ[:,i] ./ norm  # DOESN'T WORK !!!
        # end  # DOESN'T WORK !!!

        # Orthogonalize Method 
        # ψ = qr(ψ).Q * Matrix(I, size(ψ)...) # new vector
        push!(ψ_list, ψ)
        
    end

    # BerryMatrix = ψ' * ψ_first
    # BerryEnergies, BerryStates = eigen(Berry_Matrix)

    #ϕ_tot = -imag(log(det(ψ' * ψ_first)))
    
    return ψ, ψ_list, ψ_first, ψ_tilde_list, A_inv_list, A_list, E_list
end

E_Split(N) = E_list[N][2] - E_list[N][1] 
E_Gap(N) = E_list[N][3] - E_list[N][2] 
δ(N) = E_Split(N) / E_Gap(N)

function LaughlinDegeneracyBreaking(STEP)
    δ_List = []
    E_Split_List = []
    for N in 1:length(STEP)
        push!(E_Split_List, E_Split(N))
        push!(δ_List, δ(N))
    end
    return E_Split_List, δ_List
end

#   Berry Matrix
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# 
#   GÜNCELLEME!
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡
# 
#   22/04/24: HİÇ BİR YERDE KULLANMIYORUM!
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

# GÜNCELLEME!
# 22/04/24: BUNU HİÇ BİR YERDE KULLANMIYORUM!
function gram_schmidt(matrix)
    # orthogonalises the columns of the input matrix
    num_vectors = size(matrix)[2]
    orth_matrix = zeros(Complex{Float64},size(matrix))
    for vec_idx = 1:num_vectors
        orth_matrix[:, vec_idx] = matrix[:, vec_idx]
        sum = zeros(size(orth_matrix[:, 1]))
        for span_base_idx = 1:(vec_idx-1)
            # compute sum
            sum += dot(orth_matrix[:, span_base_idx], orth_matrix[:, vec_idx])*orth_matrix[:, span_base_idx]
        end
        orth_matrix[:, vec_idx] -= sum
        # normalise vector
        orth_matrix[:, vec_idx] = orth_matrix[:, vec_idx]/norm(orth_matrix[:, vec_idx])
    end
    return orth_matrix
end

# GÜNCELLEME!
# 22/04/24: BUNU HİÇ BİR YERDE KULLANMIYORUM!
function get_final_state(rec_path_1, rec_path_2, Degeneracy, delta_t)
    N_Pin = 4

    Imp_Site = [rec_path_1[1], rec_path_1[2], rec_path_2[1], rec_path_2[2]]
    V0 = [V1, 0, V2, 0]
    Impurity_Data = Impurity(V0, Imp_Site)
    
    Impurity_H = Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data)
    Impurity_H = dense((Impurity_H+Impurity_H')/2)
    
    E_Imp_0, U_Imp_0 = eigenstates(Impurity_H, Degeneracy)
    
    # Ground state enerjileri hesapla sadece!
    
    U_Imp_0 = hcat([U_Imp_0[i].data for i in 1:Degeneracy] ...)
    
    U_first = copy(U_Imp_0)
    
    STEP = 0:delta_t:1
    eig_list = []
    for (idx,imp) in ProgressBar(enumerate(rec_path_1[1:end-1]))
        Imp_Site = [imp, rec_path_1[idx+1], rec_path_2[idx], rec_path_2[idx+1]]
        for step in STEP
            V0 = [V1*(1-step), V1*step, V2*(1-step), V2*step]
            Impurity_Data = Impurity(V0, Imp_Site)
            Impurity_H = Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data)
            Impurity_H = dense((Impurity_H+Impurity_H')/2)
            E_Imp, U_Imp = eigenstates(Impurity_H, Degeneracy)
            push!(eig_list, U_Imp)
            
            U_Imp = hcat([U_Imp[i].data for i in 1:Degeneracy] ...) #convert to matrix
            A = U_Imp_0'*U_Imp #? Operator objesi ile matrisin çarpımı #?
            A_inv = inv(A)
            #U_Imp_0 = U_Imp*transpose(A_inv)
            U_Imp_0 = U_Imp*A_inv
            for i in 1:Degeneracy
                Norm_0 = sqrt(U_Imp_0[:,i]'*U_Imp_0[:,i])
                U_Imp_0[:,i] = U_Imp_0[:,i] ./ Norm_0
            end
        end
    end

    return U_Imp_0, U_first, eig_list
end

# GÜNCELLEME!
# 22/04/24: BUNU HİÇ BİR YERDE KULLANMIYORUM!
function Berry_Matrix(rec_path_1, rec_path_2, Degeneracy, delta_t)
    return get_final_state(rec_path_1, rec_path_2, Degeneracy, delta_t)[1]'*get_final_state(rec_path_1, rec_path_2, Degeneracy, delta_t)[2]
end

#   Simulation
#   ≡≡≡≡≡≡≡≡≡≡≡≡

function movie(delta_t, factor, Nx, Ny, Degeneracy, N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, co)
    println("!Nx and Ny must be equalt to each other!")
    eig_list = get_final_state(rec_path_1, rec_path_2, Degeneracy, delta_t)[3]
    STEP = 0:delta_t:1
    @gif for i in 1:length(STEP)*Nx # Bunu Nx'ten bağımsız yapmak mümkün!
        data = Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, eig_list[i])'
        Plt = Plots.heatmap(Interp(data, factor), aspect_ratio=:equal)
        # co = vcat( ( [y x] for x in 0:Nx-1 for y in 0:Ny-1 ) ... )
        # scatter!(co[:,1].+1,co[:,2].+1, series_annotations = text.([i for i in 1:N_Site], :bottom), legend=false)
    end
    return nothing
end