function braiding_path(moving_point, Nx, Ny, co)
    First_Path = Int64[]
    mod_list = [mod(moving_point,Nx),mod(moving_point,Ny)]
    bottom_site_in_the_Ny_direction = mod_list[mod_list .> 0][1]
    last_site_in_the_Ny_direction = (bottom_site_in_the_Ny_direction+(Nx*(Ny-1)))
    site_number_from_upward = (last_site_in_the_Ny_direction-moving_point)/Nx
    site_number_from_downward = (moving_point-bottom_site_in_the_Ny_direction)/Nx
    for i in 0:site_number_from_upward
        step_site = moving_point+i*Nx
        push!(First_Path, step_site)
        if step_site == last_site_in_the_Ny_direction
            for j in 0:site_number_from_downward
                push!(First_Path, bottom_site_in_the_Ny_direction+j*Nx)
            end
        end
    end
    x_co_mov_point = co[:,2][moving_point]
    A = filter!(ϵ->ϵ ∉ moving_point , findall(x->x==x_co_mov_point,co[:,2])) # 13, 14, 16, 17, 18
    AA = findall(x->x==x_co_mov_point,co[:,2]) # 13, 14, 15, 16, 17, 18
    B = findall(x->x==moving_point, AA)[1] # 3
    Second_Path = vcat(reverse(A[1:B-1]), reverse(A[B:end]))
    Third_Path = reverse(First_Path)
    Fourth_Path = push!(reverse(Second_Path), moving_point)
    rec_path_braiding = vcat(First_Path,Second_Path,Third_Path,Fourth_Path)
    return rec_path_braiding
end
 
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

function double_exchange_path(lens, start_point_1, start_point_2)
    dirs_1 = [1,-Ny, -1, Ny]
    rec_path_1 = rectangular_path(start_point_1,lens,dirs_1)
    dirs_2 = [-1, Ny, 1, -Ny]
    rec_path_2 = rectangular_path(start_point_2,lens,dirs_2)
    return rec_path_1, rec_path_2
end

function exchange_path(lens, start_point_1, start_point_2)
    dirs_1 = [1,-Ny]
    rec_path_1 = rectangular_path(start_point_1,lens,dirs_1)
    dirs_2 = [-1, Ny]
    rec_path_2 = rectangular_path(start_point_2,lens,dirs_2)
    return rec_path_1, rec_path_2
end

function braiding_path(moving_point, fixed_point, Nx, Ny, co)
    braid_path = braiding_path(moving_point, Nx, Ny, co)
    path_1 = braid_path
    path_2 = repeat([fixed_point], length(braid_path))
    return path_1, path_2
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
    p1 = scatter(co[:,1],co[:,2], series_annotations = text.([i for i in 1:Nx*Ny], :bottom), aspect_ratio = :equal, label="Lattice Sites", title="QH's Paths", xlabel=L"x",ylabel=L"y")
    p2 = plot!(M[:, 1], M[:, 2], linewidth=2, c=:blue, label="Path for first QH")
    p3 = plot!(N[:, 1], N[:, 2], linewidth=2, c=:red, label="Path for second QH")
    return savefig(p3, "Hofstadter-KM-Kagome/Many-Body/Kapit-Mueller(KM)/Braiding_Data/qh_moving_paths.png")
end

function th_AB_phase(pn, p, q, N_Pin, N_mov, number_of_plaq)
    NPhi = Int( Nx * Ny * (p/q) )
    charge = pn/(NPhi-N_Pin)
    # exch path: N_mov = 1
    # double exch path: N_mov = 2
    θ_AB = N_mov * (p/q) * charge * number_of_plaq
    ex = exp(2*im*pi*θ_AB)
    return θ_AB, ex, charge
end

function get_phases(Impurity_Data, rec_path_1, rec_path_2, basis_cut_mb, STEP, Total_H, Sub_Number_MB_Operator_List, Degeneracy)

    # INITIAL POSITION 
    V1 = Impurity_Data.V0[1]
    V2 = Impurity_Data.V0[2]
    Imp_initial = [rec_path_1[1], rec_path_1[2], rec_path_2[1], rec_path_2[2]]
    V_initial = [V1, 0, V2, 0] 
    Impurity_Data_initial = Impurity(V_initial, Imp_initial)

    Impurity_H = Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data_initial)
    E_initial, psi_initial = eigenstates(Impurity_H, Degeneracy)
    
    ψ_mat = hcat([psi_initial[i].data for i in 1:Degeneracy] ...) # matrix form
    ψ_first = copy(ψ_mat)

    Converge = []
    ψ_op = []
    ψ_mat_list = []
    imp_data = []
    ϵ_list = []

    step_break = STEP[2]-STEP[1]

    @showprogress for (idx,imp) in (enumerate(rec_path_1[1:end-1]))
        
        #Imp_Site_step = [imp, rec_path_1[idx+1], rec_path_1[mod(idx+1,length(rec_path_1))+1],rec_path_2[idx], rec_path_2[idx+1], rec_path_2[mod(idx+1,length(rec_path_2))+1]]
        Imp_Site_step = [imp, rec_path_1[idx+1], rec_path_2[idx], rec_path_2[idx+1]]
        
        for step in STEP
                
            #V0_step = [V1*(1-step-(step_break/10)), V1*(step+step_break/10), V1/100,V2*(1-step), V2*step, 0*V2/100]
            V0_step = [round(V1*(1-step),digits=2), round(V1*step, digits=2), round(V2*(1-step), digits=2), round(V2*step, digits=2)]

            Impurity_Data_step = [Impurity(V0_step, Imp_Site_step)]
            push!(imp_data, Impurity_Data_step)
           
            H = Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data_step[1])
           
            ϵ, ψ_tilde = eigenstates(H, Degeneracy) 
            push!(ψ_op, ψ_tilde)
            push!(ϵ_list, ϵ)
            
            ψ_tilde = hcat([ψ_tilde[i].data for i in 1:Degeneracy] ...)
            
            # KM Algorithm:

            A = ψ_mat'*ψ_tilde
            push!(Converge, abs(det(A)))
            A_inv = inv(A) 
            ψ_mat = ψ_tilde*A_inv
            ψ_mat = qr(ψ_mat).Q * Matrix(I, size(ψ_mat)...) # new vector
            
            # Vanderbilt:
            
            #= A = ψ_mat' * ψ_tilde
            push!(Converge, abs(det(A)))
            V, Σ, W = svd(A)
            M = V * W'
            ψ_mat = ψ_tilde * M'      
            push!(ψ_mat_list, ψ_mat) =#
        end
    end
    return ψ_first, ψ_mat, Converge, ψ_op, imp_data, ψ_mat_list, ϵ_list
end

function get_phases2(Impurity_Data, rec_path_1, rec_path_2, basis_cut_mb, STEP, Total_H, Sub_Number_MB_Operator_List, Degeneracy)

    V1 = Impurity_Data.V0[1]
    V2 = Impurity_Data.V0[2]
    Imp_Site = [rec_path_1[1], rec_path_1[2], rec_path_2[1], rec_path_2[2]]
    V0 = [V1, 0, V2, 0] 

    Impurity_Data = Impurity(V0, Imp_Site)
    Impurity_H = Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data)
    E0, ψ = eigenstates(Impurity_H, Degeneracy)
    
    ψ = hcat([ψ[i].data for i in 1:Degeneracy] ...) # matrix form
    ψ_first = copy(ψ)
    
    ψ_tilde_list = []
    ψ_tilde_op_list = []
    ψ_list = []
    ψ_updated = []
    A_inv_list = []
    A_list = []
    E_list = []
    imp_data_list = []
    Conv_list = []

    @showprogress for (idx,imp) in (enumerate(rec_path_1[1:end-1]))
        
        Imp_Site = [imp, rec_path_1[idx+1], rec_path_2[idx], rec_path_2[idx+1]]

        #= 
        push!(imp_site_list, Imp_Site)
        =#

        index = 0; Δ_step = 0.1; step = 0
        while index < 100
            index += 1
            step += Δ_step
        #for (index,step) in enumerate(STEP)
                

            #= 
            if index >= 401 && index <= 405
                continue
            end ; push!(step_list, step) 
            =#

            V0 = [V1*(1-step), V1*step, V2*(1-step), V2*step]

            Impurity_Data = [Impurity(V0, Imp_Site)]
           

            H = Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data[1])
           
            ϵ, ψ_tilde = eigenstates(H, Degeneracy) 
            
            
            ψ_tilde = hcat([ψ_tilde[i].data for i in 1:Degeneracy] ...)
            

            # KM Algorithm
            #= A = ψ'*ψ_tilde
            push!(A_list, A)
            A_inv = inv(A) 
            push!(A_inv_list, A_inv)
            ψ = ψ_tilde*A_inv
            push!(ψ_updated, ψ)  
            ψ = qr(ψ).Q * Matrix(I, size(ψ)...) # new vector =#

            # Vanderbilt:
            A = ψ' * ψ_tilde
            
            Conv = abs(det(A))
            
            if Conv > 0.9
                push!(A_list, A)    
                push!(Conv_list, Conv)
                push!(imp_data_list, Impurity_Data)
                push!(E_list, ϵ)
                push!(ψ_tilde_op_list, ψ_tilde)
                push!(ψ_tilde_list, ψ_tilde)

                V, Σ, W = svd(A)
                M = V * W'
                ψ = ψ_tilde * M'

                push!(ψ_list, ψ) 
                step += 1.5*Δ_step
            else
                step -= Δ_step
                Δ_step = Δ_step/2
            end
            
        end
    end

    # BerryMatrix = ψ' * ψ_first
    # BerryEnergies, BerryStates = eigen(Berry_Matrix)

    #ϕ_tot = -imag(log(det(ψ' * ψ_first)))
    
    return ψ, ψ_updated, ψ_list, ψ_first, ψ_tilde_list, ψ_tilde_op_list, A_inv_list, A_list, E_list, imp_data_list, Conv_list
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