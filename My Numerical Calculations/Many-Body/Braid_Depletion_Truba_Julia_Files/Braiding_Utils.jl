#   Impurity
#   ≡≡≡≡≡≡≡≡≡≡

struct Impurity
    V0::Vector{Float64}
    Imp_Site::Vector{Int64}
end

function Imp_H(Total_H, Sub_Number_MB_Operator_List, Impurity_Data)
    for imp in 1:length(Impurity_Data.V0)
        Total_H += Impurity_Data.V0[imp] * Sub_Number_MB_Operator_List[Impurity_Data.Imp_Site[imp]]
    end
    return Total_H
end

#   Moving Quasiholes
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

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