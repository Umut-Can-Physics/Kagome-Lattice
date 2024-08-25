# using QuantumOptics
using Einsum 

"""
Convert Sp Matrix to Sp Operator

#### Arguments
- `sp_basis::NLevelBasis`: Single-Particle basis.
- `N::Integer`: The total site number.
- `sp_matrix::Matrix`: Hopping phases matrix from any model.
"""
function get_sp_op(sp_basis, sp_matrix)

    #sp_basis = NLevelBasis(N)
    
    H = SparseOperator(sp_basis)

    N, = size(sp_matrix)
    
    for m in 1:N
        for n in 1:N
            H += sp_matrix[m,n] * transition(sp_basis, m, n)
        end
    end
    
    return H
end

"""
Compute the first eigen-states in a given number.

#### Arguments
- `sp_op::get_sp_op`: Single-particle operator.
- `cut_off::Int`: Limit value of eigen-state index.
"""
function get_sub_states(sp_op, cut_off)
    
    E0, states0 = eigenstates(dense(sp_op))
    states = states0[1:cut_off]
    
    return states
end

"""
Compute sub-space basis, projection and comlex conjugate of projection operator.

#### Arguments
- `states::get_sub_states`: Eigen-states of the sub-space.
- `basis::NLevelBasis`: Single-particle basis.
"""
function get_projector_op(sp_basis, states)

    sub_basis = SubspaceBasis(sp_basis,states)
    P = projector(sub_basis, sp_basis)
    Pt = dagger(P)
    
    return sub_basis, P, Pt
end

#   <font size="5"> \hat{O}_{sub}=P\hat{O}P^\dagger </font>

"""
Compute the corresponding operator in the sub-space.

#### Arguments
- `sp_op::Operator`: Single-particle operator from single-particle matrix.
- `P::get_projector_op[2]`: Projection operator.
- `Pt::get_projector_op[3]`: Complex conjugate of projection operator.
"""
function get_subspace_op(sp_op, P, Pt)
    return P*sp_op*Pt
end

#   <font size="5"> n_i=a_i^\dagger a_i </font>

"""
Compute the single-particle number operator for each lattice sites.

#### Arguments
- `N::Integer`: The total site number.
- `sp_basis::NLevelBasis`: Single-Particle basis.
- `P::get_projector_op[2]`: Projection operator.
- `Pt::get_projector_op[3]`: Complex conjugate of projection operator.
"""
function get_num_sub_list(N, P, Pt)
    sp_basis = NLevelBasis(N)
    num_sub_list = []
    for m in 1:N
        num_op = transition(sp_basis, m, m)
        num_sub_op = get_subspace_op(num_op, P, Pt)
        push!(num_sub_list, num_sub_op)
    end
    return num_sub_list
end

function get_num_list(N)
    sp_basis = NLevelBasis(N)
    num_list = []
    for m in 1:N
        num_op = transition(sp_basis, m, m)
        push!(num_list, num_op)
    end
    return num_list
end

#   <font size="5"> \hat{O}=\sum_{ij} a^\dagger_i a_j <u_i|\hat{o}|u_j> </font>

"""
Compute the many-body operator for boson particles from single-particle operator.

#### Arguments
- `mb_basis`: Many-body basis.
- `sp_op::Operator`: Single-particle operator.
"""
function get_mb_op(mb_basis, sp_op)
    
    mb_op = SparseOperator(mb_basis)
    
    #N = sp_op.basis_l.N
    N = sp_op.basis_l.shape[1]
    
    for i in 1:N
        for j in 1:N
            mb_op += sp_op.data[i,j] * transition(mb_basis, i, j)
        end
    end
    
    return mb_op
end

#! get_mb_op and get_mb_op2 have to be the same return value !#
function get_mb_op2(mb_basis, sp_op)
    return manybodyoperator(mb_basis, sp_op)
end

#   <font size="5"> \hat{V}=\sum_{ijkl}a^\dagger_ia^\dagger_ja_ka_l
#   <u_i|<u_j|\hat{v}|u_k>|u_l> </font>

#= function create_tensor(basis, dimension, size)
    Createe = zeros(ComplexF64, dimension, dimension, size)
    for i in 1:size
        Createe[:,:,i] = dense(create(basis, i)).data
    end
    return Createe
end

function destroy_tensor(basis, dimension, size)
    Destroyy = zeros(ComplexF64, dimension, dimension, size)
    for i in 1:size
        Destroyy[:,:,i] = dense(destroy(basis, i)).data
    end
    return Destroyy
end =#

function Hubbard_Interaction_op(P, Pt, cut_mb_basis, cut_off, U)
    
    P1 = P.data
    P1t = Pt.data

    @einsum coefficient[k,l,m,n] := P1[k,i] * P1[l,i] * P1t[i,m] * P1t[i,n]

    Vint_mb_cut = SparseOperator(cut_mb_basis)

    A = [destroy(cut_mb_basis, k) for k in 1:cut_off]
    At = [create(cut_mb_basis, k) for k in 1:cut_off]
    
    @showprogress for k in 1:cut_off
        for l in 1:cut_off
            for m in 1:cut_off
                for n in 1:cut_off
                    Vint_mb_cut += U/2 * coefficient[k,l,m,n] * At[k] * At[l] * A[m] * A[n]
                end
            end
        end
    end
    
    return Vint_mb_cut
end

function Hubbard_Interaction_Full(N, sp_basis, mb_basis, U)
    basis2 = sp_basis ⊗ sp_basis
    
    Vint2 = SparseOperator(basis2)

    for n in 1:N
        Vint2 += U/2*transition(sp_basis,n,n)⊗transition(sp_basis,n,n)
    end

    Vint_mb = manybodyoperator(mb_basis, Vint2)

    return Vint_mb
end

function Hubbard_Interaction_fixed_prtcl(sp_basis)

    basis2 = sp_basis ⊗ sp_basis
    
    Vint2 = SparseOperator(basis2)

    for n in 1:N
        Vint2 += U/2*transition(sp_basis,n,n)⊗transition(sp_basis,n,n)
    end

    #Vint_mb = manybodyoperator(mb_basis, Vint2)

    return Vint2
end

function Hubbard_Int_fixed_prtc_sub(Vint_mb, P, Pt, basis_sub, basis_mb_sub)
    V_int = Vint_mb.data
    P2 = (P⊗P).data # two body projection operator
    P2t = (Pt⊗Pt).data
    @einsum V2_sub[i2,j2] :=  P2[i2,k2] * V_int[k2,l2] * P2t[l2,j2]
    basis_sub_2 = basis_sub ⊗ basis_sub # two body sub basis
    V2_sub_op = Operator(basis_sub_2, V2_sub)
    V2_Sub_Op = manybodyoperator(basis_mb_sub, V2_sub_op)

    return V2_Sub_Op
end

function get_num_mb_op(N, cut_sp_basis, num_sub_list, cut_mb_basis)

    number_sp_list = [Operator(cut_sp_basis, num_sub_list[i].data) for i in 1:N]
    number_mb_list = [get_mb_op(cut_mb_basis, number_sp_list[i]) for i in 1:N]
    
    return number_mb_list
end

#= # find energies of specific particle number 
function get_filtered_energies(pn, E, V, basis)
    PN_Energies = Array{Float64}(undef, length(E), 2)
    for i in 1:length(E)
        PN_Energies[i] = round(expect(number(basis), V[i])) 
        PN_Energies[i,2] = E[i] 
    end
    
    # filter
    df = DataFrame(PN_Energies, :auto)
    df = filter(row -> (row.x1 == pn),  df)
    
    return Matrix(df)[:,2]
end

#   Allta ki fonksiyonun çalışması için, dizide ki filtre edilmiş parçacık
#   sayısı her zaman en büyük değer de olmalıdır. Örneğin, PN=[0,1,2,3,4] iken
#   filtre edilen parçacık sayısı pn=4 olmalıdır!

# Eigenstates of filtered particles
function Restricted_Hubbard_States(states, filtered_energies)
    number_of_states = length(filtered_energies)
    return states[1:number_of_states];
end =#

#= function Get_Density_Profile(N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States, index)
    Expectation_List = []
    for site in 1:N_Site
        push!(Expectation_List, expect(Sub_Number_MB_Operator_List[site], Fil_States[index]))
    end
    return real(Expectation_List)
end =#

#= function Get_Avg_Density(Nx, Ny, Degeneracy, N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States)
    Avg_Density = spzeros(Nx,Ny)
    for index in 1:Degeneracy
        Avg_Density += reshape(Get_Density_Profile(N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States, index), Nx, Ny)
    end    
    return Avg_Density / Degeneracy
end =#

#= function Get_Avg_Density_List(Nx, Ny, Degeneracy, N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States)
    Avg_Density = zeros(N_Site)
    for index in 1:Degeneracy
        Avg_Density += Get_Density_Profile(N_Site, Sub_Number_MB_Operator_List, Basis_Cut_MB, Fil_States, index)
    end    
    return Avg_Density / Degeneracy
end =#

function Interp(data, factor)
    IC = CubicSplineInterpolation((axes(data,1), axes(data,2)), data)
    finerx = LinRange(firstindex(data,1), lastindex(data,1), size(data,1) * factor)
    finery = LinRange(firstindex(data,2), lastindex(data,2), size(data,2) * factor)
    nx = length(finerx)
    ny = length(finery)
    data_interp = Array{Float64}(undef,nx,ny)
    for i ∈ 1:nx, j ∈ 1:ny
        data_interp[i,j] = IC(finerx[i],finery[j])
    end
    return finery, finerx, data_interp
end

# Find x and y coordinates from given site index

function exp_list0(Nx, Ny, a1_vec, a2_vec, Basis, site_indx, avg_density)
    
    x_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[4], 1:N)
    y_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[5], 1:N)
    
    x = hcat(x_co, y_co)[site_indx, 1]
    y = hcat(x_co, y_co)[site_indx, 2] 
    
    #!!! 
    # Burada beklenen değerlerin sıralamasının site bazında olduğunu varsaydım!!!!
    #!!!
    exp_val = real(avg_density)[site_indx] 
    
    return x, y, exp_val
end

# Find site_index from given x and y coordinates

function exp_list1(Xx, Yy)
    co_list = hcat(x_co, y_co)
    site_indx = intersect(findall(x->x==Xx, co_list[:,1]), findall(x->x==Yy, co_list[:,2]))
    return real(avg_density)[site_indx] 
end