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

#= function get_mb_op(mb_basis, sp_op)
    # Initialize thread-local SparseOperators
    thread_results = [SparseOperator(mb_basis) for _ in 1:nthreads()]
    
    N = sp_op.basis_l.shape[1]

    @threads for i in 1:N
        thread_id = threadid()  # Identify thread
        for j in 1:N
            local_mb = deepcopy(mb_basis)  # Ensure a thread-local copy
            thread_results[thread_id] += sp_op.data[i, j] * transition(mb_basis, i, j)
        end
    end
    
    # Combine thread-local results
    combined_result = SparseOperator(mb_basis)
    for result in thread_results
        combined_result += result
    end
    
    return combined_result
end =#

function get_mb_op(mb_basis, sp_op)
    
    mb_op = SparseOperator(mb_basis)
    
    N = sp_op.basis_l.shape[1]
    
    for i in 1:N
        for j in 1:N
            mb_op += sp_op.data[i,j] * transition(mb_basis, i, j)
        end
    end
    
    return mb_op
end

function get_mb_hopping(mb_basis, sp_op)
    
    mb_op = SparseOperator(mb_basis)
    
    N = sp_op.basis_l.shape[1]

    neighbor_list = [[1, 0], [-1,0], [0,1], [0,-1]]
    
    for j in 1:N
        for neighbor in neighbor_list
	    jx = mod(j-1,Nx) + 1
	    jy = Int( (j-jx)/Nx ) + 1
	    ix = mod( jx + neighbor[1] - 1, Nx ) + 1
	    iy = mod( jy + neighbor[2] - 1, Ny ) + 1
	    i = (iy-1)*Nx + ix 
	    #println(j,jx,jy," ",i,ix,iy)
            mb_op += sp_op.data[i,j] * transition(mb_basis, i, j)
        end
    end
    
    return mb_op
end



#! get_mb_op and get_mb_op2 have to be the same return !#
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

function Hubbard_Interaction_fixed_prtcl(sp_basis, U)

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
    #V2_Sub_Op = get_mb_op(basis_mb_sub, V2_sub_op) # optimized!

    return V2_Sub_Op
end

function get_num_mb_op(N, cut_sp_basis, num_sub_list, cut_mb_basis)

    number_sp_list = [Operator(cut_sp_basis, num_sub_list[i].data) for i in 1:N]
    number_mb_list = [get_mb_op(cut_mb_basis, number_sp_list[i]) for i in 1:N]
    
    return number_mb_list
end
