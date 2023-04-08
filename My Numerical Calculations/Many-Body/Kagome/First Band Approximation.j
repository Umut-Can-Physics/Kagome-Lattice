using NBInclude
@nbinclude("Kagome SP.ipynb")
using QuantumOptics

# https://juliapackages.com/p/einsum
using Einsum

using BenchmarkTools

#using LinearAlgebra
Nx = 3; Ny = 3; N=Nx*Ny*3; cut_off = 5
PN = [0,1,2]
U = 0.001
# t1=-1;L1=0 ;t2=L2=0
# t1=L1=0;t2=L2=-1
t1 = -1;L1 = 0.28;t2 = 0.3;L2 = 0.2;
# t1=t2=-1;L1=L2=0

sp_basis = NLevelBasis(N)
sp_matrix = get_SP_H(Nx, Ny, t1, L1, t2, L2);
function get_sp_op(sp_basis, sp_matrix)
  
    H = SparseOperator(sp_basis)

    for m in 1:N
        for n in 1:N
            H += sp_matrix.data[m,n] * transition(sp_basis, m, n)
        end
    end
    
    return H
end
H1 = get_sp_op(sp_basis, sp_matrix);
#check operator form
eigenenergies(dense(H1)) == eigenenergies(dense(sp_matrix))
function get_sub_states(sp_op, cut_off)
    
    E0, states0 = eigenstates(dense(sp_op))
    states = states0[1:cut_off]
    
    return states
end
sub_states = get_sub_states(H1, cut_off);
function get_projector_op(states, basis)
    
    b_sub = SubspaceBasis(basis,states)
    P = projector(b_sub, basis)
    Pt = dagger(P)
    
    return b_sub, P, Pt
end
b_sub, P, Pt = get_projector_op(sub_states, sp_basis);
function get_subspace_op(sp_op, P, Pt)
    return P*sp_op*Pt
end
H1_sub = get_subspace_op(H1, P, Pt);
states_mb = bosonstates(b_sub, PN) 
basis_mb = ManyBodyBasis(b_sub, states_mb);
function get_mb_op(basis_mb, basis, sp_op)
    
    Op_MB = SparseOperator(basis_mb)
    
    for i in 1:length(basis)
        for j in 1:length(basis)
            Op_MB += sp_op.data[i,j] * transition(basis_mb, i, j)
        end
    end
    
    return Op_MB
end
H1_MB = get_mb_op(basis_mb, b_sub, H1_sub);
@nbinclude("Kagome MB .ipynb"; regex=r"#.*executeme")
function get_hubbard_int(P, Pt, b_sub, cut_off)
    
    bcut_mb = get_Bosonic_MB_Basis(cut_off, PN)
    
    @time begin
    P1 = P.data
    P1t = Pt.data;

    @einsum P4[k,l,m,n] := P1[k,i] * P1[l,i] * P1t[i,m] * P1t[i,n]

    Vint_mb_cut = SparseOperator(bcut_mb)
        
    for k in 1:cut_off
        for l in 1:cut_off
            for m in 1:cut_off
                for n in 1:cut_off
                    a1t = create(bcut_mb, k)
                    a2t = create(bcut_mb, l)
                    a2  = destroy(bcut_mb, m)      
                    a1  = destroy(bcut_mb, n)      
                    Vint_mb_cut += U/2*P4[k,l,m,n]*a1t*a2t*a2*a1
                end
            end
        end
    end
    end
    
    return Vint_mb_cut
end
H1_Int = get_hubbard_int(P, Pt, b_sub, cut_off);
bcut_mb = get_Bosonic_MB_Basis(cut_off,PN)
H1cut = SparseOperator(bcut_mb)
H1cut.data = H1_MB.data
H_MB = H1cut + H1_Int;
using DataFrames

function get_energies(pn, Energies, states, basis)
    PN_Energies = Array{Float64}(undef, length(states[1]), 2)
    for i in 1:length(Energies)
        PN_Energies[i] = round(expect(number(basis), states[i])) #expected values (first column)
        PN_Energies[i,2] = E[i] #eigen-values (second column)
    end
    df = DataFrame(PN_Energies, :auto)
    df = filter(row -> (row.x1 == pn),  df)
    
    return df
end
pn = 3
E, S = eigenstates(dense(dense((H_MB+dagger(H_MB))/2)))
get_energies(pn, E, S, bcut_mb);
function get_num_sub_list(sp_basis, P, Pt)
    num_sub_list = []
    for m in 1:N
        NM = transition(sp_basis, m, m)
        NMP = get_subspace_op(NM, P, Pt)
        push!(num_sub_list, NMP)
    end
    return num_sub_list
end
num_sub_list = get_num_sub_list(sp_basis,P,Pt);
function get_mb_op(basis_mb, basis, sp_op)
    
    Op_MB = SparseOperator(basis_mb)
    
    for i in 1:length(basis)
        for j in 1:length(basis)
            Op_MB += sp_op.data[i,j] * transition(basis_mb, i, j)
        end
    end
    
    return Op_MB
end
function get_mb_op(basis_mb, basis, sp_op)
    
    Op_MB = SparseOperator(basis_mb)
    
    for i in 1:length(basis)
        for j in 1:length(basis)
            Op_MB += sp_op.data[i,j] * transition(basis_mb, i, j)
        end
    end
    
    return Op_MB
end
function get_num_mb_list(basis_mb, basis, num_sub_list)
    
    num_mb_list = []
    
    for m in 1:N
        NMP = get_mb_op(basis_mb, basis, num_sub_list[m])
        push!(num_mb_list, NMP)
    end
    
    return num_mb_list
end
num_mb_list = get_num_mb_list(basis_mb, b_sub, num_sub_list);
NM_MB_Array_Storage = zeros(Complex{Float64},length(bcut_mb),length(bcut_mb),N);
NM_MB_Matrix = zeros(Complex{Float64},length(bcut_mb),length(bcut_mb));
for m in 1:N
    for i in 1:length(bcut_mb)
        for j in 1:length(bcut_mb)
            NM_MB_Matrix[i,j] = num_mb_list[m].data[i,j]
        end
    end
    NM_MB_Array_Storage[:,:,m] = NM_MB_Matrix
end
BL = BR = bcut_mb
index_number_op = 6
T = NM_MB_Array_Storage[:,:,index_number_op]
Op = Operator(BL,BR,T);
expect(Op, S[1]);
index_eig_states = 3
Sum = 0
expect_list=[]
for i in 1:N
    T = NM_MB_Array_Storage[:,:,i]
    Op = Operator(BL,BR,T)
    Sum += expect(Op, S[index_eig_states])
    println(i,"\t",expect(Op, S[index_eig_states]))
    push!(expect_list,expect(Op, S[index_eig_states]))
end
plot_kagome(Nx,Ny)
# Find x and y coordinates from gicen site index

function exp_list0(site_indx)
    
    x_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[4], 1:Nx*Ny*3)
    y_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[5], 1:Nx*Ny*3)
    
    x = hcat(x_co, y_co)[site_indx, 1]
    y = hcat(x_co, y_co)[site_indx, 2] 
    exp_val = real(expect_list)[site_indx] 
    
    return x, y, exp_val
end
exp_list0(1);
# Find site_index from given x and y coordinates

function exp_list1(Xx, Yy)
    co_list = hcat(x_co, y_co)
    site_indx = intersect(findall(x->x==Xx, co_list[:,1]), findall(x->x==Yy, co_list[:,2]))
    return real(expect_list)[site_indx] 
end
x_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[4], 1:Nx*Ny*3)
y_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[5], 1:Nx*Ny*3)
site_index = 5
Xx = x_co[site_index]
Yy = y_co[site_index]
Xx = AA[:,1];Yy = AA[:,2]