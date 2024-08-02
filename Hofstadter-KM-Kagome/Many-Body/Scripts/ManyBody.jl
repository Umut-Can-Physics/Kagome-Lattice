using Revise
includet("FirstBandApproximation.jl")

function plot_square_lattice(N, Nx, Ny)
    co = vcat( ( [y x] for x in 0:Ny-1 for y in 0:Nx-1 ) ... )
    p = scatter(co[:,1],co[:,2], series_annotations = text.([i for i in 1:N], :bottom), legend=false, aspect_ratio = :equal)
    return co, display(p)
end

function H_sp(N, Nx, Ny, p, q)
    sp_basis = NLevelBasis(N)
    periodicity = 0 #periodic (select 1 for hard-wall conditions)
    sp_matrix = Hofstadter_SP(Nx, Ny, p/q, periodicity)
    H1 = get_sp_op(N, sp_matrix)
    return H1, sp_basis
end

function Sp_Op(N, matrix)
    H = get_sp_op(N, matrix)
    return dense((H'+H)/2)
end

function H_sub(N, H1, cut_off)
    sub_states = get_sub_states(H1, cut_off)
    basis_sub, P, Pt = get_projector_op(N, sub_states)
    H1_sub = get_subspace_op(H1, P, Pt)
    return H1_sub, basis_sub, P, Pt
end

function H_Kin_MB(basis_sub, PN, H1_sub)
    states_mb = bosonstates(basis_sub, PN) 
    basis_mb = ManyBodyBasis(basis_sub, states_mb)
    H1_MB = get_mb_op(basis_mb, H1_sub)
    return H1_MB
end

function get_Bosonic_MB_Basis(N, PN)
    sp_basis = NLevelBasis(N)
    N_States = bosonstates(sp_basis, PN)
    N_Basis_MB = ManyBodyBasis(sp_basis, N_States)
    return N_Basis_MB, sp_basis
end

function H_Total_Sub(PN, P, Pt, basis_cut_mb, cut_off, U, H1_MB)
    H_Kin = SparseOperator(basis_cut_mb)
    H_Kin.data = H1_MB.data
    H_Int = Hubbard_Interaction_op(P, Pt, basis_cut_mb, cut_off, U)
    H_Total = H_Kin + H_Int
    return dense((H_Total'+H_Total)/2)
end

# find energies of specific particle number 
function fixed_pn_sector(pn, E, V, b)

    PN = []

    for i in 1:length(E)
        push!(PN, round(expect(number(b), V[i])))
    end
    
    ind = findall(x->x == 1, PN.==pn)

    ϵ_fixed = []

    for i in ind
        push!(ϵ_fixed, E[i])
    end

    λ_fixed = []

    for i in ind
        push!(λ_fixed, V[i])
    end

    return ϵ_fixed, λ_fixed
end