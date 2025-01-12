using Revise
includet("FirstBandApproximation.jl")

function plot_square_lattice(N, Nx, Ny)
    co = vcat( ( [y x] for x in 0:Ny-1 for y in 0:Nx-1 ) ... )
    p = scatter(co[:,1],co[:,2], series_annotations = text.([i for i in 1:N], :bottom), legend=false, aspect_ratio = :equal)
    return co, p
end

function find_middle_site_x(Nx, co)
    x_co = co[:,1]
    if iseven(Nx)
        idx = Int((Nx / 2))
        return findall(x->x==0, x_co)[idx]
    else
        idx = Int((Nx / 2)+1/2)
        return findall(x->x==0, x_co)[idx]
    end        
end

function find_middle_site_y(Ny, co)
    y_co = co[:,2]
    if iseven(Ny)
        idx = Int((Ny / 2))+1
        return findall(x->x==Ny-1, y_co)[idx]
    else
        idx = Int((Ny / 2)+1/2)
        return findall(x->x==Ny-1, y_co)[idx]
    end        
end

function Sp_Op(sp_basis, matrix)
    H = get_sp_op(sp_basis, matrix)
    return dense((H'+H)/2)
end

function H_sub(sp_basis, H1, cut_off)
    sub_states = get_sub_states(H1, cut_off)
    basis_sub, P, Pt = get_projector_op(sp_basis, sub_states)
    H1_sub = get_subspace_op(H1, P, Pt)
    return H1_sub, basis_sub, P, Pt
end

function get_Bosonic_MB_Basis(sp_basis, pn, HardCore)
    if HardCore==false
        N_States = bosonstates(sp_basis, pn)
        N_Basis_MB = ManyBodyBasis(sp_basis, N_States)
    elseif HardCore==true
        N_States = fermionstates(sp_basis, pn)
        N_Basis_MB = ManyBodyBasis(sp_basis, N_States)
    end
    return N_Basis_MB
end

function boson_mb_basis(Sub_Basis, pn, HardCore)
    if HardCore==false
        sub_boson_states = bosonstates(Sub_Basis, pn)
        sub_mb_basis = ManyBodyBasis(Sub_Basis, sub_boson_states)
    elseif HardCore==true
        sub_boson_states = fermionstates(Sub_Basis, pn)
        sub_mb_basis = ManyBodyBasis(Sub_Basis, sub_boson_states)
    end
    return sub_mb_basis
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
    return PN, ϵ_fixed, λ_fixed
end 

"""
Get total Hamiltonian (H_kin + H_int) in sub-space (projected space)
"""
function H_Hubbard_Projection(N, pn, U, matrix, Cut_Off, HardCore)

    basis_sp = NLevelBasis(N)
    H = Sp_Op(basis_sp, matrix)
    H_sp_Sub, Sub_Basis, P, Pt = H_sub(basis_sp, H, Cut_Off)

    if HardCore==true
        basis_cut_mb = get_Bosonic_MB_Basis(Sub_Basis, pn, true)
        #basis_cut_sp = NLevelBasis(Cut_Off)
        #basis_cut_mb = get_Bosonic_MB_Basis(basis_cut_sp, pn, true)

        # type of are different but the marix are the same in contect of H_TotalSub: Sub_Basis = basis_cut_sp 

        H_MB_sub = get_mb_op(basis_cut_mb, H_sp_Sub)
        #H_MB_sub = get_mb_hopping(basis_cut_mb, H_sp_Sub)
        basis_mb_sub = boson_mb_basis(Sub_Basis, pn, true)
        H_Int = Hubbard_Interaction_fixed_prtcl(basis_sp, U)
        H_Int_Sub = Hubbard_Int_fixed_prtc_sub(H_Int, P, Pt, Sub_Basis, basis_mb_sub)
        H_TotalSub = H_MB_sub + H_Int_Sub 
        H_TotalSub = (H_TotalSub'+H_TotalSub)/2  
    
    elseif HardCore==false
        basis_cut_mb = get_Bosonic_MB_Basis(Sub_Basis, pn, false)
        H_MB_sub = get_mb_op(basis_cut_mb, H_sp_Sub)
        #H_MB_sub = get_mb_hopping(basis_cut_mb, H_sp_Sub)
        basis_mb_sub = boson_mb_basis(Sub_Basis, pn, false)
        H_Int = Hubbard_Interaction_fixed_prtcl(basis_sp, U)
        H_Int_Sub = Hubbard_Int_fixed_prtc_sub(H_Int, P, Pt, Sub_Basis, basis_mb_sub)
        H_TotalSub = H_MB_sub + H_Int_Sub 
        H_TotalSub = (H_TotalSub'+H_TotalSub)/2  
    end

    return H_TotalSub, P, Pt, basis_mb_sub, Sub_Basis
end

"""
U is necessary for finite-interaction potential (Hard-Core==false)
"""
function H_Hubbard(N, U, pn, matrix, HardCore)
    basis_sp = NLevelBasis(N)
    H = Sp_Op(basis_sp, matrix)

    if HardCore==true
        basis_mb = get_Bosonic_MB_Basis(basis_sp, pn, true)
        H_MB = get_mb_op(basis_mb, H)
        H_Total_full = H_MB 
        H_Total_full = (H_Total_full'+H_Total_full)/2
    elseif HardCore==false
        basis_mb = get_Bosonic_MB_Basis(basis_sp, pn, false)
        H_MB = get_mb_op(basis_mb, H)
        H_Int = Hubbard_Interaction_Full(N, basis_sp, basis_mb, U)
        H_Total_full = H_MB + H_Int
        H_Total_full = (H_Total_full'+H_Total_full)/2
    end

    return H_Total_full, basis_mb, basis_sp
end