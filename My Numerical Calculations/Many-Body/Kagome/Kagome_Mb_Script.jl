function get_Fermionic_MB_basis(PN)

    N = Nx*Ny*3
    NBasis = NLevelBasis(N)
    # HC=Hardcore
    # Transition operator is bosonic operator due to fermionic states doesn't include -1 coefficient in transition matrix. So, we can construct hard core many body basis by using fermionic states. Note that we work only bosons on lattice.
    HC_States = fermionstates(NBasis, [PN]) 
    HC_Basis_MB = ManyBodyBasis(NBasis, HC_States);
    
    return HC_Basis_MB
end

function Kagome_Hard_Core(mb_basis, Nx, Ny, sp_op)
    
    N = Nx*Ny*3
    # Many-Body Hamiltonian using METHOD #1 (MANUEL CONSTRUCTING)
    HC_Hamiltonian_MB = SparseOperator(mb_basis)
    # Number of states for Hard-core bosonic system is equal to number of states for fermionic system
    for m in 1:N
        for n in 1:N
            # Phase Factors= H_NN[m,n] and H_NNN[m,n]
            # Hopping Terms= transition(HC_NBasis_MB, m, n)
            # The Neighbors Condition Hides Here: H[m,n]
            HC_Hamiltonian_MB += sp_op.data[m,n] * transition(mb_basis, m, n)
        end
    end
    return HC_Hamiltonian_MB
end

function get_Bosonic_MB_Basis(N, PN)

    NBasis = NLevelBasis(N)
    States = bosonstates(NBasis, PN) 
    Basis_MB = ManyBodyBasis(NBasis, States)
    
    return Basis_MB, NBasis
end

function Kagome_Finite_U(Nx,Ny,t1,L1,t2,L2,PN,U)
    
    N = Nx*Ny*3

    sp_op = get_SP_H(Nx, Ny, t1, L1, t2, L2)
    Basis_MB, NBasis =  get_Bosonic_MB_Basis(PN)
    
    # Kinetic term
    KT = SparseOperator(Basis_MB)

    # Interaction term
    IT = SparseOperator(Basis_MB)

    for m in 1:N
        # Occupation (total particle) Operator: number()
        IT = IT + U/2 * number(Basis_MB, m) * ( number(Basis_MB, m) - identityoperator(Basis_MB) ) 

        for n in 1:N
            KT = KT + sp_op.data[m,n] * transition(Basis_MB, m, n)
        end
        
    end
    
    MB_Hamiltonian = KT + IT
    
    return MB_Hamiltonian
end
