import numpy as np

def calc_states_hofstadter(Nkx,Nky,p,q,kkx,kky,H):
    EEA=np.empty([Nkx,Nky,q])
    UUA=np.empty([Nkx,Nky,q,q],dtype=complex)
    for ikx, kx in enumerate(kkx):
        for iky, ky in enumerate(kky):
            EEA[ikx,iky,:],UUA[ikx,iky,:,:] = np.linalg.eigh(H(p,q,kx,ky))
    return UUA

def calc_states_twist_hofstadter(L_x, L_y, Nt1, Nt2, p, q, theta_x, theta_y, HMat_Theta):
    EEA=np.empty([Nt1,Nt2,L_x*L_y])
    UUA=np.empty([Nt1,Nt2,L_x*L_y,L_x*L_y],dtype=complex)
    for it1 in range(Nt1):
        for it2 in range(Nt2):
            EEA[it1,it2,:],UUA[it1,it2,:,:] = np.linalg.eigh(HMat_Theta(L_x, L_y, p, q, theta_x[it1], theta_y[it2]))
    return EEA, UUA

def calc_states_twist_kapit_mueller(L_x, L_y, Nt1, Nt2, p, q, Tx, Ty, Kapit_Mueller_Hamiltonian):
    EEA=np.empty([Nt1, Nt2, L_x*L_y])
    UUA=np.empty([Nt1,Nt2,L_x*L_y,L_x*L_y],dtype=complex)
    for it1 in range(Nt1):
        for it2 in range(Nt2):
            EEA[it1,it2,:],UUA[it1,it2,:,:] = np.linalg.eigh(Kapit_Mueller_Hamiltonian(L_x, L_y, p, q, Tx[it1], Ty[it2]))
    return EEA, UUA

def calc_states_twist_kagome(l1, l2, t1,L1,t2,L2,Nt1, Nt2, theta_1, theta_2, H):
    N = l1*l2*3
    EEA=np.zeros([Nt1,Nt2,N])
    UUA=np.zeros([Nt1,Nt2,N,N],dtype=complex)
    for it1 in range(Nt1):
        for it2 in range(Nt2):
             EEA[it1,it2,:],UUA[it1,it2,:,:] = np.linalg.eigh(H(l1,l2,t1,L1,t2,L2,theta_1[it1],theta_2[it2]))
    return UUA

def calc_states_kagome(Nkx,Nky,b1,b2,a1_b,a2_b,t1,L1,t2,L2,Hamiltonian):
    EEA=np.zeros([Nkx,Nky,3])
    UUA=np.zeros([Nkx,Nky,3,3],dtype=complex)
    for ikx in range(Nkx):
        for iky in range(Nky):
            kx,ky=ikx/Nkx*b1+iky/Nky*b2
            k1=np.dot([kx,ky],a1_b);k2=np.dot([kx,ky],a2_b);k3=k2-k1 # Kagome Convention
            EEA[ikx,iky,:],UUA[ikx,iky,:,:] = np.linalg.eigh(Hamiltonian(t1,L1,t2,L2,k1,k2,k3))
    return UUA

def calc_link_var(ψk1, ψk2):
    #return np.exp( 1j*np.angle( np.dot(np.conj(ψk1),ψk2) ) )
    s1 = np.dot(np.conj(ψk1),ψk2)
    return s1/np.abs(s1)

def calc_link_var_twist(ψk1, ψk2):
    s = np.linalg.det(np.matmul(np.conjugate(np.transpose(ψk1)), ψk2))
    return s/np.abs(s)
    
def calc_link_vars_BZ(UUA,bi):
    N1, N2 = UUA.shape[0:2]
    UU=np.zeros([N1,N2,2],dtype=complex)
    dirs = [[1,0], [0,1]]
    for i1 in range(N1):
        for i2 in range(N2):
            for idir, vdir in enumerate(dirs):
                    UU[i1,i2,idir] = calc_link_var(UUA[i1                   ,i2               ,:,bi],
                                                   UUA[np.mod(i1+vdir[0],N1),np.mod(i2+vdir[1],N2),:,bi])
    return UU   

def calc_link_vars_twist(UUA):
    N1, N2 = UUA.shape[0:2]
    UU=np.zeros([N1,N2,2],dtype=complex)
    dirs = [[1,0], [0,1]]
    for i1 in range(N1):
        for i2 in range(N2):
            for idir, vdir in enumerate(dirs):
                    UU[i1,i2,idir] = calc_link_var_twist(UUA[i1                   ,i2               ,:,:],
                                                   UUA[np.mod(i1+vdir[0],N1),np.mod(i2+vdir[1],N2),:,:])
    return UU 

def calc_F12_BZ(UU):
    return np.log( UU[:,:,0]*
                   np.roll(UU[:,:,1],-1,axis=0)*
                   np.conj(np.roll(UU[:,:,0],-1,axis=1)*
                   UU[:,:,1]) )

def calc_F12_BZ_2(UU):
    return np.log( UU[:,:,0]*
                   np.roll(UU[:,:,1],-1,axis=0)/(
                   np.roll(UU[:,:,0],-1,axis=1)*
                   UU[:,:,1] ) )