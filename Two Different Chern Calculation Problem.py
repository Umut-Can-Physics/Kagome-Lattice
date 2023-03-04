import numpy as np

#https://journals.jps.jp/doi/10.1143/JPSJ.74.1674

# Bravais Vectors
a1=np.array([2,0])
a2=np.array([1,np.sqrt(3)])
# Basis Vectors (Hopping vectors)
a1_b=a1/2
a2_b=a2/2
# Reciprocal Vectors
b1x = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * a2[1]
b1y = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * -a2[0]
b2x = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * -a1[1]
b2y = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * a1[0]
b1=np.array([b1x,b1y]) 
b2=np.array([b2x,b2y]) 

# Brilliun Zone Define and Allowed Pairs
N1 = N2 = 25
q1_list = np.arange(0,N1) / N1 
q2_list =  np.arange(0,N2) / N2 

# Hopping Parameters
# t1=-1;L1=0;t2=0;L2=0 # NN [1,0,-1]
t1=0;L1=0;t2=-1;L2=-1 # Complex NNN [3,0,-3] 
# Following results of parameters consistent with Review article
# t1=-1;L1=1;t2=0;L2=0 # Complex NN (Well-Defined Berry Curvature) [-1,0,1]
# t1=-1;L1=0.28;t2=0.3;L2=0.2 # Complex NNN (Well-Defined Berry Curvature) [-1,0,1]
def H_k(k_vec,**kwargs):
    kx = k_vec[0]
    ky = k_vec[1]
    k1=np.dot([kx,ky],a1_b);k2=np.dot([kx,ky],a2_b);k3=k2-k1 # Convention
    H = 2*t1*np.array([
    [0, np.cos(k1), np.cos(k2)],
    [np.cos(k1), 0, np.cos(k3)],        
    [np.cos(k2), np.cos(k3), 0]
    ])+2*1j*L1*np.array([
    [0, np.cos(k1), -np.cos(k2)],
    [-np.cos(k1), 0, np.cos(k3)],        
    [np.cos(k2), -np.cos(k3), 0]
    ])+2*t2*np.array([
    [0, np.cos(k2+k3), np.cos(k3-k1)],
    [np.cos(k2+k3), 0, np.cos(k1+k2)],
    [np.cos(k3-k1), np.cos(k1+k2), 0]
    ])+2*1j*L2*np.array([
    [0, -np.cos(k2+k3), np.cos(k3-k1)],
    [np.cos(k2+k3), 0, -np.cos(k1+k2)],
    [-np.cos(k3-k1), np.cos(k1+k2), 0]
    ])
    return H

def build_U(vec1, vec2):
    inner_product = np.dot(vec1, vec2.conj())
    U = inner_product / np.abs(inner_product)
    return U

def latF(k_vec, Dk):

    k =  np.array([k_vec[0], k_vec[1]])
    E, aux = np.linalg.eig(H_k(k))
    idx = E.argsort()
    psi = aux[:,idx] # column vectors

    k = np.array([k_vec[0]+Dk[0], k_vec[1]], float)
    E, aux = np.linalg.eig(H_k(k))
    idx = E.argsort()
    psiDx = aux[:,idx]

    k = np.array([k_vec[0], k_vec[1]+Dk[1]], float)
    E, aux = np.linalg.eig(H_k(k))
    idx = E.argsort()
    psiDy = aux[:,idx]

    k = np.array([k_vec[0]+Dk[0], k_vec[1]+Dk[1]], float)
    E, aux = np.linalg.eig(H_k(k))
    idx = E.argsort()
    psiDxDy = aux[:,idx]

    size=3
    U1x = np.zeros(size, dtype=complex)
    U2y = np.zeros(size, dtype=complex)
    U1y = np.zeros(size, dtype=complex)
    U2x = np.zeros(size, dtype=complex)

    for i in range(size):
        
        U1x[i] = build_U(psi[:,i], psiDx[:,i])
        
        U2y[i] = build_U(psi[:,i], psiDy[:,i])
        
        U1y[i] = build_U(psiDy[:,i], psiDxDy[:,i])
        
        U2x[i] = build_U(psiDx[:,i], psiDxDy[:,i])
        
    F12 = np.zeros(size, dtype=complex)
    F12 = np.log( U1x*U2x*1/U1y*1/U2y)
    # column vector
    return F12

# IMPORTANT: Doesn't important that whether it confirms periodic-k condition.
# In any case, it produces consistent Chern numbers.

#Dx=Delta k_x, Dy=Delta k_y (N1=N2=25)
Dx = 0.1256637061435919;Dy = 0.14510394913873714
Dk = np.array([Dx,Dy], float)

size=3
LF = np.zeros(size, dtype=complex)
Sum = np.zeros(size, dtype=complex)
Chern = np.zeros(size, dtype=complex)
Q = np.zeros([N1,N2,2])
F21=np.zeros([N1,N2],dtype=complex) 
for iq1, q1 in enumerate(q1_list):
    for iq2, q2 in enumerate(q2_list):
        Q[iq1,iq2,:] = q1*b1+q2*b2
        kx,ky=Q[iq1,iq2,:]
        k_vec = np.array([kx,ky], float)
        LF = latF(k_vec, Dk)
        Sum += LF
        F21[iq1,iq2] = LF[0]
Chern = Sum.imag/(2*np.pi)
print(Chern)

from topo_lib import hatsugai, hamiltonians

# Eigenvalues and Eigenenergies of Matrix
EEA, UUA = hatsugai.get_all_eigen(H_k, Q[:,0,0],Q[0,:,1],3,q=3)

for bi in range(3):
    bj=bi

    U1, U2 = hatsugai.get_U1_U2(bi,bj,UUA)

    F12 = hatsugai.get_F12(U1,U2)
    print(np.sum(F12.imag)/(2*np.pi))
