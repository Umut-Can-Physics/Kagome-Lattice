#################################
# MANY-BODY CHERN NUMBER KAGOME #
#################################

# Problem: Hamiltonian Fonksiyonunda: ValueError: values in indx falls outside of system

import matplotlib.pyplot as plt
import numpy as np
from quspin.operators import hamiltonian 
from quspin.basis import boson_basis_1d 

#####################
# KAGOME REAL-SPACE #
####################

sqrt = np.sqrt
#Bravais vectors
a1_vec = np.array([2,0])
a2_vec = np.array([1,sqrt(3)])
# number of unit-cell in each axes
l1 = l2 = 2
L = l1*l2 # number of unit-cells
N = 3*L # number of sites
#Basis vectors
b1=np.array([0,0]) #A atoms
b2=a1_vec/2 #B atoms
b3=a2_vec/2 #C atoms
Basis = [b1,b2,b3]
basis_colors=['red','blue','green']
coordinates = []
sites = []
atom_dict={}
for i1 in range(l1):
    for i2 in range(l2):
        Lp = a1_vec * i1 + a2_vec * i2
        center = a1_vec * i1 + a2_vec * i2 + b1
        shift = (-b2-b3)/sqrt(3)/2
        P1=np.array([center+shift,center+a1_vec+shift,center+a2_vec+shift])
        for ib, b in enumerate(Basis):
            atom_vec = Lp + b
            atom_no = 3 * i1 * l2 + 3 * i2 + ib
            site = [i1,i2,ib]
            sites.append(site)
            coordinates.append(atom_vec)
            atom_dict[tuple(site)]=atom_vec
NN = [[(0,0,1), (0,0,2), (-1,0,1),(0,-1,2),   (-1,0,2), (-1,1,1), (0,-1,1),(1,-1,2)], #A (red)
      [(0,0,-1),(0,0,1), (1,0,-1),(1,-1,1),   (0,-1,1), (1,-1,-1),(0,1,-1),(1,0,1)], #B (blue)
      [(0,0,-1),(0,0,-2),(0,1,-2),(-1,1,-1),  (-1,1,-2),(-1,0,-1),(0,1,-1),(1,0,-2)] #C (green)
      ]
# t1=L1=-1;t2=L2=0
# t1=L1=0;t2=L2=-1
t1 = -1;L1 = 0.28;t2 = -0.3;L2 = 0.2
# t1=t2=-1;L1=L2=0
hopps = [[t1+1j*L1,t1-1j*L1,t1+1j*L1,t1-1j*L1,  t2+1j*L2,t2-1j*L2,t2-1j*L2,t2+1j*L2], # from A
         [t1-1j*L1,t1+1j*L1,t1-1j*L1,t1+1j*L1,  t2-1j*L2,t2+1j*L2,t2+1j*L2,t2-1j*L2], # from B
         [t1-1j*L1,t1+1j*L1,t1+1j*L1,t1-1j*L1,  t2-1j*L2,t2+1j*L2,t2+1j*L2,t2-1j*L2]] # from C
H = np.zeros([N,N],dtype=complex)
for atom_no in range(N):
    atom_site=sites[atom_no]
    for i_delta, delta in enumerate(NN[atom_site[2]]):
        neighbor_site = np.array(atom_site)+np.array(delta)
        neighbor_site[0] = neighbor_site[0]%l1
        neighbor_site[1] = neighbor_site[1]%l2    
        neighbor_no=3*neighbor_site[0]*l2+3*neighbor_site[1]+neighbor_site[2]
        H[neighbor_no,atom_no]=hopps[atom_site[2]][i_delta]  
        
###############################
# DISCERETE TWIST ANGLE SPACE #
###############################

a1=np.array([2,0])
a2=np.array([1,np.sqrt(3)])
b1x = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * a2[1]
b1y = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * -a2[0]
b2x = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * -a1[1]
b2y = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * a1[0]
b1=np.array([b1x,b1y]) 
b2=np.array([b2x,b2y]) 
theta_size = 5
N1 = N2 = theta_size
q1_list = np.arange(0,N1) / N1
q2_list =  np.arange(0,N2) / N2
Q = []
for q1 in q1_list:
    for q2 in q2_list:
        Q.append(q1*b1+q2*b2)
Q = np.array(Q)
theta_1 = Q[:,0];theta_2=Q[:,1]
d1 = 0.62831853
d2 = 0.72551975

#########################
# MANY_BODY HAMILTONIAN #
#########################

U = 2 
# Lattice Potantiel Part ( -U/2 * n_i )
pot=[]
for m in range(N):
    pot.append([-U/2,m])
#On-Site Interaction Part ( U/2 * n_i^2 )
interact=[]
for m in range(N):
    for n in range(N):
        if m==n:
            interact.append([U/2,m,n])

basis = boson_basis_1d(L,Nb=[2]) 

def Hamiltonian(theta_1, theta_2):
    hop = []
    for atom_no in range(N):
        atom_site=sites[atom_no]
        for i_delta, delta in enumerate(NN[atom_site[2]]):
            neighbor_site = np.array(atom_site)+np.array(delta)
            neighbor_site[0] = neighbor_site[0]%l1
            neighbor_site[1] = neighbor_site[1]%l2    
            neighbor_no=3*neighbor_site[0]*l2+3*neighbor_site[1]+neighbor_site[2]
            if atom_site[0]+(l1-1)==neighbor_site[0]:
                twist_1=np.exp(1j*theta_1)
            elif atom_site[0]-(l1-1)==neighbor_site[0]:
                twist_1=np.exp(-1j*theta_1)  
            else:
                twist_1=1
            if atom_site[1]+(l2-1)==neighbor_site[1]:
                twist_2=np.exp(1j*theta_2)
            elif atom_site[1]-(l2-1)==neighbor_site[1]:
                twist_2=np.exp(-1j*theta_2)
            else:
                twist_2=1
            # Single-Particle Hamiltonian Elements
            hop.append([ twist_1*twist_2*hopps[atom_site[2]][i_delta],neighbor_no,atom_no ])
    # static=[['+-',H],['n',pot],['nn',interact]]
    static=[['+-',hop]]
    dynamic=[]
    # MB Hamiltonian
    H=hamiltonian(static,dynamic,basis=basis,dtype=np.complex64,check_symm=False,check_herm=False)
    return H

#############################
# CALCULATING CHERN NUMBERS #
#############################

# chern_array = []
# q = 3
# #0-4,4-8,8-12 (Slicing) (l1=l2=2)
# n1=0;n2=4
# S=0
# for t_1 in range(0, len(theta_1)):
#     for t_2 in range(0, len(theta_2)):
#         # Hamiltonyenin ilk 9 (l1=l2=3) enerjilerinin dejenere çıkması lazım. Hamiltonian yanlış !
#         #theta_1[0]=theta_2[0]=0 iken bile w1 D-Fold dejenere olmuyor.
#         w1, v1 = np.linalg.eig(Hamiltonian(a1_vec, a2_vec, l1, l2, N, Basis, sites, NN, hopps, theta_1[t_1], theta_2[t_2], H))
#         idx1 = np.argsort(w1)
#         v1_sorted = v1[:,idx1]
#         v11 = v1_sorted[:,n1:n2]
#         w2, v2 = np.linalg.eig(Hamiltonian(a1_vec, a2_vec, l1, l2, N, Basis, sites, NN, hopps, theta_1[t_1]+d1, theta_2[t_2], H))
#         idx2 = np.argsort(w2)
#         v2_sorted = v2[:,idx2]
#         v22 = v2_sorted[:,n1:n2]
#         w3, v3 = np.linalg.eig(Hamiltonian(a1_vec, a2_vec, l1, l2, N, Basis, sites, NN, hopps, theta_1[t_1], theta_2[t_2]+d2, H))
#         idx3 = np.argsort(w3)
#         v3_sorted = v3[:,idx3]
#         v33 = v3_sorted[:,n1:n2]
#         w4, v4 = np.linalg.eig(Hamiltonian(a1_vec, a2_vec, l1, l2, N, Basis, sites, NN, hopps, theta_1[t_1]+d1, theta_2[t_2]+d2, H))
#         idx4 = np.argsort(w4)
#         v4_sorted = v4[:,idx4]
#         v44 = v4_sorted[:,n1:n2]
#         # There is no D-Fold Degeneracy ! 
#         # Diklik Bağıntısı Sağlanıyor !
#         U1 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v11)), v22)) #U_x
#         U1 = U1 / np.absolute(U1)
#         U2 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v22)), v44)) #U_y dx
#         U2 = U2 / np.absolute(U2)
#         U3 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v33)), v44)) #U_x dy
#         U3 = U3 / np.absolute(U3)
#         U4 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v11)), v33)) #U_y
#         U4 = U4 / np.absolute(U4)
#         F = np.log(U1*U2*1/U3*1/U4)
#         S = S+F

# C = 1/(2*np.pi*1j)*S
# chern_array.append(C.real)
# print(C)
