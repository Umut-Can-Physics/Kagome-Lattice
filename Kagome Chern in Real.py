#################################################
# Kagome Chern Number Calculation in Real Space #
#################################################

import matplotlib.pyplot as plt
import numpy as np

#Bravais vectors
a1_vec = np.array([2,0])
a2_vec = np.array([1,np.sqrt(3)])

# Inputs of Lattice
l1 = l2 = 3 # number of unit-cell in each axes
L = l1*l2 # number of unit-cells
N = 3*L # number of sites

#Basis vectors
b1=np.array([0,0]) #A atoms
b2=a1_vec/2 #B atoms
b3=a2_vec/2 #C atoms
Basis = [b1,b2,b3]
basis_colors=['red','blue','green']

# Plot Kagome Latttice with Unit-Cells
fig, ax = plt.subplots()
coordinates = []
sites = []
atom_dict={}
for i1 in range(l1):
    for i2 in range(l2):
        Lp = a1_vec * i1 + a2_vec * i2
        ax.plot(Lp[0], Lp[1], marker='o', color='black', alpha=0.2)
        center = a1_vec * i1 + a2_vec * i2 + b1
        shift = (-b2-b3)/np.sqrt(3)/2
        P1=np.array([center+shift,center+a1_vec+shift,center+a2_vec+shift])
        T=plt.Polygon(P1, fill=False)
        ax.add_patch(T)        
        for ib, b in enumerate(Basis):
            atom_vec = Lp + b
            atom_no = 3 * i1 * l2 + 3 * i2 + ib
            #i1,i2 unit-cell index and ib kind of atom (A,B or C)
            site = [i1,i2,ib]
            sites.append(site)
            coordinates.append(atom_vec)
            atom_dict[tuple(site)]=atom_vec
            ax.plot(atom_vec[0], atom_vec[1], marker='o', color=basis_colors[ib], alpha=1)
            ax.annotate(str(atom_no), atom_vec)

# NN and NNN Hopping Sites Matrix          
NN = [[(0,0,1), (0,0,2), (-1,0,1),(0,-1,2),   (-1,0,2), (-1,1,1), (0,-1,1),(1,-1,2)], #A (red)
      [(0,0,-1),(0,0,1), (1,0,-1),(1,-1,1),   (0,-1,1), (1,-1,-1),(0,1,-1),(1,0,1)], #B (blue)
      [(0,0,-1),(0,0,-2),(0,1,-2),(-1,1,-1),  (-1,1,-2),(-1,0,-1),(0,1,-1),(1,0,-2)] #C (green)
      ]

#Hopping Parameters
# t1 = -1;L1 = 0;t2 = 0;L2 = 0
t1 = -1;L1 = 0.28;t2 = -0.3;L2 = 0.2

# NN and NNN Hopping Phases Matrix
hopps = [[t1+1j*L1,t1-1j*L1,t1+1j*L1,t1-1j*L1,  t2+1j*L2,t2-1j*L2,t2-1j*L2,t2+1j*L2], # from A
         [t1-1j*L1,t1+1j*L1,t1-1j*L1,t1+1j*L1,  t2-1j*L2,t2+1j*L2,t2+1j*L2,t2-1j*L2], # from B
         [t1-1j*L1,t1+1j*L1,t1+1j*L1,t1-1j*L1,  t2-1j*L2,t2+1j*L2,t2+1j*L2,t2-1j*L2]] # from C

# DISCERETE TWIST ANGLE SPACE FOR CHERN NUMBERS (KAGOME)
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

# DISCERETE TWIST ANGLE SPACE FOR CHERN NUMBERS (SQUARE)
# theta_size=10
# d1 = 2*np.pi/theta_size
# d2 = d1
# theta_1 = np.arange(0, 2*np.pi, d1)
# theta_2 = np.arange(0, 2*np.pi, d2)

# atom_index=14
# atom_site=sites[atom_index]
# atom_vec=coordinates[atom_index]
# print(atom_index,". Atom Site:",atom_site)
# for i_delta, delta in enumerate(NN[atom_site[2]]):
#     neighbor_site = np.array(atom_site)+np.array(delta)
#     neighbor_site[0] = neighbor_site[0]%l1
#     neighbor_site[1] = neighbor_site[1]%l2    
#     neighbor_no=3*neighbor_site[0]*l2+3*neighbor_site[1]+neighbor_site[2]
#     neighbor_vec=a1_vec*neighbor_site[0]+a2_vec*neighbor_site[1]+Basis[neighbor_site[2]]
#     # Tüm periyodik komşu atlamaları
#     #print("Delta:",delta, ", Neighbor Site",neighbor_site, ", Neighbor No:",neighbor_no)
#     # Sınırdan sınıra olan atlamalar
#     # i1 ve i2 eksenine (l1-1) eklenmiş yada çıkartılmış tüm siteler sınırdan sınıra olanlardır
#     # Yönleri ise eklediğime yada çıkarttığıma göre karar veriyorum
#     if atom_site[0]+(l1-1)==neighbor_site[0] or atom_site[0]-(l1-1)==neighbor_site[0] or atom_site[1]+(l2-1)==neighbor_site[1] or atom_site[1]-(l2-1)==neighbor_site[1]: 
#         print("Delta:",delta, ", Neighbor Site",neighbor_site, ", Neighbor No:",neighbor_no)

# CONSTRUCT HAMILTONIAN WITH TWIST ANGLE PHASES
H = np.zeros([N,N],dtype=complex)
def Hamiltonian(a1_vec, a2_vec, l1, l2, N, Basis, sites, NN, hopps, theta_1, theta_2, H):
    for atom_no in range(N):
        atom_site=sites[atom_no]
        for i_delta, delta in enumerate(NN[atom_site[2]]):
            neighbor_site = np.array(atom_site)+np.array(delta)
            neighbor_site[0] = neighbor_site[0]%l1
            neighbor_site[1] = neighbor_site[1]%l2    
            neighbor_no=3*neighbor_site[0]*l2+3*neighbor_site[1]+neighbor_site[2]
            # Long Range Hopping (Edge to Edge) Condition
            #Yön tanımı: i1 ekseni: sağa doğru + sola doğru -, i2 ekseni: yukarı doğru +, aşağı doğru -
            if atom_site[0]+(l1-1)==neighbor_site[0]:
                twist_1=np.exp(1j*theta_1)
            else:
                twist_1=1
            if atom_site[0]-(l1-1)==neighbor_site[0]:
                twist_1=np.exp(-1j*theta_1)  
            else:
                twist_1=1
            if atom_site[1]+(l2-1)==neighbor_site[1]:
                twist_2=np.exp(1j*theta_2)
            else: twist_2=1
            if atom_site[1]-(l2-1)==neighbor_site[1]:
                twist_2=np.exp(-1j*theta_2)
            else:
                twist_2=1
            H[neighbor_no,atom_no]=twist_1*twist_2*hopps[atom_site[2]][i_delta]
    return H

# CALCULATING CHERN NUMBERS
chern_array = []
q = 3
#0-9,9-18,9-27 (Slicing) (l1=l2=3)
n1=0;n2=9
S=0
for t_1 in range(0, len(theta_1)):
    for t_2 in range(0, len(theta_2)):
        w1, v1 = np.linalg.eig(Hamiltonian(a1_vec, a2_vec, l1, l2, N, Basis, sites, NN, hopps, theta_1[t_1], theta_2[t_2], H))
        idx1 = np.argsort(w1)
        v1_sorted = v1[:,idx1]
        v11 = v1_sorted[:,n1:n2]
        w2, v2 = np.linalg.eig(Hamiltonian(a1_vec, a2_vec, l1, l2, N, Basis, sites, NN, hopps, theta_1[t_1]+d1, theta_2[t_2], H))
        idx2 = np.argsort(w2)
        v2_sorted = v2[:,idx2]
        v22 = v2_sorted[:,n1:n2]
        w3, v3 = np.linalg.eig(Hamiltonian(a1_vec, a2_vec, l1, l2, N, Basis, sites, NN, hopps, theta_1[t_1], theta_2[t_2]+d2, H))
        idx3 = np.argsort(w3)
        v3_sorted = v3[:,idx3]
        v33 = v3_sorted[:,n1:n2]
        w4, v4 = np.linalg.eig(Hamiltonian(a1_vec, a2_vec, l1, l2, N, Basis, sites, NN, hopps, theta_1[t_1]+d1, theta_2[t_2]+d2, H))
        idx4 = np.argsort(w4)
        v4_sorted = v4[:,idx4]
        v44 = v4_sorted[:,n1:n2]
        U1 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v11)), v22)) #U_x
        U1 = U1 / np.absolute(U1)
        U2 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v22)), v44)) #U_y dx
        U2 = U2 / np.absolute(U2)
        U3 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v33)), v44)) #U_x dy
        U3 = U3 / np.absolute(U3)
        U4 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v11)), v33)) #U_y
        U4 = U4 / np.absolute(U4)
        F = np.log(U1*U2*1/U3*1/U4)
        S = S+F

C = 1/(2*np.pi*1j)*S
chern_array.append(C.real)
print(C)
