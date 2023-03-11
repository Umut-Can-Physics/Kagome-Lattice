##########################################################################
# SINGLE PARTICLE NN and NNN REAL SPACE PERIODIC B.C. KAGOME MODEL (v.2) #
##########################################################################

import matplotlib.pyplot as plt
import numpy as np

#Bravais vectors
a1_vec = np.array([2,0])
a2_vec = np.array([1,np.sqrt(3)])

l1 = l2 = 3 # number of unit-cell in each axes
L = l1*l2 # number of unit-cells
N = 3*L # number of sites

#Basis vectors
b1=np.array([0,0]) #A atoms
b2=a1_vec/2 #B atoms
b3=a2_vec/2 #C atoms
Basis = [b1,b2,b3]
basis_colors=['red','blue','green']

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
            site = [i1,i2,ib]
            sites.append(site)
            coordinates.append(atom_vec)
            atom_dict[tuple(site)]=atom_vec
            ax.plot(atom_vec[0], atom_vec[1], marker='o', color=basis_colors[ib], alpha=1)
            ax.annotate(str(atom_no), atom_vec)
ax.set_title("Real-Space Kagome \n"+str(l1)+r'$\times$'+str(l2))
ax.set_xlabel(r'$n_1$')
ax.set_ylabel(r'$n_2$')           
# i1=i2=1
# center = i1*a1_vec + i2*a2_vec + b1
# circle_1 = plt.Circle(center, np.linalg.norm(b2-b1), fill=False, color='blue')
# ax.add_patch(circle_1)

# circle_2 = plt.Circle(center, np.linalg.norm(b3-a1_vec-b1), fill=False, color='red')
# ax.add_patch(circle_2)

#NN and NNN Hopping Sites
NN = [[(0,0,1), (0,0,2), (-1,0,1),(0,-1,2),   (-1,0,2), (-1,1,1), (0,-1,1),(1,-1,2)], #A (red)
      [(0,0,-1),(0,0,1), (1,0,-1),(1,-1,1),   (0,-1,1), (1,-1,-1),(0,1,-1),(1,0,1)], #B (blue)
      [(0,0,-1),(0,0,-2),(0,1,-2),(-1,1,-1),  (-1,1,-2),(-1,0,-1),(0,1,-1),(1,0,-2)] #C (green)
      ]

# t1 = -1
# L1 = 0
# t2 = 0
# L2 = 0
t1 = -1
L1 = 0.28
t2 = 0.3
L2 = 0.2
hopps = [[t1+1j*L1,t1-1j*L1,t1+1j*L1,t1-1j*L1,  t2+1j*L2,t2-1j*L2,t2-1j*L2,t2+1j*L2], # from A
         [t1-1j*L1,t1+1j*L1,t1-1j*L1,t1+1j*L1,  t2-1j*L2,t2+1j*L2,t2+1j*L2,t2-1j*L2], # from B
         [t1-1j*L1,t1+1j*L1,t1+1j*L1,t1-1j*L1,  t2-1j*L2,t2+1j*L2,t2+1j*L2,t2-1j*L2]] # from C

H = np.zeros([N,N],dtype=complex)

# atom_index=18
# atom_site=sites[atom_index]
# # atom_vec=coordinates[atom_index]
# print(atom_index,". Atom Site:",atom_site)
# for i_delta, delta in enumerate(NN[atom_site[2]]):
#     neighbor_site = np.array(atom_site)+np.array(delta)
#     neighbor_site[0] = neighbor_site[0]%l1
#     neighbor_site[1] = neighbor_site[1]%l2    
#     neighbor_no=3*neighbor_site[0]*l2+3*neighbor_site[1]+neighbor_site[2]
#     # neighbor_vec=a1_vec*neighbor_site[0]+a2_vec*neighbor_site[1]+Basis[neighbor_site[2]]
#     print("Delta:",delta, ", Neighbor Site",neighbor_site, ", Neighbor No:",neighbor_no)

for atom_no in range(N):
    atom_site=sites[atom_no]
    # atom_vec=coordinates[atom_no]
    for i_delta, delta in enumerate(NN[atom_site[2]]):
        neighbor_site = np.array(atom_site)+np.array(delta)
        neighbor_site[0] = neighbor_site[0]%l1
        neighbor_site[1] = neighbor_site[1]%l2    
        neighbor_no=3*neighbor_site[0]*l2+3*neighbor_site[1]+neighbor_site[2]
        neighbor_vec=a1_vec*neighbor_site[0]+a2_vec*neighbor_site[1]+Basis[neighbor_site[2]]
        H[neighbor_no,atom_no]=hopps[atom_site[2]][i_delta]