####################
# Kagome Many-Body #
####################

# Python == Julia 

import matplotlib.pyplot as plt
import numpy as np
from quspin.operators import hamiltonian 
from quspin.basis import boson_basis_1d, spinless_fermion_basis_1d, boson_basis_general

#####################
# KAGOME REAL-SPACE #
####################

sqrt = np.sqrt
a1_vec = np.array([2,0])
a2_vec = np.array([1,sqrt(3)])
l1 = l2 = 2
L = l1*l2;N = 3*L 
b1=np.array([0,0]) 
b2=a1_vec/2
b3=a2_vec/2 
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
NN = [[(0,0,1), (0,0,2), (-1,0,1),(0,-1,2),   (-1,0,2), (-1,1,1), (0,-1,1),(1,-1,2)], 
      [(0,0,-1),(0,0,1), (1,0,-1),(1,-1,1),   (0,-1,1), (1,-1,-1),(0,1,-1),(1,0,1)], 
      [(0,0,-1),(0,0,-2),(0,1,-2),(-1,1,-1),  (-1,1,-2),(-1,0,-1),(0,1,-1),(1,0,-2)] 
      ]

# t1=L1=-1;t2=L2=0
# t1=L1=0;t2=L2=-1
# t1 = -1;L1 = 0.28;t2 = -0.3;L2 = 0.2
t1=t2=-1;L1=L2=0

hopps = [[t1+1j*L1,t1-1j*L1,t1+1j*L1,t1-1j*L1,  t2+1j*L2,t2-1j*L2,t2-1j*L2,t2+1j*L2], 
         [t1-1j*L1,t1+1j*L1,t1-1j*L1,t1+1j*L1,  t2-1j*L2,t2+1j*L2,t2+1j*L2,t2-1j*L2], 
         [t1-1j*L1,t1+1j*L1,t1+1j*L1,t1-1j*L1,  t2-1j*L2,t2+1j*L2,t2+1j*L2,t2-1j*L2]] 
H = np.zeros([N,N],dtype=complex)
for atom_no in range(N):
    atom_site=sites[atom_no]
    for i_delta, delta in enumerate(NN[atom_site[2]]):
        neighbor_site = np.array(atom_site)+np.array(delta)
        neighbor_site[0] = neighbor_site[0]%l1
        neighbor_site[1] = neighbor_site[1]%l2    
        neighbor_no=3*neighbor_site[0]*l2+3*neighbor_site[1]+neighbor_site[2]
        H[neighbor_no,atom_no]=hopps[atom_site[2]][i_delta]      
        
################
# KINETIC TERM #
################
        
Hoppings=[]
for i in range(N):
    for j in range(N):
        Hoppings.append([H[i,j],i,j])
print(Hoppings)

####################
# INTERACTION PART #
####################

U = 2 # U=0 iken tüm parametreler için Julia kodu ile aynı enerjileri üretiyor 
# Lattice Potantiel Part ( -U/2 * n_i )
pot=[]
for m in range(N):
    pot.append([-U/2,m])
print(pot)

#On-Site Interaction Part ( U/2 * n_i^2 )
interact=[]
for m in range(N):
    for n in range(N):
        if m==n:
            interact.append([U/2,m,n])
print(interact)

#########################
# MANY-BODY HAMILTONIAN #
#########################

PN = 2 
print("Site number=",N," | Particle number=",PN)
basis = boson_basis_1d(N,Nb=[PN]) # Boson basis
#basis_2 = spinless_fermion_basis_1d(L) # Fermion basis
print(basis)
static=[['+-',Hoppings],['n',pot],['nn',interact]]
dynamic=[]
H=hamiltonian(static,dynamic,basis=basis,dtype=np.complex64)
print(H.todense())
E,V=H.eigh()
print(np.sort(E,axis=None))
