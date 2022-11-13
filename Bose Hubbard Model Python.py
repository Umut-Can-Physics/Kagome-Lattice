#########################################
# BOSE-HUBBARD MODEL FOR SQUARE LATTICE #
#########################################

from tight_binding_approximation import * 
from quspin.operators import hamiltonian 
from quspin.basis import boson_basis_1d, spinless_fermion_basis_1d
import numpy as np 

L=L_x*L_y 
alpha=1/5

print("Site number=",L)
basis = boson_basis_1d(L,Nb=[2]) # Boson basis
#basis_2 = spinless_fermion_basis_1d(L) # Fermion basis
print(basis)
#print(basis_2)

################
# KINETIC PART #
################

# Elements of Single Particle Hamiltonian
hop=[]
for m in range(L_x*L_y):
    for n in range(L_x*L_y):
        # PBC
        if m in PerBCLat[n]:
            if np.absolute(xy[m][0]-xy[n][0])==L_x-1:
                if xy[m][0] > xy[n][0]:
                    hop.append([-np.exp(-1j*2*np.pi*alpha*xy[m][1]),m,n])
                elif xy[m][0] < xy[n][0]:
                    hop.append([-np.exp(1j*2*np.pi*alpha*xy[m][1]),m,n])
            elif m==n:
                hop.append([0,m,n])
            else:
                if xy[m][0] > xy[n][0]:
                    hop.append([-np.exp(1j*2*np.pi*alpha*xy[m][1]),m,n])
                elif xy[m][0] < xy[n][0]:
                    hop.append([-np.exp(-1j*2*np.pi*alpha*xy[m][1]),m,n])
                else:
                    hop.append([-np.exp(0),m,n])
print(hop)

####################
# INTERACTION PART #
####################

U = 2
# Lattice Potantiel Part ( -U/2 * n_i )
pot=[]
for m in range(L):
    pot.append([-U/2,m])
print(pot)

#On-Site Interaction Part ( U/2 * n_i^2 )
interact=[]
for m in range(L):
    for n in range(L):
        if m==n:
            interact.append([U/2,m,n])
print(interact)

#########################
# MANY-BODY HAMILTONIAN #
#########################

# Burada ['+-',hop] terimi tüm yönlerdeki atlama fazlarını içeriyor, yani '-+' ekleme lüzumsuz.
static=[['+-',hop],['n',pot],['nn',interact]]
dynamic=[]
H=hamiltonian(static,dynamic,basis=basis,dtype=np.complex64)
print(H.todense())
E,V=H.eigh()
print("EigenValues:\n",E)
