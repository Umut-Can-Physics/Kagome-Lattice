import sys
sys.path.insert(0, 'Hofstadter')
from square_tight_binding import *

#Periodic hamiltonian matrix
def PerHMat(L_x,L_y,alpha):
    lattice, arr, xy = lattice_2d(L_x,L_y)
    PerBCLat = PerBC(L_x, L_y)
    H = np.zeros((L_x*L_y, L_x*L_y), dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
            if m in PerBCLat[n]:
                if np.absolute(xy[m][0]-xy[n][0])==L_x-1:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*alpha*xy[m][1])
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*alpha*xy[m][1])
                else:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*alpha*xy[m][1])
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*alpha*xy[m][1])
                    else:
                        H[m][n] = -np.exp(0)
    return H

#Periodic Hamiltonian is a hermitian 
# V = PerHMat(1/5)
# VC = (np.conjugate(V)).T
# V - VC