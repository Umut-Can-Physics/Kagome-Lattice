import sys
sys.path.insert(0, 'C:/Users/Can/Dropbox/PC (2)/Desktop/My Numerical Calculations/Single Particle/Hofstadter')
from square_tight_binding import *

# hard b.c. hamiltonian matrix
def HardHMat(alfa):
    H = np.zeros((L_x*L_y, L_x*L_y), dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
            if m in HardBCLat[n]:
                if xy[m][0] > xy[n][0]:
                    H[m][n] = -np.exp(1j*2*np.pi*alfa*xy[m][1])
                elif xy[m][0] < xy[n][0]:
                    H[m][n] = -np.exp(-1j*2*np.pi*alfa*xy[m][1])
                else:
                    H[m][n]=-1
    return H

# hard hamiltonian matrix is a hermitian 
# M = HardHMat(1/5)
# MC = (np.conjugate(M)).T
# M - MC

