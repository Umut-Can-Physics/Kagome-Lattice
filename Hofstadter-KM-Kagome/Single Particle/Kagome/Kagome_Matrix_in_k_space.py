import numpy as np

def Hamiltonian(t1,L1,t2,L2,k1,k2,k3):
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