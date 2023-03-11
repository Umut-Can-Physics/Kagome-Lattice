import numpy as np

# Primitive (Bravais) vectors
a1=np.array([2,0])
a2=np.array([1,np.sqrt(3)])

# Basis Vectors (Hopping vectors)
a1_b=a1/2
a2_b=a2/2

# Reciprocal vectors
b1x = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * a2[1]
b1y = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * -a2[0]
b2x = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * -a1[1]
b2y = (2*np.pi/(a1[0]*a2[1]-a1[1]*a2[0])) * a1[0]
b1=np.array([b1x,b1y]) #(np.pi,-np.pi/np.sqrt(3))
b2=np.array([b2x,b2y]) #(0,2*np.pi/np.sqrt(3))

# t1=-1;L1=0;t2=0;L2=0 
# t1=1;L1=0.6;t2=0.3;L2=0 #NNN
# t1=-1;L1=1;t2=0;L2=0 #NN (Well-Defined Berry-Curvature)
t1=-1;L1=0.28;t2=0.3;L2=0.2 #NNN (Well-Defined Berry-Curvature)
def Hamiltonian(k1,k2,k3):
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
print("t1=",t1,"L1=",L1,"t2=",t2,"L2=",L2)