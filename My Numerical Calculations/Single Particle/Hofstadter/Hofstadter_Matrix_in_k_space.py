import numpy as np

# Initial value of q
q=5

#Hamiltonian Matrix for Hofstadter Model (doesn't include alpha=1/2 matrix)
def HMatrix(alfa, k_x, k_y):
    M = np.zeros((q,q), dtype=complex) 
    for i in range (0, q):
        M[i,i]=2*np.cos(k_y-2*np.pi*alfa*i) 
        if i==q-1: 
            M[i,i-1]=1
        elif i==0: 
            M[i,i+1]=1
        else: 
            M[i,i-1]=1
            M[i,i+1]=1
        M[0,q-1]=np.exp(-q*1.j*k_x)
        M[q-1,0]=np.exp(q*1.j*k_x)
    return M

#Special Hofstadter Matrix Just alpha=1/2
s = 2
alpha = 1/s
def HMatrix2(alpha, k_x, k_y):
    M2 = np.zeros((s,s), dtype=complex) 
    for i in range (0, s):
        M2[i,i]=2*np.cos(k_y-2*np.pi*alpha*i) 
        if i==q-1: 
            M2[i,i-1]=1
        M2[0,s-1]=1+np.exp(-s*1.j*k_x)
        M2[s-1,0]=1+np.exp(s*1.j*k_x)
    return M2

#Set Rational to Alpha Values 
def gcd(a, b): 
    if b == 0: return a
    return gcd(b, a % b)