import numpy as np

def H(l1,l2,t1,L1,t2,L2,theta_1, theta_2):
    #Bravais vectors
    a1_vec = np.array([2,0])
    a2_vec = np.array([1,np.sqrt(3)])
    
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
            shift = (-b2-b3)/np.sqrt(3)/2       
            for ib, b in enumerate(Basis):
                atom_vec = Lp + b
                atom_no = 3 * i1 * l2 + 3 * i2 + ib
                site = [i1,i2,ib]
                sites.append(site)
                coordinates.append(atom_vec)
                atom_dict[tuple(site)]=atom_vec      
    
    #NN and NNN Hopping Sites
    NN = [[(0,0,1), (0,0,2), (-1,0,1),(0,-1,2),   (-1,0,2), (-1,1,1), (0,-1,1),(1,-1,2)], #A (red)
          [(0,0,-1),(0,0,1), (1,0,-1),(1,-1,1),   (0,-1,1), (1,-1,-1),(0,1,-1),(1,0,1)], #B (blue)
          [(0,0,-1),(0,0,-2),(0,1,-2),(-1,1,-1),  (-1,1,-2),(-1,0,-1),(0,1,-1),(1,0,-2)] #C (green)
          ]

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
            Bwrap = [False, False]
            if abs( neighbor_site[0] - atom_site[0] ) == l1-1 : Bwrap[0] = True
            if abs( neighbor_site[1] - atom_site[1] ) == l2-1 : Bwrap[1] = True
            twist_1 = np.exp(1j * theta_1 / (l1-2) * (abs( neighbor_site[0] - atom_site[0] )-1 ) * np.sign(neighbor_site[0] - atom_site[0]) )
            twist_2 = np.exp(1j * theta_2 / (l2-2) * (abs( neighbor_site[1] - atom_site[1] )-1 ) * np.sign(neighbor_site[1] - atom_site[1]) )
            H[neighbor_no,atom_no] = twist_1 * twist_2 * hopps[atom_site[2]][i_delta]
    return H