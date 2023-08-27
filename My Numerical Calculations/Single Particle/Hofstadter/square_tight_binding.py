import numpy as np

def lattice_2d(L_x,L_y):
    #Create 2-D Lattice. Each numbers represents lattice sites
    lattice = np.arange(L_x*L_y).reshape(L_x,L_y)
    
    #Site coordinates on lattice in order, according to coordinate system
    x_co = np.arange(L_y)
    y_co = np.arange(L_x)
    arr = []
    for j in range(len(y_co)):
        for i in range(len(x_co)):
            arr = np.append(arr, [[x_co[i],x_co[j]]])
    xy = arr.reshape((int(arr.size/2),2))
    return lattice, arr, xy

#Find nearest-neighbors for hard-wall B.C.
def HardBC(L_x, L_y):
    lattice, arr, xy = lattice_2d(L_x,L_y)
    neighbors = {}
    for i in range(len(lattice)):
        for j, value in enumerate(lattice[i]):
            if i == 0 or i == len(lattice) - 1 or j == 0 or j == len(lattice[i]) - 1:
                new_neighbors = []
                if i != 0:
                    new_neighbors.append(lattice[i - 1][j])  
                if j != L_y - 1:
                    new_neighbors.append(lattice[i][j + 1]) 
                if i != L_x - 1:
                    new_neighbors.append(lattice[i + 1][j])  
                if j != 0:
                    new_neighbors.append(lattice[i][j - 1])
            else:
                new_neighbors = [
                    lattice[i - 1][j],  
                    lattice[i][j + 1],  
                    lattice[i + 1][j],  
                    lattice[i][j - 1]   
                ]
            neighbors[value] = new_neighbors
    return neighbors

#Find nearest-neighbors for periodic B.C.
def PerBC(L_x, L_y):
    lattice, arr, xy = lattice_2d(L_x,L_y)
    neighbors = {}
    for i in range(len(lattice)):
        for j, value in enumerate(lattice[i]):
            new_neighbors = [
                lattice[(i - 1)%L_x][j%L_y],  
                lattice[i%L_x][(j + 1)%L_y],  
                lattice[(i + 1)%L_x][j%L_y],  
                lattice[i%L_x][(j - 1)%L_y]   
            ]
            neighbors[value] = new_neighbors
    return neighbors
    
#Real Space Hamiltonian with Twisted Angle Phases
def HMat_Theta(L_x, L_y, p, q, theta_x, theta_y):
    lattice, arr, xy = lattice_2d(L_x,L_y)
    PerBCLat = PerBC(L_x, L_y)
    alpha = p/q
    H = np.zeros((L_x*L_y, L_x*L_y), dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
            #NN
            if m in PerBCLat[n]:
                #X ekseninde sınırdan sınıra atlamalar:
                if np.absolute(xy[m][0]-xy[n][0])==L_x-1:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*(alpha)*xy[m][1])*np.exp(-1j*theta_x)
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*(alpha)*xy[m][1])*np.exp(1j*theta_x)

                #Y ekseninde sınırdan sınıra atlamalar
                elif np.absolute(xy[m][1]-xy[n][1])==L_y-1:
                    if xy[m][1] > xy[n][1]:
                        H[m][n] = -np.exp(1j*theta_y)
                    elif xy[m][1] < xy[n][1]:
                        H[m][n] = -np.exp(-1j*theta_y)
                        
                #Lattice içi atlamalar
                else:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*alpha*xy[m][1])
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*alpha*xy[m][1])
                    else:
                        H[m][n] = -np.exp(0)
    return H

# import matplotlib.pyplot as plt

# #Lattice visualization 
# # fig1, ax1 = plt.subplots()
# # ax2 = ax1.twiny()
# # for i in range(L_x):
# #     for j in range(L_y):
# #         plt.plot(i, j, 'ro', markersize=7)
# # ax1.axes.get_xaxis().set_visible(False)        
# # plt.xticks(x_co)
# # plt.yticks(y_co)
# # ax1.invert_yaxis()  
# # ax1.set_ylabel(r'$X$')  
# # ax2.set_xlabel(r'$Y$')
# # plt.title(str(L_x)+r'$\times$'+str(L_y)+' '+'Lattice Structure')
# # ax1.grid()
# # ax2.grid()   

# #YAPILACAK: Komşuluklar yukarıdaki fonksiyonlardan çekilmeli!
# #YAPILACAK: Legend yeri değiştirilmeli!
# #YAPILACAK: Aşağıdaki grafikler her zaman yukarıdaki ile gelmeli.
# #Yani yukarıdaki kodu çalıştırmadan da aşağıdakiler örgü yapılarla gelmeli.
# #YAPILACAK:Hhard-wall  ve periodic B.C. görselleştirmeleri farklı grafikte olmalı.

# #Next-neighbors visualization for periodic B.C.
# def PlotPerBC(choose_x, choose_y):
#     plt.plot(choose_x, choose_y, 'bo', markersize=7, label='Selected Site')
#     plt.plot((choose_x+1)%L_x, choose_y%L_y, 'co', markersize=7, label='Periodic Neighbors')
#     plt.plot((choose_x-1)%L_x, choose_y%L_y, 'co', markersize=7)
#     plt.plot(choose_x%L_x, (choose_y+1)%L_y, 'co', markersize=7)
#     plt.plot(choose_x%L_x, (choose_y-1)%L_y, 'co', markersize=7)
#     plt.legend(loc="upper left")
#     plt.show()
    
# #PlotPerBC(1, 1)
    
# #Next-neighbors visualization for hard-wall B.C.
# def PlotHardBC(choose_x, choose_y):
#     plt.plot(choose_x, choose_y, 'bo', markersize=7, label='Selected Site')
#     if choose_x==0 or choose_x==L_x-1 or choose_y==0 or choose_y==L_y-1:
#         if choose_x!=0:
#             plt.plot(choose_x-1, choose_y, 'co', markersize=7)
#         if choose_y!=L_y-1:
#             plt.plot(choose_x, choose_y+1, 'co', markersize=7)
#         if choose_x!=L_x-1:
#             plt.plot(choose_x+1, choose_y, 'co', markersize=7)
#         if choose_y!=0:
#             plt.plot(choose_x, choose_y-1, 'co', markersize=7)
#     else:
#         plt.plot(choose_x+1, choose_y, 'co', markersize=7)
#         plt.plot(choose_x-1, choose_y, 'co', markersize=7)
#         plt.plot(choose_x, choose_y+1, 'co', markersize=7)
#         plt.plot(choose_x, choose_y-1, 'co', markersize=7)

# #PlotHardBC(1, 1)