import numpy as np

#Lattice size
L_x=L_y=5
print("lx=",L_x,", ly=",L_y)

#Create 2-D Lattice. Each numbers represents lattice sites
lattice = np.arange(L_x*L_y).reshape(L_x,L_y)

#Site coordinates on lattice in order, according to coordinate system
x_co = np.arange(L_x)
y_co = np.arange(L_y)
arr = np.empty(shape=[0,2])
for j in range(len(x_co)):
    for i in range(len(y_co)):
        arr = np.append(arr, [[x_co[i],x_co[j]]], axis=0)
xy = arr

#Find nearest-neighbors for hard-wall B.C.
def HardBC(arr):
    neighbors = {}
    for i in range(len(arr)):
        for j, value in enumerate(arr[i]):
            if i == 0 or i == len(arr) - 1 or j == 0 or j == len(arr[i]) - 1:
                new_neighbors = []
                if i != 0:
                    new_neighbors.append(arr[i - 1][j])  
                if j != len(arr[i]) - 1:
                    new_neighbors.append(arr[i][j + 1]) 
                if i != len(arr) - 1:
                    new_neighbors.append(arr[i + 1][j])  
                if j != 0:
                    new_neighbors.append(arr[i][j - 1])
            else:
                new_neighbors = [
                    arr[i - 1][j],  
                    arr[i][j + 1],  
                    arr[i + 1][j],  
                    arr[i][j - 1]   
                ]
            neighbors[value] = new_neighbors
    return neighbors

#Find nearest-neighbors for periodic B.C.
def PerBC(arr):
    neighbors = {}
    for i in range(len(arr)):
        for j, value in enumerate(arr[i]):
            new_neighbors = [
                arr[(i - 1)%L_x][j%L_y],  
                arr[i%L_x][(j + 1)%L_y],  
                arr[(i + 1)%L_x][j%L_y],  
                arr[i%L_x][(j - 1)%L_y]   
            ]
            neighbors[value] = new_neighbors
    return neighbors

#Show nearest-neighbors
HardBCLat = HardBC(lattice)
PerBCLat = PerBC(lattice)

# import matplotlib.pyplot as plt

# #YAPILACAK: Site numaraları gösterilmeli!
# #YAPILACAK: Görselleştirmede L_x, L_y yerleri farklı!

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