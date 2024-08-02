# function square_lattice(Nx,Ny)    
#     sites = collect(1:Nx*Ny)
#     coordinates = []
#     for iy in 1:Ny
#       for ix in 1:Nx
#             push!(coordinates, (ix,iy))
#         end
#     end
#     return sites, coordinates
# end

# Nx = 3
# Ny = 3
# sites, coords = square_lattice(Nx,Ny)

# reshape(sites,(Nx,Ny))|>transpose

# coords

# # executeme

# function get_square_site_no(coord,Nx)
#     return coord[1]+(coord[2]-1)*Nx
# end

# get_square_site_no(coords[1], Nx)

# function square_neighbors(sites, coords, pbc=false)
#     Nx, Ny = coords[end]
#     dirs = [ [1, 0], [-1,0], [0,1], [0,-1] ]
#     neighbors=[]
#     for site in sites
#         neighbor=[]
#         #println(site, coords[site])
#         for dir in dirs
#             neighborx = coords[site][1]+dir[1] 
#             neighbory = coords[site][2]+dir[2]
#             #print(neighborx,' ',neighbory,' ',)
#             if pbc==false && ( neighborx == 0 || neighborx == Nx+1 )
#                 continue
#             elseif pbc==false && ( neighbory == 0 || neighbory == Ny+1 )
#                 continue
#             else
#                 neighborx = mod(neighborx-1,Nx)+1
#                 neighbory = mod(neighbory-1,Ny)+1
#                 push!(neighbor,[neighborx,neighbory])
#                 #println(neighborx,' ',neighbory,' ',)
#                 neighbor_site = get_square_site_no([neighborx,neighbory],Nx)
#                 #println(neighbor_site)
#             end
#         end
#         push!(neighbors,neighbor)
#     end
#     return neighbors
# end

# pbc=false
# neighbor_list=square_neighbors(sites, coords, pbc);

# for neighbor in neighbor_list[2]
#     println(get_square_site_no(neighbor,Nx))
# end

# function Hofstadter1(Nx, Ny, alpha, periodicity)
    
#     neig = neighbors(Nx, Ny, periodicity)
#     coordinates = square_lattice(Nx, Ny)[2]
    
#     N = Nx*Ny
#     t=-1
#     H = zeros(Complex{Float64},Nx*Ny,Nx*Ny)
    
#     for m in 1:N
#         for n in 1:N
#             if m in neig[n] 
#                 if abs(coordinates[m,1]-coordinates[n,1])==Nx-1
#                     if coordinates[m,1] > coordinates[n,1]
#                         H[m,n] = t*exp(-1im*2*pi*alpha*coordinates[m,2])
#                     elseif coordinates[m,1] < coordinates[n,1]
#                         H[m,n] = t*exp(1im*2*pi*alpha*coordinates[m,2])
#                     end
#                 else
#                     if coordinates[m,1] > coordinates[n,1]
#                         H[m,n] = t*exp(1im*2*pi*alpha*coordinates[m,2])
#                     elseif coordinates[m,1] < coordinates[n,1]
#                         H[m,n] = t*exp(-1im*2*pi*alpha*coordinates[m,2])
#                     else
#                         H[m,n] = t*exp(0)
#                     end
#                 end
#             else
#                 H[m,n] = 0
#             end
#         end
#     end
    
#     return H
# end

# executeme

using OffsetArrays

# executeme

function square_lattice(Nx,Ny)    
    site_idx = range(1,Nx*Ny) 
    lattice = OffsetArray(reshape(site_idx, (Nx,Ny)), 0:Nx-1, 0:Ny-1) |> transpose
    coordinates = []
    for y in 0:Ny-1
        for x in 0:Nx-1
            coordinates = [coordinates; x; y]
        end
    end
    coordinates = reshape(coordinates, (2, Nx*Ny)) |> transpose
    
    return lattice, coordinates
end

# executeme

function neighbors(Nx, Ny, periodicity)
    
    lattice = square_lattice(Nx,Ny)[1]
    Neighbors = []

    # Periodicity On
    if periodicity == 0
           
        for j in 0:Ny-1
            for i in 0:Nx-1
                x = [lattice[mod(j,Ny),mod(i-1,Nx)],lattice[mod(j+1,Ny),mod(i,Nx)],lattice[mod(j,Ny),mod(i+1,Nx)],lattice[mod(j-1,Ny),mod(i,Nx)]]
                x = unique(x)
                push!(Neighbors,x)
            end
            
        end
    # Periodicity Off (Hard-Wall)
    elseif periodicity == 1
        
        for j in 0:Ny-1
            for i in 0:Nx-1
                if j == 0 || i == 0 || j == Ny-1 || i == Nx-1 
                    new_neighbors = []
                    if j != 0
                        push!(new_neighbors, lattice[j-1,i])  
                    end
                    if i != 0
                        push!(new_neighbors, lattice[j,i-1])  
                    end
                    if j != Ny-1
                        push!(new_neighbors, lattice[j+1,i])  
                    end
                    if i != Nx-1
                        push!(new_neighbors, lattice[j,i+1])  
                    end
                else
                    new_neighbors = [
                        lattice[j,i-1],
                        lattice[j+1,i],
                        lattice[j,i+1],
                        lattice[j-1,i]
                        ]
                    push!(Neighbors,new_neighbors)
                end
            Neighbors = push!(Neighbors,new_neighbors)
            Neighbors = unique(Neighbors)
            end
        end
        
    end
    
    return Neighbors
end

# executeme

function Hofstadter_SP(Nx, Ny, alpha ,periodicity)
    
    neig = neighbors(Nx, Ny, periodicity)
    coordinates = square_lattice(Nx, Ny)[2]
    
    N = Nx*Ny
    t = -1
    H = zeros(Complex{Float64},Nx*Ny,Nx*Ny)
    
    for m in 1:N
        for n in 1:N
            if m in neig[n] 

                if abs(coordinates[m,1]-coordinates[n,1])==Nx-1
                    if coordinates[m,1] > coordinates[n,1]
                        H[m,n] = t*exp(-1im*2*pi*alpha*coordinates[m,2])
                    elseif coordinates[m,1] < coordinates[n,1]
                        H[m,n] = t*exp(1im*2*pi*alpha*coordinates[m,2])
                    end
                    
                elseif abs(coordinates[m,2]-coordinates[n,2])==Ny-1 #Magneto Periodic BC
                    if coordinates[m,2] > coordinates[n,2]
                        H[m,n] = t*exp(1im*2*pi*alpha*coordinates[m,1]*Ny)
                    elseif coordinates[m,2] < coordinates[n,2]
                        H[m,n] = t*exp(-1im*2*pi*alpha*coordinates[m,1]*Ny)
                    end
                    
                else
                    if coordinates[m,1] > coordinates[n,1]
                        H[m,n] = t*exp(1im*2*pi*alpha*coordinates[m,2])
                    elseif coordinates[m,1] < coordinates[n,1]
                        H[m,n] = t*exp(-1im*2*pi*alpha*coordinates[m,2])
                    else
                        H[m,n] = t*exp(0)
                    end
                    
                end
            else
                
                H[m,n] = 0
            end
        end
    end
    
    return H
end

#   ToplamPlaketFaz.pdf
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

# Nx=Ny=4
# alpha = 1/2

# (2*pi)*alpha

# angle(Hofstadter_SP(Nx, Ny, alpha, 0)[13,16]*Hofstadter_SP(Nx, Ny, alpha, 0)[1,13]*Hofstadter_SP(Nx, Ny, alpha, 0)[4,1]*Hofstadter_SP(Nx, Ny, alpha, 0)[16,4])

# angle(Hofstadter_SP(Nx, Ny, alpha, 0)[5,8]*Hofstadter_SP(Nx, Ny, alpha, 0)[9,5]*Hofstadter_SP(Nx, Ny, alpha, 0)[12,9]*Hofstadter_SP(Nx, Ny, alpha, 0)[8,12])

# angle(Hofstadter_SP(Nx, Ny, alpha, 0)[7,6]*Hofstadter_SP(Nx, Ny, alpha, 0)[11,7]*Hofstadter_SP(Nx, Ny, alpha, 0)[10,11]*Hofstadter_SP(Nx, Ny, alpha, 0)[6,10])

# angle(Hofstadter_SP(Nx, Ny, alpha, 0)[15,14]*Hofstadter_SP(Nx, Ny, alpha, 0)[3,15]*Hofstadter_SP(Nx, Ny, alpha, 0)[2,3]*Hofstadter_SP(Nx, Ny, alpha, 0)[14,2])

# Hofstadter_SP(Nx, Ny, alpha, 0)[13,16]

# mod(2*pi*alpha*(1-Ny*Nx),2*pi)

# angle(Hofstadter_SP(3, 3, 1/5, 0)[8,7])+angle(Hofstadter_SP(3, 3, 1/5, 0)[2,8])+angle(Hofstadter_SP(3, 3, 1/5, 0)[1,2])+angle(Hofstadter_SP(3, 3, 1/5, 0)[7,1])

# angle(Hofstadter_SP(3, 3, 1/5, 0)[2,1])+angle(Hofstadter_SP(3, 3, 1/5, 0)[5,2])+angle(Hofstadter_SP(3, 3, 1/5, 0)[4,5])+angle(Hofstadter_SP(3, 3, 1/5, 0)[1,4])

# 5.026548245743669 - 2*pi

# using LinearAlgebra
# E,U = eigen(Hofstadter_SP(3,3, 1/3 ,0))
# E