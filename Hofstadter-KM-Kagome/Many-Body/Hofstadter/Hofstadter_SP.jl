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