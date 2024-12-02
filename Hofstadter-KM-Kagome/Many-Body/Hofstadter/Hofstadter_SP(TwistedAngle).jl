"""
DEFAULT: Periodicity == Periodic
"""
function HSP_T(Nx, Ny, ϕ, Tx, Ty)

    neig = neighbors(Nx, Ny, 0)
    coordinates = square_lattice(Nx, Ny)[2]
    
    H = zeros(Complex{Float64},Nx*Ny,Nx*Ny)
    
    for m in 1:Nx*Ny
        for n in 1:Nx*Ny
            if m in neig[n] 
                
                if abs(coordinates[m,1]-coordinates[n,1])==Nx-1
                    if coordinates[m,1] > coordinates[n,1]
                        H[m,n] = -exp(-1im*2*pi*ϕ*coordinates[m,2])*exp(-1im*Tx)
                    elseif coordinates[m,1] < coordinates[n,1]
                        H[m,n] = -exp(1im*2*pi*ϕ*coordinates[m,2])*exp(1im*Tx)
                    end
                    
                elseif abs(coordinates[m,2]-coordinates[n,2])==Ny-1 # Twist + Magneto Periodic B.C.
                    if coordinates[m,2] > coordinates[n,2]
                        H[m,n] = -exp(1im*Ty)
                    elseif coordinates[m,2] < coordinates[n,2]
                        H[m,n] = -exp(-1im*Ty)
                    end
                    
                else
                    if coordinates[m,1] > coordinates[n,1]
                        H[m,n] = -exp(1im*2*pi*ϕ*coordinates[m,2])
                    elseif coordinates[m,1] < coordinates[n,1]
                        H[m,n] = -exp(-1im*2*pi*ϕ*coordinates[m,2])
                    else
                        H[m,n] = -exp(0)
                    end
                end
                
            else
                H[m,n] = 0
            end
            
        end
    end
    return H
end