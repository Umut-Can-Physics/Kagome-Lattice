using OffsetArrays

a1_vec = [2 0]; a2_vec = [1 sqrt(3)]
b1=[0 0];b2=a1_vec/2;b3=a2_vec/2
Basis = [b1,b2,b3];

function get_sites(Nx, Ny, a1_vec, a2_vec, Basis)
    x_co = []; y_co = []
    coordinates = []; sites = []; atom_dict=Dict{}()
    for i1 in 0:Nx-1
        for i2 in 0:Ny-1
            Lp = a1_vec * i1 + a2_vec * i2
            center = a1_vec * i1 + a2_vec * i2 + b1
            shift = (-b2-b3)/sqrt(3)/2
            P1=[center+shift,center+a1_vec+shift,center+a2_vec+shift]
            for (ib, b) in zip(Iterators.countfrom(0), Basis)
                atom_vec = Lp + b
                atom_no = 3 * i1 * Ny + 3 * i2 + ib
                site = [i1,i2,ib]
                sites = push!(sites, site)
                coordinates = push!(coordinates, atom_vec)
                atom_dict[tuple(site)]=atom_vec
                x_co = push!(x_co, atom_vec[1]); y_co = push!(y_co, atom_vec[2])
            end
        end
    end
    sites = OffsetArray(sites, 0:Nx*Ny*3-1);
    x_co = OffsetArray(x_co, 0:Nx*Ny*3-1)
    y_co = OffsetArray(y_co, 0:Nx*Ny*3-1)
    
    return coordinates, sites, atom_dict, x_co, y_co
end

function plot_kagome(Nx,Ny)
    x_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[4], 1:Nx*Ny*3)
    y_co = OffsetArray(get_sites(Nx, Ny, a1_vec, a2_vec, Basis)[5], 1:Nx*Ny*3)
    N = Nx*Ny*3
    p = scatter([x_co],[y_co],series_annotations = text.(1:N, :bottom),grid=false, legend = false, aspect_ratio = :equal)
    return p
end

function get_NN_NNN(t1,L1,t2,L2)
    
    NN = [
          [[0,0,1], [0,0,2], [-1,0,1],[0,-1,2],   [-1,0,2], [-1,1,1], [0,-1,1],[1,-1,2]], 
          [[0,0,-1],[0,0,1], [1,0,-1],[1,-1,1],   [0,-1,1], [1,-1,-1],[0,1,-1],[1,0,1]], 
          [[0,0,-1],[0,0,-2],[0,1,-2],[-1,1,-1],  [-1,1,-2],[-1,0,-1],[0,1,-1],[1,0,-2]]
         ]

    NN = OffsetArray(NN, 0:2)

    hopps = [
             [t1+1im*L1,t1-1im*L1,t1+1im*L1,t1-1im*L1,  t2+1im*L2,t2-1im*L2,t2-1im*L2,t2+1im*L2], 
             [t1-1im*L1,t1+1im*L1,t1-1im*L1,t1+1im*L1,  t2-1im*L2,t2+1im*L2,t2+1im*L2,t2-1im*L2], 
             [t1-1im*L1,t1+1im*L1,t1+1im*L1,t1-1im*L1,  t2-1im*L2,t2+1im*L2,t2+1im*L2,t2-1im*L2]
            ] 

    hopps = OffsetArray(hopps, 0:2)
    
    return NN ,hopps
end

function get_H(Nx, Ny, t1, L1, t2, L2)  
    
    N = Nx*Ny*3
    
    HH = OffsetArray(zeros(Complex{Float64},N,N), 0:N-1, 0:N-1)
    
    coordinates, sites, atom_dict, x_co, y_co = get_sites(Nx, Ny, a1_vec, a2_vec, Basis)
    NN, hopps = get_NN_NNN(t1,L1,t2,L2)
    
    for atom_no in 0:N-1
        atom_site=sites[atom_no]
        for (i_delta, delta) in enumerate(NN[atom_site[3]])
            neighbor_site = atom_site+delta
            neighbor_site = OffsetArray(neighbor_site, 0:2)
            neighbor_site[0] = mod(neighbor_site[0],Nx)
            neighbor_site[1] = mod(neighbor_site[1],Ny)
            neighbor_no=3*neighbor_site[0]*Ny+3*neighbor_site[1]+neighbor_site[2]
            HH[neighbor_no,atom_no]=hopps[atom_site[3]][i_delta] 
        end
    end
    
    return HH
end

# offset array kullandığımdan linear algebra kütüphanesiyle eigen hesaplayamıyorum mecbur quantumoptics kütüphanesi ile operatöre çevirip hesapladım.
using QuantumOptics

function get_SP_H(Nx, Ny, t1, L1, t2, L2)
    
    N=Nx*Ny*3
    HH = get_H(Nx, Ny, t1, L1, t2, L2)  
    
    Basis = NLevelBasis(N)
    H = SparseOperator(Basis)
    for m in 0:N-1
        for n in 0:N-1
            H += HH[m,n] * transition(Basis, m+1, n+1)
        end
    end
    
    return H
end