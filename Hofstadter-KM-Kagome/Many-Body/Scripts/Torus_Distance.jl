function get_coords_square(Nx, Ny)
    Coords = []
    for i in 1:Nx*Ny
        push!(Coords, [i-Nx*div(i-1,Nx)-1,div(i-1,Nx)])
    end
    return Coords
end

# Find coordinate of a site in a specific coordinate distrincts
function find_co(coordinates, p) 
    coordinates = hcat(coordinates...)'
    # coordinates: Distrincts of Periodic Site
    x_co = coordinates[:,1]
    y_co = coordinates[:,2]
    x = hcat(x_co, y_co)[p, 1]
    y = hcat(x_co, y_co)[p, 2]
    return x, y
end

function distance(coordinates_p1, p1, coordinates_p2, p2)
    x2 = find_co(coordinates_p2, p2)[1]
    x1 = find_co(coordinates_p1, p1)[1]
    y2 = find_co(coordinates_p2, p2)[2]
    y1 = find_co(coordinates_p1, p1)[2]
    return sqrt( (x2 - x1)^2 + (y2 - y1)^2 )
end

# Minimum Distances of Equivalent Sites from Reference Site
function get_radii(coordinates, ref_site, N)
    distances = []
    for site_idx in 1:N 
        push!(distances, distance(coordinates, site_idx, coordinates, ref_site))
    end
    radius_list = sort(unique(distances));
    return radius_list
end

function Inner_Sites(Radius,ref_site,coords)
    
    inner_sites = []
    
    for i in 1:Nx*Ny
        
        Δx = abs(coords[i][1] - coords[ref_site][1])
        Δx = min(Δx,mod(-Δx,Nx-1)+1)

        Δy = abs(coords[i][2] - coords[ref_site][2])
        Δy = min(Δy,mod(-Δy,Ny-1)+1)
        
        if Δx^2 + Δy^2 <= (Radius)^2
            push!(inner_sites, i)
        end
        
    end
    
    return inner_sites
end

function charge_depletion_prof(radius_list,ref_site,coords,ref_par_density,avg_density)
    S_list = []
    for Radius in radius_list
        S = 0
        for i in Inner_Sites(Radius,ref_site,coords)
            S+=ref_par_density.-avg_density[(coords[i].+1)...]
        end
        push!(S_list, S)
    end
    return S_list
end