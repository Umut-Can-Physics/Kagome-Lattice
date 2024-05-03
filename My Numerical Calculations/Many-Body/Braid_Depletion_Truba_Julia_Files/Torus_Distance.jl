#   Torus Distance
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

using OffsetArrays

# Number operatörünün sıralamasına uygun olan lattice indislemesi
function get_square_lattice(Nx, Ny)
    
    site_idx = range(1,Nx*Ny) 
    latticee = reverse(OffsetArray(reshape(site_idx, (Ny,Nx)), 1:Ny, 1:Nx), dims=1)
    coordinates = []
    for x in 0:Nx-1
        for y in 0:Ny-1
            coordinates = [coordinates; y; x]
        end
    end
    coordinates = reshape(coordinates, (2, Nx*Ny)) |> transpose
    
    return coordinates, latticee
end

# Compatible with site number of plotted lattice
# Same as coordinates of get_square_lattice(Nx, Ny) function
function get_coords_square(Nx, Ny)
    Coords = []
    for i in 1:Nx*Ny
        push!(Coords, [i-Nx*div(i-1,Nx)-1,div(i-1,Nx)])
    end
    return Coords
end

# Ghost Sites
function get_ghost_sites(Nx, Ny, coordinates)
    coordinates_top_left = hcat(coordinates[:,1].-Nx,coordinates[:,2].+Ny)
    coordinates_top = hcat(coordinates[:,1],coordinates[:,2].+Ny)
    coordinates_top_right = hcat(coordinates[:,1].+Nx,coordinates[:,2].+Ny)
    coordinates_left = hcat(coordinates[:,1].-Nx,coordinates[:,2])
    coordinates_right = hcat(coordinates[:,1].+Nx,coordinates[:,2])
    coordinates_bottom_left = hcat(coordinates[:,1].-Nx,coordinates[:,2].-Ny)
    coordinates_bottom = hcat(coordinates[:,1],coordinates[:,2].-Ny)
    coordinates_bottom_right = hcat(coordinates[:,1].+Nx,coordinates[:,2].-Ny);
    co_districts = [coordinates_top_left,coordinates_top,coordinates_top_right,coordinates_left,coordinates_right,coordinates,coordinates_bottom_left,coordinates_bottom,coordinates_bottom_right];
    return coordinates_top_left, coordinates_top, coordinates_top_right, coordinates_left, coordinates_right, coordinates_bottom_left, coordinates_bottom, coordinates_bottom_right, co_districts
end

function plot_lat(Nx,Ny)
    center_x = (Nx-1)/2
    center_y = (Ny-1)/2
    scatter([coordinates_top_left[:,1]],[coordinates_top_left[:,2]],legend=false, alpha=0.3);annotate!( center_x-Nx, center_y+Ny, "top-left")
    scatter!([coordinates_top[:,1]],[coordinates_top[:,2]], alpha=0.3);annotate!(center_x, center_y+Ny, "top")
    scatter!([coordinates_top_right[:,1]],[coordinates_top_right[:,2]], alpha=0.3);annotate!(center_x+Nx, center_y+Ny, "top-right")
    scatter!([coordinates_left[:,1]],[coordinates_left[:,2]], alpha=0.3);annotate!(center_x-Nx, center_y, "left")
    scatter!([coordinates[:,1]], [coordinates[:,2]]);annotate!(center_x,center_y, "center")
    scatter!([coordinates_right[:,1]],[coordinates_right[:,2]], alpha=0.3);annotate!(center_x+Nx, center_y, "right")
    scatter!([coordinates_bottom_left[:,1]],[coordinates_bottom_left[:,2]], alpha=0.3);annotate!(center_x-Nx, center_y-Ny, "bottom-left")
    scatter!([coordinates_bottom[:,1]],[coordinates_bottom[:,2]], alpha=0.3);annotate!(center_x, center_y-Ny, "bottom")
    scatter!([coordinates_bottom_right[:,1]],[coordinates_bottom_right[:,2]], alpha=0.3);annotate!(center_x+Nx, center_y-Ny, "bottom-right")
end

# Find coordinate of a site in a specific coordinate distrincts
function find_co(coordinates, p) 
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

#   Equivalent (Periodic) Site Condition
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function get_all_sites(Nx, Ny, coordinates)
    All_Coordinates = vcat(coordinates_top_left,coordinates_top,coordinates_top_right,coordinates_left,coordinates,coordinates_right,coordinates_bottom_left,coordinates_bottom,coordinates_bottom_right)
    X_Coordinates = All_Coordinates[:,1]
    Y_Coordinates =  All_Coordinates[:,2]
    # Index and Coordinates of All Sites
    idx_idx = []
    for i in 1:9 
        for j in 1:Nx*Ny
            push!(idx_idx, j)
        end
    end
    # All_Coordinates_2, her bir Nx*Ny'lik sütun, sırasıyla bölge koordinatlarını soldan sağa ve aşağıdan yukarıya gösterir.
    # Her bir bölgede 1'den Nx*Ny'ye site indekslemesi tüm siteler için yapılır.
    All_Coordinates_2 = hcat(idx_idx,All_Coordinates)
    return All_Coordinates, X_Coordinates, Y_Coordinates, All_Coordinates_2
end

function get_eq_site(XX, YY)
    Intersect = intersect(findall(x->x==XX, All_Coordinates[:,1]), findall(x->x==YY, All_Coordinates[:,2]))
    Equivalent_Site = All_Coordinates_2[:,1][Intersect][1]
    return Equivalent_Site
end