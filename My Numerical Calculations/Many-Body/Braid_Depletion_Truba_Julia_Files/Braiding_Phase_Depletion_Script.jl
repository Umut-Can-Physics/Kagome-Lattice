#   Script of Braiding Phase \varphi_{br} by depletion profile
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function impurity_control(V_Strength, V0, Imp_Site, N_Pin)
    Impurity_Data = Impurity(V0, Imp_Site)
    factor = 9 # particle density interpolation factor
    E, Sub_Number_MB_Operator_List, basis_cut_mb, Degeneracy, Total_H, Sub_Number_MB_Operator_List, r_hubbard_states,pn, NPhi0 = Get_MB(Nx, Ny, p, q, cut_off, PN, U, Impurity_Data, factor, N_Pin)
    avg_density = Get_Avg_Density(Nx, Ny, Degeneracy, N, Sub_Number_MB_Operator_List, basis_cut_mb, r_hubbard_states)
    return E, Degeneracy, avg_density
end

function get_ref_prtcl_density(par_num,p,q,NPhi,N_Pin)
    return (par_num*p/q)/((NPhi-N_Pin))
end

# function impurity_control2(Nx, Ny, p, q, cut_off, PN, U, V_Strength, V0, Imp_Site, N_Pin)
#     Impurity_Data = Impurity(V0, Imp_Site)
#     factor = 9 # particle density interpolation factor
#     N = Nx*Ny
#     E, Sub_Number_MB_Operator_List, basis_cut_mb, Degeneracy, Total_H, Sub_Number_MB_Operator_List, r_hubbard_states,pn, NPhi0 = Get_MB(Nx, Ny, p, q, cut_off, PN, U, Impurity_Data, factor, N_Pin)
#     avg_density = Get_Avg_Density(Nx, Ny, Degeneracy, N, Sub_Number_MB_Operator_List, basis_cut_mb, r_hubbard_states)
#     return avg_density
# end

function get_radius_list(Nx, Ny, ref_site)
    N = Nx*Ny
    coordinates, latticee = get_square_lattice(Nx, Ny)
    coordinates_top_left, coordinates_top, coordinates_top_right, coordinates_left, coordinates_right, coordinates_bottom_left, coordinates_bottom, coordinates_bottom_right, co_districts = get_ghost_sites(Nx, Ny, coordinates)
    radius_list = get_radii(coordinates, ref_site, N);
    return radius_list
end

function torus_distance_func(Coords,ref_site, site_idx)
    x1 = ref_site[1]
    y1 = ref_site[2]
    x2 = site_idx[1]
    y2 = site_idx[2]
    return sqrt((mod(x2-x1,Nx))^2 + (mod(y2-y1,Ny)^2))
end

function distance_func(Coords,ref_site, site_idx)
    x1 = ref_site[1]
    y1 = ref_site[2]
    x2 = site_idx[1]
    y2 = site_idx[2]
    return sqrt( (x2-x1)^2 + (y2-y1)^2)
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

# function Inner_Sites_copy(Nx, Ny, Radius, ref_site, coords)
    
#     inner_sites = []
    
#     for i in 1:Nx*Ny
        
#         Δx = abs(coords[i][1] - coords[ref_site][1])
#         Δx = min(Δx,mod(-Δx,Nx-1)+1)

#         Δy = abs(coords[i][2] - coords[ref_site][2])
#         Δy = min(Δy,mod(-Δy,Ny-1)+1)
        
#         if Δx^2 + Δy^2 <= (Radius)^2
#             push!(inner_sites, i)
#         end
        
#     end
    
#     return inner_sites
# end

avg_density_func(j) = avg_density[(Coords[j].+1)...]
d_k_qh(j) = ref_par_density .- avg_density_func(j)

function fractional_charge(Nx, Ny, ref_site)
    d_qh_list = []
    radius_list = get_radius_list(Nx, Ny, ref_site)
    for ρ in radius_list
        d_qh_value = 0
        for j in Inner_Sites(ρ,ref_site,Coords)
            d_qh_value += d_k_qh(j) 
        end
        push!(d_qh_list, d_qh_value) 
    end
    return radius_list, d_qh_list
end

function fractional_charge_2(Nx, Ny, ref_site)
    d_qh_list = []
    radius_list = get_radius_list(Nx, Ny, ref_site)
    for ρ in radius_list
        d_qh_value = 0
        for j in Inner_Sites(ρ,ref_site,Coords)
            d_qh_value += d_k_qh(j) * distance_func(Coords, Coords[ref_site], Coords[j])^2
        end
        push!(d_qh_list, d_qh_value) 
    end
    return d_qh_list
end

function get_density_list(ref_site)
    d_qh_list = []
    radius_list_ = get_radius_list(ref_site)
    for ρ in radius_list
        for j in Inner_Sites(ρ,ref_site,Coords)
            push!(d_qh_list, d_k_qh(j)) 
        end
    end
    return d_qh_list
end

function get_braid_phase(radius_list, d_2_qh_list, d_1_qh_list,p,q,n)
    a = 1 # lattica constant
    α = p/q # flux per plaquette
    l_b = a/sqrt(2*pi*α)
    φ_br_list = []
    for ρ in radius_list[1:n]
        φ_br = 0
        for j in Inner_Sites(ρ,ref_site,Coords)
            φ_br += ( (d_2_qh_list[1:n][j] .- 2*d_1_qh_list[1:n][j] ) * distance_func(Coords, Coords[ref_site], Coords[j])^2)
        end
        push!(φ_br_list,(2*pi/(2*l_b^2)*φ_br))
    end
    return mod.(φ_br_list,2*pi)
end

# function get_braid_phase2(ρ_list, D_2_QH, D_1_QH,p,q,INNER_SITES)
#     a = 1 # lattica constant
#     α = p/q # flux per plaquette
#     l_b = a/sqrt(2*pi*α) # magnetic length
#     φ_br_list = []
#     for ρ in ρ_list
#         φ_br = 0
#         for j in Inner_Sites_copy(10, 8, ρ, ref_site, Coords)
#             φ_br += ( (D_2_QH(10,8,1,10,3,2,4,j) - 2*D_1_QH(10,7,1,10,3,2,4,j) ) * ρ^2)
#         end
#         push!(φ_br_list,(2*pi/(2*l_b^2)*φ_br))
#     end
#     return mod.(φ_br_list,2*pi)
# end