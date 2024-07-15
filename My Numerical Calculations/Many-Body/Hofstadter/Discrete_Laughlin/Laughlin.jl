ParCoord(b, i, SiteCoords) = SiteCoords[findall(x->x≠0,b[i])]

ComplexCoords(b, i, SiteCoords) = Complex.(getindex.(ParCoord(b, i, SiteCoords), 1),getindex.(ParCoord(b, i, SiteCoords), 2))

Comb(b, i, SiteCoords) = collect(combinations(ComplexCoords(b, i, SiteCoords),2))

e(b, i, SiteCoords) = exp(-sum(abs.(ComplexCoords(b, i, SiteCoords)).^2)/(4*lb))

ComplexCoordsDiff(b, i, SiteCoords) = getindex.(Comb(b, i, SiteCoords),1)-getindex.(Comb(b, i, SiteCoords),2)

"""
Wave Functions
"""
Ψ_L(b, i, m, SiteCoords) =  prod( (ComplexCoordsDiff(b, i, SiteCoords).^m) ) * e(b, i, SiteCoords) 
Ψ_1QH(b, i, SiteCoords) = Ψ_L(b, i, m, SiteCoords) * prod(ComplexCoords(b, i, SiteCoords).-QhCoord_1)
Ψ_2QH(b, i, SiteCoords) = Ψ_L(b, i, m, SiteCoords) * prod(ComplexCoords(b, i, SiteCoords).-QhCoord_1) * prod(ComplexCoords(b, i, SiteCoords).-QhCoord_2)

"""
Return coefficients of basis
"""
function get_Coefficients(b, m, SiteCoords)
    Coeff_L = [ Ψ_L(b, i, m, SiteCoords) for i in 1:length(b)]
    Coeff_1QH = [ Ψ_1QH(b, i, SiteCoords) for i in 1:length(b)]
    Coeff_2QH = [ Ψ_2QH(b, i, SiteCoords) for i in 1:length(b)]
    return Coeff_L, Coeff_1QH, Coeff_2QH
end

function density(Nx, Ny, N, mb, Coeff)
    Laughlin_Vec = normalize(Ket(mb,Coeff))
    DensArray = zeros(N)
    for j in 1:length(b)
        DensArray += abs( Laughlin_Vec.data[j] )^2 .*b[j]
    end
    DensArray = reshape(DensArray,Nx,Ny)
    return DensArray
end