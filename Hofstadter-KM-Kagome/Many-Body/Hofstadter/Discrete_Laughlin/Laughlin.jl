function BosonLabels(A::Vector{<:Real})
    result = Int[] # Initialize an empty vector to store the expanded indices

    # Iterate through the array, getting both the index (idx) and the value (val)
    for (idx, val) in enumerate(A)
        # We only care about positive values for repetition count
        if val > 0
            # Ensure we're using an integer count for repetition
            # If val is a float, like 3.5, `Int(val)` truncates to 3.
            # If you want to round, use `round(Int, val)`.
            num_repetitions = Int(val)

            # Append the current index 'num_repetitions' times to the result
            for _ in 1:num_repetitions
                push!(result, idx)
            end
        end
    end
    return result
end

function ParCoord(b, i, SiteCoords, type)
    if type == "boson"
        result = SiteCoords[BosonLabels(b[i])]
    elseif type=="fermion"
        result = SiteCoords[findall(x->x≠0, b[i])]
    end
    return result
end

ComplexCoords(b, i, SiteCoords, type) = Complex.(getindex.(ParCoord(b, i, SiteCoords, type), 1),getindex.(ParCoord(b, i, SiteCoords, type), 2))

e(b, i, SiteCoords) = exp(-sum(abs.(ComplexCoords(b, i, SiteCoords, type)).^2)/(4*lb))

# combination function automatically use the i<j condition
Comb(b, i, SiteCoords, type) = collect(combinations(ComplexCoords(b, i, SiteCoords, type), 2))

# z_i - z_j
ComplexCoordsDiff(b, i, SiteCoords, type) = getindex.(Comb(b, i, SiteCoords, type), 1) - getindex.(Comb(b, i, SiteCoords, type), 2)

# The definition of the Jacobi-Theta function
function v(a, b, z, τ, NN)
    Sum = zeros(length(z))
    for n in 1:NN
        Sum = Sum .+ exp.(im*pi*τ*(n+a)^2 .+ 2*pi*im*(n+a).*(z.+b))
    end
    return Sum
end

# Relative part of the Laughlin wave function
function Relative(basis, i, SiteCoords, Lx, Ly, NN, type)
    a = b = 1/2
    z = ComplexCoordsDiff(basis, i, SiteCoords, type)./Lx # [list]
    τ = im*Ly/Lx
    return  prod( v(a, b, z, τ, NN).^2 )
end

function Z(basis, i, SiteCoords, type)
    return sum(ComplexCoords(basis, i, SiteCoords, type))
end

# Center of mass part of the Laughlin wave function
function CenterOfMass(basis, i, SiteCoords, Lx, Ly, l, Nϕ, type)
    a = l/2 + (Nϕ-2)/4
    b = -(Nϕ-2)/2
    z = 2 * Z(basis, i, SiteCoords, type)/Lx
    τ = 2*im*Ly/Lx
    return v(a, b, z, τ, NN)
end

function GeneralizedLaughlin(basis, Lx, Ly, Nϕ, NN, type)
    #Note that Nx and Ny should be odd, because coordinates start from the center of the lattice.
    SiteCoords = [ [x,y] for x in -(Lx-1)/2:(Lx-1)/2 for y in -(Ly-1)/2:(Ly-1)/2 ]

    ψ_rel = [Relative(basis, i, SiteCoords, Lx, Ly, NN, type) for i in 1:length(basis)]

    # where l=0 and l=1 refers to the two degenerate ground state
    l = 0
    ψ_CM0 = [CenterOfMass(basis, i, SiteCoords, Lx, Ly, l, Nϕ, type) for i in 1:length(basis)]
    ψ0 = normalize(ψ_rel.*only.(ψ_CM0))

    l = 1
    ψ_CM1 = [CenterOfMass(basis, i, SiteCoords, Lx, Ly, l, Nϕ, type) for i in 1:length(basis)]
    ψ1 = normalize(ψ_rel.*only.(ψ_CM1))
    return ψ0, ψ1
end

"""
Wave Functions
"""
Ψ_L(b, i, m, SiteCoords, type) =  prod( (ComplexCoordsDiff(b, i, SiteCoords, type).^m) ) * e(b, i, SiteCoords) 
Ψ_1QH(b, i, SiteCoords, type) = Ψ_L(b, i, m, SiteCoords, type) * prod(ComplexCoords(b, i, SiteCoords, type).-QhCoord_1)
Ψ_2QH(b, i, SiteCoords, type) = Ψ_L(b, i, m, SiteCoords, type) * prod(ComplexCoords(b, i, SiteCoords, type).-QhCoord_1) * prod(ComplexCoords(b, i, SiteCoords,type).-QhCoord_2)

"""
Return coefficients of basis
"""
function get_Coefficients(b, m, SiteCoords)
    Coeff_L = [ Ψ_L(b, i, m, SiteCoords,type) for i in 1:length(b)]
    #Coeff_1QH = [ Ψ_1QH(b, i, SiteCoords, type) for i in 1:length(b)]
    #Coeff_2QH = [ Ψ_2QH(b, i, SiteCoords, type) for i in 1:length(b)]
    return Coeff_L
end

function density(Nx, Ny, N, mb, Coeff)
    Laughlin_Vec = normalize(Ket(mb, Coeff)) 
    DensArray = zeros(N)
    for j in 1:length(b)
        DensArray += abs( Laughlin_Vec.data[j] )^2 .*b[j]
    end
    DensArray = reshape(DensArray,Nx,Ny)
    return DensArray
end