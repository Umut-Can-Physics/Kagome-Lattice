#   \begin{align} \PsiL \left(\xi1,...,\xiN \right) & \propto \Pi{i<j}
#   \left(\xii - \xij \right)^2 e^{-\sum{i=1}^N |\xii|^2/4lb^2} \ \Psi{QH}
#   \left(\xi1,...,\xiN,\mathcal{Q} \right) & \propto \Pi{i=1}^N \left(\xii -
#   \mathcal{Q} \right) \PsiL \left(\xi1,...,\xi_N \right) \end{align}

using QuantumOptics
using Combinatorics
using Plots
using LaTeXStrings
using NBInclude
nbexport("Laughlin.jl","Discrete_Laughlin.ipynb")

Nx = 5
Ny = 5
N = Nx*Ny
p = 1
q = 5
PN = 5

α = p/q
lb = 1/sqrt(2*pi*α)

QhCoord = [3+im*2]

NLevel = NLevelBasis(N)

C(n,r) = Int(factorial(big(n))/(factorial(big(n-r))*factorial(big(r))))

C(N,PN)

St = fermionstates(NLevel,PN)

mb = ManyBodyBasis(NLevel,St)

SiteCoords = [ [x,y] for x in -(Nx-1)/2:(Nx-1)/2 for y in -(Ny-1)/2:(Ny-1)/2 ]

ParCoord(i) = [ SiteCoords[site] for site in findall(x->x≠0,St[i])]

# ParCoord(i, SiteCoords) = filter(>([0,0]),St[i].*SiteCoords)

i = 2

ParCoord(i, SiteCoords)

ParCoord(84, SiteCoords)

ComplexCoords(i, SiteCoords) = Complex.(getindex.(ParCoord(i, SiteCoords), 1),getindex.(ParCoord(i, SiteCoords), 2))

ComplexCoords(i, SiteCoords)

ComplexCoords(125, SiteCoords)

Comb(i, SiteCoords) = collect(combinations(ComplexCoords(i, SiteCoords),2))

Comb(125, SiteCoords)

e(i, SiteCoords) = exp(-sum(abs.(ComplexCoords(i, SiteCoords)).^2)/(4*lb))

e(125, SiteCoords)

ComplexCoordsDiff(i, SiteCoords) = getindex.(Comb(i, SiteCoords),1)-getindex.(Comb(i, SiteCoords),2)

ComplexCoordsDiff(125, SiteCoords)

Ψ_L(i, SiteCoords) =  prod( (ComplexCoordsDiff(i, SiteCoords).^2) )*e(i, SiteCoords) 
Ψ_1QH(i, SiteCoords) = Ψ_L(i, SiteCoords) 
#* prod(ComplexCoords(i, SiteCoords).-QhCoord)

SiteCoords

AllComplexCoords = Complex.(getindex.(SiteCoords, 1),getindex.(SiteCoords, 2))

AllComplexCoords.- Ref(QhCoord[1])

# ΨΨ_QH = prod(filter(x->x≠(0+0im),(AllComplexCoords .- Ref(QhCoord[1]))))

diag = collect(1:length(St))
coeff = [ Ψ_1QH(i, SiteCoords) for i in 1:length(St)]
ΨΨ = sparse(diag,diag,coeff)

∑ψψ = ΨΨ.nzval

Kett = normalize(Ket(mb,∑ψψ))

#(abs.(Kett.data).^2)

number(mb,1);

expectation = [expect(number(mb,n),Kett) for n in 1:N]

sum(expectation)

den = reshape(expectation, Nx, Ny)

heatmap(real(den),
    title=L"N_x= %$(Nx), N_y= %$(Ny), PN= %$(PN)"*"\n"*L"|Ψ_{L,QH}> ∝ ∑ Ψ_{L,QH}(z_1,z_2,...,z_n) c^{\dagger}_{z_1}...c^{\dagger}_{z_n}|vac>"*"\n"
)

zeros(25)

DensArray = zeros(N)
for j in 1:length(St)
    DensArray += abs( Ψ(j, SiteCoords)*ΨΨ_QH )^2 .*St[j]
end

heatmap(reshape(DensArray,Nx,Ny))