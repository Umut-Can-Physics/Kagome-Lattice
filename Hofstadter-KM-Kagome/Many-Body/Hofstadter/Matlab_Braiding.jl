# Started on 13.12.2024

### PARAMETERS ###
p = 1
q = 8
Lx = 8
Ly = 8
Np = 2
extra = 1  # # of extra states above the degenerate manifold (to avoid possible convergence issues to find the right number of degenerate states)
Nmax = 1  # maximum site occupation
U = 0  # interaction strength (doesn't matter for the Nmax = 1 hardcore case)
V_imp = 0.5
V_randImp = 1e-5

# Compute dependent parameters
Ns = Lx * Ly
alpha = p / q
Nphi = Int(Ns * alpha)
Nd = (Nphi - 2) - 2 * Np  # in the presence of TWO impurities !!!
degeneracy = factorial(Nd - 1 + Np) / (factorial(Nd - 1) * factorial(Np)) * (Nphi - 2) / Nd  # in the presence of TWO impurities !!!
nu_eff = Np / (Nphi - 2)  # in the presence of TWO impurities !!!

# Define model parameters and other constants
model = "Hofs"
BC = "MPBC"
# model = "KM"
gauge = 0  # 0 for symmetric, 1 for Landau

nstep = 10  # # of steps between two nearest-neighbor sites
simultaneous = 0  # path for simultaneously (1) or separately (0) moving quasiholes?

# Basis creation
firstlist = 0:Nmax
alist = collect(0:Nmax)
blist2 = []
for v1 in 1:(Ns-1)
    dummy = repeat(alist, 1, length(firstlist))
    dummy2 = Matrix{Int}(undef, 1, 0)  # Initialize as an empty matrix with one row and zero columns

    for v2 in 1:length(firstlist)
        dummy2 = hcat(dummy2, repeat([firstlist[v2]], 1, size(dummy, 2) ÷ length(firstlist)))
    end

    alist = hcat(dummy, dummy2)

end