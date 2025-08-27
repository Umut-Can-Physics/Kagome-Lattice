using LinearAlgebra
using OffsetArrays
using Plots
using ProgressMeter
using QuantumOptics
using Revise
using SparseArrays
using LaTeXStrings
using Combinatorics
includet("../Scripts/FirstBandApproximation.jl")
includet("../Scripts/ManyBody.jl")
includet("Hofstadter_SP.jl")
includet("Discrete_Laughlin/Laughlin.jl")

Nx = 4
Ny = 4
p = 1
q = 4
α = p/q
N = Nx*Ny
Nϕ = N*α
pn = 2

method = "NoProj"
HardCore = true

# fermionstates for hardcore is true
basis = basis = fermionstates(NLevelBasis(Nx*Ny),pn)
NN = 10
type = "fermion"
ψL0, ψL1 = GeneralizedLaughlin(basis, Nx, Ny, Nϕ, NN, type)

δU = 0.1
rangeU = -1:δU:1
Nev = 3
s = length(rangeU)

function TuneU(pn, Nx, Ny, p, q, N, HardCore, method, Nev)
    ψ_gs = [];ψ_fes = []; E = []
    for U in rangeU
        println("U=$U")
        ϵ, ψ = SolveModel(pn, U, Nx, Ny, p, q, N, HardCore, method, Nev)
        push!(ψ_gs, ψ[1])
        push!(ψ_fes, ψ[2])
        push!(E, ϵ)
    end
    return ψ_gs, ψ_fes, E
end

Overlaps_1 = zeros(Complex, s)
for i in 1:s
    ψ1 = ψ_gs[i].data
    ψ2 = ψL0
    Overlaps_1[i] = Overlap(ψ1, ψ2)
end

Overlaps_2 = zeros(Complex, s)
for i in 1:s
    ψ1 = ψ_gs[i].data
    ψ2 = ψL1
    Overlaps_2[i] = Overlap(ψ1, ψ2)
end

Overlaps_3 = zeros(Complex, s)
for i in 1:s
    ψ1 = ψ_fes[i].data
    ψ2 = ψL0
    Overlaps_3[i] = Overlap(ψ1, ψ2)
end

Overlaps_4 = zeros(Complex, s)
for i in 1:s
    ψ1 = ψ_fes[i].data
    ψ2 = ψL1
    Overlaps_4[i] = Overlap(ψ1, ψ2)
end

P1 = scatter(rangeU, OverlapResult_1, xlabel=L"U", ylabel="Ov")
P2 = scatter(rangeU, OverlapResult_2, xlabel=L"U", ylabel="Ov")
P3 = scatter(rangeU, OverlapResult_3, xlabel=L"U", ylabel="Ov")
P4 = scatter(rangeU, OverlapResult_4, xlabel=L"U", ylabel="Ov")
plot(P1, P2,P3, P4, layout=(2,2), size=(800,400), legend=false)

HardCore = true
ϵ, ψ = SolveModel(pn, U, Nx, Ny, p, q, N, HardCore, method, Nev)
ψ_gs = ψ[1].data
ψ_fes = ψ[2].data

ψ1 = ψ_gs
ψ2 = ψL0
O1 = Overlap(ψ1, ψ2)

ψ1 = ψ_gs
ψ2 = ψL1
O2 = Overlap(ψ1, ψ2)

ψ1 = ψ_fes
ψ2 = ψL0
O3 = Overlap(ψ1, ψ2)

ψ1 = ψ_fes
ψ2 = ψL1
O4 = Overlap(ψ1, ψ2)

scatter([1,2,3,4], [O1,O2,O3,O4],
xticks=([1, 2, 3, 4], [L"l=0", L"l=1", L"l=0", L"l=1"]),
yticks=([O1,O2,O3,O4], [L"\psi_0",L"\psi_0", L"\psi_1", L"\psi_1"]),
title="Hardcore=$(HardCore)",
legend=false,
annotations = [
        (1, O1, text(string(O1), :black, :left, 8)), # x_coord, y_coord, text_object
        (2, O2, text(string(O2), :black, :left, 8)),
        (3, O3, text(string(O3), :black, :left, 8)),
        (4, O4, text(string(O4), :black, :bottom, 8))
    ]
)
