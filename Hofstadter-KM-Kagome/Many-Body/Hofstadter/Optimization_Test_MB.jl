cd(@__DIR__)

using Pkg
println("\n",Pkg.status(), VERSION)

using QuantumOptics, SparseArrays, Plots, LinearAlgebra, ProgressMeter, Revise, OffsetArrays, LaTeXStrings, BenchmarkTools, DelimitedFiles, Serialization
includet("../Scripts/FirstBandApproximation.jl")
includet("../Scripts/ManyBody.jl")
includet("Hofstadter_SP.jl")
includet("../Scripts/Impurity.jl")
includet("../Scripts/Braiding.jl")

function get_mb_op_optimized(mb_basis, sp_op)

    mb_op = SparseOperator(mb_basis)
    N = size(sp_op.data, 1) # More robust way to get dimension

    # Pre-allocate lists to build the sparse matrix components directly
    I = Vector{Int}()
    J = Vector{Int}()
    V = Vector{ComplexF64}()

    println("Building many-body operator...")
    @showprogress for j in 1:N
        for i in 1:N
            if !iszero(sp_op.data[i, j])
                transition_op = transition(mb_basis, i, j)
                rows, cols, vals = findnz(transition_op.data)
                scale_factor = sp_op.data[i, j]
                for k in eachindex(vals)
                    push!(I, rows[k])
                    push!(J, cols[k])
                    push!(V, vals[k] * scale_factor)
                end
            end
        end
    end

    # Construct the sparse matrix directly
    sparse_data = sparse(I, J, V, only(mb_basis.shape), only(mb_basis.shape))
    mb_op.data = sparse_data

    return mb_op
end

pn = 3
Nx = Ny = 6
N = Nx*Ny
p = 1; q = Ny
matrix = Hofstadter_SP(Nx, Ny, p / q, 0);
basis_sp = NLevelBasis(N);
H1_op = Sp_Op(basis_sp, matrix);
mb_basis = get_Bosonic_MB_Basis(basis_sp, pn, false);

optimized_mb_op = get_mb_op_optimized(mb_basis, H1_op)

EE, UU = eigenstates(optimized_mb_op)

scatter(EE[1:10])