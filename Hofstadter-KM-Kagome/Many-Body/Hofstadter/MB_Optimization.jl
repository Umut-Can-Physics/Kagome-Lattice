using BenchmarkTools
using QuantumOptics
using LinearAlgebra
using Base.Threads

# Define the system and operator
N = 9*9
b = NLevelBasis(N)
op = randoperator(b)
M = 3
s = bosonstates(b, [M])
mb = ManyBodyBasis(b, s)

function MB_1(mb, op)
    # Prepare a SparseOperator for each thread
    n_threads = Threads.nthreads()
    thread_results = [SparseOperator(mb) for _ in 1:n_threads]

    # Parallel loop with independent storage for each thread
    @threads for i in 1:N
        local_op = thread_results[threadid()]
        for j in 1:N
            if op.data[i, j] != 0  # Skip zero entries for efficiency
                local_op += op.data[i, j] * transition(mb, i, j)
            end
        end
        thread_results[threadid()] = local_op  # Ensure the result is stored in thread_results
    end

    # Combine the results from all threads
    sparse_mb_op = sum(thread_results)

    return sparse_mb_op
end

MB_op_1 = @btime MB_1(mb, op); # 24 sec

function MB_2(mb, op)
    sparse_mb_op_2 = SparseOperator(mb) 
    for i in 1:N 
        for j in 1:N 
            sparse_mb_op_2 += op.data[i,j] * transition(mb, i, j)  
        end 
    end 
    return sparse_mb_op_2
end

MB_op_2 = @btime MB_2(mb, op); # 1 min 23 sec

function MB_3(mb, op)
    # Prepare transition operators once for all (i, j)
    transitions = [transition(mb, i, j) for i in 1:N, j in 1:N]

    # Initialize sparse operator
    sparse_mb_op = SparseOperator(mb)

    # Parallelized loop to populate sparse_mb_op directly
    @threads for i in 1:N
        for j in 1:N
            value = op.data[i, j]
            if value != 0  # Skip zero entries for efficiency
                @inbounds sparse_mb_op += value * transitions[i, j]
            end
        end
    end
    return sparse_mb_op
end

MB_op_3 = @btime MB_3(mb, op);

function MB_4(mb, op)
    # Prepare transition operators once for all (i, j)
    transitions = [transition(mb, i, j) for i in 1:N, j in 1:N]

    # Collect non-zero entries to construct sparse matrix in a single pass
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]

    # Populate row, column, and value arrays
    @threads for i in 1:N
        for j in 1:N
            value = op.data[i, j]
            if value != 0
                push!(row_indices, i)
                push!(col_indices, j)
                push!(values, value)
            end
        end
    end

    # Construct the sparse operator in one go
    sparse_mb_op = SparseOperator(mb, row_indices, col_indices, values)

    return sparse_mb_op
end

MB_op_4 = @btime MB_4(mb, op);

println(MB_op_1 == MB_op_2)