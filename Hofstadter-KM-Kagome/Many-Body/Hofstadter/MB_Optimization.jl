using BenchmarkTools
using QuantumOptics
using LinearAlgebra
using Base.Threads
using SparseArrays

# Define the system and operator
N = 9*9
b = NLevelBasis(N)
op = randoperator(b)
M = 3
s = bosonstates(b, [M])
mb = ManyBodyBasis(b, s)
MB_OP = manybodyoperator(op, mb)

function MB_2(N, mb, op)
    sparse_mb_op_2 = SparseOperator(mb) 
    for i in 1:N 
        for j in 1:N 
            sparse_mb_op_2 += op.data[i,j] * transition(mb, i, j)  
        end 
    end 
    return sparse_mb_op_2
end

function process_transitions_threaded(N, mb, op)
    sparse_mb_op_2 = SparseOperator(mb) 
    @threads for i in 1:N
        for j in 1:N
            local_mb = deepcopy(mb)  # Ensure a thread-local copy
            sparse_mb_op_2 += op.data[i,j] * transition(mb, i, j)  
        end
    end
    return sparse_mb_op_2
end

MB_2(N, mb, op) == process_transitions_threaded(N, mb, op) # FALSE !

function process_transitions_threaded_2(N, mb, op)
    # Initialize thread-local SparseOperators
    thread_results = [SparseOperator(mb) for _ in 1:nthreads()]
    
    @threads for i in 1:N
        thread_id = threadid()  # Identify thread
        for j in 1:N
            local_mb = deepcopy(mb)  # Ensure a thread-local copy
            thread_results[thread_id] += op.data[i, j] * transition(mb, i, j)
        end
    end
    
    # Combine thread-local results
    combined_result = SparseOperator(mb)
    for result in thread_results
        combined_result += result
    end
    
    return combined_result
end

MB_2(N, mb, op) == process_transitions_threaded_2(N, mb, op) # TRUE !