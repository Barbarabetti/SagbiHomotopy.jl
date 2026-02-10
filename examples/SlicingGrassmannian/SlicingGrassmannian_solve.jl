using SagbiHomotopy
using HomotopyContinuation
using LinearAlgebra
using Combinatorics

function get_sagbi_grassmannian(k::Int64, m::Int64)
    @var x[1:k*(m-k)] 

    M = hcat(I,transpose(reshape(x,m-k,k)))
    
    indx = collect(powerset([i for i=1:m], k, k))
    sagbi = [det(M[:,ind]) for ind in indx] 
    w = [i*(m-j+1) for i=0:k-1 for j=1:m-k ]
    return x, sagbi, w
end

#run the function for the first time without storing the time
(x1, sagbi1, w1) = get_sagbi_grassmannian(2,4);
A1 = randn(ComplexF64,4,6);
S1 = System(A1*sagbi1);
solve(S1);
 
j = 7

for m=4:j
    for i=2:floor(Int, m/2)
        println(i)
        k = i
        (x, sagbi, w) = get_sagbi_grassmannian(k,m)
                
        # Get a linear space of dimension dim(G(k,m)) in the same ambient space of the Grassmannian
        A = randn(ComplexF64,k*(m-k),binomial(m,k))
        S = System(A*sagbi)

        result, elapsed_time, allocations, _, _ = @timed solve(S)
        #open("SlicingGrassmannian/Grassmannians_timing_log_solve.txt", "a") do io
        #    println(io, "k = ", k, ", ", "m = ", m)
        #    println(io, "degree: ", length(solutions(result)))
        #    println(io, "Execution time: ", elapsed_time, " seconds")
        #    println(io, "Memory allocated: ", allocations, " bytes")
        #end 
    end
end