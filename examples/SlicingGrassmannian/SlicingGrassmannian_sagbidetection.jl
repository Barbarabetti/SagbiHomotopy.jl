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
@unique_var z[1:binomial(4,2)];
F1 = A1*z;
sagbi_homotopy(F1, sagbi1);


j = 6

for m=4:j
    for i=2:floor(Int, m/2)
        k = i
        (x, sagbi, w) = get_sagbi_grassmannian(k,m)
                
        # Get a linear space of dimension dim(G(k,m)) in the same ambient space of the Grassmannian
        A = randn(ComplexF64,k*(m-k),binomial(m,k))
        @var z[1:binomial(m,k)]
        F = A*z

        result, elapsed_time, allocations, _, _ = @timed sagbi_homotopy(F, sagbi);
        #open("SlicingGrassmannian/Grassmannians_timing_log_sagbidetection.txt", "a") do io
        #    println(io, "k = ", k, ", ", "m = ", m)
        #    println(io, "degree: ", length(result[2]))
        #    println(io, "Execution time: ", elapsed_time, " seconds")
        #    println(io, "Memory allocated: ", allocations, " bytes")
        #end 
    end
end