using SagbiHomotopy
using HomotopyContinuation
using LinearAlgebra
using Combinatorics

function get_sagbi_grassmannian(k::Int64, m::Int64)
    @var x[1:(k * (m - k))]

    M = hcat(I, transpose(reshape(x, m - k, k)))

    indx = collect(powerset([i for i in 1:m], k, k))
    sagbi = [det(M[:, ind]) for ind in indx]
    w = [i * (m - j + 1) for i in 0:(k - 1) for j in 1:(m - k) ]
    return x, sagbi, w
end

#run the function for the first time without storing the time
(x1, sagbi1, w1) = get_sagbi_grassmannian(2, 4);
A1 = randn(ComplexF64, 5, 6);
@unique_var l;
S1 = A1 * sagbi1 - l .* vcat(1, x1);
solve(S1)


j = 8

for m in 4:j
    for i in 2:floor(Int, m / 2)
        k = i
        (x, sagbi, w) = get_sagbi_grassmannian(k, m)

        A = randn(ComplexF64, k * (m - k) + 1, binomial(m, k))
        @var l
        S = A * sagbi - l .* vcat(1, x)

        result, elapsed_time, allocations, _, _ = @timed solve(S)
        #open("CCEquations/CC_timing_log_solve.txt", "a") do io
        #    println(io, "k = ", k, ", ", "m = ", m)
        #    println(io, "degree: ", length(solutions(result)))
        #    println(io, "Execution time: ", elapsed_time, " seconds")
        #    println(io, "Memory allocated: ", allocations, " bytes")
        #end
    end
end
