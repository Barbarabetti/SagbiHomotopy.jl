using SagbiHomotopy
using HomotopyContinuation

function oscillators(N, M)
    @var u[1:M, 1:N], v[1:M, 1:N]
    @var z[1:(2 * M * N), 1:(M + 3 + N - 1)]

    polynomials = Matrix{Expression}(undef, 2 * M * N, M + 3 + N - 1)

    for k in 1:N
        for i in 1:M

            polynomials[2 * i - 1 + 2 * M * (k - 1), 1] = (u[1, 1])^0
            polynomials[2 * i - 1 + 2 * M * (k - 1), 2] = u[i, k]
            polynomials[2 * i - 1 + 2 * M * (k - 1), 3] = v[i, k]

            polynomials[2 * i + 2 * M * (k - 1), 1] = (u[1, 1])^0
            polynomials[2 * i + 2 * M * (k - 1), 2] = u[i, k]
            polynomials[2 * i + 2 * M * (k - 1), 3] = v[i, k]

            for j in 1:M
                polynomials[2 * i - 1 + 2 * M * (k - 1), j + 3] = u[i, k] * (u[j, k]^2 + v[j, k]^2)
            end

            for j in 1:(N - 1)
                lst = filter(j -> j != k, collect(1:N))
                polynomials[2 * i - 1 + 2 * M * (k - 1), j + M + 3] = v[i, lst[j]]
            end

            for j in 1:M
                polynomials[2 * i + 2 * M * (k - 1), j + 3] = v[i, k] * (u[j, k]^2 + v[j, k]^2)
            end

            for j in 1:(N - 1)
                lst = filter(j -> j != k, collect(1:N))
                polynomials[2 * i + 2 * M * (k - 1), j + M + 3] = u[i, lst[j]]
            end
        end
    end
    sagbi = [polynomials[i, :] for i in 1:(2 * M * N)]

    c = rand(Float64, 2 * M * N, M + 3 + N - 1)
    lin_sys = [[sum(z[i, :] .* c[i, :])] for i in 1:(2 * M * N)]

    return (lin_sys, sagbi)
end

#run the function for the first time without storing the time
F1, sagbi1 = oscillators(2, 1);
S1 = vcat([subs(F1[i], variables(F1[i]) => sagbi1[i]) for i in 1:length(sagbi1)]...);
solve(S1);


N = 2
M = 3


for s in 1:N
    for r in 1:M
        F, sagbi = oscillators(s, r)
        S = vcat([subs(F[i], variables(F[i]) => sagbi[i]) for i in 1:length(sagbi)]...)

        result, elapsed_time, allocations, _, _ = @timed solve(S)
        #open("Oscillators/oscillators_timing_log_solve.txt", "a") do io
        #    println(io, "N = ", s, ", ", "M = ", r)
        #    println(io, "degree: ", length(solutions(result)))
        #    println(io, "Execution time: ", elapsed_time, " seconds")
        #    println(io, "Memory allocated: ", allocations, " bytes")
        #end
    end
end
