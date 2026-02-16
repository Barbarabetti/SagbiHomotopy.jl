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
F1, poss_sagbi = oscillators(1, 3);
sagbi_homotopy(F1, poss_sagbi);


N = 2
M = 3
@time get_weight(sagbi1)

for s in 2:N
    for r in 1:M
        F, sagbi = oscillators(s, r)

        result, elapsed_time, allocations, _, _ = @timed sagbi_homotopy(F, sagbi, degreeCheck = false, varyLinearPart = true)
        #open("Oscillators/oscillators_timing_log_sagbidetection.txt", "a") do io
        #    println(io, "N = ", s, ", ", "M = ", r)
        #    println(io, "degree: ", length(result[2]))
        #    println(io, "Execution time: ", elapsed_time, " seconds")
        #    println(io, "Memory allocated: ", allocations, " bytes")
        #end
    end
end
