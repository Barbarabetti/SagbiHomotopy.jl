using BenchmarkTools
using Random
using SagbiHomotopy
using HomotopyContinuation

function oscillator_instance(N::Int, M::Int; seed::Int = 20260210)
    rng = MersenneTwister(seed)

    @var u[1:M, 1:N] v[1:M, 1:N]
    L = M + 3 + N - 1
    @var z[1:(2 * M * N), 1:L]

    polynomials = Matrix{HomotopyContinuation.Expression}(undef, 2 * M * N, L)
    for k in 1:N
        for i in 1:M
            r1 = 2 * i - 1 + 2 * M * (k - 1)
            r2 = 2 * i + 2 * M * (k - 1)

            polynomials[r1, 1] = (u[1, 1])^0
            polynomials[r1, 2] = u[i, k]
            polynomials[r1, 3] = v[i, k]

            polynomials[r2, 1] = (u[1, 1])^0
            polynomials[r2, 2] = u[i, k]
            polynomials[r2, 3] = v[i, k]

            for j in 1:M
                polynomials[r1, j + 3] = u[i, k] * (u[j, k]^2 + v[j, k]^2)
                polynomials[r2, j + 3] = v[i, k] * (u[j, k]^2 + v[j, k]^2)
            end

            # Coupling terms (present when N>1).
            for j in 1:(N - 1)
                lst = filter(jj -> jj != k, collect(1:N))
                polynomials[r1, j + M + 3] = v[i, lst[j]]
                polynomials[r2, j + M + 3] = u[i, lst[j]]
            end
        end
    end

    sagbi = [polynomials[i, :] for i in 1:(2 * M * N)]
    c = rand(rng, Float64, 2 * M * N, L)
    lin_sys = [[sum(z[i, :] .* c[i, :])] for i in 1:(2 * M * N)]

    # Weight pattern from the oscillator example: u-variables get weight 2, v-variables get weight 1.
    w = vcat(fill(2, M * N), fill(1, M * N))

    full_eqs = [subs(lin_sys[i][1], collect(z[i, :]) => sagbi[i]) for i in 1:length(sagbi)]

    return (lin_sys = lin_sys, sagbi = sagbi, w = w, full_eqs = full_eqs)
end

function benchmark_oscillators!(suite::BenchmarkTools.BenchmarkGroup)
    grp = BenchmarkTools.BenchmarkGroup()
    suite["Oscillators"] = grp

    # Section 5.2 / tutorial-style nonlinear resonator system.
    # Case 1: smallest nontrivial instance (N=1, M=2).
    inst_12 = oscillator_instance(1, 2; seed = 20260210)
    case_12 = BenchmarkTools.BenchmarkGroup()
    grp["N=1 M=2"] = case_12
    case_12["sagbi_homotopy"] = @benchmarkable sagbi_homotopy(
        $(inst_12.lin_sys),
        $(inst_12.sagbi);
        weight = $(inst_12.w),
        degreeCheck = false,
    ) seconds = 20
    case_12["direct solve"] = @benchmarkable HomotopyContinuation.solve(
        $(inst_12.full_eqs);
        show_progress = false,
    ) seconds = 20

    # Case 2: heavier instance (N=2, M=2).
    # This introduces cross-oscillator coupling terms (the extra "N-1" columns).
    inst_22 = oscillator_instance(2, 2; seed = 20260211)
    case_22 = BenchmarkTools.BenchmarkGroup()
    grp["N=2 M=2"] = case_22
    case_22["sagbi_homotopy"] = @benchmarkable sagbi_homotopy(
        $(inst_22.lin_sys),
        $(inst_22.sagbi);
        weight = $(inst_22.w),
        degreeCheck = false,
    ) seconds = 20
    case_22["direct solve"] = @benchmarkable HomotopyContinuation.solve(
        $(inst_22.full_eqs);
        show_progress = false,
    ) seconds = 20

    return nothing
end
