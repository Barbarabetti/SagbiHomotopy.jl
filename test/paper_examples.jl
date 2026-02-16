@testset "Paper examples" begin
    Random.seed!(20260210)


    @testset "Weight deformation sanity" begin
        Random.seed!(20260210)
        @var x y z
        sagbi = [[x, y, (x^2 + y^2), 1], [y, z, (x^2 + y^2), (x^3 + z^3)]]
        w = get_weight(sagbi)
        @test w !== nothing

        vars = HomotopyContinuation.variables(SagbiHomotopy.flatten_expressions(sagbi))
        @unique_var t

        for poly in SagbiHomotopy.flatten_expressions(sagbi)
            deformed = weight_deformation_for_poly(poly, vars, w; t = t)
            @test subs(deformed, t => 1) == poly
            @test subs(deformed, t => 0) == get_leading_monomial(poly, vars, w)
        end
    end

    @testset "Example 4.1 (P3)" begin
        Random.seed!(20260210)
        @var x y z
        sagbi = [x^2 + 1, y^2 + 1, x * y + z^2, 1]
        # Use the weight reported in the paper (Example 4.1).
        w = Int[-2, -2, -3]
        @test degree_map(sagbi) == 8
        @test degree_monomial_map(sagbi, w) == 8

        # The fiber cardinality for a *generic* linear section is 8.
        # To avoid rare unlucky random choices, try a small fixed list of RNG seeds.
        @var p[1:4]
        got = nothing
        for seed in (20260210, 20260211, 20260212)
            rng = MersenneTwister(seed)
            lin_sys = [rand(rng, -100:100, 3, 4) * p]
            # Skip the internal randomized degree check here; we test degrees above.
            res, sols = sagbi_homotopy(lin_sys, sagbi; weight = w, degreeCheck = false)
            if length(sols) == 8
                got = (res, sols)
                break
            end
        end
        @test got !== nothing
        res, sols = got
        @test length(sols) == 8
        @test length(solutions(res)) == 8
    end

    @testset "Example 4.2 (semimixed)" begin
        Random.seed!(20260210)
        @var x y z
        sagbi = [[x, y, (x^2 + y^2), 1], [y, z, (x^2 + y^2), (x^3 + z^3)]]
        w = get_weight(sagbi)
        @test w !== nothing
        @test degree_map(sagbi) == 1
        @test degree_monomial_map(sagbi, w) == 1

        rng = MersenneTwister(20260210)
        @var p[1:4]
        @var q[1:4]
        lin_sys = [
            rand(rng, -100:100, 2, 4) * p,
            rand(rng, -100:100, 1, 4) * q,
        ]

        res, sols = sagbi_homotopy(lin_sys, sagbi; weight = w)
        @test length(sols) == 6
        @test length(solutions(res)) == 6
    end

    @testset "Example 4.3 (base locus)" begin
        Random.seed!(20260210)
        @var x y
        sagbi = [[x * (x^2 + y^2 - 2 * x), x * (5 - 4 * y), y * (x^2 + y^2 - 2 * x), y * (5 - 4 * y)]]
        # Use the weight reported in the paper (Example 4.3).
        w = Int[-2, -1]
        @test degree_map(sagbi) == 2
        @test degree_monomial_map(sagbi, w) == 1

        rng = MersenneTwister(20260210)
        @var p[1:4]
        lin_sys = [rand(rng, ComplexF64, 2, 4) * p]

        res, sols = sagbi_homotopy(lin_sys, sagbi; weight = w)
        @test length(sols) == 2

        res2, sols2 = sagbi_homotopy(lin_sys, sagbi; weight = w, getBaseLocus = true)
        @test length(sols2) == 4
        @test length(solutions(res2)) == 2
    end

    @testset "Example 3.5 (degree drop)" begin
        Random.seed!(20260210)
        @var x y z
        sagbi = [x^3, y^3, x^6 * y^3 + x * z^2 + y^3]
        # Weight reported in the paper (Example 3.5).
        w = Int[-3, -3, -8]

        @test degree_map(sagbi) == 3
        @test degree_monomial_map(sagbi, w) == 18
    end

    @testset "Section 5.1 (Gr(2,6) linear slice)" begin
        # From the snippet in Section 5.1: slicing Gr(2,6) with 8 hyperplanes gives 14 solutions.
        rng = MersenneTwister(20260210)
        (_, sagbi, w) = get_sagbi_grassmannian(2, 6)

        A = randn(rng, ComplexF64, 8, 15)
        @var z[1:15]
        F = A * z

        res, sols = sagbi_homotopy(F, sagbi; weight = w, degreeCheck = false)
        @test length(sols) == 14
        @test length(solutions(res)) == 14
    end

    @testset "Section 5.1 (Gr(2,5) linear slice)" begin
        # Gr(2,5) has dimension 6 and degree 5 in the Plücker embedding.
        rng = MersenneTwister(20260210)
        (_, sagbi, w) = get_sagbi_grassmannian(2, 5)

        A = randn(rng, ComplexF64, 6, 10)
        @var z[1:10]
        F = A * z

        res, sols = sagbi_homotopy(F, sagbi; weight = w, degreeCheck = false)
        @test length(sols) == 5
        @test length(solutions(res)) == 5
    end

    @testset "Section 5.1 (Gr(3,6) linear slice)" begin
        # Gr(3,6) has dimension 9 and degree 42 in the Plücker embedding.
        rng = MersenneTwister(20260210)
        (_, sagbi, w) = get_sagbi_grassmannian(3, 6)

        A = randn(rng, ComplexF64, 9, 20)
        @var z[1:20]
        F = A * z

        res, sols = sagbi_homotopy(F, sagbi; weight = w, degreeCheck = false)
        @test length(sols) == 42
        @test length(solutions(res)) == 42
    end

    @testset "Section 5.2 (nonlinear resonator / oscillators, N=1 M=2)" begin
        # Derived from the MathRepo tutorial / oscillators scripts (Section 5.2 context).
        # For N=1, M=2 and weight \omega = (2,2,1,1), we get 25 solutions.
        rng = MersenneTwister(20260210)
        N = 1
        M = 2

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
            end
        end
        sagbi = [polynomials[i, :] for i in 1:(2 * M * N)]

        # Deterministic linear equations: each equation is a linear combination of its SAGBI row.
        c = rand(rng, Float64, 2 * M * N, L)
        lin_sys = [[sum(z[i, :] .* c[i, :])] for i in 1:(2 * M * N)]

        w = vcat(fill(2, M), fill(1, M))
        res, sols = sagbi_homotopy(lin_sys, sagbi; weight = w, degreeCheck = false)
        @test length(sols) == 25
        @test length(solutions(res)) == 25

        # Cross-check with a direct solve of the expanded polynomial system.
        full_eqs = [subs(lin_sys[i][1], collect(z[i, :]) => sagbi[i]) for i in 1:length(sagbi)]
        res_full = HomotopyContinuation.solve(full_eqs; show_progress = false)
        @test length(HomotopyContinuation.solutions(res_full)) == 25
    end
end
