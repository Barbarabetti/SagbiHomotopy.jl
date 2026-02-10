using BenchmarkTools
using Random
using SagbiHomotopy
using HomotopyContinuation

function benchmark_grassmannians!(suite::BenchmarkTools.BenchmarkGroup)
    # Section 5.1: linear slice of Gr(2,6).
    rng = MersenneTwister(20260210)
    (_, sagbi, w) = get_sagbi_grassmannian(2, 6)

    A = randn(rng, ComplexF64, 8, 15)
    @var z[1:15]
    F = A * z

    SUITE["Grassmannians"]["Gr(2,6) slice sagbi_homotopy"] = @benchmarkable sagbi_homotopy(
        $F,
        $sagbi;
        weight = $w,
        degreeCheck = false,
    ) seconds = 20

    return nothing
end
