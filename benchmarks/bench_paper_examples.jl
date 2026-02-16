using BenchmarkTools
using Random
using SagbiHomotopy
using HomotopyContinuation

function benchmark_paper_examples!(suite::BenchmarkTools.BenchmarkGroup)
    # Example 4.2 (semimixed) from the paper/tests.
    @var x y z
    sagbi = [[x, y, (x^2 + y^2), 1], [y, z, (x^2 + y^2), (x^3 + z^3)]]
    w = get_weight(sagbi)

    rng = MersenneTwister(20260210)
    @var p[1:4]
    @var q[1:4]
    lin_sys = [rand(rng, -100:100, 2, 4) * p,
        rand(rng, -100:100, 1, 4) * q]

    SUITE["semimixed"]["get_weight"] = @benchmarkable get_weight($sagbi) seconds = 10
    SUITE["semimixed"]["sagbi_homotopy"] = @benchmarkable sagbi_homotopy(
        $lin_sys,
        $sagbi;
        weight = $w,
        degreeCheck = false,
    ) seconds = 10

    return nothing
end
