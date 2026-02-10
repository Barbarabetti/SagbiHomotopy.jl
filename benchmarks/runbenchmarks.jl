using BenchmarkTools
using SagbiHomotopy
using Oscar
using HomotopyContinuation
using Random

const SUITE = BenchmarkGroup()

include("bench_paper_examples.jl")
include("bench_grassmannians.jl")
include("bench_oscillators.jl")

benchmark_paper_examples!(SUITE)
benchmark_grassmannians!(SUITE)
benchmark_oscillators!(SUITE)

BenchmarkTools.tune!(SUITE)

results = BenchmarkTools.run(SUITE; verbose = true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
