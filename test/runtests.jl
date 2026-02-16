using SagbiHomotopy
using Test
using Oscar
using HomotopyContinuation
using Random


@testset "code quality" begin
    include("code_quality.jl")
end

@testset "paper examples" begin
    include("paper_examples.jl")
end
