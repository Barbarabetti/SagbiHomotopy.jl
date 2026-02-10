using SagbiHomotopy
using Test

@testset "Concretely typed" begin
    using CheckConcreteStructs
    all_concrete(SagbiHomotopy)

end

@testset "ExplicitImports" begin
    using ExplicitImports

    @test check_no_implicit_imports(SagbiHomotopy) == nothing
    @test check_all_explicit_imports_via_owners(SagbiHomotopy) == nothing
    @test check_all_explicit_imports_are_public(SagbiHomotopy) == nothing
    @test check_no_stale_explicit_imports(SagbiHomotopy) == nothing
    @test check_all_qualified_accesses_via_owners(SagbiHomotopy) == nothing
    @test check_all_qualified_accesses_are_public(SagbiHomotopy) == nothing
    @test check_no_self_qualified_accesses(SagbiHomotopy) == nothing
end

@testset "Code quality" begin
    using Aqua
    Aqua.test_all(SagbiHomotopy)
end

@testset "Code linting" begin
    using JET
    JET.test_package(SagbiHomotopy; target_defined_modules = true)
end

@testset "Sagbi Homotopy Computation" begin
    using Oscar
    using HomotopyContinuation
    using SagbiHomotopy

    HomotopyContinuation.ModelKit.@var x, y, z

    sagbi = [[x, y, (x^2 + y^2), 1], [y, z, (x^2 + y^2), (x^3 + z^3)]]
    w = get_weight(sagbi)

    @var p[1:4]
    @var q[1:4]

    lin_sys = [rand(ComplexF64, 2, 4) * p, rand(ComplexF64, 1, 4) * q]

    sagbi_homotopy(lin_sys, sagbi; weight = w)
end
