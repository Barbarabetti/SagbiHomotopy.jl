
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
    JET.test_package(SagbiHomotopy)
end
