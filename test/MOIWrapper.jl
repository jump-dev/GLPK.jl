using MathOptInterface

const MOI  = MathOptInterface
const MOIT = MathOptInterface.Test

@testset "Unit Tests" begin
    config = MOIT.TestConfig()
    solver = GLPK.Optimizer()

    # MOIT.basic_constraint_tests(solver, config)

    MOIT.unittest(solver, config, [
        # These are excluded because GLPK does not support quadratics.
        "solve_qp_edge_cases",
        "solve_qcp_edge_cases"
    ])

    MOIT.modificationtest(solver, config, [
        # This is excluded because LQOI does not support setting the constraint
        # function.
        "solve_func_scalaraffine_lessthan"
    ])
end

@testset "Linear tests" begin
    solver = GLPK.Optimizer()
    MOIT.contlineartest(solver, MOIT.TestConfig(), [
        # GLPK returns InfeasibleOrUnbounded
        "linear8a",
        # Requires infeasiblity certificate for variable bounds
        "linear12"
    ])
end

@testset "Linear Conic tests" begin
    MOIT.lintest(GLPK.Optimizer(), MOIT.TestConfig(infeas_certificates=false))
end

@testset "Integer Linear tests" begin
    MOIT.intlineartest(GLPK.Optimizer(), MOIT.TestConfig(), [
        # int2 is excluded because SOS constraints are not supported.
        "int2"
    ])
end

@testset "ModelLike tests" begin
    solver = GLPK.Optimizer()
    @testset "nametest" begin
        MOIT.nametest(solver)
    end
    @testset "validtest" begin
        MOIT.validtest(solver)
    end
    @testset "emptytest" begin
        MOIT.emptytest(solver)
    end
    @testset "orderedindicestest" begin
        # MOIT.orderedindicestest(solver)
    end
    @testset "canaddconstrainttest" begin
        MOIT.canaddconstrainttest(solver, Float64, Complex{Float64})
    end
    @testset "copytest" begin
        MOIT.copytest(solver, GLPK.Optimizer())
    end
end

@testset "Parameter setting" begin
    solver = GLPK.Optimizer(tm_lim=1, ord_alg=2, alien=3)
    @test solver.simplex.tm_lim == 1
    @test solver.intopt.tm_lim == 1
    @test solver.interior.ord_alg == 2
    @test solver.intopt.alien == 3
end
