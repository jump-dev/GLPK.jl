using MathOptInterface

const MOI  = MathOptInterface
const MOIT = MathOptInterface.Test

@testset "Unit Tests" begin
    config = MOIT.TestConfig()
    solver = GLPK.Optimizer()

    # MOIT.basic_constraint_tests(solver, config)

    MOIT.unittest(solver, config, [
        "solve_qp_edge_cases",
        "solve_qcp_edge_cases"
    ])

    MOIT.modificationtest(solver, config, [
        "solve_func_scalaraffine_lessthan"
    ])
end

@testset "Linear tests" begin
    linconfig_nocertificate = MOIT.TestConfig()
    solver = GLPK.Optimizer()
    MOIT.contlineartest(solver, linconfig_nocertificate, [
        # GLPK returns InfeasibleOrUnbounded
        "linear8a",
        # Need to query an interval
        "linear10",
        # Requires infeasiblity certificate for variable bounds
        "linear12"
    ])
    MOIT.linear10test(solver, MOIT.TestConfig(query=false))
end

@testset "Linear Conic tests" begin
    linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
    solver = GLPK.Optimizer()
    MOIT.lintest(solver, linconfig_nocertificate)
end

@testset "Integer Linear tests" begin
    intconfig = MOIT.TestConfig()
    solver_mip = GLPK.Optimizer()
    MOIT.intlineartest(solver_mip, intconfig, ["int2","int1"])
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
        solver2 = GLPK.Optimizer()
        MOIT.copytest(solver,solver2)
    end
end
