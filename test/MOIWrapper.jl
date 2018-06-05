using GLPK, Base.Test, MathOptInterface, MathOptInterface.Test

const MOI  = MathOptInterface
const MOIT = MathOptInterface.Test

@testset "MathOptInterfaceGLPK" begin
    @testset "Unit Tests" begin
        @testset "LP solver" begin
            config = MOIT.TestConfig()
            solver = GLPKOptimizerLP()

            MOIT.basic_constraint_tests(solver, config;
                exclude = [
                    (MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64})
                ]
            )

            MOIT.basic_constraint_tests(solver, config;
                get_constraint_function = false,
                get_constraint_set      = false,
                include = [
                    (MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64})
                ]
            )

            MOIT.unittest(solver, config, [
                "solve_qp_edge_cases",
                "solve_qcp_edge_cases"
            ])
        end
        @testset "MIP solver" begin
            solver = GLPKOptimizerMIP()
            config = MOIT.TestConfig(duals=false)

            MOIT.basic_constraint_tests(solver, config;
                exclude = [
                    (MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64})
                ]
            )

            MOIT.basic_constraint_tests(solver, config;
                get_constraint_function = false,
                get_constraint_set      = false,
                include = [
                    (MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64})
                ]
            )

            MOIT.unittest(solver, config, [
                "solve_qp_edge_cases",
                "solve_qcp_edge_cases"
            ])
        end
    end

    @testset "Linear tests" begin
        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
        solver = GLPKOptimizerLP()
        MOIT.contlineartest(solver, linconfig_nocertificate, ["linear10"])

        linconfig_nocertificate_noduals = MOIT.TestConfig(duals=false,infeas_certificates=false)
        solver_mip = GLPKOptimizerMIP()
        MOIT.contlineartest(solver_mip, linconfig_nocertificate_noduals, ["linear10","linear8b","linear8c"])

        # Intervals
        linconfig_noquery = MOIT.TestConfig(query=false)
        MOIT.linear10test(solver, linconfig_noquery)
    end

    @testset "Linear Conic tests" begin
        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
        solver = GLPKOptimizerLP()
        MOIT.lintest(solver, linconfig_nocertificate)

        linconfig_nocertificate_noduals = MOIT.TestConfig(duals=false,infeas_certificates=false)
        solver_mip = GLPKOptimizerMIP()
        MOIT.lintest(solver_mip, linconfig_nocertificate_noduals)
    end

    @testset "Integer Linear tests" begin
        intconfig = MOIT.TestConfig()
        solver_mip = GLPKOptimizerMIP()
        MOIT.intlineartest(solver_mip, intconfig, ["int2","int1"])
    end

    @testset "ModelLike tests - MIP" begin
        intconfig = MOIT.TestConfig()
        solver = GLPKOptimizerMIP()
        MOIT.validtest(solver)
        MOIT.emptytest(solver)
        solver2 = GLPKOptimizerMIP()
        MOIT.copytest(solver,solver2)
    end
    @testset "ModelLike tests - LP" begin
        intconfig = MOIT.TestConfig()
        solver = GLPKOptimizerLP()
        MOIT.validtest(solver)
        MOIT.emptytest(solver)
        solver2 = GLPKOptimizerLP()
        MOIT.copytest(solver,solver2)
    end
end
