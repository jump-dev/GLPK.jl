using GLPK, Test

const MOI  = GLPK.MathOptInterface
const MOIT = MOI.Test

const OPTIMIZER = MOI.Bridges.full_bridge_optimizer(
    GLPK.Optimizer(), Float64
)
const CONFIG = MOIT.TestConfig()

@testset "Unit Tests" begin
    MOIT.basic_constraint_tests(OPTIMIZER, CONFIG)
    MOIT.unittest(OPTIMIZER, CONFIG, [
        # FIXME `NumberOfThreads` not supported
        "number_threads",
        # These are excluded because GLPK does not support quadratics.
        "solve_qcp_edge_cases",
        "solve_qp_edge_cases",
        "delete_soc_variables",

        # Tested below because the termination status is different.
        "solve_zero_one_with_bounds_3",

        # TODO(odow): not implemented.
        "number_threads",
    ])
    @testset "solve_zero_one_with_bounds_3" begin
        MOI.empty!(OPTIMIZER)
        MOI.Utilities.loadfromstring!(OPTIMIZER,"""
            variables: x
            maxobjective: 2.0x
            c1: x in ZeroOne()
            c2: x >= 0.2
            c3: x <= 0.5
        """)
        MOI.optimize!(OPTIMIZER)
        # We test this here because the TerminationStatus is INVALID_MODEL not
        # INFEASIBLE.
        @test MOI.get(OPTIMIZER, MOI.TerminationStatus()) == MOI.INVALID_MODEL
    end
    MOIT.modificationtest(OPTIMIZER, CONFIG)
end

@testset "Linear tests" begin
@testset "Default Solver"  begin
        MOIT.contlineartest(OPTIMIZER, MOIT.TestConfig(basis = true), [
            # This requires an infeasiblity certificate for a variable bound.
            "linear12",
            # VariablePrimalStart not supported.
            "partial_start"
        ])
    end
    @testset "No certificate" begin
        MOIT.linear12test(OPTIMIZER, MOIT.TestConfig(infeas_certificates=false))
    end
end

@testset "Linear Conic tests" begin
    MOIT.lintest(OPTIMIZER, CONFIG)
end

@testset "Integer Linear tests" begin
    MOIT.intlineartest(OPTIMIZER, CONFIG, [
        "int2", "indicator1", "indicator2", "indicator3", "indicator4"
    ])
end

@testset "ModelLike tests" begin
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "GLPK"

    @testset "default_objective_test" begin
        MOIT.default_objective_test(OPTIMIZER)
    end

    @testset "default_status_test" begin
        MOIT.default_status_test(OPTIMIZER)
    end

    @testset "nametest" begin
        MOIT.nametest(OPTIMIZER)
    end

    @testset "validtest" begin
        MOIT.validtest(OPTIMIZER)
    end

    @testset "emptytest" begin
        MOIT.emptytest(OPTIMIZER)
    end

    @testset "orderedindicestest" begin
        MOIT.orderedindicestest(OPTIMIZER)
    end

    @testset "copytest" begin
        MOIT.copytest(
            OPTIMIZER,
            MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)
        )
    end

    @testset "scalar_function_constant_not_zero" begin
        MOIT.scalar_function_constant_not_zero(OPTIMIZER)
    end

    @testset "start_values_test" begin
        # We don't support ConstraintDualStart or ConstraintPrimalStart yet.
        # @test_broken MOIT.start_values_test(GLPK.Optimizer(), OPTIMIZER)
    end

    @testset "supports_constrainttest" begin
        # supports_constrainttest needs VectorOfVariables-in-Zeros,
        # MOIT.supports_constrainttest(GLPK.Optimizer(), Float64, Float32)
        # but supports_constrainttest is broken via bridges:
        MOI.empty!(OPTIMIZER)
        MOI.add_variable(OPTIMIZER)
        @test  MOI.supports_constraint(OPTIMIZER, MOI.SingleVariable, MOI.EqualTo{Float64})
        @test  MOI.supports_constraint(OPTIMIZER, MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64})
        # This test is broken for some reason:
        @test_broken !MOI.supports_constraint(OPTIMIZER, MOI.ScalarAffineFunction{Int}, MOI.EqualTo{Float64})
        @test !MOI.supports_constraint(OPTIMIZER, MOI.ScalarAffineFunction{Int}, MOI.EqualTo{Int})
        @test !MOI.supports_constraint(OPTIMIZER, MOI.SingleVariable, MOI.EqualTo{Int})
        @test  MOI.supports_constraint(OPTIMIZER, MOI.VectorOfVariables, MOI.Zeros)
        @test !MOI.supports_constraint(OPTIMIZER, MOI.VectorOfVariables, MOI.EqualTo{Float64})
        @test !MOI.supports_constraint(OPTIMIZER, MOI.SingleVariable, MOI.Zeros)
        @test !MOI.supports_constraint(OPTIMIZER, MOI.VectorOfVariables, MOIT.UnknownVectorSet)
    end

    @testset "set_lower_bound_twice" begin
        MOIT.set_lower_bound_twice(OPTIMIZER, Float64)
    end

    @testset "set_upper_bound_twice" begin
        MOIT.set_upper_bound_twice(OPTIMIZER, Float64)
    end
end

@testset "Parameter setting" begin
    solver = GLPK.Optimizer(tm_lim=1, ord_alg=2, alien=3)
    @test solver.simplex_param.tm_lim == 1
    @test solver.intopt_param.tm_lim == 1
    @test solver.interior_param.ord_alg == 2
    @test solver.intopt_param.alien == 3
end

@testset "Issue #79" begin
    @testset "An unbounded integer model" begin
        model = GLPK.Optimizer()
        MOI.Utilities.loadfromstring!(model, """
            variables: x, y
            minobjective: -5.0x + y
            c1: x in Integer()
            c2: x in LessThan(1.0)
        """)
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
    end

    @testset "An infeasible integer model" begin
        model = GLPK.Optimizer()
        MOI.Utilities.loadfromstring!(model, """
            variables: x
            minobjective: -5.0x
            c1: x in Integer()
            c2: x in LessThan(1.0)
            c3: 1.0x in GreaterThan(2.0)
        """)
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
    end
end

@testset "Issue #70" begin
    model = GLPK.Optimizer()
    x = MOI.add_variable(model)
    f =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0)
    s = MOI.LessThan(2.0)
    c = MOI.add_constraint(model, f, s)
    row = GLPK._info(model, c).row
    @test GLPK.get_row_type(model.inner, row) == GLPK.UP
    @test GLPK.get_row_lb(model.inner, row) == -GLPK.DBL_MAX
    @test GLPK.get_row_ub(model.inner, row) == 2.0
    # Modify the constraint set and verify that the internal constraint
    # has the correct bounds afterwards
    MOI.set(model, MOI.ConstraintSet(), c, MOI.LessThan(1.0))
    @test GLPK.get_row_type(model.inner, row) == GLPK.UP
    @test GLPK.get_row_lb(model.inner, row) == -GLPK.DBL_MAX
    @test GLPK.get_row_ub(model.inner, row) == 1.0
end

@testset "Infeasible bounds" begin
    model = GLPK.Optimizer()
    x = MOI.add_variable(model)
    MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Interval(1.0, -1.0))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INVALID_MODEL
end

@testset "RawParameter" begin
    model = GLPK.Optimizer(method = GLPK.SIMPLEX)
    exception = ErrorException("Invalid option: cb_func. Use the MOI attribute `GLPK.CallbackFunction` instead.")
    @test_throws exception MOI.set(model, MOI.RawParameter("cb_func"), (cb) -> nothing)
    MOI.set(model, MOI.RawParameter("tm_lim"), 100)
    @test MOI.get(model, MOI.RawParameter("tm_lim")) == 100
    @test_throws ErrorException MOI.get(model, MOI.RawParameter(:tm_lim))
    @test_throws ErrorException MOI.set(model, MOI.RawParameter(:tm_lim), 120)
    param = MOI.RawParameter("bad")
    @test_throws MOI.UnsupportedAttribute(param) MOI.set(model, param, 1)
    @test_throws MOI.UnsupportedAttribute(param) MOI.get(model, param)

    model = GLPK.Optimizer(method = GLPK.INTERIOR)
    exception = ErrorException("Invalid option: cb_func. Use the MOI attribute `GLPK.CallbackFunction` instead.")
    @test_throws exception MOI.set(model, MOI.RawParameter("cb_func"), (cb) -> nothing)
    MOI.set(model, MOI.RawParameter("tm_lim"), 100)
    @test MOI.get(model, MOI.RawParameter("tm_lim")) == 100
    @test_throws MOI.UnsupportedAttribute(param) MOI.set(model, MOI.RawParameter("bad"), 1)
    @test_throws MOI.UnsupportedAttribute(param) MOI.get(model, MOI.RawParameter("bad"))

    model = GLPK.Optimizer(method = GLPK.EXACT)
    exception = ErrorException("Invalid option: cb_func. Use the MOI attribute `GLPK.CallbackFunction` instead.")
    @test_throws exception MOI.set(model, MOI.RawParameter("cb_func"), (cb) -> nothing)
    MOI.set(model, MOI.RawParameter("tm_lim"), 100)
    @test MOI.get(model, MOI.RawParameter("tm_lim")) == 100
    @test_throws MOI.UnsupportedAttribute(param) MOI.set(model, MOI.RawParameter("bad"), 1)
    @test_throws MOI.UnsupportedAttribute(param) MOI.get(model, MOI.RawParameter("bad"))

    model = GLPK.Optimizer()
    MOI.set(model, MOI.RawParameter("mip_gap"), 0.001)
    @test MOI.get(model, MOI.RawParameter("mip_gap")) == 0.001
end

@testset "TimeLimitSec issue #110" begin
    model = GLPK.Optimizer(method = GLPK.SIMPLEX)
    MOI.set(model, MOI.TimeLimitSec(), nothing)
    @test MOI.get(model, MOI.RawParameter("tm_lim")) == typemax(Int32)
    MOI.set(model, MOI.TimeLimitSec(), 100)
    @test MOI.get(model, MOI.RawParameter("tm_lim")) == 100000
    @test MOI.get(model, MOI.TimeLimitSec()) == 100
    # conversion between ms and sec
    MOI.set(model, MOI.RawParameter("tm_lim"), 100)
    @test isapprox(MOI.get(model, MOI.TimeLimitSec()), 0.1)
end

@testset "RelativeGap" begin
    model = GLPK.Optimizer()
    MOI.Utilities.loadfromstring!(model, """
        variables: x
        minobjective: 1.0x
        c1: x in Integer()
        c2: x in GreaterThan(1.5)
    """)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.RelativeGap()) == 0.0

    model = GLPK.Optimizer()
    MOI.Utilities.loadfromstring!(model, """
        variables: x
        minobjective: 1.0x
        c1: x in GreaterThan(1.5)
    """)
    MOI.optimize!(model)
    @test_throws ErrorException MOI.get(model, MOI.RelativeGap())
end

@testset "Extra name tests" begin
    model = GLPK.Optimizer()
    @testset "Variables" begin
        MOI.empty!(model)
        x = MOI.add_variables(model, 2)
        MOI.set(model, MOI.VariableName(), x[1], "x1")
        @test MOI.get(model, MOI.VariableIndex, "x1") == x[1]
        MOI.set(model, MOI.VariableName(), x[1], "x2")
        @test MOI.get(model, MOI.VariableIndex, "x1") === nothing
        @test MOI.get(model, MOI.VariableIndex, "x2") == x[1]
        MOI.set(model, MOI.VariableName(), x[2], "x1")
        @test MOI.get(model, MOI.VariableIndex, "x1") == x[2]
        MOI.set(model, MOI.VariableName(), x[1], "x1")
        @test_throws ErrorException MOI.get(model, MOI.VariableIndex, "x1")
    end

    @testset "Variable bounds" begin
        MOI.empty!(model)
        x = MOI.add_variable(model)
        c1 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(0.0))
        c2 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.LessThan(1.0))
        MOI.set(model, MOI.ConstraintName(), c1, "c1")
        @test MOI.get(model, MOI.ConstraintIndex, "c1") == c1
        MOI.set(model, MOI.ConstraintName(), c1, "c2")
        @test MOI.get(model, MOI.ConstraintIndex, "c1") === nothing
        @test MOI.get(model, MOI.ConstraintIndex, "c2") == c1
        MOI.set(model, MOI.ConstraintName(), c2, "c1")
        @test MOI.get(model, MOI.ConstraintIndex, "c1") == c2
        MOI.set(model, MOI.ConstraintName(), c1, "c1")
        @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "c1")
    end

    @testset "Affine constraints" begin
        MOI.empty!(model)
        x = MOI.add_variable(model)
        f = MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0)
        c1 = MOI.add_constraint(model, f, MOI.GreaterThan(0.0))
        c2 = MOI.add_constraint(model, f, MOI.LessThan(1.0))
        MOI.set(model, MOI.ConstraintName(), c1, "c1")
        @test MOI.get(model, MOI.ConstraintIndex, "c1") == c1
        MOI.set(model, MOI.ConstraintName(), c1, "c2")
        @test MOI.get(model, MOI.ConstraintIndex, "c1") === nothing
        @test MOI.get(model, MOI.ConstraintIndex, "c2") == c1
        MOI.set(model, MOI.ConstraintName(), c2, "c1")
        @test MOI.get(model, MOI.ConstraintIndex, "c1") == c2
        MOI.set(model, MOI.ConstraintName(), c1, "c1")
        @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "c1")
    end
end

@testset "Issue #102" begin
    model = GLPK.Optimizer()
    x = MOI.add_variable(model)
    MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(0.0))
    MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Integer())
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 3.0)
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) == 3.0
    @test MOI.get(model, MOI.ObjectiveBound()) == 3.0
end

@testset "Issue #116" begin
    model = GLPK.Optimizer(method = GLPK.EXACT)
    x = MOI.add_variables(model, 2)
    c1 = MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(1.0, x[1]), MOI.ScalarAffineTerm(1.0, x[2])],
            0.0
        ),
        MOI.EqualTo(1.0)
    )
    MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.EqualTo(1.0))
    MOI.add_constraint(model, MOI.SingleVariable(x[2]), MOI.EqualTo(1.0))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
    @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    @test MOI.get(model, MOI.ConstraintDual(), c1) == -1
end

@testset "Default parameters" begin
    model = GLPK.Optimizer()
    @test MOI.get(model, MOI.RawParameter("msg_lev")) == GLPK.MSG_ERR
    @test MOI.get(model, MOI.RawParameter("presolve")) == GLPK.OFF
    model = GLPK.Optimizer(msg_lev = GLPK.MSG_ALL, presolve = true)
    @test MOI.get(model, MOI.RawParameter("msg_lev")) == GLPK.MSG_ALL
    @test MOI.get(model, MOI.RawParameter("presolve")) == GLPK.ON
end

@testset "Duplicate names" begin
    @testset "Variables" begin
        model = GLPK.Optimizer()
        (x, y, z) = MOI.add_variables(model, 3)
        MOI.set(model, MOI.VariableName(), x, "x")
        MOI.set(model, MOI.VariableName(), y, "x")
        MOI.set(model, MOI.VariableName(), z, "z")
        @test MOI.get(model, MOI.VariableIndex, "z") == z
        @test_throws ErrorException MOI.get(model, MOI.VariableIndex, "x")
        MOI.set(model, MOI.VariableName(), y, "y")
        @test MOI.get(model, MOI.VariableIndex, "x") == x
        @test MOI.get(model, MOI.VariableIndex, "y") == y
        MOI.set(model, MOI.VariableName(), z, "x")
        @test_throws ErrorException MOI.get(model, MOI.VariableIndex, "x")
        MOI.delete(model, x)
        @test MOI.get(model, MOI.VariableIndex, "x") == z
    end
    @testset "SingleVariable" begin
        model = GLPK.Optimizer()
        x = MOI.add_variables(model, 3)
        c = MOI.add_constraints(model, MOI.SingleVariable.(x), MOI.GreaterThan(0.0))
        MOI.set(model, MOI.ConstraintName(), c[1], "x")
        MOI.set(model, MOI.ConstraintName(), c[2], "x")
        MOI.set(model, MOI.ConstraintName(), c[3], "z")
        @test MOI.get(model, MOI.ConstraintIndex, "z") == c[3]
        @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
        MOI.set(model, MOI.ConstraintName(), c[2], "y")
        @test MOI.get(model, MOI.ConstraintIndex, "x") == c[1]
        @test MOI.get(model, MOI.ConstraintIndex, "y") == c[2]
        MOI.set(model, MOI.ConstraintName(), c[3], "x")
        @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
        MOI.delete(model, c[1])
        @test MOI.get(model, MOI.ConstraintIndex, "x") == c[3]
        MOI.set(model, MOI.ConstraintName(), c[2], "x")
        @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
        MOI.delete(model, x[3])
        @test MOI.get(model, MOI.ConstraintIndex, "x") == c[2]
    end
    @testset "ScalarAffineFunction" begin
        model = GLPK.Optimizer()
        x = MOI.add_variables(model, 3)
        fs = [
            MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, xi)], 0.0)
            for xi in x
        ]
        c = MOI.add_constraints(model, fs, MOI.GreaterThan(0.0))
        MOI.set(model, MOI.ConstraintName(), c[1], "x")
        MOI.set(model, MOI.ConstraintName(), c[2], "x")
        MOI.set(model, MOI.ConstraintName(), c[3], "z")
        @test MOI.get(model, MOI.ConstraintIndex, "z") == c[3]
        @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
        MOI.set(model, MOI.ConstraintName(), c[2], "y")
        @test MOI.get(model, MOI.ConstraintIndex, "x") == c[1]
        @test MOI.get(model, MOI.ConstraintIndex, "y") == c[2]
        MOI.set(model, MOI.ConstraintName(), c[3], "x")
        @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
        MOI.delete(model, c[1])
        @test MOI.get(model, MOI.ConstraintIndex, "x") == c[3]
    end
end

@testset "Duals with equal bounds" begin
    model = GLPK.Optimizer()
    x = MOI.add_variable(model)
    xl = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(1.0))
    xu = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.LessThan(1.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(x))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ConstraintDual(), xl) == 1.0
    @test MOI.get(model, MOI.ConstraintDual(), xu) == 0.0
end

# TODO move to MOI
@testset "PR #121" begin
    model = GLPK.Optimizer()
    ci = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}(1)
    @test_throws MOI.InvalidIndex(ci) MOI.get(model, MOI.ConstraintSet(), ci)
    @test_throws MOI.InvalidIndex(ci) MOI.get(model, MOI.ConstraintFunction(), ci)
    @test_throws MOI.InvalidIndex(ci) MOI.delete(model, ci)
end

@testset "Non-ascii names" begin
    model = GLPK.Optimizer()
    x = MOI.add_variable(model)
    MOI.set(model, MOI.VariableName(), x, "ω")
    @test MOI.get(model, MOI.VariableName(), x) == "ω"
    c = MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
        MOI.GreaterThan(0.0)
    )
    MOI.set(model, MOI.ConstraintName(), c, "ω")
    @test MOI.get(model, MOI.ConstraintName(), c) == "ω"
end
