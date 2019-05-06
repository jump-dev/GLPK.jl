using LinQuadOptInterface

const MOI  = LinQuadOptInterface.MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

@testset "Unit Tests" begin
    config = MOIT.TestConfig()
    solver = GLPK.Optimizer()

    MOIT.basic_constraint_tests(solver, config)

    MOIT.unittest(solver, config, [
        # These are excluded because GLPK does not support quadratics.
        "solve_qp_edge_cases",
        "solve_qcp_edge_cases"
    ])

    MOIT.modificationtest(solver, config)
end

@testset "Linear tests" begin
    solver = GLPK.Optimizer()
    MOIT.contlineartest(solver, MOIT.TestConfig(), [
        # GLPK returns InfeasibleOrUnbounded
        "linear8a",
        # Requires infeasiblity certificate for variable bounds
        "linear12",
        # FIXME
        "partial_start"
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
    @test MOI.get(solver, MOI.SolverName()) == "GLPK"
    @testset "default_objective_test" begin
         MOIT.default_objective_test(solver)
     end
     @testset "default_status_test" begin
         MOIT.default_status_test(solver)
     end
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
        MOIT.orderedindicestest(solver)
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
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE_OR_UNBOUNDED
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

@testset "Callbacks" begin
    @testset "Lazy cut" begin
        model = GLPK.Optimizer()
        MOI.Utilities.loadfromstring!(model, """
            variables: x, y
            maxobjective: y
            c1: x in Integer()
            c2: y in Integer()
            c3: x in Interval(0.0, 2.0)
            c4: y in Interval(0.0, 2.0)
        """)
        x = MOI.get(model, MOI.VariableIndex, "x")
        y = MOI.get(model, MOI.VariableIndex, "y")

        # We now define our callback function that takes the callback handle.
        # Note that we can access model, x, and y because this function is
        # defined inside the same scope.
        cb_calls = Int32[]
        function callback_function(cb_data::GLPK.CallbackData)
            reason = GLPK.ios_reason(cb_data.tree)
            push!(cb_calls, reason)
            if reason == GLPK.IROWGEN
                GLPK.load_variable_primal!(cb_data)
                x_val = MOI.get(model, MOI.VariablePrimal(), x)
                y_val = MOI.get(model, MOI.VariablePrimal(), y)
                # We have two constraints, one cutting off the top
                # left corner and one cutting off the top right corner, e.g.
                # (0,2) +---+---+ (2,2)
                #       |xx/ \xx|
                #       |x/   \x|
                #       |/     \|
                # (0,1) +   +   + (2,1)
                #       |       |
                # (0,0) +---+---+ (2,0)
                TOL = 1e-6  # Allow for some impreciseness in the solution
                if y_val - x_val > 1 + TOL
                    GLPK.add_lazy_constraint!(cb_data,
                        MOI.ScalarAffineFunction{Float64}(
                            MOI.ScalarAffineTerm.([-1.0, 1.0], [x, y]),
                            0.0
                        ),
                        MOI.LessThan{Float64}(1.0)
                    )
                elseif y_val + x_val > 3 + TOL
                    GLPK.add_lazy_constraint!(cb_data,
                        MOI.ScalarAffineFunction{Float64}(
                            MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]),
                            0.0
                        ),
                        MOI.LessThan{Float64}(3.0)
                    )
                end
            end
        end
        MOI.set(model, GLPK.CallbackFunction(), callback_function)
        MOI.optimize!(model)
        @test MOI.get(model, MOI.VariablePrimal(), x) == 1
        @test MOI.get(model, MOI.VariablePrimal(), y) == 2
        @test length(cb_calls) > 0
        @test GLPK.ISELECT in cb_calls
        @test GLPK.IPREPRO in cb_calls
        @test GLPK.IROWGEN in cb_calls
        @test GLPK.IBINGO in cb_calls
        @test !(GLPK.IHEUR in cb_calls)
    end
end

@testset "Constant objectives" begin
    # Test that setting the objective constant actually propagates
    # to the solver, rather than being cached in the LinQuadOptInterface
    # layer.
    model = GLPK.Optimizer()
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 1.5))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) == 1.5
    @test GLPK.get_obj_val(model.inner) == 1.5

    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], -2.0))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) == -2.0
    @test GLPK.get_obj_val(model.inner) == -2.0
end

@testset "Issue #70" begin
    model = GLPK.Optimizer()
    x = MOI.add_variable(model)
    f =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0)
    s = MOI.LessThan(2.0)
    c = MOI.add_constraint(model, f, s)
    row = model[c]
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
    @testset "SingleVariable" begin
        model = GLPK.Optimizer()
        x = MOI.add_variable(model)
        MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Interval(1.0, -1.0))
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.INVALID_MODEL
    end
    @testset "ScalarAffine" begin
        model = GLPK.Optimizer()
        x = MOI.add_variable(model)
        MOI.add_constraint(model,
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0),
            MOI.Interval(1.0, -1.0))
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.INVALID_MODEL
    end
end
