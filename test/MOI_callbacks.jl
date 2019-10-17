using GLPK, Test, Random

const MOI = GLPK.MOI

function callback_simple_model()
    model = GLPK.Optimizer()
    MOI.Utilities.loadfromstring!(model, """
        variables: x, y
        maxobjective: y
        c1: x in Integer()
        c2: y in Integer()
        c3: x in Interval(0.0, 2.5)
        c4: y in Interval(0.0, 2.5)
    """)
    x = MOI.get(model, MOI.VariableIndex, "x")
    y = MOI.get(model, MOI.VariableIndex, "y")
    return model, x, y
end

function callback_knapsack_model()
    model = GLPK.Optimizer()
    N = 30
    x = MOI.add_variables(model, N)
    MOI.add_constraints(model, MOI.SingleVariable.(x), MOI.ZeroOne())
    Random.seed!(1)
    item_weights, item_values = rand(N), rand(N)
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0)
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_values, x), 0.0)
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    return model, x, item_weights
end

@testset "All Callbacks" begin
    @testset "LazyConstraintCallback" begin
        @testset "LazyConstraint" begin
            model, x, y = callback_simple_model()
            lazy_called = false
            MOI.set(model, MOI.LazyConstraintCallback(), (cb_data) -> begin
                lazy_called = true
                x_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), x)
                y_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), y)
                if y_val - x_val > 1 + 1e-6
                    MOI.submit(
                        model,
                        MOI.LazyConstraint(cb_data),
                        MOI.ScalarAffineFunction{Float64}(
                            MOI.ScalarAffineTerm.([-1.0, 1.0], [x, y]),
                            0.0
                        ),
                        MOI.LessThan{Float64}(1.0)
                    )
                elseif y_val + x_val > 3 + 1e-6
                    MOI.submit(
                        model,
                        MOI.LazyConstraint(cb_data),
                        MOI.ScalarAffineFunction{Float64}(
                            MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]),
                            0.0
                        ), MOI.LessThan{Float64}(3.0)
                    )
                end
            end)
            MOI.optimize!(model)
            @test lazy_called
            @test MOI.get(model, MOI.VariablePrimal(), x) == 1
            @test MOI.get(model, MOI.VariablePrimal(), y) == 2
        end
        @testset "OptimizeInProgress" begin
            model, x, y = callback_simple_model()
            MOI.set(model, MOI.LazyConstraintCallback(), (cb_data) -> begin
                @test_throws(
                    MOI.OptimizeInProgress(MOI.VariablePrimal()),
                    MOI.get(model, MOI.VariablePrimal(), x)
                )
                @test_throws(
                    MOI.OptimizeInProgress(MOI.ObjectiveValue()),
                    MOI.get(model, MOI.ObjectiveValue())
                )
                @test_throws(
                    MOI.OptimizeInProgress(MOI.ObjectiveBound()),
                    MOI.get(model, MOI.ObjectiveBound())
                )
            end)
            MOI.optimize!(model)
        end
        @testset "UserCut" begin
            model, x, y = callback_simple_model()
            MOI.set(model, MOI.LazyConstraintCallback(), (cb_data) -> begin
                MOI.submit(
                    model,
                    MOI.UserCut(cb_data),
                    MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
                    MOI.LessThan(2.0)
                )
            end)
            @test_throws(
                MOI.InvalidCallbackUsage(
                    MOI.LazyConstraintCallback(),
                    MOI.UserCut(model.callback_data)
                ),
                MOI.optimize!(model)
            )
        end
        @testset "HeuristicSolution" begin
            model, x, y = callback_simple_model()
            MOI.set(model, MOI.LazyConstraintCallback(), (cb_data) -> begin
                MOI.submit(
                    model,
                    MOI.HeuristicSolution(cb_data),
                    [x, y],
                    [1.0, 2.0]
                )
            end)
            @test_throws(
                MOI.InvalidCallbackUsage(
                    MOI.LazyConstraintCallback(),
                    MOI.HeuristicSolution(model.callback_data)
                ),
                MOI.optimize!(model)
            )
        end
    end

    @testset "UserCutCallback" begin
        @testset "UserCut" begin
            model, x, item_weights = callback_knapsack_model()
            user_cut_submitted = false
            MOI.set(model, MOI.UserCutCallback(), (cb_data) -> begin
                terms = MOI.ScalarAffineTerm{Float64}[]
                accumulated = 0.0
                for (i, xi) in enumerate(x)
                    if MOI.get(model, MOI.CallbackVariablePrimal(cb_data), xi) > 0.0
                        push!(terms, MOI.ScalarAffineTerm(1.0, xi))
                        accumulated += item_weights[i]
                    end
                end
                if accumulated > 10.0
                    MOI.submit(
                        model,
                        MOI.UserCut(cb_data),
                        MOI.ScalarAffineFunction{Float64}(terms, 0.0),
                        MOI.LessThan{Float64}(length(terms) - 1)
                    )
                    user_cut_submitted = true
                end
            end)
            MOI.optimize!(model)
            @test user_cut_submitted
        end
        @testset "LazyConstraint" begin
            model, x, item_weights = callback_knapsack_model()
            MOI.set(model, MOI.UserCutCallback(), (cb_data) -> begin
                MOI.submit(
                    model,
                    MOI.LazyConstraint(cb_data),
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0),
                    MOI.LessThan(5.0)
                )
            end)
            @test_throws(
                MOI.InvalidCallbackUsage(
                    MOI.UserCutCallback(),
                    MOI.LazyConstraint(model.callback_data)
                ),
                MOI.optimize!(model)
            )
        end
        @testset "HeuristicSolution" begin
            model, x, item_weights = callback_knapsack_model()
            MOI.set(model, MOI.UserCutCallback(), (cb_data) -> begin
                MOI.submit(
                    model,
                    MOI.HeuristicSolution(cb_data),
                    x,
                    fill(0.0, length(x))
                )
            end)
            @test_throws(
                MOI.InvalidCallbackUsage(
                    MOI.UserCutCallback(),
                    MOI.HeuristicSolution(model.callback_data)
                ),
                MOI.optimize!(model)
            )
        end
    end

    @testset "HeuristicCallback" begin
        @testset "HeuristicSolution" begin
            model, x, item_weights = callback_knapsack_model()
            solution_accepted = false
            solution_rejected = false
            MOI.set(model, MOI.HeuristicCallback(), (cb_data) -> begin
                x_vals = MOI.get.(model, MOI.CallbackVariablePrimal(cb_data), x)
                if MOI.submit(
                    model,
                    MOI.HeuristicSolution(cb_data),
                    x,
                    floor.(x_vals)
                ) == MOI.HEURISTIC_SOLUTION_ACCEPTED
                    solution_accepted = true
                end
                if MOI.submit(
                    model,
                    MOI.HeuristicSolution(cb_data),
                    [x[1]],
                    [1.0]
                ) == MOI.HEURISTIC_SOLUTION_REJECTED
                    solution_rejected = true
                end
            end)
            MOI.optimize!(model)
            @test solution_accepted
            @test solution_rejected
        end
        @testset "LazyConstraint" begin
            model, x, item_weights = callback_knapsack_model()
            MOI.set(model, MOI.HeuristicCallback(), (cb_data) -> begin
                MOI.submit(
                    model,
                    MOI.LazyConstraint(cb_data),
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0),
                    MOI.LessThan(5.0)
                )
            end)
            @test_throws(
                MOI.InvalidCallbackUsage(
                    MOI.HeuristicCallback(),
                    MOI.LazyConstraint(model.callback_data)
                ),
                MOI.optimize!(model)
            )
        end
        @testset "UserCut" begin
            model, x, item_weights = callback_knapsack_model()
            MOI.set(model, MOI.HeuristicCallback(), (cb_data) -> begin
                MOI.submit(
                    model,
                    MOI.UserCut(cb_data),
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0),
                    MOI.LessThan(5.0)
                )
            end)
            @test_throws(
                MOI.InvalidCallbackUsage(
                    MOI.HeuristicCallback(),
                    MOI.UserCut(model.callback_data)
                ),
                MOI.optimize!(model)
            )
        end
    end

    @testset "GLPK.CallbackFunction" begin
        @testset "OptimizeInProgress" begin
            model, x, y = callback_simple_model()
            MOI.set(model, GLPK.CallbackFunction(), (cb_data) -> begin
                @test_throws(
                    MOI.OptimizeInProgress(MOI.VariablePrimal()),
                    MOI.get(model, MOI.VariablePrimal(), x)
                )
                @test_throws(
                    MOI.OptimizeInProgress(MOI.ObjectiveValue()),
                    MOI.get(model, MOI.ObjectiveValue())
                )
                @test_throws(
                    MOI.OptimizeInProgress(MOI.ObjectiveBound()),
                    MOI.get(model, MOI.ObjectiveBound())
                )
            end)
            MOI.optimize!(model)
        end
        @testset "LazyConstraint" begin
            model, x, y = callback_simple_model()
            cb_calls = Int32[]
            function callback_function(cb_data)
                reason = GLPK.ios_reason(cb_data.tree)
                push!(cb_calls, reason)
                if reason != GLPK.IROWGEN
                    return
                end
                x_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), x)
                y_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), y)
                if y_val - x_val > 1 + 1e-6
                    MOI.submit(model, MOI.LazyConstraint(cb_data),
                        MOI.ScalarAffineFunction{Float64}(
                            MOI.ScalarAffineTerm.([-1.0, 1.0], [x, y]),
                            0.0
                        ),
                        MOI.LessThan{Float64}(1.0)
                    )
                elseif y_val + x_val > 3 + 1e-6
                    MOI.submit(model, MOI.LazyConstraint(cb_data),
                        MOI.ScalarAffineFunction{Float64}(
                            MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]),
                            0.0
                        ),
                        MOI.LessThan{Float64}(3.0)
                    )
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
        @testset "UserCut" begin
            model, x, item_weights = callback_knapsack_model()
            user_cut_submitted = false
            cb_calls = Int32[]
            MOI.set(model, GLPK.CallbackFunction(), (cb_data) -> begin
                reason = GLPK.ios_reason(cb_data.tree)
                push!(cb_calls, reason)
                if reason != GLPK.ICUTGEN
                    return
                end
                terms = MOI.ScalarAffineTerm{Float64}[]
                accumulated = 0.0
                for (i, xi) in enumerate(x)
                    if MOI.get(model, MOI.CallbackVariablePrimal(cb_data), xi) > 0.0
                        push!(terms, MOI.ScalarAffineTerm(1.0, xi))
                        accumulated += item_weights[i]
                    end
                end
                if accumulated > 10.0
                    MOI.submit(
                        model,
                        MOI.UserCut(cb_data),
                        MOI.ScalarAffineFunction{Float64}(terms, 0.0),
                        MOI.LessThan{Float64}(length(terms) - 1)
                    )
                    user_cut_submitted = true
                end
            end)
            MOI.optimize!(model)
            @test user_cut_submitted
            @test GLPK.ICUTGEN in cb_calls
        end
        @testset "HeuristicSolution" begin
            model, x, item_weights = callback_knapsack_model()
            solution_accepted = false
            solution_rejected = false
            cb_calls = Int32[]
            MOI.set(model, GLPK.CallbackFunction(), (cb_data) -> begin
                reason = GLPK.ios_reason(cb_data.tree)
                push!(cb_calls, reason)
                if reason != GLPK.IHEUR
                    return
                end
                x_vals = MOI.get.(model, MOI.CallbackVariablePrimal(cb_data), x)
                if MOI.submit(
                    model,
                    MOI.HeuristicSolution(cb_data),
                    x,
                    floor.(x_vals)
                ) == MOI.HEURISTIC_SOLUTION_ACCEPTED
                    solution_accepted = true
                end
                if MOI.submit(
                    model,
                    MOI.HeuristicSolution(cb_data),
                    [x[1]],
                    [1.0]
                ) == MOI.HEURISTIC_SOLUTION_REJECTED
                    solution_rejected = true
                end
            end)
            MOI.optimize!(model)
            @test solution_accepted
            @test solution_rejected
            @test GLPK.IHEUR in cb_calls
        end
    end
end
