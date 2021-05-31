module TestCallbacks

using GLPK, Test, Random

const MOI = GLPK.MOI

function runtests()
    for name in names(@__MODULE__; all = true)
        if !startswith("$(name)", "test_")
            continue
        elseif startswith("$(name)", "test_no_cache_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        else
            @testset "$(name)" begin
                getfield(@__MODULE__, name)(true)
                getfield(@__MODULE__, name)(false)
            end
        end
    end
    return
end

function _callback_simple_model(cache)
    model = if cache
        m = MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            GLPK.Optimizer(),
        )
        MOI.Utilities.reset_optimizer(m)
        m
    else
        GLPK.Optimizer()
    end
    MOI.set(model, MOI.Silent(), true)
    MOI.Utilities.loadfromstring!(
        model,
        """
    variables: x, y
    maxobjective: 1.0 * y
    c1: x in Integer()
    c2: y in Integer()
    c3: x in Interval(0.0, 2.5)
    c4: y in Interval(0.0, 2.5)
""",
    )
    x = MOI.get(model, MOI.VariableIndex, "x")
    y = MOI.get(model, MOI.VariableIndex, "y")
    return model, x, y
end

function _callback_knapsack_model(cache)
    model = if cache
        m = MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            GLPK.Optimizer(),
        )
        MOI.Utilities.reset_optimizer(m)
        m
    else
        GLPK.Optimizer()
    end
    MOI.set(model, MOI.Silent(), true)
    N = 30
    x = MOI.add_variables(model, N)
    MOI.add_constraints(model, MOI.SingleVariable.(x), MOI.ZeroOne())
    Random.seed!(1)
    item_weights, item_values = rand(N), rand(N)
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0),
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_values, x), 0.0),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    return model, x, item_weights
end

function test_lazy_constraint(cache)
    model, x, y = _callback_simple_model(cache)
    lazy_called = false
    MOI.set(
        model,
        MOI.LazyConstraintCallback(),
        (cb_data) -> begin
            lazy_called = true
            x_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), x)
            y_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), y)
            @test MOI.supports(model, MOI.LazyConstraint(cb_data))
            status = MOI.get(model, MOI.CallbackNodeStatus(cb_data))
            if [x_val, y_val] ≈ round.(Int, [x_val, y_val])
                atol = 1e-7
                @test status == MOI.CALLBACK_NODE_STATUS_INTEGER
            else
                @test status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            end
            if y_val - x_val > 1 + 1e-6
                MOI.submit(
                    model,
                    MOI.LazyConstraint(cb_data),
                    MOI.ScalarAffineFunction{Float64}(
                        MOI.ScalarAffineTerm.([-1.0, 1.0], [x, y]),
                        0.0,
                    ),
                    MOI.LessThan{Float64}(1.0),
                )
            elseif y_val + x_val > 3 + 1e-6
                MOI.submit(
                    model,
                    MOI.LazyConstraint(cb_data),
                    MOI.ScalarAffineFunction{Float64}(
                        MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]),
                        0.0,
                    ),
                    MOI.LessThan{Float64}(3.0),
                )
            end
        end,
    )
    @test MOI.supports(model, MOI.LazyConstraintCallback())
    MOI.optimize!(model)
    @test lazy_called
    @test MOI.get(model, MOI.VariablePrimal(), x) == 1
    @test MOI.get(model, MOI.VariablePrimal(), y) == 2
end

function test_lazy_constraint_optimize_in_progress(cache)
    model, x, y = _callback_simple_model(cache)
    MOI.set(
        model,
        MOI.LazyConstraintCallback(),
        (cb_data) -> begin
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
        end,
    )
    return MOI.optimize!(model)
end

function test_no_cache_LazyConstraint_UserCut()
    model, x, y = _callback_simple_model(false)
    MOI.set(
        model,
        MOI.LazyConstraintCallback(),
        (cb_data) -> begin
            MOI.submit(
                model,
                MOI.UserCut(cb_data),
                MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
                MOI.LessThan(2.0),
            )
        end,
    )
    @test_throws(
        MOI.InvalidCallbackUsage(
            MOI.LazyConstraintCallback(),
            MOI.UserCut(model.callback_data),
        ),
        MOI.optimize!(model)
    )
end

function test_no_cache_LazyConstraint_HeuristicSolution()
    model, x, y = _callback_simple_model(false)
    MOI.set(
        model,
        MOI.LazyConstraintCallback(),
        (cb_data) -> begin
            MOI.submit(model, MOI.HeuristicSolution(cb_data), [x, y], [1.0, 2.0])
        end,
    )
    @test_throws(
        MOI.InvalidCallbackUsage(
            MOI.LazyConstraintCallback(),
            MOI.HeuristicSolution(model.callback_data),
        ),
        MOI.optimize!(model)
    )
end

function test_UserCut(cache)
    model, x, item_weights = _callback_knapsack_model(cache)
    user_cut_submitted = false
    MOI.set(
        model,
        MOI.UserCutCallback(),
        (cb_data) -> begin
            terms = MOI.ScalarAffineTerm{Float64}[]
            accumulated = 0.0
            for (i, xi) in enumerate(x)
                if MOI.get(model, MOI.CallbackVariablePrimal(cb_data), xi) > 0.0
                    push!(terms, MOI.ScalarAffineTerm(1.0, xi))
                    accumulated += item_weights[i]
                end
            end
            @test MOI.supports(model, MOI.UserCut(cb_data))
            status = MOI.get(model, MOI.CallbackNodeStatus(cb_data))
            @test status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            if accumulated > 10.0
                MOI.submit(
                    model,
                    MOI.UserCut(cb_data),
                    MOI.ScalarAffineFunction{Float64}(terms, 0.0),
                    MOI.LessThan{Float64}(length(terms) - 1),
                )
                user_cut_submitted = true
            end
        end,
    )
    @test MOI.supports(model, MOI.UserCutCallback())
    MOI.optimize!(model)
    @test user_cut_submitted
end

function test_no_cache_UserCut_LazyConstraint()
    model, x, item_weights = _callback_knapsack_model(false)
    MOI.set(
        model,
        MOI.UserCutCallback(),
        (cb_data) -> begin
            MOI.submit(
                model,
                MOI.LazyConstraint(cb_data),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0),
                MOI.LessThan(5.0),
            )
        end,
    )
    @test_throws(
        MOI.InvalidCallbackUsage(
            MOI.UserCutCallback(),
            MOI.LazyConstraint(model.callback_data),
        ),
        MOI.optimize!(model)
    )
end

function test_no_cache_UserCut_HeuristicSolution()
    model, x, item_weights = _callback_knapsack_model(false)
    MOI.set(
        model,
        MOI.UserCutCallback(),
        (cb_data) -> begin
            MOI.submit(
                model,
                MOI.HeuristicSolution(cb_data),
                x,
                fill(0.0, length(x)),
            )
        end,
    )
    @test_throws(
        MOI.InvalidCallbackUsage(
            MOI.UserCutCallback(),
            MOI.HeuristicSolution(model.callback_data),
        ),
        MOI.optimize!(model)
    )
end

function test_Heuristic(cache)
    model, x, item_weights = _callback_knapsack_model(cache)
    solution_accepted = false
    solution_rejected = false
    MOI.set(
        model,
        MOI.HeuristicCallback(),
        (cb_data) -> begin
            x_vals = MOI.get.(model, MOI.CallbackVariablePrimal(cb_data), x)
            @test MOI.supports(model, MOI.HeuristicSolution(cb_data))
            if MOI.submit(
                model,
                MOI.HeuristicSolution(cb_data),
                x,
                floor.(x_vals),
            ) == MOI.HEURISTIC_SOLUTION_ACCEPTED
                solution_accepted = true
            end
            if MOI.submit(
                model,
                MOI.HeuristicSolution(cb_data),
                [x[1]],
                [1.0],
            ) == MOI.HEURISTIC_SOLUTION_REJECTED
                solution_rejected = true
            end
            status = MOI.get(model, MOI.CallbackNodeStatus(cb_data))
            @test status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        end,
    )
    @test MOI.supports(model, MOI.HeuristicCallback())
    MOI.optimize!(model)
    @test solution_accepted
    @test solution_rejected
end

function test_no_cache_Heuristic_LazyConstraint()
    model, x, item_weights = _callback_knapsack_model(false)
    MOI.set(
        model,
        MOI.HeuristicCallback(),
        (cb_data) -> begin
            MOI.submit(
                model,
                MOI.LazyConstraint(cb_data),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0),
                MOI.LessThan(5.0),
            )
        end,
    )
    @test_throws(
        MOI.InvalidCallbackUsage(
            MOI.HeuristicCallback(),
            MOI.LazyConstraint(model.callback_data),
        ),
        MOI.optimize!(model)
    )
end

function test_no_cache_Heuristic_UserCut()
    model, x, item_weights = _callback_knapsack_model(false)
    MOI.set(
        model,
        MOI.HeuristicCallback(),
        (cb_data) -> begin
            MOI.submit(
                model,
                MOI.UserCut(cb_data),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0),
                MOI.LessThan(5.0),
            )
        end,
    )
    @test_throws(
        MOI.InvalidCallbackUsage(
            MOI.HeuristicCallback(),
            MOI.UserCut(model.callback_data),
        ),
        MOI.optimize!(model)
    )
end

function test_CallbackFunction_OptimizeInProgress(cache)
    model, x, y = _callback_simple_model(cache)
    MOI.set(
        model,
        GLPK.CallbackFunction(),
        (cb_data) -> begin
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
        end,
    )
    @test MOI.supports(model, GLPK.CallbackFunction())
    return MOI.optimize!(model)
end

function test_CallbackFunction_LazyConstraint(cache)
    model, x, y = _callback_simple_model(cache)
    cb_calls = Int32[]
    function callback_function(cb_data)
        reason = GLPK.glp_ios_reason(cb_data.tree)
        push!(cb_calls, reason)
        if reason != GLPK.GLP_IROWGEN
            return
        end
        x_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), x)
        y_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), y)
        status = MOI.get(model, MOI.CallbackNodeStatus(cb_data))
        if [x_val, y_val] ≈ round.(Int, [x_val, y_val])
            atol = 1e-7
            @test status == MOI.CALLBACK_NODE_STATUS_INTEGER
        else
            @test status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        end
        if y_val - x_val > 1 + 1e-6
            MOI.submit(
                model,
                MOI.LazyConstraint(cb_data),
                MOI.ScalarAffineFunction{Float64}(
                    MOI.ScalarAffineTerm.([-1.0, 1.0], [x, y]),
                    0.0,
                ),
                MOI.LessThan{Float64}(1.0),
            )
        elseif y_val + x_val > 3 + 1e-6
            MOI.submit(
                model,
                MOI.LazyConstraint(cb_data),
                MOI.ScalarAffineFunction{Float64}(
                    MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]),
                    0.0,
                ),
                MOI.LessThan{Float64}(3.0),
            )
        end
    end
    MOI.set(model, GLPK.CallbackFunction(), callback_function)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), x) == 1
    @test MOI.get(model, MOI.VariablePrimal(), y) == 2
    @test length(cb_calls) > 0
    @test GLPK.GLP_ISELECT in cb_calls
    @test GLPK.GLP_IPREPRO in cb_calls
    @test GLPK.GLP_IROWGEN in cb_calls
    @test GLPK.GLP_IBINGO in cb_calls
    @test !(GLPK.GLP_IHEUR in cb_calls)
end

function test_CallbackFunction_UserCut(cache)
    model, x, item_weights = _callback_knapsack_model(cache)
    user_cut_submitted = false
    cb_calls = Int32[]
    MOI.set(
        model,
        GLPK.CallbackFunction(),
        (cb_data) -> begin
            reason = GLPK.glp_ios_reason(cb_data.tree)
            push!(cb_calls, reason)
            if reason != GLPK.GLP_ICUTGEN
                status = MOI.get(model, MOI.CallbackNodeStatus(cb_data))
                @test status !== nothing
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
                    MOI.LessThan{Float64}(length(terms) - 1),
                )
                user_cut_submitted = true
            end
            status = MOI.get(model, MOI.CallbackNodeStatus(cb_data))
            @test status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        end,
    )
    MOI.optimize!(model)
    @test user_cut_submitted
    @test GLPK.GLP_ICUTGEN in cb_calls
end

function test_CallbackFunction_HeuristicSolution(cache)
    model, x, item_weights = _callback_knapsack_model(cache)
    solution_accepted = false
    solution_rejected = false
    cb_calls = Int32[]
    MOI.set(
        model,
        GLPK.CallbackFunction(),
        (cb_data) -> begin
            reason = GLPK.glp_ios_reason(cb_data.tree)
            push!(cb_calls, reason)
            if reason != GLPK.GLP_IHEUR
                return
            end
            x_vals =
                MOI.get.(model, MOI.CallbackVariablePrimal(cb_data), x)
            if MOI.submit(
                model,
                MOI.HeuristicSolution(cb_data),
                x,
                floor.(x_vals),
            ) == MOI.HEURISTIC_SOLUTION_ACCEPTED
                solution_accepted = true
            end
            if MOI.submit(
                model,
                MOI.HeuristicSolution(cb_data),
                [x[1]],
                [1.0],
            ) == MOI.HEURISTIC_SOLUTION_REJECTED
                solution_rejected = true
            end
            status = MOI.get(model, MOI.CallbackNodeStatus(cb_data))
            @test status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        end,
    )
    MOI.optimize!(model)
    @test solution_accepted
    @test solution_rejected
    @test GLPK.GLP_IHEUR in cb_calls
end

function test_broadcasting(cache)
    model, x, _ = _callback_knapsack_model(cache)
    f(cb_data, x) = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), x)
    solutions = Vector{Float64}[]
    MOI.set(
        model,
        GLPK.CallbackFunction(),
        (cb_data) -> begin
            if GLPK.glp_ios_reason(cb_data.tree) == GLPK.GLP_IHEUR
                push!(solutions, f.(cb_data, x))
            end
        end,
    )
    MOI.optimize!(model)
    @test length(solutions) > 0
    @test length(solutions[1]) == length(x)
end

end

TestCallbacks.runtests()
