using JuMP, GLPK, Test

# This special macro loads a module called GLPKExtensions, which exports the
# `@lazy_constraint` macro.
GLPK.@load_extensions

model = JuMP.direct_model(GLPK.Optimizer())

# Create the basic model. It is set up in such a way that the root relaxation
# is fractional.
@variable(model, 0 <= x <= 2, Int)
@variable(model, 0 <= y <= 4, Int)
@constraint(model, y <= 3.5 + x)
@constraint(model, y <= 4.1 - 0.2x)
@objective(model, Max, y)

# This vector is going to cache the value of `cb_where` everytime our callback
# gets called.
cb_calls = Int32[]

# Here is the callback function. We set the solver-specific attribute
# `Gurobi.CallbackFunction()`, passing a function (in this case, anonymous) as
# the third argument.
MOI.set(model, GLPK.CallbackFunction(), (cb_data) -> begin
    reason = GLPK.ios_reason(cb_data.tree)
    # Cache the value of `cb_where` in `cb_calls`. Note how we can access
    # variables in the outer scope.
    push!(cb_calls, reason)
    if reason == GLPK.IROWGEN
        GLPK.get_col_prim(cb_data)
        # Double check for sanity's sake that we have a feasible (given the
        # current constraints) point.
        @assert JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        # Get the values of x and y.
        x_val, y_val = JuMP.value(x), JuMP.value(y)
        # Add the lazy constraints using the `@lazy_constraint` macro that was
        # loaded by `GLPK.@load_extensions`.
        if y_val - x_val > 1.1 + 1e-6
            @lazy_constraint(cb_data, y <= 1.1 + x)
        elseif y_val + x_val > 3 + 1e-6
            @lazy_constraint(cb_data, y <= 3 - x)
        end
    elseif reason == GLPK.IHEUR
        # GLPK has a fractional solution to the current problem.
        # Load the solution into the model.
        GLPK.get_col_prim(cb_data)
        # Double check for sanity's sake that we have a feasible (given the
        # current constraints) point.
        @assert JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        # Get the values of x and y.
        x_val, y_val = JuMP.value(x), JuMP.value(y)
        # Provide a heuristic solution. We don't need to provide a value for all
        # variables.
        GLPK.ios_heur_sol(cb_data, Dict(x => 1.0, y => 2.0))
    end
    return
end)

# Solve the model.
JuMP.optimize!(model)

# Check the solution.
@test JuMP.value(x) == 1
@test JuMP.value(y) == 2

# Check that our callback has been used.
@test length(cb_calls) > 0
@test GLPK.ISELECT in cb_calls
@test GLPK.IPREPRO in cb_calls
@test GLPK.IROWGEN in cb_calls
@test GLPK.IHEUR in cb_calls
