module GLPKExtensions

using GLPK, JuMP

export @lazy_constraint
export set_callbacks, add_lazy_constraint, add_heuristic_solution

"""
    @lazy_constraint(model, cb_data, expr)

Add a lazy constraint (where the constraint is given by `expr`) to `model`.

This can be called *only* from a lazy callback set by `JuMP.set_callbacks`.

### Examples

    @lazy_constraint(model, cb_data, 2x + y <+ 1)
"""
macro lazy_constraint(model, cb_data, expr)
    code = quote
        lazy_con = @build_constraint $expr
        add_lazy_constraint(
            $model,
            $cb_data,
            JuMP.moi_function(lazy_con.func),
            lazy_con.set
        )

    end
    quote
        let
            $(esc(code))
        end
    end
end

"""
    set_callbacks(model::Model; lazy = nothing, heuristic = nothing)

Set an (optional) lazy and heuristic callback. This can be used only when
`JuMP.mode(mode) == DIRECT`.

Each callback must be a function that takes two arguments: the backend
`optimizer` object (returned by `JuMP.backend(model)`), and a solver-specific
type `cb_data`.

### Example

    set_callbacks(model;
        lazy = (optimizer, cb_data) -> println("Called from lazy callback.")
    )
"""
function set_callbacks(model::Model; lazy = nothing, heuristic = nothing)
    if JuMP.mode(model) != JuMP.DIRECT
        error("You must use a solver in DIRECT mode to use callbacks in JuMP.")
    end
    _set_callbacks(model, JuMP.backend(model), lazy, heuristic)
end

function _set_callbacks(model, optimizer, lazy, heuristic)
    error("The model $(typeof(optimizer)) does not support callbacks.")
end

"""
    add_lazy_constraint(model, cb_data, func, set)

Add a lazy constraint `func`-in-`set` to `model`.

This can be called only from a lazy callback set by `set_callbacks`.
"""
function add_lazy_constraint end


"""
    add_heuristic_solution(model, cb_data, sol::Dict{JuMP.VariableRef, Float64})

Provide the heuristic solution given by the variable-value mapping of `sol` to
`model`.

This can be called only from a heuristic callback set by `set_callbacks`.
"""
function add_heuristic_solution end

###
### GLPK specific functionality.
###

function _set_callbacks(model::Model, optimizer::GLPK.Optimizer,
                        lazy_callback, heuristic_callback)
    MOI.set(optimizer, GLPK.CallbackFunction(), (cb_data) -> begin
        reason = GLPK.ios_reason(cb_data.tree)
        if lazy_callback !== nothing && reason == GLPK.IROWGEN
            GLPK.get_col_prim(cb_data)
            lazy_callback(model, cb_data)
        elseif heuristic_callback !== nothing && reason == GLPK.IHEUR
            GLPK.get_col_prim(cb_data)
            heuristic_callback(model, cb_data)
        end
        return
    end)
end

function add_lazy_constraint(model::Model, cb_data::GLPK.CallbackData, func, set)
    GLPK.addrow!(cb_data, func, set)
end

function add_heuristic_solution(model::Model, cb_data::GLPK.CallbackData,
                                sol::Dict{JuMP.VariableRef, Float64})
    GLPK.ios_heur_sol(cb_data, Dict(
        JuMP.index(variable) => value for (variable, value) in sol))
end

end
