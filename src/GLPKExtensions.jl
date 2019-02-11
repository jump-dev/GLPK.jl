module GLPKJuMPExtensions

using GLPK, JuMP

function JuMP._set_callbacks(model::Model, optimizer::GLPK.Optimizer,
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

function JuMP.add_lazy_constraint(
        model::Model, cb_data::GLPK.CallbackData, func, set)
    GLPK.addrow!(cb_data, func, set)
end

function JuMP.add_heuristic_solution(model::Model, cb_data::GLPK.CallbackData,
                                     sol::Dict{JuMP.VariableRef, Float64})
    GLPK.ios_heur_sol(cb_data, Dict(
        JuMP.index(variable) => value for (variable, value) in sol))
end

end
