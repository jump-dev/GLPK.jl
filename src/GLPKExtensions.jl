module GLPKExtensions

using GLPK, JuMP

export @lazy_constraint

macro lazy_constraint(cb_data, expr)
    code = quote
        lazy_con = @build_constraint $expr
        GLPK.addrow!(
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

function GLPK.ios_heur_sol(cb_data::GLPK.CallbackData,
                           sol::Dict{JuMP.VariableRef, Float64})
    GLPK.ios_heur_sol(
        cb_data,
        Dict(JuMP.index(variable) => value for (variable, value) in sol)
    )
    return
end


end
