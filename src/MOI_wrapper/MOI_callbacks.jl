# ==============================================================================
#    Generic Callbacks in GLPK
# ==============================================================================

"""
    CallbackFunction()

Set a generic GLPK callback function.
"""
struct CallbackFunction <: MOI.AbstractCallback end

MOI.supports(::Optimizer, ::CallbackFunction) = true

function MOI.set(model::Optimizer, ::CallbackFunction, callback::Function)
    model.has_generic_callback = true
    _set_callback(model, (cb_data) -> begin
        model.callback_state = _CB_GENERIC
        callback(cb_data)
        model.callback_state = _CB_NONE
    end)
    return
end

# ==============================================================================
#    MOI callbacks
# ==============================================================================

function _default_moi_callback(model::Optimizer)
    return (cb_data) -> begin
        reason = glp_ios_reason(cb_data.tree)
        if reason == GLP_IROWGEN && model.lazy_callback !== nothing
            model.callback_state = _CB_LAZY
            model.lazy_callback(cb_data)
        elseif reason == GLP_ICUTGEN && model.user_cut_callback !== nothing
            model.callback_state = _CB_USER_CUT
            model.user_cut_callback(cb_data)
        elseif reason == GLP_IHEUR && model.heuristic_callback !== nothing
            model.callback_state = _CB_HEURISTIC
            model.heuristic_callback(cb_data)
        end
        return
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.CallbackVariablePrimal{CallbackData},
    x::MOI.VariableIndex,
)
    subproblem = glp_ios_get_prob(attr.callback_data.tree)
    return glp_get_col_prim(subproblem, _info(model, x).column)
end

function MOI.get(model::Optimizer, attr::MOI.CallbackNodeStatus{CallbackData})
    reason = glp_ios_reason(attr.callback_data.tree)
    if reason == GLP_ISELECT
        return MOI.CALLBACK_NODE_STATUS_UNKNOWN
    elseif reason == GLP_IPREPRO
        return MOI.CALLBACK_NODE_STATUS_UNKNOWN
    elseif reason == GLP_IROWGEN
        # From the the GLPK documentation:
        #
        # The callback routine is called with the reason code GLP_IROWGEN if LP
        # relaxation of the current subproblem has just been solved to
        # optimality and its objective value is better than the best known
        # integer feasible solution.
        #
        # This can mean the solution is integer _or_ fractional, so we need to
        # check.
        subproblem = glp_ios_get_prob(attr.callback_data.tree)
        for info in values(model.variable_info)
            if info.type == _CONTINUOUS
                continue
            end
            x = glp_get_col_prim(subproblem, info.column)
            if abs(x - round(Int, x)) > 1e-7
                return MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            end
        end
        return MOI.CALLBACK_NODE_STATUS_INTEGER
    elseif reason == GLP_IHEUR
        return MOI.CALLBACK_NODE_STATUS_FRACTIONAL
    elseif reason == GLP_ICUTGEN
        return MOI.CALLBACK_NODE_STATUS_FRACTIONAL
    elseif reason == GLP_IBRANCH
        return MOI.CALLBACK_NODE_STATUS_FRACTIONAL
    else
        @assert reason == GLP_IBINGO
        return MOI.CALLBACK_NODE_STATUS_INTEGER
    end
end

# ==============================================================================
#    MOI.LazyConstraint
# ==============================================================================

MOI.supports(::Optimizer, ::MOI.LazyConstraintCallback) = true

MOI.supports(::Optimizer, ::MOI.LazyConstraint{CallbackData}) = true

function MOI.set(model::Optimizer, ::MOI.LazyConstraintCallback, cb::Function)
    model.lazy_callback = cb
    return
end

function MOI.submit(
    model::Optimizer,
    cb::MOI.LazyConstraint{CallbackData},
    f::MOI.ScalarAffineFunction{Float64},
    s::Union{
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
        MOI.EqualTo{Float64},
    },
)
    if model.callback_state == _CB_USER_CUT
        throw(MOI.InvalidCallbackUsage(MOI.UserCutCallback(), cb))
    elseif model.callback_state == _CB_HEURISTIC
        throw(MOI.InvalidCallbackUsage(MOI.HeuristicCallback(), cb))
    end
    key = CleverDicts.add_item(model.affine_constraint_info, _ConstraintInfo(s))
    model.affine_constraint_info[key].row = length(model.affine_constraint_info)
    indices, coefficients = _indices_and_coefficients(model, f)
    sense, rhs = _sense_and_rhs(s)
    inner = glp_ios_get_prob(cb.callback_data.tree)
    _add_affine_constraint(
        inner,
        indices,
        coefficients,
        sense,
        rhs - f.constant,
    )
    return
end

# ==============================================================================
#    MOI.UserCutCallback
# ==============================================================================

MOI.supports(::Optimizer, ::MOI.UserCutCallback) = true
MOI.supports(::Optimizer, ::MOI.UserCut{CallbackData}) = true

function MOI.set(model::Optimizer, ::MOI.UserCutCallback, cb::Function)
    model.user_cut_callback = cb
    return
end

function MOI.submit(
    model::Optimizer,
    cb::MOI.UserCut{CallbackData},
    f::MOI.ScalarAffineFunction{Float64},
    s::Union{MOI.LessThan{Float64},MOI.GreaterThan{Float64}},
)
    if model.callback_state == _CB_LAZY
        throw(MOI.InvalidCallbackUsage(MOI.LazyConstraintCallback(), cb))
    elseif model.callback_state == _CB_HEURISTIC
        throw(MOI.InvalidCallbackUsage(MOI.HeuristicCallback(), cb))
    end
    indices, coefficients = _indices_and_coefficients(model, f)
    sense, rhs = _sense_and_rhs(s)
    bound_type = sense == Cchar('G') ? GLP_LO : GLP_UP
    glp_ios_add_row(
        cb.callback_data.tree,
        "",
        101,
        0,
        Cint(length(indices)),
        offset(indices),
        offset(coefficients),
        bound_type,
        rhs,
    )
    return
end

# ==============================================================================
#    MOI.HeuristicCallback
# ==============================================================================

MOI.supports(::Optimizer, ::MOI.HeuristicCallback) = true
MOI.supports(::Optimizer, ::MOI.HeuristicSolution{CallbackData}) = true

function MOI.set(model::Optimizer, ::MOI.HeuristicCallback, cb::Function)
    model.heuristic_callback = cb
    return
end

function MOI.submit(
    model::Optimizer,
    cb::MOI.HeuristicSolution{CallbackData},
    variables::Vector{MOI.VariableIndex},
    values::MOI.Vector{Float64},
)
    if model.callback_state == _CB_LAZY
        throw(MOI.InvalidCallbackUsage(MOI.LazyConstraintCallback(), cb))
    elseif model.callback_state == _CB_USER_CUT
        throw(MOI.InvalidCallbackUsage(MOI.UserCutCallback(), cb))
    end
    solution = fill(NaN, MOI.get(model, MOI.NumberOfVariables()))
    for (var, value) in zip(variables, values)
        solution[_info(model, var).column] = value
    end
    ret = glp_ios_heur_sol(cb.callback_data.tree, offset(solution))
    return ret == 0 ? MOI.HEURISTIC_SOLUTION_ACCEPTED :
           MOI.HEURISTIC_SOLUTION_REJECTED
end
