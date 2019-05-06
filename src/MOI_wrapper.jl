using LinQuadOptInterface

const LQOI = LinQuadOptInterface
const MOI  = LQOI.MOI

const SUPPORTED_OBJECTIVES = [
    LQOI.SinVar,
    LQOI.Linear
]

const SUPPORTED_CONSTRAINTS = [
    (LQOI.Linear, LQOI.EQ),
    (LQOI.Linear, LQOI.LE),
    (LQOI.Linear, LQOI.GE),
    (LQOI.Linear, LQOI.IV),
    (LQOI.SinVar, LQOI.EQ),
    (LQOI.SinVar, LQOI.LE),
    (LQOI.SinVar, LQOI.GE),
    (LQOI.SinVar, LQOI.IV),
    (LQOI.VecVar, MOI.Nonnegatives),
    (LQOI.VecVar, MOI.Nonpositives),
    (LQOI.VecVar, MOI.Zeros),
    (LQOI.VecLin, MOI.Nonnegatives),
    (LQOI.VecLin, MOI.Nonpositives),
    (LQOI.VecLin, MOI.Zeros),
    (LQOI.SinVar, MOI.ZeroOne),
    (LQOI.SinVar, MOI.Integer)
]

"""
    AbstractCallbackData

An abstract type to prevent recursive type definition of Optimizer and
CallbackData, each of which need the other type in a field.
"""
abstract type AbstractCallbackData end

mutable struct Optimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    presolve::Bool
    method::Symbol
    interior::GLPK.InteriorParam
    intopt::GLPK.IntoptParam
    simplex::SimplexParam
    solver_status::Int32
    last_solved_by_mip::Bool
    # See note associated with abstractcallbackdata. When using, make sure to
    # add a type assertion, since this will always be the concrete type
    # CallbackData.
    callback_data::AbstractCallbackData
    objective_bound::Float64
    callback_function::Function
    # See https://github.com/JuliaOpt/GLPKMathProgInterface.jl/pull/15
    # for why this is necesary. GLPK interacts weirdly with binary variables and
    # bound modification. So lets set binary variables as "Integer" with [0,1]
    # bounds that we enforce just before solve.
    binaries::Vector{Int}
    Optimizer(::Nothing) = new()
end

"""
    CallbackData
"""
mutable struct CallbackData <: AbstractCallbackData
    model::Optimizer
    tree::Ptr{Cvoid}
    CallbackData(model::Optimizer) = new(model, C_NULL)
end

"""
    __internal_callback__(tree::Ptr{Cvoid}, info::Ptr{Cvoid})

Dummy callback function for internal use only. Responsible for updating the
objective bound.
"""
function __internal_callback__(tree::Ptr{Cvoid}, info::Ptr{Cvoid})
    callback_data = unsafe_pointer_to_objref(info)::CallbackData
    model = callback_data.model
    callback_data.tree = tree
    node = GLPK.ios_best_node(tree)
    if node != 0
        model.objective_bound = GLPK.ios_node_bound(tree, node)
    end
    model.callback_function(callback_data)
    return nothing
end

function Optimizer(;presolve=false, method=:Simplex, kwargs...)
    optimizer = Optimizer(nothing)
    MOI.empty!(optimizer)
    optimizer.presolve = presolve
    optimizer.method   = method
    optimizer.interior = GLPK.InteriorParam()
    optimizer.intopt   = GLPK.IntoptParam()
    optimizer.simplex  = GLPK.SimplexParam()
    solver_status = Int32(0)
    optimizer.last_solved_by_mip = false
    optimizer.callback_data = CallbackData(optimizer)
    optimizer.intopt.cb_func = @cfunction(__internal_callback__, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))
    optimizer.intopt.cb_info = pointer_from_objref(optimizer.callback_data)
    optimizer.objective_bound = NaN
    optimizer.callback_function = (cb_data::CallbackData) -> nothing
    optimizer.binaries = Int[]
    # Parameters
    if presolve
        optimizer.simplex.presolve = GLPK.ON
        optimizer.intopt.presolve  = GLPK.ON
    end
    optimizer.interior.msg_lev = GLPK.MSG_ERR
    optimizer.intopt.msg_lev   = GLPK.MSG_ERR
    optimizer.simplex.msg_lev  = GLPK.MSG_ERR
    for (key, value) in kwargs
        set_interior = set_parameter(optimizer.interior, key, value)
        set_intopt   = set_parameter(optimizer.intopt, key, value)
        set_simplex  = set_parameter(optimizer.simplex, key, value)
        if !set_interior && !set_intopt && !set_simplex
            @warn("Ignoring option: $(key) => $(value)")
        end
    end
    return optimizer
end

MOI.get(::Optimizer, ::MOI.SolverName) = "GLPK"

function LQOI.get_objective_bound(model::Optimizer)
    if !model.last_solved_by_mip
        return LQOI.get_objectivesense(model) == MOI.MIN_SENSE ? -Inf : Inf
    end
    constant = LQOI.get_constant_objective(model)
    # @mlubin and @ccoey observed some cases where mip_status == OPT and objval
    # and objbound didn't match. In that case, they return mip_obj_val, but
    # objbound may still be incorrect in cases where GLPK terminates early.
    if GLPK.mip_status(model.inner) == GLPK.OPT
        return GLPK.mip_obj_val(model.inner) + constant
    else
        return model.objective_bound + constant
    end
end

"""
    set_parameter(param_store, key::Symbol, value)

Set the field name `key` in a `param_store` type (that is one of `InteriorParam`,
`IntoptParam`, or `SimplexParam`) to `value`.
"""
function set_parameter(param_store, key::Symbol, value)
    if key in [:cb_func, :cb_info]
        @warn("Ignored option: $(string(k)). Use the MOI attribute " *
              "`GLPK.CallbackFunction` instead.")
        return true
    end

    if key in fieldnames(typeof(param_store))
        field_type = typeof(getfield(param_store, key))
        setfield!(param_store, key, convert(field_type, value))
        return true
    else
        return false
    end
end

LQOI.LinearQuadraticModel(::Type{Optimizer}, env) = GLPK.Prob()

LQOI.supported_objectives(model::Optimizer) = SUPPORTED_OBJECTIVES
LQOI.supported_constraints(model::Optimizer) = SUPPORTED_CONSTRAINTS

function LQOI.set_constant_objective!(model::Optimizer, value)
    GLPK.set_obj_coef(model.inner, 0, value)
end
LQOI.get_constant_objective(model::Optimizer) = GLPK.get_obj_coef(model.inner, 0)

"""
    get_col_bound_type(lower::Float64, upper::Float64)

Return the GLPK type of the variable bound given a lower bound of `lower` and an
upper bound of `upper`.
"""
function get_col_bound_type(lower::Float64, upper::Float64)
    GLPK_INFINITY = VERSION >= v"0.7-" ? floatmax(Float64) : realmax(Float64)
    if lower == upper
        return GLPK.FX
    elseif lower <= -GLPK_INFINITY
        return upper >= GLPK_INFINITY ? GLPK.FR : GLPK.UP
    else
        return upper >= GLPK_INFINITY ? GLPK.LO : GLPK.DB
    end
end

"""
    set_variable_bound(model::Optimizer, column::Int, lower::Float64, upper::Float64)

Set the bounds of the variable in column `column` to `[lower, upper]`.
"""
function set_variable_bound(model::Optimizer, column::Int, lower::Float64,
                            upper::Float64)
    bound_type = get_col_bound_type(lower, upper)
    # Disable preemptive checking of variable bounds for the case when lower
    # > upper. If you solve a model with lower > upper, the
    # TerminationStatus will be InvalidModel.
    prev_preemptive_check = GLPK.jl_get_preemptive_check()
    GLPK.jl_set_preemptive_check(false)
    GLPK.set_col_bnds(model.inner, column, bound_type, lower, upper)
    # Reset the preemptive check.
    GLPK.jl_set_preemptive_check(prev_preemptive_check)
end

function LQOI.change_variable_bounds!(model::Optimizer, columns::Vector{Int},
        new_bounds::Vector{Float64}, senses::Vector{Cchar})
    for (column, bound, sense) in zip(columns, new_bounds, senses)
        if sense == Cchar('L')
            lower_bound = bound
            upper_bound = GLPK.get_col_ub(model.inner, column)
        elseif sense == Cchar('U')
            lower_bound = GLPK.get_col_lb(model.inner, column)
            upper_bound = bound
        else
            error("Invalid variable bound sense: $(sense)")
        end
        set_variable_bound(model, column, lower_bound, upper_bound)
    end
end

function LQOI.get_variable_lowerbound(model::Optimizer, col)
    return GLPK.get_col_lb(model.inner, col)
end

function LQOI.get_variable_upperbound(model::Optimizer, col)
    return GLPK.get_col_ub(model.inner, col)
end

function LQOI.get_number_linear_constraints(model::Optimizer)
    return GLPK.get_num_rows(model.inner)
end

function LQOI.add_linear_constraints!(model::Optimizer,
        A::LQOI.CSRMatrix{Float64}, senses::Vector{Cchar}, rhs::Vector{Float64})
    nrows = length(rhs)
    if nrows <= 0
        error("Number of rows must be more than zero.")
    elseif nrows == 1
        add_row!(model.inner, A.columns, A.coefficients, senses[1], rhs[1])
    else
        push!(A.row_pointers, length(A.columns)+1)
        for i in 1:nrows
            indices = A.row_pointers[i]:A.row_pointers[i+1]-1
            add_row!(model.inner, A.columns[indices], A.coefficients[indices],
                    senses[i], rhs[i])
        end
        pop!(A.row_pointers)
    end
end

function LQOI.add_ranged_constraints!(model::Optimizer,
        A::LQOI.CSRMatrix{Float64}, lowerbound::Vector{Float64},
        upperbound::Vector{Float64})
    row1 = GLPK.get_num_rows(model.inner)
    LQOI.add_linear_constraints!(model, A,
                                 fill(Cchar('R'), length(lowerbound)),
                                 lowerbound)
    row2 = GLPK.get_num_rows(model.inner)
    for (row, lower, upper) in zip(row1+1:row2, lowerbound, upperbound)
        # Disable preemptive checking of variable bounds for the case when lower
        # > upper. If you solve a model with lower > upper, the
        # TerminationStatus will be InvalidModel.
        prev_preemptive_check = GLPK.jl_get_preemptive_check()
        GLPK.jl_set_preemptive_check(false)
        GLPK.set_row_bnds(model.inner, row, GLPK.DB, lower, upper)
        # Reset the preemptive check.
        GLPK.jl_set_preemptive_check(prev_preemptive_check)
    end
end

function LQOI.modify_ranged_constraints!(model::Optimizer,
        rows::Vector{Int}, lowerbounds::Vector{Float64},
        upperbounds::Vector{Float64})
    for (row, lower, upper) in zip(rows, lowerbounds, upperbounds)
        LQOI.change_rhs_coefficient!(model, row, lower)
        # Disable preemptive checking of variable bounds for the case when lower
        # > upper. If you solve a model with lower > upper, the
        # TerminationStatus will be InvalidModel.
        prev_preemptive_check = GLPK.jl_get_preemptive_check()
        GLPK.jl_set_preemptive_check(false)
        GLPK.set_row_bnds(model.inner, row, GLPK.DB, lower, upper)
        # Reset the preemptive check.
        GLPK.jl_set_preemptive_check(prev_preemptive_check)
    end
end

function LQOI.get_range(model::Optimizer, row::Int)
    return GLPK.get_row_lb(model.inner, row), GLPK.get_row_ub(model.inner, row)
end

"""
    add_row!(problem::GLPK.Prob, columns::Vector{Int},
             coefficients::Vector{Float64}, sense::Cchar, rhs::Real)

Helper function to add a row to the problem. Sense must be one of `'E'` (ax == b),
`'G'` (ax >= b), `'L'` (ax <= b) , or `'R'` (b <= ax).

If the sense is `'R'` the `rhs` should be the lower bound, and the bounds should
be set in a new API call to enforce the upper bound.
"""
function add_row!(problem::GLPK.Prob, columns::Vector{Int},
                 coefficients::Vector{Float64}, sense::Cchar, rhs::Real)
    if length(columns) != length(coefficients)
        error("columns and coefficients have different lengths.")
    end
    GLPK.add_rows(problem, 1)
    num_rows = GLPK.get_num_rows(problem)
    GLPK.set_mat_row(problem, num_rows, columns, coefficients)
    # According to http://most.ccib.rutgers.edu/glpk.pdf page 22,
    # the `lb` argument is ignored for constraint types with no
    # lower bound (GLPK.UP) and the `ub` argument is ignored for
    # constraint types with no upper bound (GLPK.LO). We pass
    # ±DBL_MAX for those unused bounds since (a) we have to pass
    # something, and (b) it is consistent with the other usages of
    # ±DBL_MAX to represent infinite bounds in the rest of the
    # GLPK interface.
    if sense == Cchar('E')
        GLPK.set_row_bnds(problem, num_rows, GLPK.FX, rhs, rhs)
    elseif sense == Cchar('G')
        GLPK.set_row_bnds(problem, num_rows, GLPK.LO, rhs, GLPK.DBL_MAX)
    elseif sense == Cchar('L')
        GLPK.set_row_bnds(problem, num_rows, GLPK.UP, -GLPK.DBL_MAX, rhs)
    elseif sense == Cchar('R')
        GLPK.set_row_bnds(problem, num_rows, GLPK.DB, rhs, GLPK.DBL_MAX)
    else
        error("Invalid row sense: $(sense)")
    end
end

function LQOI.get_rhs(model::Optimizer, row)
    sense = GLPK.get_row_type(model.inner, row)
    if sense == GLPK.LO || sense == GLPK.FX || sense == GLPK.DB
        return GLPK.get_row_lb(model.inner, row)
    else
        return GLPK.get_row_ub(model.inner, row)
    end
end

function LQOI.get_linear_constraint(model::Optimizer, row::Int)
    # note: we return 1-indexed columns here
    return GLPK.get_mat_row(model.inner, row)
end

function LQOI.change_rhs_coefficient!(model::Optimizer, row::Int,
                                      rhs::Real)
    current_lower = GLPK.get_row_lb(model.inner, row)
    current_upper = GLPK.get_row_ub(model.inner, row)
    # `get_row_lb` and `get_row_ub` return ±DBL_MAX for rows with no
    # lower or upper bound. See page 30 of the GLPK user manual
    # http://most.ccib.rutgers.edu/glpk.pdf
    if current_lower == current_upper
        GLPK.set_row_bnds(model.inner, row,  GLPK.FX, rhs, rhs)
    elseif current_lower > -GLPK.DBL_MAX && current_upper < GLPK.DBL_MAX
        GLPK.set_row_bnds(model.inner, row,  GLPK.FX, rhs, rhs)
    elseif current_lower > -GLPK.DBL_MAX
        GLPK.set_row_bnds(model.inner, row,  GLPK.LO, rhs, GLPK.DBL_MAX)
    elseif current_upper < GLPK.DBL_MAX
        GLPK.set_row_bnds(model.inner, row,  GLPK.UP, -GLPK.DBL_MAX, rhs)
    else
        error("Cannot set right-hand side of a free constraint.")
    end
end

function LQOI.change_objective_coefficient!(model::Optimizer, col, coef)
    GLPK.set_obj_coef(model.inner, col, coef)
end

function LQOI.change_matrix_coefficient!(model::Optimizer, row, col, coef)
    columns, coefficients = GLPK.get_mat_row(model.inner, row)
    index = something(findfirst(isequal(col), columns), 0)
    if index > 0
        coefficients[index] = coef
    else
        push!(columns, col)
        push!(coefficients, coef)
    end
    GLPK.set_mat_row(model.inner, row, columns, coefficients)
end

function LQOI.delete_linear_constraints!(model::Optimizer, start_row, stop_row)
    GLPK.std_basis(model.inner)
    indices = collect(start_row:stop_row)
    GLPK.del_rows(model.inner, length(indices), indices)
end

function LQOI.change_variable_types!(model::Optimizer,
        columns::Vector{Int}, variable_types::Vector)
    model.binaries = Int[]
    for (column, variable_type) in zip(columns, variable_types)
        if variable_type == Cint('I')
            GLPK.set_col_kind(model.inner, column, GLPK.IV)
        elseif variable_type == Cint('C')
            GLPK.set_col_kind(model.inner, column, GLPK.CV)
        elseif variable_type == Cint('B')
            # note: we lie to GLPK here and set it as an integer variable. See
            # the comment in the definition of Optimizer.
            GLPK.set_col_kind(model.inner, column, GLPK.IV)
            push!(model.binaries, column)
        else
            error("Invalid variable type: $(vtype).")
        end
    end
end

function LQOI.change_linear_constraint_sense!(model::Optimizer, rows, senses)
    for (row, sense) in zip(rows, senses)
        change_row_sense!(model, row, sense)
    end
end

"""
    change_row_sense!(model::Optimizer, row, sense)

Convert a linear constraint into another type of linear constraint by changing
the comparison sense.

Constraint types supported are 'E' (equality), 'L' (less-than), and 'G'
(greater-than).

For example, `ax <= b` can become `ax >= b` or `ax == b`.
"""
function change_row_sense!(model::Optimizer, row::Int, sense)
    old_sense = GLPK.get_row_type(model.inner, row)
    if old_sense == GLPK.DB
        error("Cannot transform sense of an interval constraint.")
    elseif old_sense == GLPK.FR
        error("Cannot transform sense of a free constraint.")
    end
    new_sense = ROW_SENSE_MAP[sense]
    if old_sense == new_sense
        error("Cannot transform constraint with same sense.")
    elseif new_sense == GLPK.DB
        error("Cannot transform constraint to ranged constraint.")
    end
    if old_sense == GLPK.FX || old_sense == GLPK.LO
        right_hand_side = GLPK.get_row_lb(model.inner, row)
    else
        right_hand_side = GLPK.get_row_ub(model.inner, row)
    end
    # According to http://most.ccib.rutgers.edu/glpk.pdf page 22,
    # the `lb` argument is ignored for constraint types with no
    # lower bound (GLPK.UP) and the `ub` argument is ignored for
    # constraint types with no upper bound (GLPK.LO). We pass
    # ±DBL_MAX for those unused bounds since (a) we have to pass
    # something, and (b) it is consistent with the other usages of
    # ±DBL_MAX to represent infinite bounds in the rest of the
    # GLPK interface.
    if new_sense == GLPK.FX
        GLPK.set_row_bnds(model.inner, row, new_sense, right_hand_side, right_hand_side)
    elseif new_sense == GLPK.LO
        GLPK.set_row_bnds(model.inner, row, new_sense, right_hand_side, GLPK.DBL_MAX)
    elseif new_sense == GLPK.UP
        GLPK.set_row_bnds(model.inner, row, new_sense, -GLPK.DBL_MAX, right_hand_side)
    end
end

const ROW_SENSE_MAP = Dict(
    Cchar('E') => GLPK.FX,
    Cchar('R') => GLPK.DB,
    Cchar('L') => GLPK.UP,
    Cchar('G') => GLPK.LO
)

function LQOI.add_sos_constraint!(model::Optimizer, columns, weights, sos_type)
    GLPK.add_sos!(instance.inner, sos_type, columns, weights)
end

function LQOI.set_linear_objective!(model::Optimizer, columns, coefficients)
    ncols = GLPK.get_num_cols(model.inner)
    new_coefficients = zeros(ncols)
    for (col, coef) in zip(columns, coefficients)
        new_coefficients[col] += coef
    end
    for (col, coef) in zip(1:ncols, new_coefficients)
        GLPK.set_obj_coef(model.inner, col, coef)
    end
end

function LQOI.change_objective_sense!(model::Optimizer, sense)
    if sense == :min
        GLPK.set_obj_dir(model.inner, GLPK.MIN)
    elseif sense == :max
        GLPK.set_obj_dir(model.inner, GLPK.MAX)
    else
        error("Invalid objective sense: $(sense)")
    end
end

function LQOI.get_linear_objective!(model::Optimizer, x::Vector{Float64})
    @assert length(x) == GLPK.get_num_cols(model.inner)
    for col in 1:length(x)
        x[col] = GLPK.get_obj_coef(model.inner, col)
    end
end

function LQOI.get_objectivesense(model::Optimizer)
    sense = GLPK.get_obj_dir(model.inner)
    if sense == GLPK.MIN
        return MOI.MIN_SENSE
    elseif sense == GLPK.MAX
        return MOI.MAX_SENSE
    else
        error("Invalid objective sense: $(sense)")
    end
end

function LQOI.get_number_variables(model::Optimizer)
    GLPK.get_num_cols(model.inner)
end

function LQOI.add_variables!(model::Optimizer, number_to_add::Int)
    num_variables = GLPK.get_num_cols(model.inner)
    GLPK.add_cols(model.inner, number_to_add)
    for i in 1:number_to_add
        GLPK.set_col_bnds(model.inner, num_variables+i, GLPK.FR, 0.0, 0.0)
    end
end

function LQOI.delete_variables!(model::Optimizer, col, col2)
    GLPK.std_basis(model.inner)
    columns = collect(col:col2)
    GLPK.del_cols(model.inner, length(columns), columns)
end

function MOI.write_to_file(model::Optimizer, lp_file_name::String)
    GLPK.write_lp(model.inner, lp_file_name)
end

"""
    _certificates_potentially_available(model::Optimizer)

Return true if an infeasiblity certificate or an unbounded ray is potentially
available (i.e., the model has been solved using either the Simplex or Exact
methods).
"""
function _certificates_potentially_available(model::Optimizer)
    !model.last_solved_by_mip && (model.method == :Simplex || model.method == :Exact)
end

function LQOI.get_termination_status(model::Optimizer)
    if model.solver_status == GLPK.EMIPGAP  # MIP Gap
        return MOI.OPTIMAL
    elseif model.solver_status == GLPK.EBOUND  # Invalid bounds
        return MOI.INVALID_MODEL
    elseif model.solver_status == GLPK.ETMLIM  # Time limit
        return MOI.TIME_LIMIT
    elseif model.solver_status == GLPK.ENODFS  # No feasible dual
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif model.solver_status == GLPK.ENOPFS  # No feasible primal
        return MOI.INFEASIBLE
    elseif model.solver_status == GLPK.EFAIL  # Solver fail
        return MOI.NUMERICAL_ERROR
    elseif model.solver_status == GLPK.ESTOP  # Callback
        return MOI.INTERRUPTED
    end
    status = get_status(model)
    if status == GLPK.OPT
        return MOI.OPTIMAL
    elseif status == GLPK.INFEAS
        return MOI.INFEASIBLE
    elseif status == GLPK.UNBND
        return MOI.DUAL_INFEASIBLE
    elseif status == GLPK.FEAS
        return MOI.SLOW_PROGRESS
    elseif status == GLPK.NOFEAS
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif status == GLPK.UNDEF
        return MOI.OTHER_ERROR
    else
        error("Invalid termination status: $(status)")
    end
end

"""
    get_status(model::Optimizer)

Get the status from GLPK depending on which method was used to solve the model.
"""
function get_status(model::Optimizer)
    if model.last_solved_by_mip
        return GLPK.mip_status(model.inner)
    else
        if model.method == :Simplex || model.method == :Exact
            return GLPK.get_status(model.inner)
        elseif model.method == :InteriorPoint
            return GLPK.ipt_status(model.inner)
        end
        _throw_invalid_method(model)
    end
end

function LQOI.get_primal_status(model::Optimizer)
    status = get_status(model)
    if status == GLPK.OPT
        return MOI.FEASIBLE_POINT
    elseif status == GLPK.UNBND && _certificates_potentially_available(model)
        return MOI.INFEASIBILITY_CERTIFICATE
    end
    return MOI.NO_SOLUTION
end

function LQOI.get_dual_status(model::Optimizer)
    if !model.last_solved_by_mip
        status = get_status(model.inner)
        if status == GLPK.OPT
            return MOI.FEASIBLE_POINT
        elseif status == GLPK.INFEAS && _certificates_potentially_available(model)
            return MOI.INFEASIBILITY_CERTIFICATE
        end
    end
    return MOI.NO_SOLUTION
end

"""
    copy_function_result!(dest::Vector, foo, model::GLPK.Prob)

A helper function that loops through the indices in `dest` and stores the result
of `foo(model, i)` for the `i`th index.
"""
function copy_function_result!(dest::Vector, foo, model::GLPK.Prob)
    for i in eachindex(dest)
        dest[i] = foo(model, i)
    end
end
function copy_function_result!(dest::Vector, foo, model::Optimizer)
    copy_function_result!(dest, foo, model.inner)
end

"""
    _throw_invalid_method(instance::Optimizer)

A helper function to throw an error when the method is set incorrectly. Mainly
used to enforce type-stability in functions that have a run-time switch on the
method.
"""
function _throw_invalid_method(instance::Optimizer)
    error("Method is $(instance.method), but it must be one of :Simplex, " *
          ":Exact, or :InteriorPoint.")
end

function LQOI.get_variable_primal_solution!(model::Optimizer, place)
    if model.last_solved_by_mip
        copy_function_result!(place, GLPK.mip_col_val, model)
    else
        if model.method == :Simplex || model.method == :Exact
            copy_function_result!(place, GLPK.get_col_prim, model)
        elseif model.method == :InteriorPoint
            copy_function_result!(place, GLPK.ipt_col_prim, model)
        else
            _throw_invalid_method(model)
        end
    end
end

function LQOI.get_linear_primal_solution!(model::Optimizer, place)
    if model.last_solved_by_mip
        copy_function_result!(place, GLPK.mip_row_val, model)
    else
        if model.method == :Simplex || model.method == :Exact
            copy_function_result!(place, GLPK.get_row_prim, model)
        elseif model.method == :InteriorPoint
            copy_function_result!(place, GLPK.ipt_row_prim, model)
        else
            _throw_invalid_method(model)
        end
    end
end

function LQOI.get_variable_dual_solution!(model::Optimizer, place)
    @assert !model.last_solved_by_mip
    if model.method == :Simplex || model.method == :Exact
        copy_function_result!(place, GLPK.get_col_dual, model)
    elseif model.method == :InteriorPoint
        copy_function_result!(place, GLPK.ipt_col_dual, model)
    else
        _throw_invalid_method(model)
    end
end

function LQOI.get_linear_dual_solution!(model::Optimizer, place)
    @assert !model.last_solved_by_mip
    if model.method == :Simplex || model.method == :Exact
        copy_function_result!(place, GLPK.get_row_dual, model)
    elseif model.method == :InteriorPoint
        copy_function_result!(place, GLPK.ipt_row_dual, model)
    else
        _throw_invalid_method(model)
    end
end

function LQOI.get_objective_value(model::Optimizer)
    if model.last_solved_by_mip
        return GLPK.mip_obj_val(model.inner)
    else
        if model.method == :Simplex || model.method == :Exact
            return GLPK.get_obj_val(model.inner)
        elseif model.method == :InteriorPoint
            return GLPK.ipt_obj_val(model.inner)
        end
        _throw_invalid_method(model)
    end
end

function LQOI.solve_linear_problem!(model::Optimizer)
    model.last_solved_by_mip = false
    if model.method == :Simplex
        model.solver_status = GLPK.simplex(model.inner, model.simplex)
    elseif model.method == :Exact
        model.solver_status = GLPK.exact(model.inner, model.simplex)
    elseif model.method == :InteriorPoint
        model.solver_status = GLPK.interior(model.inner, model.interior)
    else
        _throw_invalid_method(model)
    end
    return
end

"""
    round_bounds_to_integer(model)::Tuple{Bool, Vector{Float64}, Vector{Float64}}

GLPK does not allow integer variables with fractional bounds. Therefore, we
round the bounds of binary and integer variables to integer values prior to
solving.

Returns a tuple of the original bounds, along with a Boolean flag indicating if
they need to be reset after solve.
"""
function round_bounds_to_integer(model::Optimizer)
    num_variables = GLPK.get_num_cols(model.inner)
    lower_bounds = map(i->GLPK.get_col_lb(model.inner, i), 1:num_variables)
    upper_bounds = map(i->GLPK.get_col_ub(model.inner, i), 1:num_variables)
    variable_types = get_variable_types(model)
    bounds_modified = false
    for (col, variable_type) in enumerate(variable_types)
        new_lower = ceil(lower_bounds[col])
        new_upper = floor(upper_bounds[col])
        if variable_type == :Bin
            new_lower = max(0.0, ceil(lower_bounds[col]))
            new_upper = min(1.0, floor(upper_bounds[col]))
        elseif variable_type != :Int
            continue
        end
        if lower_bounds[col] != new_lower || upper_bounds[col] != new_upper
            set_variable_bound(model, col, new_lower, new_upper)
            bounds_modified = true
        end
    end
    return bounds_modified, lower_bounds, upper_bounds
end

function LQOI.solve_mip_problem!(model::Optimizer)
    bounds_modified, lower_bounds, upper_bounds = round_bounds_to_integer(model)
    # Because we're muddling with the presolve in this function, cache the
    # original setting so that it can be reset.
    presolve_cache = model.intopt.presolve
    try
        # GLPK.intopt requires a starting basis for the LP relaxation. There are
        # two ways to get this. If presolve=GLPK.ON, then the presolve will find
        # a basis. If presolve=GLPK.OFF, then we should solve the problem via
        # GLPK.simplex first.
        if model.intopt.presolve == GLPK.OFF
            GLPK.simplex(model.inner, model.simplex)
            if GLPK.get_status(model.inner) != GLPK.OPT
                # We didn't find an optimal solution to the LP relaxation, so
                # let's turn presolve on and let intopt figure out what the
                # problem is.
                model.intopt.presolve = GLPK.ON
            end
        end
        model.solver_status = GLPK.intopt(model.inner, model.intopt)
        model.last_solved_by_mip = true
    finally
        if bounds_modified
            for (col, (lower, upper)) in enumerate(zip(lower_bounds, upper_bounds))
                set_variable_bound(model, col, lower, upper)
            end
        end
        # Restore the original presolve setting.
        model.intopt.presolve = presolve_cache
    end
end

const VARIABLE_TYPE_MAP = Dict(
    GLPK.CV => :Cont,
    GLPK.IV => :Int,
    GLPK.BV => :Bin
)

"""
    get_variable_types(model::Optimizer)

Return a vector of symbols (one element for each variable) of the variable type.
The symbols are given by the key-value mapping in `GLPK.VARIABLE_TYPE_MAP`.
"""
function get_variable_types(model::Optimizer)
    ncol = GLPK.get_num_cols(model.inner)
    col_types = fill(:Cont, ncol)
    for i in 1:ncol
        col_type = GLPK.get_col_kind(model.inner, i)
        col_types[i] = VARIABLE_TYPE_MAP[col_type]
        if i in model.binaries
            # See the note in the definition of Optimizer about this.
            col_types[i] = :Bin
        elseif col_types[i] == :Bin
            # We never set a variable as binary. GLPK must have made a mistake
            # and inferred so based on bounds?
            col_types[i] = :Int
        end
    end
    return col_types
end

include("infeasibility_certificates.jl")

function LQOI.get_farkas_dual!(model::Optimizer, place)
    get_infeasibility_ray(model, place)
end

function LQOI.get_unbounded_ray!(model::Optimizer, place)
    get_unbounded_ray(model, place)
end

# ==============================================================================
#    Callbacks in GLPK
# ==============================================================================
"""
    CallbackFunction

The attribute to set the callback function in GLPK. The function takes a single
argument of type `CallbackData`.
"""
struct CallbackFunction <: MOI.AbstractOptimizerAttribute end
function MOI.set(model::Optimizer, ::CallbackFunction, foo::Function)
    model.callback_function = foo
end

"""
    load_variable_primal!(cb_data::CallbackData)

Load the variable primal solution in a callback.

This can only be called in a callback from `GLPK.IROWGEN`. After it is called,
you can access the `VariablePrimal` attribute as usual.
"""
function load_variable_primal!(cb_data::CallbackData)
    model = cb_data.model
    if GLPK.ios_reason(cb_data.tree) != GLPK.IROWGEN
        error("load_variable_primal! can only be called when reason is GLPK.IROWGEN.")
    end
    subproblem = GLPK.ios_get_prob(cb_data.tree)
    copy_function_result!(model.variable_primal_solution, GLPK.get_col_prim, subproblem)
end

"""
    add_lazy_constraint!(cb_data::GLPK.CallbackData, func::LQOI.Linear, set::S) where S <: Union{LQOI.LE, LQOI.GE, LQOI.EQ}

Add a lazy constraint to the model `cb_data.model`. This can only be called in a
callback from `GLPK.IROWGEN`.
"""
function add_lazy_constraint!(cb_data::CallbackData, func::LQOI.Linear, set::S) where S <: Union{LQOI.LE, LQOI.GE, LQOI.EQ}
    model = cb_data.model
    add_row!(
        GLPK.ios_get_prob(cb_data.tree),
        [LQOI.get_column(model, term.variable_index) for term in func.terms],
        [term.coefficient for term in func.terms],
        LQOI.backend_type(model, set),
        MOI.Utilities.getconstant(set)
    )
end
