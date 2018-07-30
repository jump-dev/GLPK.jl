using LinQuadOptInterface

import Compat.LinearAlgebra
using Nullables

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

mutable struct Optimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    presolve::Bool
    method::Symbol
    interior::GLPK.InteriorParam
    intopt::GLPK.IntoptParam
    simplex::SimplexParam
    solver_status::Int32
    last_solved_by_mip::Bool
    # See https://github.com/JuliaOpt/GLPKMathProgInterface.jl/pull/15
    # for why this is necesary. GLPK interacts weirdly with binary variables and
    # bound modification. So lets set binary variables as "Integer" with [0,1]
    # bounds that we enforce just before solve.
    binaries::Vector{Int}
    Optimizer(::Nothing) = new()
end
function Optimizer(presolve=false, method=:Simplex; kwargs...)
    optimizer = Optimizer(nothing)
    MOI.empty!(optimizer)
    optimizer.presolve = presolve
    optimizer.method   = method
    optimizer.interior = GLPK.InteriorParam()
    optimizer.intopt   = GLPK.IntoptParam()
    optimizer.simplex  = GLPK.SimplexParam()
    solver_status      = Int32(0)
    optimizer.last_solved_by_mip = false
    optimizer.binaries = Int[]
    return optimizer
end

LQOI.LinearQuadraticModel(::Type{Optimizer}, env) = GLPK.Prob()

LQOI.supported_objectives(model::Optimizer) = SUPPORTED_OBJECTIVES
LQOI.supported_constraints(model::Optimizer) = SUPPORTED_CONSTRAINTS

function LQOI.change_variable_bounds!(model::Optimizer,
          columns::Vector{Int}, new_bounds::Vector{Float64},
          senses::Vector{Cchar})
    for (column, bound, sense) in zip(columns, new_bounds, senses)
        new_upper = Inf
        new_lower = -Inf
        bound_type = GLPK.DB
        if sense == Cchar('E')
            new_upper = bound
            new_lower = bound
            bound_type = GLPK.FX
        elseif sense == Cchar('L')
            current_upper_bound = GLPK.get_col_ub(model.inner, column)
            new_lower = bound
            if current_upper_bound < 1e100
                new_upper = current_upper_bound
                bound_type = GLPK.DB
            else
                new_upper = Inf
                bound_type = GLPK.LO
            end
        elseif sense == Cchar('U')
            current_lower_bound = GLPK.get_col_lb(model.inner, column)
            new_upper = bound
            if current_lower_bound > -1e100
                new_lower = current_lower_bound
                bound_type = GLPK.DB
            else
                new_lower = -Inf
                bound_type = GLPK.UP
            end
        else
            error("Variable sense $(sense) not recognised.")
        end
        if new_upper > 1e100 && new_lower < -1e100
            bound_type = GLPK.FR
            new_upper = Inf
            new_lower = -Inf
        elseif new_lower == new_upper
            bound_type = GLPK.FX
        end
        GLPK.set_col_bnds(model.inner, column, bound_type, new_lower,
                          new_upper)
    end
end

function LQOI.get_variable_lowerbound(model::Optimizer, col)
    GLPK.get_col_lb(model.inner, col)
end

function LQOI.get_variable_upperbound(model::Optimizer, col)
    GLPK.get_col_ub(model.inner, col)
end

function LQOI.get_number_linear_constraints(model::Optimizer)
    GLPK.get_num_rows(model.inner)
end

function LQOI.add_linear_constraints!(model::Optimizer,
        A::LQOI.CSRMatrix{Float64}, senses::Vector{Cchar}, rhs::Vector{Float64})
    nrows = length(rhs)
    if nrows <= 0
        error("no row to be added")
    elseif nrows == 1
        addrow!(model.inner, A.columns, A.coefficients, senses[1], rhs[1])
    else
        push!(A.row_pointers, length(A.columns)+1)
        for i in 1:nrows
            indices = A.row_pointers[i]:A.row_pointers[i+1]-1
            addrow!(model.inner, A.columns[indices], A.coefficients[indices],
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
        GLPK.set_row_bnds(model.inner, row, GLPK.DB, lower, upper)
    end
end

function LQOI.modify_ranged_constraints!(model::Optimizer,
        rows::Vector{Int}, lowerbounds::Vector{Float64},
        upperbounds::Vector{Float64})
    for (row, lower, upper) in zip(rows, lowerbounds, upperbounds)
        LQOI.change_rhs_coefficient!(model, row, lower)
        GLPK.set_row_bnds(model.inner, row, GLPK.DB, lower, upper)
    end
end

function LQOI.get_range(model::Optimizer, row::Int)
    GLPK.get_row_lb(model.inner, row), GLPK.get_row_ub(model.inner, row)
end

function addrow!(lp::GLPK.Prob, colidx::Vector, colcoef::Vector, sense::Cchar, rhs::Real)
    if length(colidx) != length(colcoef)
        error("colidx and colcoef have different legths")
    end
    GLPK.add_rows(lp, 1)
    m = GLPK.get_num_rows(lp)
    GLPK.set_mat_row(lp, m, colidx, colcoef)
    if sense == Cchar('E')
        bt = GLPK.FX
        rowlb = rhs
        rowub = rhs
    elseif sense == Cchar('G')
        bt = GLPK.LO
        rowlb = rhs
        rowub = Inf
    elseif sense == Cchar('L')
        bt = GLPK.UP
        rowlb = -Inf
        rowub = rhs
    elseif sense == Cchar('R')
        # start with lower
        bt = GLPK.DB
        rowlb = rhs
        rowub = Inf
    else
        error("row type $(sense) not valid")
        bt = GLPK.FR
    end
    GLPK.set_row_bnds(lp, m, bt, rowlb, rowub)
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
    columns, coefficients = GLPK.get_mat_row(model.inner, row)
    # note: we return 1-indexed columns here
    return columns, coefficients
end

function LQOI.change_rhs_coefficient!(model::Optimizer, row::Int,
                                      rhs::Real)
    current_lower = GLPK.get_row_lb(model.inner, row)
    current_upper = GLPK.get_row_ub(model.inner, row)
    if current_lower == current_upper
        bound_type = GLPK.FX
        new_lower_bound = rhs
        new_upper_bound = rhs
    elseif current_lower > -Inf && current_upper < Inf
        bound_type = GLPK.FX
        new_lower_bound = rhs
        new_upper_bound = rhs
    elseif current_lower > -Inf
        bound_type = GLPK.LO
        new_lower_bound = rhs
        new_upper_bound = Inf
    elseif current_upper < Inf
        bound_type = GLPK.UP
        new_lower_bound = -Inf
        new_upper_bound = rhs
    else
        error("not valid rhs")
    end
    GLPK.set_row_bnds(model.inner, row, bound_type, new_lower_bound,
                      new_upper_bound)
end

function LQOI.change_objective_coefficient!(model::Optimizer, col, coef)
    GLPK.set_obj_coef(model.inner, col, coef)
end

function LQOI.change_matrix_coefficient!(model::Optimizer, row, col, coef)
    columns, coefficients = GLPK.get_mat_row(model.inner, row)
    index = findfirst(columns, col)
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
    new_vtype = GLPK.CV
    model.binaries = Int[]
    for (col, vtype) in zip(columns, variable_types)
        if vtype == Cint('I')
            new_vtype = GLPK.IV
        elseif vtype == Cint('C')
            new_vtype = GLPK.CV
        elseif vtype == Cint('B')
            new_vtype = GLPK.IV
            push!(model.binaries, col)
        else
            error("Invalid variable type: $(vtype).")
        end
        GLPK.set_col_kind(model.inner, col, new_vtype)
    end
end

function LQOI.change_linear_constraint_sense!(model::Optimizer, rows, senses)
    for (row, sense) in zip(rows, senses)
        changesense!(model, row, sense)
    end
end

function changesense!(model::Optimizer, row, sense)
    m = model.inner
    oldsense = GLPK.get_row_type(m, row)
    newsense = sense_char_to_glpk(sense)
    if oldsense == newsense
        return nothing
    end

    if newsense == GLPK.FX
        rowub = rowlb
    elseif oldsense == GLPK.DB
        if newsense == GLPK.UP
            rowlb = -Inf
        else
            rowub = Inf
        end
    else
        rowlb = get_row_lb(m, row)
        rowub = get_row_ub(m, row)
        if newsense == GLPK.UP
            rowub = rowlb
            rowlb = -Inf
        else
            rowlb = rowub
            rowub = Inf
        end
    end

    GLPK.set_row_bnds(m, row, newsense, rowlb, rowub)

    nothing
end

function sense_char_to_glpk(sense)
    if sense == Cchar('E')
        return GLPK.FX
    elseif sense == Cchar('R')
        return GLPK.DB
    elseif sense == Cchar('L')
        return GLPK.UP
    elseif sense == Cchar('G')
        return GLPK.LO
    else
        error("Invalid constraint sense: $(sense).")
    end
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
        error("Unrecognized objective sense $sense")
    end
end

function LQOI.get_linear_objective!(model::Optimizer, x::Vector{Float64})
    @assert length(x) == GLPK.get_num_cols(model.inner)
    for col in 1:length(x)
        x[col] = GLPK.get_obj_coef(model.inner, col)
    end
end

function LQOI.get_objectivesense(model::Optimizer)
    s = GLPK.get_obj_dir(model.inner)
    if s == GLPK.MIN
        return MOI.MinSense
    elseif s == GLPK.MAX
        return MOI.MaxSense
    else
        error("Internal library error")
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

function _certificates_available(model::Optimizer)
    !model.last_solved_by_mip && (model.method == :Simplex || model.method == :Exact)
end

function LQOI.get_termination_status(model::Optimizer)
    if model.last_solved_by_mip
        if model.solver_status in [GLPK.EMIPGAP, GLPK.ETMLIM, GLPK.ESTOP]
            return MOI.OtherLimit
        end
    end
    status = get_status(model)
    if status == GLPK.OPT
        return MOI.Success
    elseif status == GLPK.INFEAS
        if _certificates_available(model)
            return MOI.Success
        else
            return MOI.InfeasibleNoResult
        end
    elseif status == GLPK.UNBND
        if _certificates_available(model)
            return MOI.Success
        else
            return MOI.UnboundedNoResult
        end
    elseif status == GLPK.FEAS
        return MOI.SlowProgress
    elseif status == GLPK.NOFEAS
        return MOI.InfeasibleOrUnbounded
    elseif status == GLPK.UNDEF
        return MOI.OtherError
    else
        error("Status $(status) not recognised.")
    end
end

function get_status(model::Optimizer)
    if model.last_solved_by_mip
        return GLPK.mip_status(model.inner)
    else
        _check_lp_method(model.method)
        if model.method == :Simplex || model.method == :Exact
            return GLPK.get_status(model.inner)
        elseif model.method == :InteriorPoint
            return GLPK.ipt_status(model.inner)
        end
    end
end

function LQOI.get_primal_status(model::Optimizer)
    status = get_status(model)
    if status == GLPK.OPT
        return MOI.FeasiblePoint
    elseif status == GLPK.UNBND && _certificates_available(model)
        return MOI.InfeasibilityCertificate
    end
    return MOI.UnknownResultStatus
end

function LQOI.get_dual_status(model::Optimizer)
    if !model.last_solved_by_mip
        status = get_status(model.inner)
        if status == GLPK.OPT
            return MOI.FeasiblePoint
        elseif status == GLPK.INFEAS && _certificates_available(model)
            return MOI.InfeasibilityCertificate
        end
    end
    return MOI.UnknownResultStatus
end

function copy_function_result!(dest, foo, model)
    for i in eachindex(dest)
        dest[i] = foo(model.inner, i)
    end
end

function _check_lp_method(method)
    if !(method in [:Simplex, :Exact, :InteriorPoint])
        error("Method must be one of :Simplex, :Exact, or :InteriorPoint.")
    end
end

function LQOI.get_variable_primal_solution!(model::Optimizer, place)
    if model.last_solved_by_mip
        copy_function_result!(place, GLPK.mip_col_val, model)
    else
        _check_lp_method(model.method)
        if model.method == :Simplex || model.method == :Exact
            copy_function_result!(place, GLPK.get_col_prim, model)
        elseif model.method == :InteriorPoint
            copy_function_result!(place, GLPK.ipt_col_prim, model)
        end
    end
end

function LQOI.get_linear_primal_solution!(model::Optimizer, place)
    if model.last_solved_by_mip
        copy_function_result!(place, GLPK.mip_row_val, model)
    else
        _check_lp_method(model.method)
        if model.method == :Simplex || model.method == :Exact
            copy_function_result!(place, GLPK.get_row_prim, model)
        elseif model.method == :InteriorPoint
            copy_function_result!(place, GLPK.ipt_row_prim, model)
        end
    end
end

function LQOI.get_variable_dual_solution!(model::Optimizer, place)
    @assert !model.last_solved_by_mip
    _check_lp_method(model.method)
    if model.method == :Simplex || model.method == :Exact
        copy_function_result!(place, GLPK.get_col_dual, model)
    elseif model.method == :InteriorPoint
        copy_function_result!(place, GLPK.ipt_col_dual, model)
    end
end

function LQOI.get_linear_dual_solution!(model::Optimizer, place)
    @assert !model.last_solved_by_mip
    _check_lp_method(model.method)
    if model.method == :Simplex || model.method == :Exact
        copy_function_result!(place, GLPK.get_row_dual, model)
    elseif model.method == :InteriorPoint
        copy_function_result!(place, GLPK.ipt_row_dual, model)
    end
end

function LQOI.get_objective_value(model::Optimizer)
    if model.last_solved_by_mip
        return GLPK.mip_obj_val(model.inner)
    else
        _check_lp_method(model.method)
        if model.method == :Simplex || model.method == :Exact
            return GLPK.get_obj_val(model.inner)
        elseif model.method == :InteriorPoint
            return GLPK.ipt_obj_val(model.inner)
        end
    end
end

function LQOI.solve_linear_problem!(model::Optimizer)
    _check_lp_method(model.method)
    model.last_solved_by_mip = false
    if model.method == :Simplex
        return GLPK.simplex(model.inner, model.simplex)
    elseif model.method == :Exact
        return GLPK.exact(model.inner, model.simplex)
    elseif model.method == :InteriorPoint
        return GLPK.interior(model.inner, model.interior)
    end
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
    try
        if model.intopt.presolve == GLPK.OFF
            status = GLPK.simplex(model.inner, model.simplex)
            if status != 0
                model.last_solved_by_mip = false
                model.solver_status = status
                return
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
    end
end

function set_variable_bound(model::Optimizer, col::Int, lower::Float64,
                            upper::Float64)
    if upper >= realmax(Float64)
        upper = Inf
    end
    if lower <= -realmax(Float64)
        lower = -Inf
    end
    if upper < Inf
        if lower > -Inf
            if lower == upper
                GLPK.set_col_bnds(model.inner, col, GLPK.FX, lower, upper)
            else
                GLPK.set_col_bnds(model.inner, col, GLPK.DB, lower, upper)
            end
        else
            GLPK.set_col_bnds(model.inner, col, GLPK.UP, 0.0, upper)
        end
    else
        if lower > -Inf
            GLPK.set_col_bnds(model.inner, col, GLPK.LO, lower, 0.0)
        else
            GLPK.set_col_bnds(model.inner, col, GLPK.FR, 0.0, 0.0)
        end
    end
end

const VARIABLE_TYPE_MAP = Dict(
    GLPK.CV => :Cont,
    GLPK.IV => :Int,
    GLPK.BV => :Bin
)

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
