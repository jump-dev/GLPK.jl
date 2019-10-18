import MathOptInterface

const MOI = MathOptInterface
const CleverDicts = MOI.Utilities.CleverDicts

@enum(TypeEnum, CONTINUOUS, BINARY, INTEGER)
@enum(BoundEnum, NONE, LESS_THAN, GREATER_THAN, LESS_AND_GREATER_THAN, INTERVAL, EQUAL_TO)
@enum(ObjectiveEnum, SINGLE_VARIABLE, SCALAR_AFFINE)
@enum(MethodEnum, SIMPLEX, INTERIOR, EXACT)
@enum(CallbackState, CB_NONE, CB_GENERIC, CB_LAZY, CB_USER_CUT, CB_HEURISTIC)

mutable struct VariableInfo
    index::MOI.VariableIndex
    column::Int
    bound::BoundEnum
    type::TypeEnum
    name::String
    # Storage for constraint names associated with variables because GLPK
    # can only store names for variables and proper constraints.
    # We can perform an optimization and only store three strings for the
    # constraint names because, at most, there can be three SingleVariable
    # constraints, e.g., LessThan, GreaterThan, and Integer.
    lessthan_name::String
    greaterthan_interval_or_equalto_name::String
    type_constraint_name::String
    function VariableInfo(index::MOI.VariableIndex, column::Int)
        return new(index, column, NONE, CONTINUOUS, "", "", "", "")
    end
end

struct ConstraintKey
    value::Int
end
CleverDicts.key_to_index(k::ConstraintKey) = k.value
CleverDicts.index_to_key(::Type{ConstraintKey}, index) = ConstraintKey(index)

mutable struct ConstraintInfo
    row::Int
    set::MOI.AbstractSet
    name::String
    ConstraintInfo(set) = new(0, set, "")
end

# Dummy callback function for internal use only. Responsible for updating the
# objective bound, saving the mip gap, and calling the user's callback.
function _internal_callback(tree::Ptr{Cvoid}, info::Ptr{Cvoid})
    callback_data = unsafe_pointer_to_objref(info)::CallbackData
    model = callback_data.model
    callback_data.tree = tree
    node = GLPK.ios_best_node(tree)
    if node != 0
        model.objective_bound = GLPK.ios_node_bound(tree, node)
        model.relative_gap = GLPK.ios_mip_gap(tree)
    end
    model.callback_function(callback_data)
    return nothing
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    # The low-level GLPK problem.
    inner::GLPK.Prob
    presolve::Bool
    method::MethodEnum

    interior_param::GLPK.InteriorParam
    intopt_param::GLPK.IntoptParam
    simplex_param::GLPK.SimplexParam
    solver_status::Int32
    last_solved_by_mip::Bool
    num_binaries::Int
    num_integers::Int

    objective_bound::Float64
    relative_gap::Float64
    solve_time::Float64
    callback_data::Any

    # A flag to keep track of MOI.Silent, which over-rides the print_level
    # parameter.
    silent::Bool

    # An enum to remember what objective is currently stored in the model.
    objective_type::ObjectiveEnum

    # A flag to keep track of MOI.FEASIBILITY_SENSE, since Gurobi only stores
    # MIN_SENSE or MAX_SENSE. This allows us to differentiate between MIN_SENSE
    # and FEASIBILITY_SENSE.
    is_feasibility::Bool

    variable_info::CleverDicts.CleverDict{MOI.VariableIndex, VariableInfo}
    affine_constraint_info::CleverDicts.CleverDict{ConstraintKey, ConstraintInfo}

    # Mappings from variable and constraint names to their indices. These are
    # lazily built on-demand, so most of the time, they are `nothing`.
    name_to_variable::Union{Nothing, Dict{String, Union{Nothing, MOI.VariableIndex}}}
    name_to_constraint_index::Union{Nothing, Dict{String, Union{Nothing, MOI.ConstraintIndex}}}

    optimize_not_called::Bool

    # These two flags allow us to distinguish between FEASIBLE_POINT and
    # INFEASIBILITY_CERTIFICATE when querying VariablePrimal and ConstraintDual.
    unbounded_ray::Union{Vector{Float64}, Nothing}
    infeasibility_cert::Union{Vector{Float64}, Nothing}

    # Callback fields.
    has_generic_callback::Bool
    callback_state::CallbackState
    lazy_callback::Union{Nothing, Function}
    user_cut_callback::Union{Nothing, Function}
    heuristic_callback::Union{Nothing, Function}

    """
        Optimizer(;kwargs...)

    Create a new Optimizer object. Common keywords include

     - `method::MethodEnum = SIMPLEX` Use the simplex method. Other options are `EXACT` and `INTERIOR`.
     - `tm_lim::Float64`              Set a time limit
     - `msg_lev::Int`                 Control the log level

    See the GLPK pdf documentation for a full list of parameters.
    """
    function Optimizer(; presolve = false, method = SIMPLEX, kwargs...)
        model = new()
        model.inner = GLPK.Prob()
        model.presolve = presolve
        model.method = method

        model.interior_param = GLPK.InteriorParam()
        model.intopt_param = GLPK.IntoptParam()
        model.simplex_param = GLPK.SimplexParam()

        MOI.set(model, MOI.RawParameter("msg_lev"), GLPK.MSG_ERR)
        if model.presolve
            MOI.set(model, MOI.RawParameter("presolve"), GLPK.ON)
        end
        for (key, val) in kwargs
            MOI.set(model, MOI.RawParameter(String(key)), val)
        end
        model.silent = false
        model.variable_info = CleverDicts.CleverDict{MOI.VariableIndex, VariableInfo}()
        model.affine_constraint_info = CleverDicts.CleverDict{ConstraintKey, ConstraintInfo}()

        MOI.empty!(model)

        return model
    end
end

mutable struct CallbackData
    c_callback::Base.CFunction
    tree::Ptr{Cvoid}
    exception::Union{Nothing, Exception}
end

# Dummy callback function for internal use only. Responsible for updating the
# objective bound, saving the mip gap, and calling the user's callback.
#
# !!! Very Important Note
#
# If Julia throws an exception from within the callback, the GLPK model does not
# gracefully exit! Instead, it throws a `glp_delete_prob: operation not
# supported` error when Julia tries to finalize `prob`. This is very annoying to
# debug.
#
# As a work-around, we catch all Julia exceptions with a try-catch, terminate
# the callback with `ios_terminate`, store the exception in `cb_data.exception`,
# and then re-throw it once we have gracefully exited from the callback.
#
# See also: the note in `_solve_mip_problem`.
function set_callback(model::Optimizer, callback_function::Function)
    internal_callback = (tree::Ptr{Cvoid}, info::Ptr{Cvoid}) -> begin
        cb_data = unsafe_pointer_to_objref(info)::CallbackData
        node = GLPK.ios_best_node(tree)
        if node != 0
            model.objective_bound = GLPK.ios_node_bound(tree, node)
            model.relative_gap = GLPK.ios_mip_gap(tree)
        end
        try
            cb_data.tree = tree
            callback_function(cb_data)
        catch ex
            GLPK.ios_terminate(tree)
            cb_data.exception = ex
        end
        return Cint(0)
    end
    c_callback = @cfunction(
        $internal_callback, Cint, (Ptr{Cvoid}, Ptr{Cvoid})
    )
    model.callback_data = CallbackData(c_callback, C_NULL, nothing)
    model.intopt_param.cb_func = c_callback.ptr
    model.intopt_param.cb_info = pointer_from_objref(model.callback_data)
    return
end

Base.show(io::IO, model::Optimizer) = print(io, "A GLPK model")

function MOI.empty!(model::Optimizer)
    GLPK.erase_prob(model.inner)
    model.solver_status = GLPK.UNDEF
    model.last_solved_by_mip = false
    model.num_binaries = 0
    model.num_integers = 0
    model.objective_bound = NaN
    model.relative_gap = NaN
    model.solve_time = NaN
    model.objective_type = SCALAR_AFFINE
    model.is_feasibility = true
    model.optimize_not_called = true
    empty!(model.variable_info)
    empty!(model.affine_constraint_info)
    model.name_to_variable = nothing
    model.name_to_constraint_index = nothing
    model.unbounded_ray = nothing
    model.infeasibility_cert = nothing
    model.has_generic_callback = false
    model.callback_state = CB_NONE
    model.lazy_callback = nothing
    model.user_cut_callback = nothing
    model.heuristic_callback = nothing
    set_callback(model, cb_data -> nothing)
    return
end

function MOI.is_empty(model::Optimizer)
    model.objective_type != SCALAR_AFFINE && return false
    model.is_feasibility == false && return false
    !isempty(model.variable_info) && return false
    !isempty(model.affine_constraint_info) && return false
    model.name_to_variable !== nothing && return false
    model.name_to_constraint_index !== nothing && return false
    model.unbounded_ray !== nothing && return false
    model.infeasibility_cert !== nothing && return false
    model.has_generic_callback && return false
    model.callback_state != CB_NONE && return false
    model.lazy_callback !== nothing && return false
    model.user_cut_callback !== nothing && return false
    model.heuristic_callback !== nothing && return false
    return true
end

MOI.get(::Optimizer, ::MOI.SolverName) = "GLPK"

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{F}
) where {F <: Union{
    MOI.SingleVariable,
    MOI.ScalarAffineFunction{Float64}
}}
    return true
end

function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.SingleVariable}, ::Type{F}
) where {F <: Union{
    MOI.EqualTo{Float64}, MOI.LessThan{Float64}, MOI.GreaterThan{Float64},
    MOI.Interval{Float64}, MOI.ZeroOne, MOI.Integer
}}
    return true
end

function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{F}
) where {F <: Union{
    MOI.EqualTo{Float64}, MOI.LessThan{Float64}, MOI.GreaterThan{Float64}
}}
    return true
end

const SCALAR_SETS = Union{
    MOI.GreaterThan{Float64}, MOI.LessThan{Float64},
    MOI.EqualTo{Float64}, MOI.Interval{Float64}
}

MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true
MOI.supports(::Optimizer, ::MOI.ConstraintName, ::Type{<:MOI.ConstraintIndex}) = true

MOI.supports(::Optimizer, ::MOI.Name) = true
MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.RawParameter) = true

"""
    set_parameter(param_store, key::Symbol, value)::Bool

Set the field name `key` in a `param_store` type (that is one of `InteriorParam`,
`IntoptParam`, or `SimplexParam`) to `value`.

Returns a `Bool` indicating if the parameter was set.
"""
function set_parameter(param_store, key::Symbol, value)
    if key == :cb_func || key == :cb_info
        error("Invalid option: $(string(key)). Use the MOI attribute " *
              "`GLPK.CallbackFunction` instead.")
    elseif key in fieldnames(typeof(param_store))
        field_type = typeof(getfield(param_store, key))
        setfield!(param_store, key, convert(field_type, value))
        return true
    end
    return false
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value)
    if typeof(param.name) != String
        error("GLPK.jl requires strings as arguments to `RawParameter`.")
    end
    key = Symbol(param.name)
    set_interior = set_parameter(model.interior_param, key, value)
    set_intopt = set_parameter(model.intopt_param, key, value)
    set_simplex = set_parameter(model.simplex_param, key, value)
    if !set_interior && !set_intopt && !set_simplex
        throw(MOI.UnsupportedAttribute(param))
    end
    return
end

function MOI.get(model::Optimizer, param::MOI.RawParameter)
    if typeof(param.name) != String
        error("GLPK.jl requires strings as arguments to `RawParameter`.")
    end
    name = Symbol(param.name)
    if (model.method == SIMPLEX || model.method == EXACT) && name in fieldnames(GLPK.SimplexParam)
        return getfield(model.simplex_param, name)
    elseif model.method == INTERIOR && name in fieldnames(GLPK.InteriorParam)
        return getfield(model.interior_param, name)
    elseif name in fieldnames(GLPK.IntoptParam)
        return getfield(model.intopt_param, name)
    end
    throw(MOI.UnsupportedAttribute(param))
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, limit::Union{Nothing,Real})
    if limit === nothing
        MOI.set(model, MOI.RawParameter("tm_lim"), typemax(Int32))
    else
        limit_ms = ceil(Int32, limit * 1000)
        MOI.set(model, MOI.RawParameter("tm_lim"), limit_ms)
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    # convert internal ms to sec
    return MOI.get(model, MOI.RawParameter("tm_lim")) / 1000
end

MOI.Utilities.supports_default_copy_to(::Optimizer, ::Bool) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kwargs...)
    return MOI.Utilities.automatic_copy_to(dest, src; kwargs...)
end

function MOI.get(model::Optimizer, ::MOI.ListOfVariableAttributesSet)
    return MOI.AbstractVariableAttribute[MOI.VariableName()]
end

function MOI.get(model::Optimizer, ::MOI.ListOfModelAttributesSet)
    obj_func_type = MOI.get(model, MOI.ObjectiveFunctionType())
    attributes = [
        MOI.ObjectiveSense(),
        MOI.ObjectiveFunction{obj_func_type}()
    ]
    if MOI.get(model, MOI.Name()) != ""
        push!(attributes, MOI.Name())
    end
    return attributes
end

function MOI.get(model::Optimizer, ::MOI.ListOfConstraintAttributesSet)
    return MOI.AbstractConstraintAttribute[MOI.ConstraintName()]
end

function _indices_and_coefficients(
    indices::AbstractVector{Int}, coefficients::AbstractVector{Float64},
    model::Optimizer, f::MOI.ScalarAffineFunction{Float64}
)
    i = 1
    for term in f.terms
        indices[i] = _info(model, term.variable_index).column
        coefficients[i] = term.coefficient
        i += 1
    end
    return indices, coefficients
end

function _indices_and_coefficients(
    model::Optimizer, f::MOI.ScalarAffineFunction{Float64}
)
    f_canon = MOI.Utilities.canonical(f)
    nnz = length(f_canon.terms)
    indices = Vector{Int}(undef, nnz)
    coefficients = Vector{Float64}(undef, nnz)
    _indices_and_coefficients(indices, coefficients, model, f_canon)
    return indices, coefficients
end

_sense_and_rhs(s::MOI.LessThan{Float64}) = (Cchar('L'), s.upper)
_sense_and_rhs(s::MOI.GreaterThan{Float64}) = (Cchar('G'), s.lower)
_sense_and_rhs(s::MOI.EqualTo{Float64}) = (Cchar('E'), s.value)

###
### Variables
###

# Short-cuts to return the VariableInfo associated with an index.
function _info(model::Optimizer, key::MOI.VariableIndex)
    if haskey(model.variable_info, key)
        return model.variable_info[key]
    end
    throw(MOI.InvalidIndex(key))
end

function MOI.add_variable(model::Optimizer)
    # Initialize `VariableInfo` with a dummy `VariableIndex` and a column,
    # because we need `add_item` to tell us what the `VariableIndex` is.
    index = CleverDicts.add_item(
        model.variable_info, VariableInfo(MOI.VariableIndex(0), 0)
    )
    info = _info(model, index)
    # Now, set `.index` and `.column`.
    info.index = index
    info.column = length(model.variable_info)
    GLPK.add_cols(model.inner, 1)
    GLPK.set_col_bnds(model.inner, info.column, GLPK.FR, 0.0, 0.0)
    return index
end

function MOI.add_variables(model::Optimizer, N::Int)
    indices = Vector{MOI.VariableIndex}(undef, N)
    num_variables = length(model.variable_info)
    GLPK.add_cols(model.inner, N)
    for i in 1:N
        # Initialize `VariableInfo` with a dummy `VariableIndex` and a column,
        # because we need `add_item` to tell us what the `VariableIndex` is.
        index = CleverDicts.add_item(
            model.variable_info, VariableInfo(MOI.VariableIndex(0), 0)
        )
        info = _info(model, index)
        # Now, set `.index` and `.column`.
        info.index = index
        info.column = num_variables + i
        GLPK.set_col_bnds(model.inner, info.column, GLPK.FR, 0.0, 0.0)
        indices[i] = index
    end
    return indices
end

function MOI.is_valid(model::Optimizer, v::MOI.VariableIndex)
    return haskey(model.variable_info, v)
end

function MOI.delete(model::Optimizer, v::MOI.VariableIndex)
    info = _info(model, v)
    GLPK.std_basis(model.inner)
    GLPK.del_cols(model.inner, 1, [info.column])
    delete!(model.variable_info, v)
    for other_info in values(model.variable_info)
        if other_info.column > info.column
            other_info.column -= 1
        end
    end
    model.name_to_variable = nothing
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(model::Optimizer, ::Type{MOI.VariableIndex}, name::String)
    if model.name_to_variable === nothing
        _rebuild_name_to_variable(model)
    end
    if haskey(model.name_to_variable, name)
        variable = model.name_to_variable[name]
        if variable === nothing
            error("Duplicate name detected: $(name)")
        end
        return variable
    end
    return nothing
end

function _rebuild_name_to_variable(model::Optimizer)
    model.name_to_variable = Dict{String, Union{Nothing, MOI.VariableIndex}}()
    for (index, info) in model.variable_info
        if isempty(info.name)
            continue
        end
        if haskey(model.name_to_variable, info.name)
            model.name_to_variable[info.name] = nothing
        else
            model.name_to_variable[info.name] = index
        end
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex)
    return _info(model, v).name
end

function MOI.set(
    model::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex, name::String
)
    info = _info(model, v)
    info.name = name
    if !isempty(name) && isascii(name)
        # Note: GLPK errors if we try to set non-ascii column names.
        GLPK.set_col_name(model.inner, info.column, name)
    end
    model.name_to_variable = nothing
    return
end

###
### Objectives
###

function MOI.set(
    model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense
)
    if sense == MOI.MIN_SENSE
        GLPK.set_obj_dir(model.inner, GLPK.MIN)
        model.is_feasibility = false
    elseif sense == MOI.MAX_SENSE
        GLPK.set_obj_dir(model.inner, GLPK.MAX)
        model.is_feasibility = false
    else
        @assert sense == MOI.FEASIBILITY_SENSE
        GLPK.set_obj_dir(model.inner, GLPK.MIN)
        model.is_feasibility = true
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    sense = GLPK.get_obj_dir(model.inner)
    if model.is_feasibility
        return MOI.FEASIBILITY_SENSE
    elseif sense == GLPK.MAX
        return MOI.MAX_SENSE
    else
        @assert sense == GLPK.MIN
        return MOI.MIN_SENSE
    end
end

function MOI.set(
    model::Optimizer, ::MOI.ObjectiveFunction{F}, f::F
) where {F <: MOI.SingleVariable}
    MOI.set(
        model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        convert(MOI.ScalarAffineFunction{Float64}, f)
    )
    model.objective_type = SINGLE_VARIABLE
    return
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunction{F}) where {F}
    obj = MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    return convert(F, obj)
end

function MOI.set(
    model::Optimizer, ::MOI.ObjectiveFunction{F}, f::F
) where {F <: MOI.ScalarAffineFunction{Float64}}
    num_vars = length(model.variable_info)
    obj = zeros(Float64, num_vars)
    for term in f.terms
        column = _info(model, term.variable_index).column
        obj[column] += term.coefficient
    end
    for (col, coef) in enumerate(obj)
        GLPK.set_obj_coef(model.inner, col, coef)
    end
    GLPK.set_obj_coef(model.inner, 0, f.constant)
    model.objective_type = SCALAR_AFFINE
end

function MOI.get(
    model::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}
)
    dest = zeros(length(model.variable_info))
    for col in 1:length(dest)
        dest[col] = GLPK.get_obj_coef(model.inner, col)
    end
    terms = MOI.ScalarAffineTerm{Float64}[]
    for (index, info) in model.variable_info
        coefficient = dest[info.column]
        iszero(coefficient) && continue
        push!(terms, MOI.ScalarAffineTerm(coefficient, index))
    end
    constant = GLPK.get_obj_coef(model.inner, 0)
    return MOI.ScalarAffineFunction(terms, constant)
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarConstantChange{Float64}
)
    GLPK.set_obj_coef(model.inner, 0, chg.new_constant)
    return
end

##
##  SingleVariable-in-Set constraints.
##

function _info(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, <:Any}
)
    var_index = MOI.VariableIndex(c.value)
    if haskey(model.variable_info, var_index)
        return _info(model, var_index)
    end
    return throw(MOI.InvalidIndex(c))
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    if haskey(model.variable_info, MOI.VariableIndex(c.value))
        info = _info(model, c)
        return info.bound == LESS_THAN || info.bound == LESS_AND_GREATER_THAN
    end
    return false
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    if haskey(model.variable_info, MOI.VariableIndex(c.value))
        info = _info(model, c)
        return info.bound == GREATER_THAN || info.bound == LESS_AND_GREATER_THAN
    end
    return false
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        _info(model, c).bound == INTERVAL
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        _info(model, c).bound == EQUAL_TO
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        _info(model, c).type == BINARY
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        _info(model, c).type == INTEGER
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable, <:Any}
)
    MOI.throw_if_not_valid(model, c)
    return MOI.SingleVariable(MOI.VariableIndex(c.value))
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable, <:Any}, ::MOI.SingleVariable
)
    return throw(MOI.SettingSingleVariableFunctionNotAllowed())
end

_bounds(s::MOI.GreaterThan{Float64}) = (s.lower, nothing)
_bounds(s::MOI.LessThan{Float64}) = (nothing, s.upper)
_bounds(s::MOI.EqualTo{Float64}) = (s.value, s.value)
_bounds(s::MOI.Interval{Float64}) = (s.lower, s.upper)

function _throw_if_existing_lower(
    bound::BoundEnum, var_type::TypeEnum, new_set::Type{<:MOI.AbstractSet},
    variable::MOI.VariableIndex
)
    existing_set = if bound == LESS_AND_GREATER_THAN || bound == GREATER_THAN
        MOI.GreaterThan{Float64}
    elseif bound == INTERVAL
        MOI.Interval{Float64}
    elseif bound == EQUAL_TO
        MOI.EqualTo{Float64}
    else
        nothing  # Also covers `NONE` and `LESS_THAN`.
    end
    if existing_set !== nothing
        throw(MOI.LowerBoundAlreadySet{existing_set, new_set}(variable))
    end
end

function _throw_if_existing_upper(
    bound::BoundEnum, var_type::TypeEnum, new_set::Type{<:MOI.AbstractSet},
    variable::MOI.VariableIndex
)
    existing_set = if bound == LESS_AND_GREATER_THAN || bound == LESS_THAN
        MOI.LessThan{Float64}
    elseif bound == INTERVAL
        MOI.Interval{Float64}
    elseif bound == EQUAL_TO
        MOI.EqualTo{Float64}
    else
        nothing  # Also covers `NONE` and `GREATER_THAN`.
    end
    if existing_set !== nothing
        throw(MOI.UpperBoundAlreadySet{existing_set, new_set}(variable))
    end
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.SingleVariable, s::S
) where {S <: SCALAR_SETS}
    info = _info(model, f.variable)
    if S <: MOI.LessThan{Float64}
        _throw_if_existing_upper(info.bound, info.type, S, f.variable)
        info.bound = info.bound == GREATER_THAN ? LESS_AND_GREATER_THAN : LESS_THAN
    elseif S <: MOI.GreaterThan{Float64}
        _throw_if_existing_lower(info.bound, info.type, S, f.variable)
        info.bound = info.bound == LESS_THAN ? LESS_AND_GREATER_THAN : GREATER_THAN
    elseif S <: MOI.EqualTo{Float64}
        _throw_if_existing_lower(info.bound, info.type, S, f.variable)
        _throw_if_existing_upper(info.bound, info.type, S, f.variable)
        info.bound = EQUAL_TO
    else
        @assert S <: MOI.Interval{Float64}
        _throw_if_existing_lower(info.bound, info.type, S, f.variable)
        _throw_if_existing_upper(info.bound, info.type, S, f.variable)
        info.bound = INTERVAL
    end
    index = MOI.ConstraintIndex{MOI.SingleVariable, typeof(s)}(f.variable.value)
    MOI.set(model, MOI.ConstraintSet(), index, s)
    return index
end

function _set_variable_bound(
    model::Optimizer, column::Int, lower::Union{Nothing, Float64},
    upper::Union{Nothing, Float64}
)
    if lower === nothing
        lower = GLPK.get_col_lb(model.inner, column)
    end
    if upper === nothing
        upper = GLPK.get_col_ub(model.inner, column)
    end
    bound_type = if lower == upper
        GLPK.FX
    elseif lower <= -GLPK.DBL_MAX
        upper >= GLPK.DBL_MAX ? GLPK.FR : GLPK.UP
    else
        upper >= GLPK.DBL_MAX ? GLPK.LO : GLPK.DB
    end
    if upper < lower
        # Here, we disable GLPK's pre-emptive checks, because otherwise GLPK
        # will through an error complaining about invalid bounds when `upper <
        # lower`. This let's us throw `INFEASIBLE` instead of erroring.
        prev_preemptive = GLPK.jl_get_preemptive_check()
        GLPK.jl_set_preemptive_check(false)
        GLPK.set_col_bnds(model.inner, column, bound_type, lower, upper)
        GLPK.jl_set_preemptive_check(prev_preemptive)
    else
        GLPK.set_col_bnds(model.inner, column, bound_type, lower, upper)
    end
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    _set_variable_bound(model, info.column, nothing, Inf)
    if info.bound == LESS_AND_GREATER_THAN
        info.bound = GREATER_THAN
    else
        info.bound = NONE
    end
    info.lessthan_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    _set_variable_bound(model, info.column, -Inf, nothing)
    info.bound = info.bound == LESS_AND_GREATER_THAN ? LESS_THAN : NONE
    info.greaterthan_interval_or_equalto_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    _set_variable_bound(model, info.column, -Inf, Inf)
    info.bound = NONE
    info.greaterthan_interval_or_equalto_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    _set_variable_bound(model, info.column, -Inf, Inf)
    info.bound = NONE
    info.greaterthan_interval_or_equalto_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    MOI.throw_if_not_valid(model, c)
    lower = GLPK.get_col_lb(model.inner, _info(model, c).column)
    return MOI.GreaterThan(lower)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    MOI.throw_if_not_valid(model, c)
    upper = GLPK.get_col_ub(model.inner, _info(model, c).column)
    return MOI.LessThan(upper)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
)
    MOI.throw_if_not_valid(model, c)
    lower = GLPK.get_col_lb(model.inner, _info(model, c).column)
    return MOI.EqualTo(lower)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    lower = GLPK.get_col_lb(model.inner, info.column)
    upper = GLPK.get_col_ub(model.inner, info.column)
    return MOI.Interval(lower, upper)
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}, s::S
) where {S<:SCALAR_SETS}
    MOI.throw_if_not_valid(model, c)
    lower, upper = _bounds(s)
    info = _info(model, c)
    _set_variable_bound(model, info.column, lower, upper)
    return
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.SingleVariable, ::MOI.ZeroOne
)
    info = _info(model, f.variable)
    # See https://github.com/JuliaOpt/GLPKMathProgInterface.jl/pull/15
    # for why this is necesary. GLPK interacts weirdly with binary variables and
    # bound modification. So let's set binary variables as "Integer" with [0,1]
    # bounds that we enforce just before solve.
    GLPK.set_col_kind(model.inner, info.column, GLPK.IV)
    info.type = BINARY
    model.num_binaries += 1
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}(f.variable.value)
end

function MOI.delete(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    GLPK.set_col_kind(model.inner, info.column, GLPK.CV)
    info.type = CONTINUOUS
    model.num_binaries -= 1
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}
)
    MOI.throw_if_not_valid(model, c)
    return MOI.ZeroOne()
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.SingleVariable, ::MOI.Integer
)
    info = _info(model, f.variable)
    GLPK.set_col_kind(model.inner, info.column, GLPK.IV)
    info.type = INTEGER
    model.num_integers += 1
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer}(f.variable.value)
end

function MOI.delete(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer}
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    GLPK.set_col_kind(model.inner, info.column, GLPK.CV)
    info.type = CONTINUOUS
    model.num_integers -= 1
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer}
)
    MOI.throw_if_not_valid(model, c)
    return MOI.Integer()
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where {S}
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    if S <: MOI.LessThan
        return info.lessthan_name
    elseif S <: Union{MOI.GreaterThan, MOI.Interval, MOI.EqualTo}
        return info.greaterthan_interval_or_equalto_name
    else
        @assert S <: Union{MOI.ZeroOne, MOI.Integer}
        return info.type_constraint_name
    end
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}, name::String
) where {S}
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    if S <: MOI.LessThan
        old_name = info.lessthan_name
        info.lessthan_name = name
    elseif S <: Union{MOI.GreaterThan, MOI.Interval, MOI.EqualTo}
        old_name = info.greaterthan_interval_or_equalto_name
        info.greaterthan_interval_or_equalto_name = name
    else
        @assert S <: Union{MOI.ZeroOne, MOI.Integer}
        old_name = info.type_constraint_name
        info.type_constraint_name = name
    end
    model.name_to_constraint_index = nothing
    return
end

###
### ScalarAffineFunction-in-Set
###

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}}
)
    key = ConstraintKey(c.value)
    if haskey(model.affine_constraint_info, key)
        return model.affine_constraint_info[key]
    end
    throw(MOI.InvalidIndex(c))
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}
) where {S}
    key = ConstraintKey(c.value)
    if haskey(model.affine_constraint_info, key)
        info = model.affine_constraint_info[key]
        return typeof(info.set) == S
    else
        return false
    end
end

"""
    _add_affine_constraint(
        problem::GLPK.Prob, columns::Vector{Int}, coefficients::Vector{Float64},
        sense::Cchar, rhs::Float64
    )

Helper function to add a row to the problem. Sense must be one of `'E'`
(ax == b), `'G'` (ax >= b), `'L'` (ax <= b).
"""
function _add_affine_constraint(
    problem::GLPK.Prob, columns::Vector{Int}, coefficients::Vector{Float64},
    sense::Cchar, rhs::Float64
)
    if length(columns) != length(coefficients)
        error("columns and coefficients have different lengths.")
    end
    GLPK.add_rows(problem, 1)
    num_rows = GLPK.get_num_rows(problem)
    GLPK.set_mat_row(problem, num_rows, columns, coefficients)
    # According to http://most.ccib.rutgers.edu/glpk.pdf page 22, the `lb`
    # argument is ignored for constraint types with no lower bound (GLPK.UP) and
    # the `ub` argument is ignored for constraint types with no upper bound
    # (GLPK.LO). We pass ±DBL_MAX for those unused bounds since (a) we have to
    # pass something, and (b) it is consistent with the other usages of ±DBL_MAX
    # to represent infinite bounds in the rest of the GLPK interface.
    if sense == Cchar('E')
        GLPK.set_row_bnds(problem, num_rows, GLPK.FX, rhs, rhs)
    elseif sense == Cchar('G')
        GLPK.set_row_bnds(problem, num_rows, GLPK.LO, rhs, GLPK.DBL_MAX)
    else
        @assert sense == Cchar('L')
        GLPK.set_row_bnds(problem, num_rows, GLPK.UP, -GLPK.DBL_MAX, rhs)
    end
    return
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.ScalarAffineFunction{Float64},
    s::Union{MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.EqualTo{Float64}}
)
    if !iszero(f.constant)
        throw(MOI.ScalarFunctionConstantNotZero{Float64, typeof(f), typeof(s)}(f.constant))
    end
    key = CleverDicts.add_item(model.affine_constraint_info, ConstraintInfo(s))
    model.affine_constraint_info[key].row = length(model.affine_constraint_info)
    indices, coefficients = _indices_and_coefficients(model, f)
    sense, rhs = _sense_and_rhs(s)
    _add_affine_constraint(model.inner, indices, coefficients, sense, rhs)
    return MOI.ConstraintIndex{typeof(f), typeof(s)}(key.value)
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    row = _info(model, c).row
    GLPK.std_basis(model.inner)
    GLPK.del_rows(model.inner, 1, [row])
    for info in values(model.affine_constraint_info)
        if info.row > row
            info.row -= 1
        end
    end
    key = ConstraintKey(c.value)
    delete!(model.affine_constraint_info, key)
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}
) where {S}
    row = _info(model, c).row
    sense = GLPK.get_row_type(model.inner, row)
    if sense == GLPK.LO || sense == GLPK.FX || sense == GLPK.DB
        return S(GLPK.get_row_lb(model.inner, row))
    else
        return S(GLPK.get_row_ub(model.inner, row))
    end
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}, s::S
) where {S <: Union{MOI.LessThan, MOI.GreaterThan, MOI.EqualTo}}
    row = _info(model, c).row
    if S <: MOI.LessThan
        GLPK.set_row_bnds(model.inner, row, GLPK.UP, -GLPK.DBL_MAX, s.upper)
    elseif S <: MOI.GreaterThan
        GLPK.set_row_bnds(model.inner, row, GLPK.LO, s.lower, GLPK.DBL_MAX)
    else
        @assert S <: MOI.EqualTo
        GLPK.set_row_bnds(model.inner, row, GLPK.FX, s.value, s.value)
    end
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}
) where {S}
    indices, values = GLPK.get_mat_row(model.inner, _info(model, c).row)
    terms = MOI.ScalarAffineTerm{Float64}[]
    for (col, val) in zip(indices, values)
        iszero(val) && continue
        push!(
            terms,
            MOI.ScalarAffineTerm(
                val,
                model.variable_info[CleverDicts.LinearIndex(col)].index
            )
        )
    end
    return MOI.ScalarAffineFunction(terms, 0.0)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any},
    name::String
)
    info = _info(model, c)
    old_name = info.name
    info.name = name
    if !isempty(name) && isascii(name)
        # Note: GLPK errors if we try to set non-ascii row names.
        GLPK.set_row_name(model.inner, info.row, name)
    end
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(model::Optimizer, ::Type{MOI.ConstraintIndex}, name::String)
    if model.name_to_constraint_index === nothing
        _rebuild_name_to_constraint_index(model)
    end
    if haskey(model.name_to_constraint_index, name)
        constr = model.name_to_constraint_index[name]
        if constr === nothing
            error("Duplicate constraint name detected: $(name)")
        end
        return constr
    end
    return nothing
end

function MOI.get(
    model::Optimizer, C::Type{MOI.ConstraintIndex{F, S}}, name::String
) where {F, S}
    index = MOI.get(model, MOI.ConstraintIndex, name)
    if typeof(index) == C
        return index::MOI.ConstraintIndex{F, S}
    end
    return nothing
end

function _rebuild_name_to_constraint_index(model::Optimizer)
    model.name_to_constraint_index = Dict{String, Union{Nothing, MOI.ConstraintIndex}}()
    for (key, info) in model.affine_constraint_info
        if isempty(info.name)
            continue
        end
        _set_name_to_constraint_index(
            model,
            info.name,
            MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, typeof(info.set)}(key.value)
        )
    end
    for (key, info) in model.variable_info
        if !isempty(info.lessthan_name)
            _set_name_to_constraint_index(
                model,
                info.lessthan_name,
                MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}(key.value)
            )
        end
        if !isempty(info.greaterthan_interval_or_equalto_name)
            S = if info.bound == GREATER_THAN || info.bound == LESS_AND_GREATER_THAN
                MOI.GreaterThan{Float64}
            elseif info.bound == EQUAL_TO
                MOI.EqualTo{Float64}
            else
                @assert info.bound == INTERVAL
                MOI.Interval{Float64}
            end
            _set_name_to_constraint_index(
                model,
                info.greaterthan_interval_or_equalto_name,
                MOI.ConstraintIndex{MOI.SingleVariable, S}(key.value)
            )
        end
        if !isempty(info.type_constraint_name)
            S = if info.type == BINARY
                MOI.ZeroOne
            else
                @assert info.type == INTEGER
                MOI.Integer
            end
            _set_name_to_constraint_index(
                model,
                info.type_constraint_name,
                MOI.ConstraintIndex{MOI.SingleVariable, S}(key.value)
            )
        end
    end
    return
end

function _set_name_to_constraint_index(
    model::Optimizer, name::String, index::MOI.ConstraintIndex
)
    if haskey(model.name_to_constraint_index, name)
        model.name_to_constraint_index[name] = nothing
    else
        model.name_to_constraint_index[name] = index
    end
    return
end

###
### Optimize methods.
###

function _solve_linear_problem(model::Optimizer)
    model.last_solved_by_mip = false
    if model.method == SIMPLEX
        model.solver_status = GLPK.simplex(model.inner, model.simplex_param)
    elseif model.method == EXACT
        model.solver_status = GLPK.exact(model.inner, model.simplex_param)
    else
        @assert model.method == INTERIOR
        model.solver_status = GLPK.interior(model.inner, model.interior_param)
    end
    return
end

"""
    _round_bounds_to_integer(model)

GLPK does not allow integer variables with fractional bounds. Therefore, we
round the bounds of binary and integer variables to integer values prior to
solving.

Returns a tuple `(column, lower, upper)` for the bounds that need to be reset.
"""
function _round_bounds_to_integer(model::Optimizer)
    bounds_to_reset = Tuple{Int, Float64, Float64}[]
    for (key, info) in model.variable_info
        if info.type == BINARY || info.type == INTEGER
            lb = GLPK.get_col_lb(model.inner, info.column)
            ub = GLPK.get_col_ub(model.inner, info.column)
            new_lb = ceil(lb)
            new_ub = floor(ub)
            if info.type == BINARY
                new_lb = max(0.0, new_lb)
                new_ub = min(1.0, new_ub)
            end
            if new_lb != lb || new_ub != ub
                push!(bounds_to_reset, (info.column, lb, ub))
                _set_variable_bound(model, info.column, new_lb, new_ub)
            end
        end
    end
    return bounds_to_reset
end

function _solve_mip_problem(model::Optimizer)
    bounds_to_reset = _round_bounds_to_integer(model)
    # Because we're muddling with the presolve in this function, cache the
    # original setting so that it can be reset.
    presolve_cache = model.intopt_param.presolve
    try
        # GLPK.intopt requires a starting basis for the LP relaxation. There are
        # two ways to get this. If presolve=GLPK.ON, then the presolve will find
        # a basis. If presolve=GLPK.OFF, then we should solve the problem via
        # GLPK.simplex first.
        if model.intopt_param.presolve == GLPK.OFF
            GLPK.simplex(model.inner, model.simplex_param)
            if GLPK.get_status(model.inner) != GLPK.OPT
                # We didn't find an optimal solution to the LP relaxation, so
                # let's turn presolve on and let intopt figure out what the
                # problem is.
                model.intopt_param.presolve = GLPK.ON
            end
        end
        model.solver_status = GLPK.intopt(model.inner, model.intopt_param)
        model.last_solved_by_mip = true

        # !!! Very Important Note
        #
        #  This next bit is _very_ important! See the note associated with
        #  set_callback.
        if model.callback_data.exception !== nothing
            throw(model.callback_data.exception)
        end

    finally
        for (column, lower, upper) in bounds_to_reset
            _set_variable_bound(model, column, lower, upper)
        end
        model.intopt_param.presolve = presolve_cache
    end
    return
end

include("infeasibility_certificates.jl")

function check_moi_callback_validity(model::Optimizer)
    has_moi_callback =
        model.lazy_callback !== nothing ||
        model.user_cut_callback !== nothing ||
        model.heuristic_callback !== nothing
    if has_moi_callback && model.has_generic_callback
        error("Cannot use `GLPK.CallbackFunction` as well as `MOI.AbstractCallbackFunction`.")
    end
    return has_moi_callback
end

function MOI.optimize!(model::Optimizer)
    start_time = time()
    model.optimize_not_called = false
    model.infeasibility_cert = nothing
    model.unbounded_ray = nothing

    # Initialize callbacks if necessary.
    if check_moi_callback_validity(model)
        MOI.set(model, CallbackFunction(), default_moi_callback(model))
        model.has_generic_callback = false
    end

    if model.num_binaries > 0 || model.num_integers > 0
        _solve_mip_problem(model)
    else
        _solve_linear_problem(model)
    end
    if MOI.get(model, MOI.ResultCount()) > 0
        if MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
            model.unbounded_ray = fill(NaN, GLPK.get_num_cols(model.inner))
            get_unbounded_ray(model, model.unbounded_ray)
        end
        if MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
            model.infeasibility_cert = fill(NaN, GLPK.get_num_rows(model.inner))
            get_infeasibility_ray(model, model.infeasibility_cert)
        end
    end
    model.solve_time = time() - start_time
    return
end

function _throw_if_optimize_in_progress(model::Optimizer, attr)
    if model.callback_state != CB_NONE
        throw(MOI.OptimizeInProgress(attr))
    end
end

# GLPK has a complicated status reporting system because it can be solved via
# multiple different solution algorithms. Regardless of the algorithm, the
# return value is stored in `model.solver_status`.
#
# Note that the first status (`Int32(0)`) should map to a `SUCCESS` status,
# because it doesn't imply anything about the solution. If `solver_status` is
# `Int32(0)`, then a solution-specific status can be queried with `_get_status`.

const RAW_SIMPLEX_STRINGS = Dict{Int32, Tuple{MOI.TerminationStatusCode, String}}(
    GLPK.EBADB  => (MOI.INVALID_MODEL,   "Unable to start the search, because the initial basis specified in the problem object is invalid—the number of basic (auxiliary and structural) variables is not the same as the number of rows in the problem object."),
    GLPK.ESING  => (MOI.NUMERICAL_ERROR, "Unable to start the search, because the basis matrix corresponding to the initial basis is singular within the working precision."),
    GLPK.ECOND  => (MOI.NUMERICAL_ERROR, "Unable to start the search, because the basis matrix corresponding to the initial basis is ill-conditioned, i.e. its condition number is too large."),
    GLPK.EBOUND => (MOI.INVALID_MODEL,   "Unable to start the search, because some double-bounded (auxiliary or structural) variables have incorrect bounds."),
    GLPK.EFAIL  => (MOI.NUMERICAL_ERROR, "The search was prematurely terminated due to the solver failure."),
    GLPK.EOBJLL => (MOI.OBJECTIVE_LIMIT, "The search was prematurely terminated, because the objective function being maximized has reached its lower limit and continues decreasing (the dual simplex only)."),
    GLPK.EOBJUL => (MOI.OBJECTIVE_LIMIT, "The search was prematurely terminated, because the objective function being minimized has reached its upper limit and continues increasing (the dual simplex only)."),
    GLPK.EITLIM => (MOI.ITERATION_LIMIT, "The search was prematurely terminated, because the simplex iteration limit has been exceeded."),
    GLPK.ETMLIM => (MOI.TIME_LIMIT,      "The search was prematurely terminated, because the time limit has been exceeded."),
    GLPK.ENOPFS => (MOI.INFEASIBLE,      "The LP problem instance has no primal feasible solution (only if the LP presolver is used)."),
    GLPK.ENODFS => (MOI.DUAL_INFEASIBLE, "The LP problem instance has no dual feasible solution (only if the LP presolver is used).")
)

const RAW_EXACT_STRINGS = Dict{Int32, Tuple{MOI.TerminationStatusCode, String}}(
    GLPK.EBADB  => (MOI.INVALID_MODEL,   "Unable to start the search, because the initial basis specified in the problem object is invalid—the number of basic (auxiliary and structural) variables is not the same as the number of rows in the problem object."),
    GLPK.ESING  => (MOI.NUMERICAL_ERROR, "Unable to start the search, because the basis matrix corresponding to the initial basis is exactly singular."),
    GLPK.EBOUND => (MOI.INVALID_MODEL,   "Unable to start the search, because some double-bounded (auxiliary or structural) variables have incorrect bounds."),
    GLPK.EFAIL  => (MOI.INVALID_MODEL,   "The problem instance has no rows/columns."),
    GLPK.EITLIM => (MOI.ITERATION_LIMIT, "The search was prematurely terminated, because the simplex iteration limit has been exceeded."),
    GLPK.ETMLIM => (MOI.TIME_LIMIT,      "The search was prematurely terminated, because the time limit has been exceeded.")
)

const RAW_INTERIOR_STRINGS = Dict{Int32, Tuple{MOI.TerminationStatusCode, String}}(
    GLPK.EFAIL   => (MOI.INVALID_MODEL,   "The problem instance has no rows/columns."),
    GLPK.ENOCVG  => (MOI.SLOW_PROGRESS,   "Very slow convergence or divergence."),
    GLPK.EITLIM  => (MOI.ITERATION_LIMIT, "Iteration limit exceeded."),
    GLPK.EINSTAB => (MOI.NUMERICAL_ERROR, "Numerical instability on solving Newtonian system.")
)

const RAW_INTOPT_STRINGS = Dict{Int32, Tuple{MOI.TerminationStatusCode, String}}(
    GLPK.EBOUND  => (MOI.INVALID_MODEL,   "Unable to start the search, because some double-bounded (auxiliary or structural) variables have incorrect bounds."),
    GLPK.ENOPFS  => (MOI.INFEASIBLE,      "Unable to start the search, because LP relaxation of the MIP problem instance has no primal feasible solution. (This code may appear only if the presolver is enabled.)"),
    GLPK.ENODFS  => (MOI.DUAL_INFEASIBLE, "Unable to start the search, because LP relaxation of the MIP problem instance has no dual feasible solution. In other word, this code means that if the LP relaxation has at least one primal feasible solution, its optimal solution is unbounded, so if the MIP problem has at least one integer feasible solution, its (integer) optimal solution is also unbounded. (This code may appear only if the presolver is enabled.)"),
    GLPK.EFAIL   => (MOI.INVALID_MODEL,   "The search was prematurely terminated due to the solver failure."),
    GLPK.EMIPGAP => (MOI.OPTIMAL,         "The search was prematurely terminated, because the relative mip gap tolerance has been reached."),
    GLPK.ETMLIM  => (MOI.TIME_LIMIT,      "The search was prematurely terminated, because the time limit has been exceeded."),
    GLPK.ESTOP   => (MOI.INTERRUPTED,     "The search was prematurely terminated by application. (This code may appear only if the advanced solver interface is used.)")
)

const RAW_SOLUTION_STATUS = Dict{Int32, Tuple{MOI.TerminationStatusCode, String}}(
    GLPK.OPT    => (MOI.OPTIMAL,            "Solution is optimal"),
    GLPK.FEAS   => (MOI.LOCALLY_SOLVED,     "Solution is feasible"),
    GLPK.INFEAS => (MOI.LOCALLY_INFEASIBLE, "Solution is infeasible"),
    GLPK.NOFEAS => (MOI.INFEASIBLE,         "No feasible primal-dual solution exists."),
    GLPK.UNBND  => (MOI.DUAL_INFEASIBLE,    "Problem has unbounded solution"),
    GLPK.UNDEF  => (MOI.OTHER_ERROR,        "Solution is undefined")
)

function MOI.get(model::Optimizer, attr::MOI.RawStatusString)
    _throw_if_optimize_in_progress(model, attr)
    if model.solver_status == Int32(0)
        (_, msg) = _get_status(model)
        return msg
    elseif model.last_solved_by_mip
        return RAW_INTOPT_STRINGS[model.solver_status][2]
    elseif model.method == SIMPLEX
        return RAW_SIMPLEX_STRINGS[model.solver_status][2]
    elseif model.method == EXACT
        return RAW_EXACT_STRINGS[model.solver_status][2]
    else
        @assert model.method == INTERIOR
        return RAW_INTERIOR_STRINGS[model.solver_status][2]
    end
end

function _get_status(model::Optimizer)
    status_code = if model.last_solved_by_mip
        GLPK.mip_status(model.inner)
    elseif model.method == SIMPLEX || model.method == EXACT
        GLPK.get_status(model.inner)
    else
        @assert model.method == INTERIOR
        GLPK.ipt_status(model.inner)
    end
    return RAW_SOLUTION_STATUS[status_code]
end

"""
    _certificates_potentially_available(model::Optimizer)

Return true if an infeasiblity certificate or an unbounded ray is potentially
available (i.e., the model has been solved using either the Simplex or Exact
methods).
"""
function _certificates_potentially_available(model::Optimizer)
    return !model.last_solved_by_mip && (model.method == SIMPLEX || model.method == EXACT)
end

function MOI.get(model::Optimizer, attr::MOI.TerminationStatus)
    _throw_if_optimize_in_progress(model, attr)
    if model.optimize_not_called
        return MOI.OPTIMIZE_NOT_CALLED
    elseif model.solver_status != Int32(0)
        # The solver did not exit successfully for some reason.
        if model.last_solved_by_mip
            return RAW_INTOPT_STRINGS[model.solver_status][1]
        elseif model.method == SIMPLEX
            return RAW_SIMPLEX_STRINGS[model.solver_status][1]
        elseif model.method == INTERIOR
            return RAW_INTERIOR_STRINGS[model.solver_status][1]
        else
            @assert model.method == EXACT
            return RAW_EXACT_STRINGS[model.solver_status][1]
        end
    else
        (status, _) = _get_status(model)
        return status
    end
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    _throw_if_optimize_in_progress(model, attr)
    if attr.N != 1
        return MOI.NO_SOLUTION
    end
    (status, _) = _get_status(model)
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        return MOI.FEASIBLE_POINT
    elseif status == MOI.LOCALLY_INFEASIBLE
        return MOI.INFEASIBLE_POINT
    elseif status == MOI.DUAL_INFEASIBLE
        if _certificates_potentially_available(model)
            return MOI.INFEASIBILITY_CERTIFICATE
        end
    else
        @assert status == MOI.INFEASIBLE || status == MOI.OTHER_ERROR
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    _throw_if_optimize_in_progress(model, attr)
    if attr.N != 1 || model.last_solved_by_mip
        return MOI.NO_SOLUTION
    end
    (status, _) = _get_status(model)
    if status == MOI.OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif status == MOI.INFEASIBLE || status == MOI.LOCALLY_INFEASIBLE
        if _certificates_potentially_available(model)
            return MOI.INFEASIBILITY_CERTIFICATE
        end
    end
    return MOI.NO_SOLUTION
end

function _get_col_dual(model::Optimizer, column::Int)
    @assert !model.last_solved_by_mip
    if model.method == SIMPLEX || model.method == EXACT
        return _dual_multiplier(model) * GLPK.get_col_dual(model.inner, column)
    else
        @assert model.method == INTERIOR
        return _dual_multiplier(model) * GLPK.ipt_col_dual(model.inner, column)
    end
end

function _get_col_primal(model::Optimizer, column::Int)
    if model.last_solved_by_mip
        return GLPK.mip_col_val(model.inner, column)
    elseif model.method == SIMPLEX || model.method == EXACT
        return GLPK.get_col_prim(model.inner, column)
    else
        @assert model.method == INTERIOR
        return GLPK.ipt_col_prim(model.inner, column)
    end
end

function _get_row_primal(model::Optimizer, row::Int)
    if model.last_solved_by_mip
        return GLPK.mip_row_val(model.inner, row)
    elseif model.method == SIMPLEX || model.method == EXACT
        return GLPK.get_row_prim(model.inner, row)
    else
        @assert model.method == INTERIOR
        return GLPK.ipt_row_prim(model.inner, row)
    end
end

function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, x::MOI.VariableIndex)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    if model.unbounded_ray !== nothing
        return model.unbounded_ray[_info(model, x).column]
    else
        return _get_col_primal(model, _info(model, x).column)
    end
end

function MOI.get(
    model::Optimizer, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.SingleVariable, <:Any}
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end

function MOI.get(
    model::Optimizer, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    return _get_row_primal(model, _info(model, c).row)
end

function _dual_multiplier(model::Optimizer)
    return MOI.get(model, MOI.ObjectiveSense()) == MOI.MIN_SENSE ? 1.0 : -1.0
end

function MOI.get(
    model::Optimizer, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    column = _info(model, c).column
    reduced_cost = if model.method == SIMPLEX || model.method == EXACT
        GLPK.get_col_dual(model.inner, column)
    else
        @assert model.method == INTERIOR
        GLPK.ipt_col_dual(model.inner, column)
    end
    sense = MOI.get(model, MOI.ObjectiveSense())
    # The following is a heuristic for determining whether the reduced cost
    # (i.e., the column dual) applies to the lower or upper bound. It can be
    # wrong by at most `tol_dj`.
    if sense == MOI.MIN_SENSE && reduced_cost < 0
        # If minimizing, the reduced cost must be negative (ignoring
        # tolerances).
        return reduced_cost
    elseif sense == MOI.MAX_SENSE && reduced_cost > 0
        # If minimizing, the reduced cost must be positive (ignoring
        # tolerances). However, because of the MOI dual convention, we return a
        # negative value.
        return -reduced_cost
    else
        # The reduced cost, if non-zero, must related to the lower bound.
        return 0.0
    end
end

function MOI.get(
    model::Optimizer, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    column = _info(model, c).column
    reduced_cost = if model.method == SIMPLEX || model.method == EXACT
        GLPK.get_col_dual(model.inner, column)
    else
        @assert model.method == INTERIOR
        GLPK.ipt_col_dual(model.inner, column)
    end
    sense = MOI.get(model, MOI.ObjectiveSense())
    # The following is a heuristic for determining whether the reduced cost
    # (i.e., the column dual) applies to the lower or upper bound. It can be
    # wrong by at most `tol_dj`.
    if sense == MOI.MIN_SENSE && reduced_cost > 0
        # If minimizing, the reduced cost must be negative (ignoring
        # tolerances).
        return reduced_cost
    elseif sense == MOI.MAX_SENSE && reduced_cost < 0
        # If minimizing, the reduced cost must be positive (ignoring
        # tolerances). However, because of the MOI dual convention, we return a
        # negative value.
        return -reduced_cost
    else
        # The reduced cost, if non-zero, must related to the lower bound.
        return 0.0
    end
end

function MOI.get(
    model::Optimizer, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where {S <: Union{MOI.EqualTo, MOI.Interval}}
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    return _get_col_dual(model, _info(model, c).column)
end

function MOI.get(
    model::Optimizer, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    row = _info(model, c).row
    if model.infeasibility_cert !== nothing
        return model.infeasibility_cert[row]
    else
        @assert !model.last_solved_by_mip
        if model.method == SIMPLEX || model.method == EXACT
            return _dual_multiplier(model) * GLPK.get_row_dual(model.inner, row)
        else
            @assert model.method == INTERIOR
            return _dual_multiplier(model) * GLPK.ipt_row_dual(model.inner, row)
        end
    end
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    if model.last_solved_by_mip
        return GLPK.mip_obj_val(model.inner)
    elseif model.method == SIMPLEX || model.method == EXACT
        return GLPK.get_obj_val(model.inner)
    else
        @assert model.method == INTERIOR
        return GLPK.ipt_obj_val(model.inner)
    end
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveBound)
    _throw_if_optimize_in_progress(model, attr)
    if !model.last_solved_by_mip
        return MOI.get(model, MOI.ObjectiveSense()) == MOI.MIN_SENSE ? -Inf : Inf
    end
    # @mlubin and @ccoey observed some cases where mip_status == OPT and objval
    # and objbound didn't match. In that case, they return mip_obj_val, but
    # objbound may still be incorrect in cases where GLPK terminates early.
    if GLPK.mip_status(model.inner) == GLPK.OPT
        return GLPK.mip_obj_val(model.inner)
    end
    return model.objective_bound
end

function MOI.get(model::Optimizer, attr::MOI.DualObjectiveValue)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    return MOI.Utilities.get_fallback(model, attr, Float64)
end

function MOI.get(model::Optimizer, attr::MOI.RelativeGap)
    _throw_if_optimize_in_progress(model, attr)
    if !model.last_solved_by_mip
        error("RelativeGap is only available for models with integer variables.")
    end
    return model.relative_gap
end

function MOI.get(model::Optimizer, attr::MOI.SolveTime)
    _throw_if_optimize_in_progress(model, attr)
    return model.solve_time
end

function MOI.get(model::Optimizer, attr::MOI.ResultCount)
    _throw_if_optimize_in_progress(model, attr)
    (status, _) = _get_status(model)
    if status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.LOCALLY_INFEASIBLE)
        return 1
    elseif status in (MOI.DUAL_INFEASIBLE, MOI. INFEASIBLE, MOI.LOCALLY_INFEASIBLE)
        if _certificates_potentially_available(model)
            return 1
        end
    end
    return 0
end

function MOI.get(model::Optimizer, ::MOI.Silent)
    return model.silent
end

function MOI.set(model::Optimizer, ::MOI.Silent, flag::Bool)
    model.silent = flag
    output_flag = flag ? GLPK.OFF : MOI.get(model, MOI.RawParameter("msg_lev"))
    MOI.set(model, MOI.RawParameter("msg_lev"), output_flag)
    return
end

function MOI.get(model::Optimizer, ::MOI.Name)
    return GLPK.get_prob_name(model.inner)
end

function MOI.set(model::Optimizer, ::MOI.Name, name::String)
    GLPK.set_prob_name(model.inner, name)
    return
end

MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = length(model.variable_info)
function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return sort!(collect(keys(model.variable_info)), by = x -> x.value)
end

MOI.get(model::Optimizer, ::MOI.RawSolver) = model.inner

function MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{F, S}) where {F, S}
    # TODO: this could be more efficient.
    return length(MOI.get(model, MOI.ListOfConstraintIndices{F, S}()))
end

_bound_enums(::Type{<:MOI.LessThan}) = (LESS_THAN, LESS_AND_GREATER_THAN)
_bound_enums(::Type{<:MOI.GreaterThan}) = (GREATER_THAN, LESS_AND_GREATER_THAN)
_bound_enums(::Type{<:MOI.Interval}) = (INTERVAL,)
_bound_enums(::Type{<:MOI.EqualTo}) = (EQUAL_TO,)
_bound_enums(::Any) = (nothing,)

_type_enums(::Type{MOI.ZeroOne}) = (BINARY,)
_type_enums(::Type{MOI.Integer}) = (INTEGER,)
_type_enums(::Any) = (nothing,)

function MOI.get(
    model::Optimizer, ::MOI.ListOfConstraintIndices{MOI.SingleVariable, S}
) where {S}
    indices = MOI.ConstraintIndex{MOI.SingleVariable, S}[]
    for (key, info) in model.variable_info
        if info.bound in _bound_enums(S) || info.type in _type_enums(S)
            push!(indices, MOI.ConstraintIndex{MOI.SingleVariable, S}(key.value))
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64}, S}
) where {S}
    indices = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}[]
    for (key, info) in model.affine_constraint_info
        if typeof(info.set) == S
            push!(indices, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}(key.value))
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(model::Optimizer, ::MOI.ListOfConstraints)
    constraints = Set{Tuple{DataType, DataType}}()
    for info in values(model.variable_info)
        if info.bound == NONE
        elseif info.bound == LESS_THAN
            push!(constraints, (MOI.SingleVariable, MOI.LessThan{Float64}))
        elseif info.bound == GREATER_THAN
            push!(constraints, (MOI.SingleVariable, MOI.GreaterThan{Float64}))
        elseif info.bound == LESS_AND_GREATER_THAN
            push!(constraints, (MOI.SingleVariable, MOI.LessThan{Float64}))
            push!(constraints, (MOI.SingleVariable, MOI.GreaterThan{Float64}))
        elseif info.bound == EQUAL_TO
            push!(constraints, (MOI.SingleVariable, MOI.EqualTo{Float64}))
        elseif info.bound == INTERVAL
            push!(constraints, (MOI.SingleVariable, MOI.Interval{Float64}))
        end
        if info.type == CONTINUOUS
        elseif info.type == BINARY
            push!(constraints, (MOI.SingleVariable, MOI.ZeroOne))
        elseif info.type == INTEGER
            push!(constraints, (MOI.SingleVariable, MOI.Integer))
        end
    end
    for info in values(model.affine_constraint_info)
        push!(constraints, (MOI.ScalarAffineFunction{Float64}, typeof(info.set)))
    end
    return collect(constraints)
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunctionType)
    if model.objective_type == SINGLE_VARIABLE
        return MOI.SingleVariable
    else
        @assert model.objective_type == SCALAR_AFFINE
        return MOI.ScalarAffineFunction{Float64}
    end
end

# TODO(odow): is there a way to modify a single element, rather than the whole
# row?
function MOI.modify(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any},
    chg::MOI.ScalarCoefficientChange{Float64}
)
    row = _info(model, c).row
    col = _info(model, chg.variable).column
    columns, coefficients = GLPK.get_mat_row(model.inner, row)
    index = something(findfirst(isequal(col), columns), 0)
    if index > 0
        coefficients[index] = chg.new_coefficient
    else
        push!(columns, col)
        push!(coefficients, chg.new_coefficient)
    end
    GLPK.set_mat_row(model.inner, row, columns, coefficients)
    return
end

function MOI.modify(
    model::Optimizer,
    c::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarCoefficientChange{Float64}
)
    GLPK.set_obj_coef(
        model.inner, _info(model, chg.variable).column, chg.new_coefficient
    )
    return
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:SCALAR_SETS},
    f::MOI.ScalarAffineFunction{Float64}
)
    if !iszero(f.constant)
        throw(MOI.ScalarFunctionConstantNotZero(f.constant))
    end
    row = _info(model, c).row
    indices, coefficients = _indices_and_coefficients(model, f)
    GLPK.set_mat_row(model.inner, row, indices, coefficients)
    return
end

function MOI.get(
    model::Optimizer, attr::MOI.ConstraintBasisStatus,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}
) where {S <: SCALAR_SETS}
    _throw_if_optimize_in_progress(model, attr)
    row = _info(model, c).row
    cbasis = GLPK.get_row_stat(model.inner, row)
    if cbasis == GLPK.BS
        return MOI.BASIC
    elseif cbasis == GLPK.NL || cbasis == GLPK.NU || cbasis == GLPK.NF || cbasis == GLPK.NS
        return MOI.NONBASIC
    else
        error("CBasis value of $(cbasis) isn't defined.")
    end
end

function MOI.get(
    model::Optimizer, attr::MOI.ConstraintBasisStatus,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where {S <: SCALAR_SETS}
    _throw_if_optimize_in_progress(model, attr)
    column = _info(model, c).column
    vbasis = GLPK.get_col_stat(model.inner, column)
    if vbasis == GLPK.BS
        return MOI.BASIC
    elseif vbasis == GLPK.NL
        if S <: MOI.LessThan
            return MOI.BASIC
        elseif !(S <: MOI.Interval)
            return MOI.NONBASIC
        else
            return MOI.NONBASIC_AT_LOWER
        end
    elseif vbasis == GLPK.NU
        MOI.NONBASIC_AT_UPPER
        if S <: MOI.GreaterThan
            return MOI.BASIC
        elseif !(S <: MOI.Interval)
            return MOI.NONBASIC
        else
            return MOI.NONBASIC_AT_UPPER
        end
    elseif vbasis == GLPK.NF
        return MOI.NONBASIC
    elseif vbasis == GLPK.NS
        return MOI.NONBASIC
    else
        error("VBasis value of $(vbasis) isn't defined.")
    end
end
