# =======================
#   `copy_to` function
# =======================

const DoubleDicts = MathOptInterface.Utilities.DoubleDicts

_add_bounds(::Vector{Float64}, ub, i, s::MOI.LessThan{Float64}) = ub[i] = s.upper
_add_bounds(lb, ::Vector{Float64}, i, s::MOI.GreaterThan{Float64}) = lb[i] = s.lower
_add_bounds(lb, ub, i, s::MOI.EqualTo{Float64}) = lb[i], ub[i] = s.value, s.value
_add_bounds(lb, ub, i, s::MOI.Interval{Float64}) = lb[i], ub[i] = s.lower, s.upper

_bound_type(::Type{MOI.Interval{Float64}}) = INTERVAL
_bound_type(::Type{MOI.EqualTo{Float64}}) = EQUAL_TO
_bound_type(::Type{S}) where S = NONE

function _extract_bound_data(src, mapping, lb, ub, bound_type, s::Type{S}) where S
    dict = DoubleDicts.with_type(mapping.conmap, MOI.SingleVariable, S)
    type = _bound_type(s)
    list = MOI.get(src, MOI.ListOfConstraintIndices{MOI.SingleVariable, S}())
    add_sizehint!(dict, length(list))
    for con_index in list
        f = MOI.get(src, MOI.ConstraintFunction(), con_index)
        s = MOI.get(src, MOI.ConstraintSet(), con_index)
        column = mapping[f.variable].value
        _add_bounds(lb, ub, column, s)
        bound_type[column] = type
        dict[con_index] = MOI.ConstraintIndex{MOI.SingleVariable, S}(column)
    end
end

_add_type(type, i, ::MOI.Integer) = type[i] = INTEGER
_add_type(type, i, ::MOI.ZeroOne) = type[i] = BINARY

function _extract_type_data(src, mapping, var_type, ::Type{S}) where S
    dict = DoubleDicts.with_type(mapping.conmap, MOI.SingleVariable, S)
    list = MOI.get(src, MOI.ListOfConstraintIndices{MOI.SingleVariable, S}())
    add_sizehint!(dict, length(list))
    for con_index in list
        f = MOI.get(src, MOI.ConstraintFunction(), con_index)
        column = mapping[f.variable].value
        _add_type(var_type, column, S())
        dict[con_index] = MOI.ConstraintIndex{MOI.SingleVariable, S}(column)
    end
end

function _init_index_map(src)
    x_src = MOI.get(src, MOI.ListOfVariableIndices())
    N = Cint(length(x_src))

    is_contiguous = true
    # assuming all indexes are different
    for x in x_src
        if !(1 <= x.value <= N)
            is_contiguous = false
        end
    end

    mapping = if is_contiguous
        MOIU.IndexMap(N)
    else
        MOIU.IndexMap()
    end
    for i = 1:N
        mapping[x_src[i]] = MOI.VariableIndex(i)
    end

    return N, mapping
end

_bounds2(s::MOI.GreaterThan{Float64}) = (s.lower, Inf)
_bounds2(s::MOI.LessThan{Float64}) = (-Inf, s.upper)
_bounds2(s::MOI.EqualTo{Float64}) = (s.value, s.value)
_bounds2(s::MOI.Interval{Float64}) = (s.lower, s.upper)

function add_sizehint!(vec, n)
    len = length(vec)
    return sizehint!(vec, len + n)
end

function _extract_row_data(src, mapping, lb, ub, I, J, V, ::Type{S}) where S
    dict = mapping.conmap[MOI.ScalarAffineFunction{Float64}, S]

    list = MOI.get(
        src, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64}, S}()
    )::Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}}

    N = length(list)
    add_sizehint!(lb, N)
    add_sizehint!(ub, N)
    add_sizehint!(dict, N)

    # first loop caches functions and counts terms to be added
    n_terms = 0
    function_cache = Array{MOI.ScalarAffineFunction{Float64}}(undef, N)
    for i in 1:N
        pre_function = MOI.get(src, MOI.ConstraintFunction(), list[i])
        f = if MOIU.is_canonical(pre_function)
            pre_function
        else
            # no duplicates are allowed in GLPK
            MOIU.canonical(pre_function)
        end
        function_cache[i] = f
        l, u = _bounds2(MOI.get(src, MOI.ConstraintSet(), list[i]))
        push!(lb, l - f.constant)
        push!(ub, u - f.constant)
        n_terms += length(f.terms)
    end

    non_zeros = length(I)
    row = non_zeros == 0 ? 1 : I[end] + 1
    # resize + setindex is faster than sizehint! + push
    # makes difference because I, J, V can be huge
    resize!(I, non_zeros + n_terms)
    resize!(J, non_zeros + n_terms)
    resize!(V, non_zeros + n_terms)
    for i in 1:N
        f = function_cache[i]
        for term in f.terms
            non_zeros += 1
            I[non_zeros] = row
            J[non_zeros] = Cint(mapping[term.variable_index].value::Int64)
            V[non_zeros] = term.coefficient
        end
        row += 1
        ind = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}(row)
        dict[list[i]] = ind
    end
    return
end

function test_data(src, dest)
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if !MOI.supports_constraint(dest, F, S)
            throw(MOI.UnsupportedConstraint{F, S}("GLPK.Optimizer does not support constraints of type $F-in-$S."))
        end
    end
    fobj_type = MOI.get(src, MOI.ObjectiveFunctionType())
    if !MOI.supports(dest, MOI.ObjectiveFunction{fobj_type}())
        throw(MOI.UnsupportedAttribute(MOI.ObjectiveFunction(fobj_type)))
    end
end

function _add_all_variables(model::Optimizer, N, lower, upper, bound_type,
    var_type)
    glp_add_cols(model, N)
    sizehint!(model.variable_info, N)

    for i in 1:N
        bound = get_moi_bound_type(lower[i], upper[i], bound_type)
        # We started from empty model.variable_info, hence we assume ordering
        index = CleverDicts.add_item(
            model.variable_info, VariableInfo(MOI.VariableIndex(i), i, bound, var_type[i])
        )
        glp_bound_type = get_glp_bound_type(lower[i], upper[i])
        glp_set_col_bnds(model, i, glp_bound_type, lower[i], upper[i])
        if var_type[i] == BINARY
            model.num_binaries += 1
        end
        if var_type[i] == INTEGER
            model.num_integers += 1
        end
    end
    return nothing
end

function _add_all_constraints(dest::Optimizer, rl, ru, I, J, V)

    n_constraints = length(rl)

    glp_add_rows(dest, n_constraints)
    glp_load_matrix(dest, length(I), offset(I), offset(J), offset(V))

    sizehint!(dest.affine_constraint_info, n_constraints)
    for i in 1:n_constraints
        # assume ordered indexing
        # assume no range constraints
        if rl[i] == ru[i]
            glp_set_row_bnds(dest, i, GLP_FX, rl[i], ru[i])
            CleverDicts.add_item(dest.affine_constraint_info,
                ConstraintInfo(i, MOI.EqualTo{Float64}(rl[i])))
        elseif ru[i] == Inf
            glp_set_row_bnds(dest, i, GLP_LO, rl[i], GLP_DBL_MAX)
            CleverDicts.add_item(dest.affine_constraint_info,
                ConstraintInfo(i, MOI.GreaterThan{Float64}(rl[i])))
        else
            glp_set_row_bnds(dest, i, GLP_UP, -GLP_DBL_MAX, ru[i])
            CleverDicts.add_item(dest.affine_constraint_info,
                ConstraintInfo(i, MOI.LessThan{Float64}(ru[i])))
        end
    end
    return
end

function MOI.copy_to(
# function _copy_to(
    dest::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false
)

    @assert MOI.is_empty(dest)
    test_data(src, dest)

    N, mapping = _init_index_map(src)
    cl, cu = fill(-Inf, N), fill(Inf, N)
    bound_type = fill(NONE, N)
    var_type = fill(CONTINUOUS, N)

    _extract_bound_data(src, mapping, cl, cu, bound_type, MOI.GreaterThan{Float64})
    _extract_bound_data(src, mapping, cl, cu, bound_type, MOI.LessThan{Float64})
    _extract_bound_data(src, mapping, cl, cu, bound_type, MOI.EqualTo{Float64})
    _extract_bound_data(src, mapping, cl, cu, bound_type, MOI.Interval{Float64})

    _extract_type_data(src, mapping, var_type, MOI.Integer)
    _extract_type_data(src, mapping, var_type, MOI.ZeroOne)

    _add_all_variables(dest, N, cl, cu, bound_type, var_type)

    rl, ru, I, J, V = Float64[], Float64[], Cint[], Cint[], Float64[]

    _extract_row_data(src, mapping, rl, ru, I, J, V, MOI.GreaterThan{Float64})
    _extract_row_data(src, mapping, rl, ru, I, J, V, MOI.LessThan{Float64})
    _extract_row_data(src, mapping, rl, ru, I, J, V, MOI.EqualTo{Float64})
    # range constraints not supported
    # _extract_row_data(src, mapping, rl, ru, I, J, V, MOI.Interval{Float64})

    _add_all_constraints(dest, rl, ru, I, J, V)

    # Copy model attributes:
    # obj function and sense are passed here
    MOIU.pass_attributes(dest, src, copy_names, mapping)
    variables = MOI.get(src, MOI.ListOfVariableIndices())
    MOIU.pass_attributes(dest, src, copy_names, mapping, variables)
    pass_constraint_attributes(dest, src, copy_names, mapping)
    return mapping
end

function pass_constraint_attributes(dest, src, copy_names, mapping)
    ctr_types = MOI.get(src, MOI.ListOfConstraints())
    for (F,S) in ctr_types
        pass_constraint_attributes(dest, src, copy_names, mapping, F, S)
    end
    return
end
function pass_constraint_attributes(dest, src, copy_names, mapping,
    ::Type{F}, ::Type{S}) where {F,S}
    indices = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    MOIU.pass_attributes(dest, src, copy_names, mapping, indices)
    return
end