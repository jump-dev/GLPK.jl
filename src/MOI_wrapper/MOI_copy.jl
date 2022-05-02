# Copyright (c) 2012 GLPK.jl contributors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the Licence, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

struct _OptimizerCache
    "Column lower bounds"
    cl::Vector{Float64}
    "Column upper bounds"
    cu::Vector{Float64}
    "Column bound types"
    bounds::Vector{_VariableBound}
    "Column types"
    types::Vector{_VariableType}
    "Row lower bounds"
    rl::Vector{Float64}
    "Row upper bounds"
    ru::Vector{Float64}
    "A matrix in 1-indexed sparse triplet form"
    I::Vector{Cint}
    J::Vector{Cint}
    V::Vector{Float64}
    function _OptimizerCache(N::Int)
        return new(
            fill(-Inf, N),
            fill(Inf, N),
            fill(_NONE, N),
            fill(_CONTINUOUS, N),
            Float64[],
            Float64[],
            Cint[],
            Cint[],
            Float64[],
        )
    end
end

"""
    _validate_constraint_types(dest::Optimizer, src::MOI.ModelLike)

Throw an error if unsupported constraint or objective types are present in
`src`.
"""
function _validate_constraint_types(dest::Optimizer, src::MOI.ModelLike)
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        if !MOI.supports_constraint(dest, F, S)
            throw(
                MOI.UnsupportedConstraint{F,S}(
                    "GLPK.Optimizer does not support constraints of type $F-in-$S.",
                ),
            )
        end
        for attr in MOI.get(src, MOI.ListOfConstraintAttributesSet{F,S}())
            if !MOI.supports(dest, attr, MOI.ConstraintIndex{F,S})
                throw(MOI.UnsupportedAttribute(attr))
            end
        end
    end
    for attr in MOI.get(src, MOI.ListOfModelAttributesSet())
        if !MOI.supports(dest, attr)
            throw(MOI.UnsupportedAttribute(attr))
        end
    end
    for attr in MOI.get(src, MOI.ListOfVariableAttributesSet())
        if !MOI.supports(dest, attr, MOI.VariableIndex)
            throw(MOI.UnsupportedAttribute(attr))
        end
    end
    return
end

function _init_index_map(src::MOI.ModelLike)
    variables = MOI.get(src, MOI.ListOfVariableIndices())
    map = MOI.Utilities.IndexMap()
    N = 0
    for x in variables
        N += 1
        map[x] = MOI.VariableIndex(N)
    end
    return variables, map
end

_add_set_data(cache, i, s::MOI.LessThan{Float64}) = (cache.cu[i] = s.upper)
_add_set_data(cache, i, s::MOI.GreaterThan{Float64}) = (cache.cl[i] = s.lower)

function _add_set_data(cache, i, s::MOI.EqualTo{Float64})
    cache.cl[i] = s.value
    cache.cu[i] = s.value
    cache.bounds[i] = _EQUAL_TO
    return
end

function _add_set_data(cache, i, s::MOI.Interval{Float64})
    cache.cl[i] = s.lower
    cache.cu[i] = s.upper
    cache.bounds[i] = _INTERVAL
    return
end

_add_set_data(cache, i, ::MOI.Integer) = (cache.types[i] = _INTEGER)

_add_set_data(cache, i, ::MOI.ZeroOne) = (cache.types[i] = _BINARY)

function _extract_variable_data(src, map, cache, ::Type{S}) where {S}
    ci_map = map.con_map[MOI.VariableIndex, S]
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{MOI.VariableIndex,S}())
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        column = map[f].value
        _add_set_data(cache, column, s)
        ci_map[ci] = MOI.ConstraintIndex{MOI.VariableIndex,S}(column)
    end
    return
end

function _extract_row_data(src, map, cache, ::Type{S}) where {S}
    F = MOI.ScalarAffineFunction{Float64}
    ci_map = map.con_map[F, S]
    nnz = length(cache.I)
    row = nnz == 0 ? 1 : cache.I[end] + 1
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        if !MOI.Utilities.is_canonical(f)
            f = MOI.Utilities.canonical(f)
        end
        l, u = _bounds(MOI.get(src, MOI.ConstraintSet(), ci))
        push!(cache.rl, l === nothing ? -Inf : l - f.constant)
        push!(cache.ru, u === nothing ? Inf : u - f.constant)
        resize!(cache.I, nnz + length(f.terms))
        resize!(cache.J, nnz + length(f.terms))
        resize!(cache.V, nnz + length(f.terms))
        for term in f.terms
            nnz += 1
            cache.I[nnz] = row
            cache.J[nnz] = Cint(map[term.variable].value::Int64)
            cache.V[nnz] = term.coefficient
        end
        ci_map[ci] = MOI.ConstraintIndex{F,S}(row)
        row += 1
    end
    return
end

function _get_moi_bound_type(lower, upper, type)
    if type == _INTERVAL
        return _INTERVAL
    elseif type == _EQUAL_TO
        return _EQUAL_TO
    elseif lower == upper
        return _LESS_AND_GREATER_THAN
    elseif lower <= -GLP_DBL_MAX
        return upper >= GLP_DBL_MAX ? _NONE : _LESS_THAN
    else
        return upper >= GLP_DBL_MAX ? _GREATER_THAN : _LESS_AND_GREATER_THAN
    end
end

function _add_all_variables(model::Optimizer, cache::_OptimizerCache)
    N = length(cache.cl)
    glp_add_cols(model, N)
    sizehint!(model.variable_info, N)
    @inbounds for i in 1:N
        bound = _get_moi_bound_type(cache.cl[i], cache.cu[i], cache.bounds[i])
        CleverDicts.add_item(
            model.variable_info,
            _VariableInfo(MOI.VariableIndex(i), i, bound, cache.types[i], ""),
        )
        glp_bound_type = _get_glp_bound_type(cache.cl[i], cache.cu[i])
        glp_set_col_bnds(model, i, glp_bound_type, cache.cl[i], cache.cu[i])
        if cache.types[i] == _BINARY
            model.num_binaries += 1
            glp_set_col_kind(model, i, GLP_IV)
        elseif cache.types[i] == _INTEGER
            model.num_integers += 1
            glp_set_col_kind(model, i, GLP_IV)
        end
    end
    return
end

function _add_all_constraints(dest::Optimizer, cache::_OptimizerCache)
    N = length(cache.rl)
    if N == 0
        return  # GLPK doesn't like it if we add 0 rows...
    end
    glp_add_rows(dest, N)
    glp_load_matrix(
        dest,
        length(cache.I),
        offset(cache.I),
        offset(cache.J),
        offset(cache.V),
    )
    sizehint!(dest.affine_constraint_info, N)
    @inbounds for i in 1:N
        l, u = cache.rl[i], cache.ru[i]
        if l == -Inf
            glp_set_row_bnds(dest, i, GLP_UP, -GLP_DBL_MAX, u)
            CleverDicts.add_item(
                dest.affine_constraint_info,
                _ConstraintInfo(i, MOI.LessThan{Float64}(u)),
            )
        elseif u == Inf
            glp_set_row_bnds(dest, i, GLP_LO, l, GLP_DBL_MAX)
            CleverDicts.add_item(
                dest.affine_constraint_info,
                _ConstraintInfo(i, MOI.GreaterThan{Float64}(l)),
            )
        else
            @assert l ≈ u  # Must be an equal-to constraint!
            glp_set_row_bnds(dest, i, GLP_FX, l, u)
            CleverDicts.add_item(
                dest.affine_constraint_info,
                _ConstraintInfo(i, MOI.EqualTo{Float64}(l)),
            )
        end
    end
    return
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    @assert MOI.is_empty(dest)
    _validate_constraint_types(dest, src)
    # Initialize the problem storage
    variables, map = _init_index_map(src)
    cache = _OptimizerCache(length(variables))
    # Extract the problem data
    #   Variable bounds:
    _extract_variable_data(src, map, cache, MOI.GreaterThan{Float64})
    _extract_variable_data(src, map, cache, MOI.LessThan{Float64})
    _extract_variable_data(src, map, cache, MOI.EqualTo{Float64})
    _extract_variable_data(src, map, cache, MOI.Interval{Float64})
    #   Variable types:
    _extract_variable_data(src, map, cache, MOI.Integer)
    _extract_variable_data(src, map, cache, MOI.ZeroOne)
    #   Affine constraints:
    _extract_row_data(src, map, cache, MOI.GreaterThan{Float64})
    _extract_row_data(src, map, cache, MOI.LessThan{Float64})
    _extract_row_data(src, map, cache, MOI.EqualTo{Float64})
    # Add the problem data
    _add_all_variables(dest, cache)
    _add_all_constraints(dest, cache)
    # Copy model attributes:
    MOI.Utilities.pass_attributes(dest, src, map)
    MOI.Utilities.pass_attributes(dest, src, map, variables)
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        indices = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        MOI.Utilities.pass_attributes(dest, src, map, indices)
    end
    return map
end
