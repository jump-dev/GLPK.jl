# The functions _get_infeasibility_ray and _get_unbounded_ray are adapted from
# code taken from the LEMON C++ optimization library. This is the copyright
# notice:
#
### Copyright (C) 2003-2010
### Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
### (Egervary Research Group on Combinatorial Optimization, EGRES).
###
### Permission to use, modify and distribute this software is granted
### provided that this copyright notice appears in all copies. For
### precise terms see the accompanying LICENSE file.
###
### This software is provided "AS IS" with no warranty of any kind,
### express or implied, and with no claim as to its suitability for any
### purpose.

"""
    _get_infeasibility_ray(model::Optimizer, ray::Vector{Float64})

Compute a Farkas certificate of primal infeasibility (unbounded dual ray) and
store in `ray`. Returns `true` if successful, and `false` if a ray cannot be
found.

Assumes `ray` has been initialized to all `0.0`s.
"""
function _get_infeasibility_ray(model::Optimizer, ray::Vector{Float64})
    m = glp_get_num_rows(model)
    n = glp_get_num_cols(model)
    @assert length(ray) == m
    # Solve with dual simplex to find unbounded ray.
    param = glp_smcp()
    glp_init_smcp(param)
    param.msg_lev = GLP_MSG_ERR
    param.meth = GLP_DUAL
    status = glp_simplex(model, param)
    if status != 0 || glp_get_status(model) != GLP_NOFEAS
        return false  # Something went wrong finding an unbounded ray.
    end
    unbounded_index = glp_get_unbnd_ray(model)
    if unbounded_index == 0
        return false  # Something went wrong finding an unbounded ray.
    end
    # If the primal value exceeds the upper bound, then the unbounded_index
    # wants to increase. Otherwise, it must want to decrease.
    scale = if unbounded_index <= m
        primal = glp_get_row_prim(model, unbounded_index)
        primal > glp_get_row_ub(model, unbounded_index) ? 1 : -1
    else
        primal = glp_get_col_prim(model, unbounded_index - m)
        primal > glp_get_col_ub(model, unbounded_index - m) ? 1 : -1
    end
    if unbounded_index <= m
        ray[unbounded_index] = scale
    end
    nnz = m + n
    vind = Vector{Cint}(undef, nnz)
    vval = Vector{Cdouble}(undef, nnz)
    len = glp_eval_tab_row(model, unbounded_index, offset(vind), offset(vval))
    for i = 1:len
        if vind[i] <= m
            ray[vind[i]] = -scale * vval[i]
        end
    end
    return true
end

"""
    _get_unbounded_ray(model::Optimizer, ray::Vector{Float64})::Bool

Compute an unbounded primal ray and store in `ray`. Returns `true` if
successful, and `false` if a ray cannot be found.

Assumes the primal has been solved with primal simplex and is proven unbounded.
Assumes `ray` has been initialized to all `0.0`s.
"""
function _get_unbounded_ray(model::Optimizer, ray::Vector{Float64})
    m = glp_get_num_rows(model)
    n = glp_get_num_cols(model)
    @assert length(ray) == n
    unbounded_index = glp_get_unbnd_ray(model)
    if unbounded_index == 0
        return false  # Something went wrong finding an unbounded ray.
    end
    dual = if unbounded_index <= m
        glp_get_row_dual(model, unbounded_index)
    else
        glp_get_col_dual(model, unbounded_index - m)
    end
    scale = xor(glp_get_obj_dir(model) == GLP_MAX, dual > 0) ? -1 : 1
    if unbounded_index > m
        ray[unbounded_index-m] = scale
    end
    nnz = m + n
    vind = Vector{Cint}(undef, nnz)
    vval = Vector{Cdouble}(undef, nnz)
    len = glp_eval_tab_col(model, unbounded_index, offset(vind), offset(vval))
    for i = 1:len
        if vind[i] > m
            ray[vind[i]-m] = scale * vval[i]
        end
    end
    return true
end
