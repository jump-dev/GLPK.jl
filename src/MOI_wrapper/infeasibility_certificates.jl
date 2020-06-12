# The functions getinfeasibilityray and getunboundedray are adapted from code
# taken from the LEMON C++ optimization library. This is the copyright notice:
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

Get the Farkas certificate of primal infeasiblity.

Can only be called when glp_simplex is used as the solver.
"""
function _get_infeasibility_ray(model::Optimizer, ray::Vector{Float64})
    num_rows = glp_get_num_rows(model.inner)
    @assert length(ray) == num_rows
    ur = glp_get_unbnd_ray(model.inner)
    if ur != 0
        if ur <= num_rows
            k = ur
            status      = glp_get_row_stat(model.inner, k)
            bind        = glp_get_row_bind(model.inner, k)
            primal      = glp_get_row_prim(model.inner, k)
            upper_bound = glp_get_row_ub(model.inner, k)
        else
            k = ur - num_rows
            status      = glp_get_col_stat(model.inner, k)
            bind        = glp_get_col_bind(model.inner, k)
            primal      = glp_get_col_prim(model.inner, k)
            upper_bound = glp_get_col_ub(model.inner, k)
        end
        if status != GLP_BS
            error("unbounded ray is primal (use getunboundedray)")
        end
        ray[bind] = (primal > upper_bound) ? -1 : 1
        glp_btran(model.inner, ray)
    else
        eps = 1e-7
        # We need to factorize here because sometimes GLPK will prove
        # infeasibility before it has a factorized basis in memory.
        glp_factorize(model.inner)
        for row in 1:num_rows
            idx = glp_get_bhead(model.inner, row)
            if idx <= num_rows
                k = idx
                primal      = glp_get_row_prim(model.inner, k)
                upper_bound = glp_get_row_ub(model.inner, k)
                lower_bound = glp_get_row_lb(model.inner, k)
            else
                k = idx - num_rows
                primal      = glp_get_col_prim(model.inner, k)
                upper_bound = glp_get_col_ub(model.inner, k)
                lower_bound = glp_get_col_lb(model.inner, k)
            end
            if primal > upper_bound + eps
                ray[row] = -1
            elseif primal < lower_bound - eps
                ray[row] = 1
            else
                continue # ray[row] == 0
            end
            if idx <= num_rows
                ray[row] *= glp_get_rii(model.inner, k)
            else
                ray[row] /= glp_get_sjj(model.inner, k)
            end
        end
        glp_btran(model.inner, ray)
        for row in 1:num_rows
            ray[row] /= glp_get_rii(model.inner, row)
        end
    end
end

"""
    _get_unbounded_ray(model::Optimizer, ray::Vector{Float64})

Get the certificate of primal unboundedness.

Can only be called when glp_simplex is used as the solver.
"""
function _get_unbounded_ray(model::Optimizer, ray::Vector{Float64})
    num_rows = glp_get_num_rows(model.inner)
    n = glp_get_num_cols(model.inner)
    @assert length(ray) == n
    ur = glp_get_unbnd_ray(model.inner)
    if ur <= num_rows
        k = ur
        status = glp_get_row_stat(model.inner, k)
        dual = glp_get_row_dual(model.inner, k)
    else
        k = ur - num_rows
        status = glp_get_col_stat(model.inner, k)
        dual = glp_get_col_dual(model.inner, k)
        ray[k] = 1
    end
    if status == GLP_BS
        error("unbounded ray is dual (use _get_infeasibility_ray)")
    end
    nnz = n + num_rows
    indices, coefficients = zeros(Cint, nnz), zeros(Cdouble, nnz)
    len = glp_eval_tab_col(
        model.inner,
        nnz,
        pointer(indices) - sizeof(Cint),
        pointer(coefficients) - sizeof(Cdouble),
    )
    splice!(indices, (len+1):nnz)
    splice!(coefficients, (len+1):nnz)
    for (row, coef) in zip(indices, coefficients)
        if row > num_rows
            ray[row - num_rows] = coef
        end
    end
    if xor(glp_get_obj_dir(model.inner) == GLP_MAX, dual > 0)
        ray *= -1
    end
end
