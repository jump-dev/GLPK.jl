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
    get_infeasibility_ray(model::Optimizer, ray::Vector{Float64})

Get the Farkas certificate of primal infeasiblity.

Can only be called when GLPK.simplex is used as the solver.
"""
function get_infeasibility_ray(model::Optimizer, ray::Vector{Float64})
    num_rows = GLPK.get_num_rows(model.inner)
    @assert length(ray) == num_rows
    ur = GLPK.get_unbnd_ray(model.inner)
    if ur != 0
        if ur <= num_rows
            k = ur
            status      = GLPK.get_row_stat(model.inner, k)
            bind        = GLPK.get_row_bind(model.inner, k)
            primal      = GLPK.get_row_prim(model.inner, k)
            upper_bound = GLPK.get_row_ub(model.inner, k)
        else
            k = ur - num_rows
            status      = GLPK.get_col_stat(model.inner, k)
            bind        = GLPK.get_col_bind(model.inner, k)
            primal      = GLPK.get_col_prim(model.inner, k)
            upper_bound = GLPK.get_col_ub(model.inner, k)
        end
        if status != GLPK.BS
            error("unbounded ray is primal (use getunboundedray)")
        end
        ray[bind] = (primal > upper_bound) ? -1 : 1
        GLPK.btran(model.inner, ray)
    else
        eps = 1e-7
        # We need to factorize here because sometimes GLPK will prove
        # infeasibility before it has a factorized basis in memory.
        GLPK.factorize(model.inner)
        for row in 1:num_rows
            idx = GLPK.get_bhead(model.inner, row)
            if idx <= num_rows
                k = idx
                primal      = GLPK.get_row_prim(model.inner, k)
                upper_bound = GLPK.get_row_ub(model.inner, k)
                lower_bound = GLPK.get_row_lb(model.inner, k)
            else
                k = idx - num_rows
                primal      = GLPK.get_col_prim(model.inner, k)
                upper_bound = GLPK.get_col_ub(model.inner, k)
                lower_bound = GLPK.get_col_lb(model.inner, k)
            end
            if primal > upper_bound + eps
                ray[row] = -1
            elseif primal < lower_bound - eps
                ray[row] = 1
            else
                continue # ray[row] == 0
            end
            if idx <= num_rows
                ray[row] *= GLPK.get_rii(model.inner, k)
            else
                ray[row] /= GLPK.get_sjj(model.inner, k)
            end
        end
        GLPK.btran(model.inner, ray)
        for row in 1:num_rows
            ray[row] /= GLPK.get_rii(model.inner, row)
        end
    end
end

"""
    get_unbounded_ray(model::Optimizer, ray::Vector{Float64})

Get the certificate of primal unboundedness.

Can only be called when GLPK.simplex is used as the solver.
"""
function get_unbounded_ray(model::Optimizer, ray::Vector{Float64})
    num_rows = GLPK.get_num_rows(model.inner)
    n = GLPK.get_num_cols(model.inner)
    @assert length(ray) == n
    ur = GLPK.get_unbnd_ray(model.inner)
    if ur <= num_rows
        k = ur
        status = GLPK.get_row_stat(model.inner, k)
        dual = GLPK.get_row_dual(model.inner, k)
    else
        k = ur - num_rows
        status = GLPK.get_col_stat(model.inner, k)
        dual = GLPK.get_col_dual(model.inner, k)
        ray[k] = 1
    end
    if status == GLPK.BS
        error("unbounded ray is dual (use getinfeasibilityray)")
    end
    indices, values = GLPK.eval_tab_col(model.inner, ur)
    for (row, coef) in zip(indices, values)
        if row > num_rows
            ray[row - num_rows] = coef
        end
    end
    if xor(GLPK.get_obj_dir(model.inner) == GLPK.MAX, dual > 0)
        ray *= -1
    end
end
