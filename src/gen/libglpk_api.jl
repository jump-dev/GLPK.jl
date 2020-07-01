# Julia wrapper for header: glpk.h
# Automatically generated using Clang.jl


function glp_create_prob()
    ccall((:glp_create_prob, libglpk), Ptr{glp_prob}, ())
end

function glp_set_prob_name(P, name)
    ccall((:glp_set_prob_name, libglpk), Cvoid, (Ptr{glp_prob}, Cstring), P, name)
end

function glp_set_obj_name(P, name)
    ccall((:glp_set_obj_name, libglpk), Cvoid, (Ptr{glp_prob}, Cstring), P, name)
end

function glp_set_obj_dir(P, dir)
    ccall((:glp_set_obj_dir, libglpk), Cvoid, (Ptr{glp_prob}, Cint), P, dir)
end

function glp_add_rows(P, nrs)
    ccall((:glp_add_rows, libglpk), Cint, (Ptr{glp_prob}, Cint), P, nrs)
end

function glp_add_cols(P, ncs)
    ccall((:glp_add_cols, libglpk), Cint, (Ptr{glp_prob}, Cint), P, ncs)
end

function glp_set_row_name(P, i, name)
    ccall((:glp_set_row_name, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cstring), P, i, name)
end

function glp_set_col_name(P, j, name)
    ccall((:glp_set_col_name, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cstring), P, j, name)
end

function glp_set_row_bnds(P, i, type, lb, ub)
    ccall((:glp_set_row_bnds, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cint, Cdouble, Cdouble), P, i, type, lb, ub)
end

function glp_set_col_bnds(P, j, type, lb, ub)
    ccall((:glp_set_col_bnds, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cint, Cdouble, Cdouble), P, j, type, lb, ub)
end

function glp_set_obj_coef(P, j, coef)
    ccall((:glp_set_obj_coef, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cdouble), P, j, coef)
end

function glp_set_mat_row(P, i, len, ind, val)
    ccall((:glp_set_mat_row, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}), P, i, len, ind, val)
end

function glp_set_mat_col(P, j, len, ind, val)
    ccall((:glp_set_mat_col, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}), P, j, len, ind, val)
end

function glp_load_matrix(P, ne, ia, ja, ar)
    ccall((:glp_load_matrix, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), P, ne, ia, ja, ar)
end

function glp_check_dup(m, n, ne, ia, ja)
    ccall((:glp_check_dup, libglpk), Cint, (Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}), m, n, ne, ia, ja)
end

function glp_sort_matrix(P)
    ccall((:glp_sort_matrix, libglpk), Cvoid, (Ptr{glp_prob},), P)
end

function glp_del_rows(P, nrs, num)
    ccall((:glp_del_rows, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Ptr{Cint}), P, nrs, num)
end

function glp_del_cols(P, ncs, num)
    ccall((:glp_del_cols, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Ptr{Cint}), P, ncs, num)
end

function glp_copy_prob(dest, prob, names)
    ccall((:glp_copy_prob, libglpk), Cvoid, (Ptr{glp_prob}, Ptr{glp_prob}, Cint), dest, prob, names)
end

function glp_erase_prob(P)
    ccall((:glp_erase_prob, libglpk), Cvoid, (Ptr{glp_prob},), P)
end

function glp_delete_prob(P)
    ccall((:glp_delete_prob, libglpk), Cvoid, (Ptr{glp_prob},), P)
end

function glp_get_prob_name(P)
    ccall((:glp_get_prob_name, libglpk), Cstring, (Ptr{glp_prob},), P)
end

function glp_get_obj_name(P)
    ccall((:glp_get_obj_name, libglpk), Cstring, (Ptr{glp_prob},), P)
end

function glp_get_obj_dir(P)
    ccall((:glp_get_obj_dir, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_num_rows(P)
    ccall((:glp_get_num_rows, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_num_cols(P)
    ccall((:glp_get_num_cols, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_row_name(P, i)
    ccall((:glp_get_row_name, libglpk), Cstring, (Ptr{glp_prob}, Cint), P, i)
end

function glp_get_col_name(P, j)
    ccall((:glp_get_col_name, libglpk), Cstring, (Ptr{glp_prob}, Cint), P, j)
end

function glp_get_row_type(P, i)
    ccall((:glp_get_row_type, libglpk), Cint, (Ptr{glp_prob}, Cint), P, i)
end

function glp_get_row_lb(P, i)
    ccall((:glp_get_row_lb, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, i)
end

function glp_get_row_ub(P, i)
    ccall((:glp_get_row_ub, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, i)
end

function glp_get_col_type(P, j)
    ccall((:glp_get_col_type, libglpk), Cint, (Ptr{glp_prob}, Cint), P, j)
end

function glp_get_col_lb(P, j)
    ccall((:glp_get_col_lb, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, j)
end

function glp_get_col_ub(P, j)
    ccall((:glp_get_col_ub, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, j)
end

function glp_get_obj_coef(P, j)
    ccall((:glp_get_obj_coef, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, j)
end

function glp_get_num_nz(P)
    ccall((:glp_get_num_nz, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_mat_row(P, i, ind, val)
    ccall((:glp_get_mat_row, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{Cint}, Ptr{Cdouble}), P, i, ind, val)
end

function glp_get_mat_col(P, j, ind, val)
    ccall((:glp_get_mat_col, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{Cint}, Ptr{Cdouble}), P, j, ind, val)
end

function glp_create_index(P)
    ccall((:glp_create_index, libglpk), Cvoid, (Ptr{glp_prob},), P)
end

function glp_find_row(P, name)
    ccall((:glp_find_row, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, name)
end

function glp_find_col(P, name)
    ccall((:glp_find_col, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, name)
end

function glp_delete_index(P)
    ccall((:glp_delete_index, libglpk), Cvoid, (Ptr{glp_prob},), P)
end

function glp_set_rii(P, i, rii)
    ccall((:glp_set_rii, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cdouble), P, i, rii)
end

function glp_set_sjj(P, j, sjj)
    ccall((:glp_set_sjj, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cdouble), P, j, sjj)
end

function glp_get_rii(P, i)
    ccall((:glp_get_rii, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, i)
end

function glp_get_sjj(P, j)
    ccall((:glp_get_sjj, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, j)
end

function glp_scale_prob(P, flags)
    ccall((:glp_scale_prob, libglpk), Cvoid, (Ptr{glp_prob}, Cint), P, flags)
end

function glp_unscale_prob(P)
    ccall((:glp_unscale_prob, libglpk), Cvoid, (Ptr{glp_prob},), P)
end

function glp_set_row_stat(P, i, stat)
    ccall((:glp_set_row_stat, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cint), P, i, stat)
end

function glp_set_col_stat(P, j, stat)
    ccall((:glp_set_col_stat, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cint), P, j, stat)
end

function glp_std_basis(P)
    ccall((:glp_std_basis, libglpk), Cvoid, (Ptr{glp_prob},), P)
end

function glp_adv_basis(P, flags)
    ccall((:glp_adv_basis, libglpk), Cvoid, (Ptr{glp_prob}, Cint), P, flags)
end

function glp_cpx_basis(P)
    ccall((:glp_cpx_basis, libglpk), Cvoid, (Ptr{glp_prob},), P)
end

function glp_simplex(P, parm)
    ccall((:glp_simplex, libglpk), Cint, (Ptr{glp_prob}, Ptr{glp_smcp}), P, parm)
end

function glp_exact(P, parm)
    ccall((:glp_exact, libglpk), Cint, (Ptr{glp_prob}, Ptr{glp_smcp}), P, parm)
end

function glp_init_smcp(parm)
    ccall((:glp_init_smcp, libglpk), Cvoid, (Ptr{glp_smcp},), parm)
end

function glp_get_status(P)
    ccall((:glp_get_status, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_prim_stat(P)
    ccall((:glp_get_prim_stat, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_dual_stat(P)
    ccall((:glp_get_dual_stat, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_obj_val(P)
    ccall((:glp_get_obj_val, libglpk), Cdouble, (Ptr{glp_prob},), P)
end

function glp_get_row_stat(P, i)
    ccall((:glp_get_row_stat, libglpk), Cint, (Ptr{glp_prob}, Cint), P, i)
end

function glp_get_row_prim(P, i)
    ccall((:glp_get_row_prim, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, i)
end

function glp_get_row_dual(P, i)
    ccall((:glp_get_row_dual, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, i)
end

function glp_get_col_stat(P, j)
    ccall((:glp_get_col_stat, libglpk), Cint, (Ptr{glp_prob}, Cint), P, j)
end

function glp_get_col_prim(P, j)
    ccall((:glp_get_col_prim, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, j)
end

function glp_get_col_dual(P, j)
    ccall((:glp_get_col_dual, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, j)
end

function glp_get_unbnd_ray(P)
    ccall((:glp_get_unbnd_ray, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_it_cnt(P)
    ccall((:glp_get_it_cnt, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_set_it_cnt(P, it_cnt)
    ccall((:glp_set_it_cnt, libglpk), Cvoid, (Ptr{glp_prob}, Cint), P, it_cnt)
end

function glp_interior(P, parm)
    ccall((:glp_interior, libglpk), Cint, (Ptr{glp_prob}, Ptr{glp_iptcp}), P, parm)
end

function glp_init_iptcp(parm)
    ccall((:glp_init_iptcp, libglpk), Cvoid, (Ptr{glp_iptcp},), parm)
end

function glp_ipt_status(P)
    ccall((:glp_ipt_status, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_ipt_obj_val(P)
    ccall((:glp_ipt_obj_val, libglpk), Cdouble, (Ptr{glp_prob},), P)
end

function glp_ipt_row_prim(P, i)
    ccall((:glp_ipt_row_prim, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, i)
end

function glp_ipt_row_dual(P, i)
    ccall((:glp_ipt_row_dual, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, i)
end

function glp_ipt_col_prim(P, j)
    ccall((:glp_ipt_col_prim, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, j)
end

function glp_ipt_col_dual(P, j)
    ccall((:glp_ipt_col_dual, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, j)
end

function glp_set_col_kind(P, j, kind)
    ccall((:glp_set_col_kind, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cint), P, j, kind)
end

function glp_get_col_kind(P, j)
    ccall((:glp_get_col_kind, libglpk), Cint, (Ptr{glp_prob}, Cint), P, j)
end

function glp_get_num_int(P)
    ccall((:glp_get_num_int, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_num_bin(P)
    ccall((:glp_get_num_bin, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_intopt(P, parm)
    ccall((:glp_intopt, libglpk), Cint, (Ptr{glp_prob}, Ptr{glp_iocp}), P, parm)
end

function glp_init_iocp(parm)
    ccall((:glp_init_iocp, libglpk), Cvoid, (Ptr{glp_iocp},), parm)
end

function glp_mip_status(P)
    ccall((:glp_mip_status, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_mip_obj_val(P)
    ccall((:glp_mip_obj_val, libglpk), Cdouble, (Ptr{glp_prob},), P)
end

function glp_mip_row_val(P, i)
    ccall((:glp_mip_row_val, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, i)
end

function glp_mip_col_val(P, j)
    ccall((:glp_mip_col_val, libglpk), Cdouble, (Ptr{glp_prob}, Cint), P, j)
end

function glp_check_kkt(P, sol, cond, ae_max, ae_ind, re_max, re_ind)
    ccall((:glp_check_kkt, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}), P, sol, cond, ae_max, ae_ind, re_max, re_ind)
end

function glp_print_sol(P, fname)
    ccall((:glp_print_sol, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_read_sol(P, fname)
    ccall((:glp_read_sol, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_write_sol(P, fname)
    ccall((:glp_write_sol, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_print_ranges(P, len, list, flags, fname)
    ccall((:glp_print_ranges, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{Cint}, Cint, Cstring), P, len, list, flags, fname)
end

function glp_print_ipt(P, fname)
    ccall((:glp_print_ipt, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_read_ipt(P, fname)
    ccall((:glp_read_ipt, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_write_ipt(P, fname)
    ccall((:glp_write_ipt, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_print_mip(P, fname)
    ccall((:glp_print_mip, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_read_mip(P, fname)
    ccall((:glp_read_mip, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_write_mip(P, fname)
    ccall((:glp_write_mip, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_bf_exists(P)
    ccall((:glp_bf_exists, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_factorize(P)
    ccall((:glp_factorize, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_bf_updated(P)
    ccall((:glp_bf_updated, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_get_bfcp(P, parm)
    ccall((:glp_get_bfcp, libglpk), Cvoid, (Ptr{glp_prob}, Ptr{glp_bfcp}), P, parm)
end

function glp_set_bfcp(P, parm)
    ccall((:glp_set_bfcp, libglpk), Cvoid, (Ptr{glp_prob}, Ptr{glp_bfcp}), P, parm)
end

function glp_get_bhead(P, k)
    ccall((:glp_get_bhead, libglpk), Cint, (Ptr{glp_prob}, Cint), P, k)
end

function glp_get_row_bind(P, i)
    ccall((:glp_get_row_bind, libglpk), Cint, (Ptr{glp_prob}, Cint), P, i)
end

function glp_get_col_bind(P, j)
    ccall((:glp_get_col_bind, libglpk), Cint, (Ptr{glp_prob}, Cint), P, j)
end

function glp_ftran(P, x)
    ccall((:glp_ftran, libglpk), Cvoid, (Ptr{glp_prob}, Ptr{Cdouble}), P, x)
end

function glp_btran(P, x)
    ccall((:glp_btran, libglpk), Cvoid, (Ptr{glp_prob}, Ptr{Cdouble}), P, x)
end

function glp_warm_up(P)
    ccall((:glp_warm_up, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_eval_tab_row(P, k, ind, val)
    ccall((:glp_eval_tab_row, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{Cint}, Ptr{Cdouble}), P, k, ind, val)
end

function glp_eval_tab_col(P, k, ind, val)
    ccall((:glp_eval_tab_col, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{Cint}, Ptr{Cdouble}), P, k, ind, val)
end

function glp_transform_row(P, len, ind, val)
    ccall((:glp_transform_row, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{Cint}, Ptr{Cdouble}), P, len, ind, val)
end

function glp_transform_col(P, len, ind, val)
    ccall((:glp_transform_col, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{Cint}, Ptr{Cdouble}), P, len, ind, val)
end

function glp_prim_rtest(P, len, ind, val, dir, eps)
    ccall((:glp_prim_rtest, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint, Cdouble), P, len, ind, val, dir, eps)
end

function glp_dual_rtest(P, len, ind, val, dir, eps)
    ccall((:glp_dual_rtest, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint, Cdouble), P, len, ind, val, dir, eps)
end

function glp_analyze_bound(P, k, value1, var1, value2, var2)
    ccall((:glp_analyze_bound, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}), P, k, value1, var1, value2, var2)
end

function glp_analyze_coef(P, k, coef1, var1, value1, coef2, var2, value2)
    ccall((:glp_analyze_coef, libglpk), Cvoid, (Ptr{glp_prob}, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}), P, k, coef1, var1, value1, coef2, var2, value2)
end

function glp_ios_reason(T)
    ccall((:glp_ios_reason, libglpk), Cint, (Ptr{glp_tree},), T)
end

function glp_ios_get_prob(T)
    ccall((:glp_ios_get_prob, libglpk), Ptr{glp_prob}, (Ptr{glp_tree},), T)
end

function glp_ios_tree_size(T, a_cnt, n_cnt, t_cnt)
    ccall((:glp_ios_tree_size, libglpk), Cvoid, (Ptr{glp_tree}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), T, a_cnt, n_cnt, t_cnt)
end

function glp_ios_curr_node(T)
    ccall((:glp_ios_curr_node, libglpk), Cint, (Ptr{glp_tree},), T)
end

function glp_ios_next_node(T, p)
    ccall((:glp_ios_next_node, libglpk), Cint, (Ptr{glp_tree}, Cint), T, p)
end

function glp_ios_prev_node(T, p)
    ccall((:glp_ios_prev_node, libglpk), Cint, (Ptr{glp_tree}, Cint), T, p)
end

function glp_ios_up_node(T, p)
    ccall((:glp_ios_up_node, libglpk), Cint, (Ptr{glp_tree}, Cint), T, p)
end

function glp_ios_node_level(T, p)
    ccall((:glp_ios_node_level, libglpk), Cint, (Ptr{glp_tree}, Cint), T, p)
end

function glp_ios_node_bound(T, p)
    ccall((:glp_ios_node_bound, libglpk), Cdouble, (Ptr{glp_tree}, Cint), T, p)
end

function glp_ios_best_node(T)
    ccall((:glp_ios_best_node, libglpk), Cint, (Ptr{glp_tree},), T)
end

function glp_ios_mip_gap(T)
    ccall((:glp_ios_mip_gap, libglpk), Cdouble, (Ptr{glp_tree},), T)
end

function glp_ios_node_data(T, p)
    ccall((:glp_ios_node_data, libglpk), Ptr{Cvoid}, (Ptr{glp_tree}, Cint), T, p)
end

function glp_ios_row_attr(T, i, attr)
    ccall((:glp_ios_row_attr, libglpk), Cvoid, (Ptr{glp_tree}, Cint, Ptr{glp_attr}), T, i, attr)
end

function glp_ios_pool_size(T)
    ccall((:glp_ios_pool_size, libglpk), Cint, (Ptr{glp_tree},), T)
end

function glp_ios_add_row(T, name, klass, flags, len, ind, val, type, rhs)
    ccall((:glp_ios_add_row, libglpk), Cint, (Ptr{glp_tree}, Cstring, Cint, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint, Cdouble), T, name, klass, flags, len, ind, val, type, rhs)
end

function glp_ios_del_row(T, i)
    ccall((:glp_ios_del_row, libglpk), Cvoid, (Ptr{glp_tree}, Cint), T, i)
end

function glp_ios_clear_pool(T)
    ccall((:glp_ios_clear_pool, libglpk), Cvoid, (Ptr{glp_tree},), T)
end

function glp_ios_can_branch(T, j)
    ccall((:glp_ios_can_branch, libglpk), Cint, (Ptr{glp_tree}, Cint), T, j)
end

function glp_ios_branch_upon(T, j, sel)
    ccall((:glp_ios_branch_upon, libglpk), Cvoid, (Ptr{glp_tree}, Cint, Cint), T, j, sel)
end

function glp_ios_select_node(T, p)
    ccall((:glp_ios_select_node, libglpk), Cvoid, (Ptr{glp_tree}, Cint), T, p)
end

function glp_ios_heur_sol(T, x)
    ccall((:glp_ios_heur_sol, libglpk), Cint, (Ptr{glp_tree}, Ptr{Cdouble}), T, x)
end

function glp_ios_terminate(T)
    ccall((:glp_ios_terminate, libglpk), Cvoid, (Ptr{glp_tree},), T)
end

function glp_init_mpscp(parm)
    ccall((:glp_init_mpscp, libglpk), Cvoid, (Ptr{glp_mpscp},), parm)
end

function glp_read_mps(P, fmt, parm, fname)
    ccall((:glp_read_mps, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{glp_mpscp}, Cstring), P, fmt, parm, fname)
end

function glp_write_mps(P, fmt, parm, fname)
    ccall((:glp_write_mps, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{glp_mpscp}, Cstring), P, fmt, parm, fname)
end

function glp_init_cpxcp(parm)
    ccall((:glp_init_cpxcp, libglpk), Cvoid, (Ptr{glp_cpxcp},), parm)
end

function glp_read_lp(P, parm, fname)
    ccall((:glp_read_lp, libglpk), Cint, (Ptr{glp_prob}, Ptr{glp_cpxcp}, Cstring), P, parm, fname)
end

function glp_write_lp(P, parm, fname)
    ccall((:glp_write_lp, libglpk), Cint, (Ptr{glp_prob}, Ptr{glp_cpxcp}, Cstring), P, parm, fname)
end

function glp_read_prob(P, flags, fname)
    ccall((:glp_read_prob, libglpk), Cint, (Ptr{glp_prob}, Cint, Cstring), P, flags, fname)
end

function glp_write_prob(P, flags, fname)
    ccall((:glp_write_prob, libglpk), Cint, (Ptr{glp_prob}, Cint, Cstring), P, flags, fname)
end

function glp_mpl_alloc_wksp()
    ccall((:glp_mpl_alloc_wksp, libglpk), Ptr{glp_tran}, ())
end

function glp_mpl_init_rand(tran, seed)
    ccall((:glp_mpl_init_rand, libglpk), Cvoid, (Ptr{glp_tran}, Cint), tran, seed)
end

function glp_mpl_read_model(tran, fname, skip)
    ccall((:glp_mpl_read_model, libglpk), Cint, (Ptr{glp_tran}, Cstring, Cint), tran, fname, skip)
end

function glp_mpl_read_data(tran, fname)
    ccall((:glp_mpl_read_data, libglpk), Cint, (Ptr{glp_tran}, Cstring), tran, fname)
end

function glp_mpl_generate(tran, fname)
    ccall((:glp_mpl_generate, libglpk), Cint, (Ptr{glp_tran}, Cstring), tran, fname)
end

function glp_mpl_build_prob(tran, prob)
    ccall((:glp_mpl_build_prob, libglpk), Cvoid, (Ptr{glp_tran}, Ptr{glp_prob}), tran, prob)
end

function glp_mpl_postsolve(tran, prob, sol)
    ccall((:glp_mpl_postsolve, libglpk), Cint, (Ptr{glp_tran}, Ptr{glp_prob}, Cint), tran, prob, sol)
end

function glp_mpl_free_wksp(tran)
    ccall((:glp_mpl_free_wksp, libglpk), Cvoid, (Ptr{glp_tran},), tran)
end

function glp_read_cnfsat(P, fname)
    ccall((:glp_read_cnfsat, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_check_cnfsat(P)
    ccall((:glp_check_cnfsat, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_write_cnfsat(P, fname)
    ccall((:glp_write_cnfsat, libglpk), Cint, (Ptr{glp_prob}, Cstring), P, fname)
end

function glp_minisat1(P)
    ccall((:glp_minisat1, libglpk), Cint, (Ptr{glp_prob},), P)
end

function glp_intfeas1(P, use_bound, obj_bound)
    ccall((:glp_intfeas1, libglpk), Cint, (Ptr{glp_prob}, Cint, Cint), P, use_bound, obj_bound)
end

function glp_init_env()
    ccall((:glp_init_env, libglpk), Cint, ())
end

function glp_version()
    ccall((:glp_version, libglpk), Cstring, ())
end

function glp_config(option)
    ccall((:glp_config, libglpk), Cstring, (Cstring,), option)
end

function glp_free_env()
    ccall((:glp_free_env, libglpk), Cint, ())
end

function glp_puts(s)
    ccall((:glp_puts, libglpk), Cvoid, (Cstring,), s)
end

function glp_term_out(flag)
    ccall((:glp_term_out, libglpk), Cint, (Cint,), flag)
end

function glp_term_hook(func, info)
    ccall((:glp_term_hook, libglpk), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), func, info)
end

function glp_open_tee(name)
    ccall((:glp_open_tee, libglpk), Cint, (Cstring,), name)
end

function glp_close_tee()
    ccall((:glp_close_tee, libglpk), Cint, ())
end

function glp_error_(file, line)
    ccall((:glp_error_, libglpk), glp_errfunc, (Cstring, Cint), file, line)
end

function glp_at_error()
    ccall((:glp_at_error, libglpk), Cint, ())
end

function glp_assert_(expr, file, line)
    ccall((:glp_assert_, libglpk), Cvoid, (Cstring, Cstring, Cint), expr, file, line)
end

function glp_error_hook(func, info)
    ccall((:glp_error_hook, libglpk), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), func, info)
end

function glp_alloc(n, size)
    ccall((:glp_alloc, libglpk), Ptr{Cvoid}, (Cint, Cint), n, size)
end

function glp_realloc(ptr, n, size)
    ccall((:glp_realloc, libglpk), Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Cint), ptr, n, size)
end

function glp_free(ptr)
    ccall((:glp_free, libglpk), Cvoid, (Ptr{Cvoid},), ptr)
end

function glp_mem_limit(limit)
    ccall((:glp_mem_limit, libglpk), Cvoid, (Cint,), limit)
end

function glp_mem_usage(count, cpeak, total, tpeak)
    ccall((:glp_mem_usage, libglpk), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Csize_t}, Ptr{Csize_t}), count, cpeak, total, tpeak)
end

function glp_time()
    ccall((:glp_time, libglpk), Cdouble, ())
end

function glp_difftime(t1, t0)
    ccall((:glp_difftime, libglpk), Cdouble, (Cdouble, Cdouble), t1, t0)
end

function glp_create_graph(v_size, a_size)
    ccall((:glp_create_graph, libglpk), Ptr{glp_graph}, (Cint, Cint), v_size, a_size)
end

function glp_set_graph_name(G, name)
    ccall((:glp_set_graph_name, libglpk), Cvoid, (Ptr{glp_graph}, Cstring), G, name)
end

function glp_add_vertices(G, nadd)
    ccall((:glp_add_vertices, libglpk), Cint, (Ptr{glp_graph}, Cint), G, nadd)
end

function glp_set_vertex_name(G, i, name)
    ccall((:glp_set_vertex_name, libglpk), Cvoid, (Ptr{glp_graph}, Cint, Cstring), G, i, name)
end

function glp_add_arc(G, i, j)
    ccall((:glp_add_arc, libglpk), Ptr{glp_arc}, (Ptr{glp_graph}, Cint, Cint), G, i, j)
end

function glp_del_vertices(G, ndel, num)
    ccall((:glp_del_vertices, libglpk), Cvoid, (Ptr{glp_graph}, Cint, Ptr{Cint}), G, ndel, num)
end

function glp_del_arc(G, a)
    ccall((:glp_del_arc, libglpk), Cvoid, (Ptr{glp_graph}, Ptr{glp_arc}), G, a)
end

function glp_erase_graph(G, v_size, a_size)
    ccall((:glp_erase_graph, libglpk), Cvoid, (Ptr{glp_graph}, Cint, Cint), G, v_size, a_size)
end

function glp_delete_graph(G)
    ccall((:glp_delete_graph, libglpk), Cvoid, (Ptr{glp_graph},), G)
end

function glp_create_v_index(G)
    ccall((:glp_create_v_index, libglpk), Cvoid, (Ptr{glp_graph},), G)
end

function glp_find_vertex(G, name)
    ccall((:glp_find_vertex, libglpk), Cint, (Ptr{glp_graph}, Cstring), G, name)
end

function glp_delete_v_index(G)
    ccall((:glp_delete_v_index, libglpk), Cvoid, (Ptr{glp_graph},), G)
end

function glp_read_graph(G, fname)
    ccall((:glp_read_graph, libglpk), Cint, (Ptr{glp_graph}, Cstring), G, fname)
end

function glp_write_graph(G, fname)
    ccall((:glp_write_graph, libglpk), Cint, (Ptr{glp_graph}, Cstring), G, fname)
end

function glp_mincost_lp(P, G, names, v_rhs, a_low, a_cap, a_cost)
    ccall((:glp_mincost_lp, libglpk), Cvoid, (Ptr{glp_prob}, Ptr{glp_graph}, Cint, Cint, Cint, Cint, Cint), P, G, names, v_rhs, a_low, a_cap, a_cost)
end

function glp_mincost_okalg(G, v_rhs, a_low, a_cap, a_cost, sol, a_x, v_pi)
    ccall((:glp_mincost_okalg, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cint, Cint, Ptr{Cdouble}, Cint, Cint), G, v_rhs, a_low, a_cap, a_cost, sol, a_x, v_pi)
end

function glp_mincost_relax4(G, v_rhs, a_low, a_cap, a_cost, crash, sol, a_x, a_rc)
    ccall((:glp_mincost_relax4, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cint, Cint, Cint, Ptr{Cdouble}, Cint, Cint), G, v_rhs, a_low, a_cap, a_cost, crash, sol, a_x, a_rc)
end

function glp_maxflow_lp(P, G, names, s, t, a_cap)
    ccall((:glp_maxflow_lp, libglpk), Cvoid, (Ptr{glp_prob}, Ptr{glp_graph}, Cint, Cint, Cint, Cint), P, G, names, s, t, a_cap)
end

function glp_maxflow_ffalg(G, s, t, a_cap, sol, a_x, v_cut)
    ccall((:glp_maxflow_ffalg, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cint, Ptr{Cdouble}, Cint, Cint), G, s, t, a_cap, sol, a_x, v_cut)
end

function glp_check_asnprob(G, v_set)
    ccall((:glp_check_asnprob, libglpk), Cint, (Ptr{glp_graph}, Cint), G, v_set)
end

function glp_asnprob_lp(P, form, G, names, v_set, a_cost)
    ccall((:glp_asnprob_lp, libglpk), Cint, (Ptr{glp_prob}, Cint, Ptr{glp_graph}, Cint, Cint, Cint), P, form, G, names, v_set, a_cost)
end

function glp_asnprob_okalg(form, G, v_set, a_cost, sol, a_x)
    ccall((:glp_asnprob_okalg, libglpk), Cint, (Cint, Ptr{glp_graph}, Cint, Cint, Ptr{Cdouble}, Cint), form, G, v_set, a_cost, sol, a_x)
end

function glp_asnprob_hall(G, v_set, a_x)
    ccall((:glp_asnprob_hall, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint), G, v_set, a_x)
end

function glp_cpp(G, v_t, v_es, v_ls)
    ccall((:glp_cpp, libglpk), Cdouble, (Ptr{glp_graph}, Cint, Cint, Cint), G, v_t, v_es, v_ls)
end

function glp_read_mincost(G, v_rhs, a_low, a_cap, a_cost, fname)
    ccall((:glp_read_mincost, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cint, Cint, Cstring), G, v_rhs, a_low, a_cap, a_cost, fname)
end

function glp_write_mincost(G, v_rhs, a_low, a_cap, a_cost, fname)
    ccall((:glp_write_mincost, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cint, Cint, Cstring), G, v_rhs, a_low, a_cap, a_cost, fname)
end

function glp_read_maxflow(G, s, t, a_cap, fname)
    ccall((:glp_read_maxflow, libglpk), Cint, (Ptr{glp_graph}, Ptr{Cint}, Ptr{Cint}, Cint, Cstring), G, s, t, a_cap, fname)
end

function glp_write_maxflow(G, s, t, a_cap, fname)
    ccall((:glp_write_maxflow, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cint, Cstring), G, s, t, a_cap, fname)
end

function glp_read_asnprob(G, v_set, a_cost, fname)
    ccall((:glp_read_asnprob, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cstring), G, v_set, a_cost, fname)
end

function glp_write_asnprob(G, v_set, a_cost, fname)
    ccall((:glp_write_asnprob, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cstring), G, v_set, a_cost, fname)
end

function glp_read_ccdata(G, v_wgt, fname)
    ccall((:glp_read_ccdata, libglpk), Cint, (Ptr{glp_graph}, Cint, Cstring), G, v_wgt, fname)
end

function glp_write_ccdata(G, v_wgt, fname)
    ccall((:glp_write_ccdata, libglpk), Cint, (Ptr{glp_graph}, Cint, Cstring), G, v_wgt, fname)
end

function glp_netgen(G, v_rhs, a_cap, a_cost, parm)
    ccall((:glp_netgen, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cint, Ptr{Cint}), G, v_rhs, a_cap, a_cost, parm)
end

function glp_netgen_prob(nprob, parm)
    ccall((:glp_netgen_prob, libglpk), Cvoid, (Cint, Ptr{Cint}), nprob, parm)
end

function glp_gridgen(G, v_rhs, a_cap, a_cost, parm)
    ccall((:glp_gridgen, libglpk), Cint, (Ptr{glp_graph}, Cint, Cint, Cint, Ptr{Cint}), G, v_rhs, a_cap, a_cost, parm)
end

function glp_rmfgen(G, s, t, a_cap, parm)
    ccall((:glp_rmfgen, libglpk), Cint, (Ptr{glp_graph}, Ptr{Cint}, Ptr{Cint}, Cint, Ptr{Cint}), G, s, t, a_cap, parm)
end

function glp_weak_comp(G, v_num)
    ccall((:glp_weak_comp, libglpk), Cint, (Ptr{glp_graph}, Cint), G, v_num)
end

function glp_strong_comp(G, v_num)
    ccall((:glp_strong_comp, libglpk), Cint, (Ptr{glp_graph}, Cint), G, v_num)
end

function glp_top_sort(G, v_num)
    ccall((:glp_top_sort, libglpk), Cint, (Ptr{glp_graph}, Cint), G, v_num)
end

function glp_wclique_exact(G, v_wgt, sol, v_set)
    ccall((:glp_wclique_exact, libglpk), Cint, (Ptr{glp_graph}, Cint, Ptr{Cdouble}, Cint), G, v_wgt, sol, v_set)
end
