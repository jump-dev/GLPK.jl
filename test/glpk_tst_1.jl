using Compat.Test
using Compat.SparseArrays
import GLPK

# Basically, same example as in the GLPK manual
# but using many more functions just to test them out

function glpk_tst_1()

    # suppress output
    prev_term_out = GLPK.term_out(GLPK.OFF)

    # define the problem as in the original GLPK manual example
    lp = GLPK.Prob()
    GLPK.set_prob_name(lp, "sample")
    GLPK.set_obj_name(lp, "OBJECTIVE")
    GLPK.set_obj_dir(lp, GLPK.MAX)
    GLPK.add_rows(lp, 3)
    GLPK.set_row_name(lp, 1, "p")
    GLPK.set_row_bnds(lp, 1, GLPK.UP, 0.0, 100.0)
    GLPK.set_row_name(lp, 2, "q")
    GLPK.set_row_bnds(lp, 2, GLPK.UP, 0.0, 600.0)
    GLPK.set_row_name(lp, 3, "r")
    GLPK.set_row_bnds(lp, 3, GLPK.UP, 0.0, 300.0)
    GLPK.add_cols(lp, 3)
    GLPK.set_col_name(lp, 1, "x1")
    GLPK.set_col_bnds(lp, 1, GLPK.LO, 0.0, 0.0)
    GLPK.set_obj_coef(lp, 1, 10.0)
    GLPK.set_col_name(lp, 2, "x2")
    GLPK.set_col_bnds(lp, 2, GLPK.LO, 0.0, 0.0)
    GLPK.set_obj_coef(lp, 2, 6.0)
    GLPK.set_col_name(lp, 3, "x3")
    GLPK.set_col_bnds(lp, 3, GLPK.LO, 0.0, 0.0)
    GLPK.set_obj_coef(lp, 3, 4.0)
    ia = zeros(Int, 9)
    ja = zeros(Int, 9)
    ar = zeros(Float64, 9)
    ia[1] = 1; ja[1] = 1; ar[1] = 1.0
    ia[2] = 1; ja[2] = 2; ar[2] = 1.0
    ia[3] = 1; ja[3] = 3; ar[3] = 1.0
    ia[4] = 2; ja[4] = 1; ar[4] = 10.0
    ia[5] = 2; ja[5] = 2; ar[5] = 4.0
    ia[6] = 2; ja[6] = 3; ar[6] = 5.0
    ia[7] = 3; ja[7] = 1; ar[7] = 2.0
    ia[8] = 3; ja[8] = 2; ar[8] = 2.0
    ia[9] = 3; ja[9] = 3; ar[9] = 6.0
    @test GLPK.check_dup(3, 3, ia, ja) == 0
    GLPK.load_matrix(lp, 9, ia, ja, ar)

    # mess up with the problem definition, test
    # alternative way to set rows and columns
    GLPK.del_rows(lp, [3])
    GLPK.add_rows(lp, 1)
    GLPK.set_row_name(lp, 3, "r")
    GLPK.set_row_bnds(lp, 3, GLPK.UP, 0.0, 300.0)

    GLPK.set_mat_row(lp, 1, nothing, nothing)
    GLPK.set_mat_row(lp, 1, 0, nothing, nothing)
    GLPK.set_mat_col(lp, 1, nothing, nothing)
    GLPK.set_mat_col(lp, 1, 0, nothing, nothing)

    GLPK.set_mat_row(lp, 1, 3, [1,2,3], [1., 1., 1.])
    GLPK.set_mat_row(lp, 1, [1,2,3], [1., 1., 1.])
    GLPK.set_mat_row(lp, 2, 3, [1,2,3], [10., 4., 5.])
    GLPK.set_mat_row(lp, 3, 3, [1,2,3], [2., 2., 6.])

    GLPK.set_mat_col(lp, 1, 3, [1,2,3], [1., 10., 2.])
    GLPK.set_mat_col(lp, 1, [1,2,3], [1., 10., 2.])
    GLPK.set_mat_col(lp, 2, 3, [1,2,3], [1., 4., 2.])
    GLPK.set_mat_col(lp, 3, 3, [1,2,3], [1., 5., 6.])

    # test setting the problem coefficients via sparse
    # matrix
    sm = sparse(ia, ja, ar)
    GLPK.load_matrix(lp, sm)

    # test random scaling of the coefficients and
    # letting the library figure out the scaling
    # and order by itself
    for i = 1:GLPK.get_num_rows(lp)
        GLPK.set_rii(lp, i, rand())
    end
    for i = 1:GLPK.get_num_cols(lp)
        GLPK.set_sjj(lp, i, rand())
    end
    GLPK.scale_prob(lp, GLPK.SF_AUTO)
    GLPK.unscale_prob(lp)
    GLPK.sort_matrix(lp)

    # test general problem properties
    @test GLPK.get_prob_name(lp) == "sample"
    @test GLPK.get_obj_name(lp) == "OBJECTIVE"
    @test GLPK.get_obj_dir(lp) == GLPK.MAX
    @test GLPK.get_num_rows(lp) == 3
    @test GLPK.get_num_cols(lp) == 3
    @test GLPK.get_obj_coef(lp, 0) == 0.0
    @test GLPK.get_num_nz(lp) == 9

    # test rows properties
    @test GLPK.get_row_name(lp, 1) == "p"
    @test GLPK.get_row_type(lp, 1) == GLPK.UP
    if VERSION >= v"0.7-"
        @test GLPK.get_row_lb(lp, 1) == -floatmax(Float64)
    else
        @test GLPK.get_row_lb(lp, 1) == -realmax(Float64)
    end
    @test GLPK.get_row_ub(lp, 1) == 100.0
    @test GLPK.get_mat_row(lp, 1) == ([1,2,3], [1.0,1.0,1.0])
    @test GLPK.get_rii(lp, 1) == 1.0

    @test GLPK.get_row_name(lp, 2) == "q"
    @test GLPK.get_row_type(lp, 2) == GLPK.UP
    if VERSION >= v"0.7-"
        @test GLPK.get_row_lb(lp, 2) == -floatmax(Float64)
    else
        @test GLPK.get_row_lb(lp, 2) == -realmax(Float64)
    end
    @test GLPK.get_row_ub(lp, 2) == 600.0
    @test GLPK.get_mat_row(lp, 2) == ([1,2,3], [10.0,4.0,5.0])
    @test GLPK.get_rii(lp, 2) == 1.0

    @test GLPK.get_row_name(lp, 3) == "r"
    @test GLPK.get_row_type(lp, 3) == GLPK.UP
    if VERSION >= v"0.7-"
        @test GLPK.get_row_lb(lp, 3) == -floatmax(Float64)
    else
        @test GLPK.get_row_lb(lp, 3) == -realmax(Float64)
    end
    @test GLPK.get_row_ub(lp, 3) == 300.0
    @test GLPK.get_mat_row(lp, 3) == ([1,2,3], [2.0,2.0,6.0])
    @test GLPK.get_rii(lp, 3) == 1.0

    # test columns properties
    @test GLPK.get_col_name(lp, 1) == "x1"
    @test GLPK.get_col_type(lp, 1) == GLPK.LO
    @test GLPK.get_col_lb(lp, 1) == 0.0
    if VERSION >= v"0.7-"
        @test GLPK.get_col_ub(lp, 1) == floatmax(Float64)
    else
        @test GLPK.get_col_ub(lp, 1) == realmax(Float64)
    end
    @test GLPK.get_obj_coef(lp, 1) == 10.0

    @test GLPK.get_col_name(lp, 2) == "x2"
    @test GLPK.get_col_type(lp, 2) == GLPK.LO
    @test GLPK.get_col_lb(lp, 2) == 0.0
    if VERSION >= v"0.7-"
        @test GLPK.get_col_ub(lp, 2) == floatmax(Float64)
    else
        @test GLPK.get_col_ub(lp, 2) == realmax(Float64)
    end
    @test GLPK.get_obj_coef(lp, 2) == 6.0

    @test GLPK.get_col_name(lp, 3) == "x3"
    @test GLPK.get_col_type(lp, 3) == GLPK.LO
    @test GLPK.get_col_lb(lp, 3) == 0.0
    if VERSION >= v"0.7-"
        @test GLPK.get_col_ub(lp, 3) == floatmax(Float64)
    else
        @test GLPK.get_col_ub(lp, 3) == realmax(Float64)
    end
    @test GLPK.get_obj_coef(lp, 3) == 4.0

    # test index creation for retreiving row/column indices
    # by name
    GLPK.create_index(lp)
    for (ri,rn) in enumerate(["p", "q", "r"])
        @test GLPK.find_row(lp, rn) == ri
    end
    @test GLPK.find_row(lp, "?") == 0
    for (ci,cn) in enumerate(["x1", "x2", "x3"])
        @test GLPK.find_col(lp, cn) == ci
    end
    @test GLPK.find_col(lp, "?") == 0
    GLPK.delete_index(lp)

    # mess up the status of rows and columns, and
    # test various basis-generation algorithms
    for i = 1:GLPK.get_num_rows(lp)
        GLPK.set_row_stat(lp, i, GLPK.BS)
    end
    for i = 1:GLPK.get_num_cols(lp)
        GLPK.set_col_stat(lp, i, GLPK.NF)
    end

    GLPK.std_basis(lp)
    GLPK.cpx_basis(lp)
    GLPK.adv_basis(lp)

    # set parameters, reset them, set again
    param = GLPK.SimplexParam()
    param.msg_lev = GLPK.MSG_ERR
    param.presolve = GLPK.ON
    GLPK.init_smcp(param)
    param["msg_lev"] = GLPK.MSG_ERR
    param["presolve"] = GLPK.ON

    # verify status before solving
    @test GLPK.get_status(lp) == GLPK.UNDEF
    @test GLPK.get_prim_stat(lp) == GLPK.UNDEF
    @test GLPK.get_dual_stat(lp) == GLPK.UNDEF

    # now solve, simplex algorithm
    flag = GLPK.simplex(lp, param)

    # verify status after solving
    @test flag == 0
    @test GLPK.get_status(lp) == GLPK.OPT
    @test GLPK.get_prim_stat(lp) == GLPK.FEAS
    @test GLPK.get_dual_stat(lp) == GLPK.FEAS
    @test GLPK.get_unbnd_ray(lp) == 0

    ae_max, ae_ind, re_max, re_ind = GLPK.check_kkt(lp, GLPK.SOL, GLPK.KKT_PE)
    @test ae_max == 0 ? (ae_ind == 0) : (1 <= ae_ind <= GLPK.get_num_rows(lp))
    @test re_max == 0 ? (re_ind == 0) : (1 <= re_ind <= GLPK.get_num_rows(lp))
    ae_max, ae_ind, re_max, re_ind = GLPK.check_kkt(lp, GLPK.SOL, GLPK.KKT_PB)
    @test ae_max == 0 ? (ae_ind == 0) : (1 <= ae_ind <= GLPK.get_num_rows(lp) + GLPK.get_num_cols(lp))
    @test re_max == 0 ? (re_ind == 0) : (1 <= re_ind <= GLPK.get_num_rows(lp) + GLPK.get_num_cols(lp))

    # solve again, using simplex with exact arithmetics
    flag = GLPK.exact(lp, param)

    # verify status after solving
    @test flag == 0
    @test GLPK.get_status(lp) == GLPK.OPT
    @test GLPK.get_prim_stat(lp) == GLPK.FEAS
    @test GLPK.get_dual_stat(lp) == GLPK.FEAS
    @test GLPK.get_unbnd_ray(lp) == 0

    ae_max, ae_ind, re_max, re_ind = GLPK.check_kkt(lp, GLPK.SOL, GLPK.KKT_PE)
    @test ae_max == 0 ? (ae_ind == 0) : (1 <= ae_ind <= GLPK.get_num_rows(lp))
    @test re_max == 0 ? (re_ind == 0) : (1 <= re_ind <= GLPK.get_num_rows(lp))
    ae_max, ae_ind, re_max, re_ind = GLPK.check_kkt(lp, GLPK.SOL, GLPK.KKT_PB)
    @test ae_max == 0 ? (ae_ind == 0) : (1 <= ae_ind <= GLPK.get_num_rows(lp) + GLPK.get_num_cols(lp))
    @test re_max == 0 ? (re_ind == 0) : (1 <= re_ind <= GLPK.get_num_rows(lp) + GLPK.get_num_cols(lp))

    # verify results
    tol = 1e-10
    @test GLPK.get_row_stat(lp, 1) == GLPK.NU
    @test GLPK.get_row_prim(lp, 1) == 100.0
    @test abs(GLPK.get_row_dual(lp, 1) - 10.0 / 3) < tol

    @test GLPK.get_row_stat(lp, 2) == GLPK.NU
    @test GLPK.get_row_prim(lp, 2) == 600.0
    @test abs(GLPK.get_row_dual(lp, 2) - 2.0 / 3) < tol

    @test GLPK.get_row_stat(lp, 3) == GLPK.BS
    @test GLPK.get_row_prim(lp, 3) == 200.0
    @test abs(GLPK.get_row_dual(lp, 3) - 0.0) < tol

    @test GLPK.get_col_stat(lp, 1) == GLPK.BS
    @test abs(GLPK.get_col_prim(lp, 1) - 100.0 / 3) < tol
    @test abs(GLPK.get_col_dual(lp, 1) - 0.0) < tol

    @test GLPK.get_col_stat(lp, 2) == GLPK.BS
    @test abs(GLPK.get_col_prim(lp, 2) - 200.0 / 3) < tol
    @test abs(GLPK.get_col_dual(lp, 2) - 0.0) < tol

    @test GLPK.get_col_stat(lp, 3) == GLPK.NL
    @test abs(GLPK.get_col_prim(lp, 3) - 0.0) < tol
    @test abs(GLPK.get_col_dual(lp, 3) - (-8.0 / 3)) < tol

    @test abs(GLPK.get_obj_val(lp) - 2200.0 / 3) < tol

    # solve again, using interior point method
    param = GLPK.InteriorParam()
    param.msg_lev = GLPK.MSG_ERR
    flag = GLPK.interior(lp, param)

    # verify status after solving
    @test flag == 0
    @test GLPK.get_status(lp) == GLPK.OPT
    @test GLPK.get_prim_stat(lp) == GLPK.FEAS
    @test GLPK.get_dual_stat(lp) == GLPK.FEAS
    @test GLPK.get_unbnd_ray(lp) == 0

    ae_max, ae_ind, re_max, re_ind = GLPK.check_kkt(lp, GLPK.IPT, GLPK.KKT_PE)
    @test ae_max == 0 ? (ae_ind == 0) : (1 <= ae_ind <= GLPK.get_num_rows(lp))
    @test re_max == 0 ? (re_ind == 0) : (1 <= re_ind <= GLPK.get_num_rows(lp))
    ae_max, ae_ind, re_max, re_ind = GLPK.check_kkt(lp, GLPK.IPT, GLPK.KKT_PB)
    @test ae_max == 0 ? (ae_ind == 0) : (1 <= ae_ind <= GLPK.get_num_rows(lp) + GLPK.get_num_cols(lp))
    @test re_max == 0 ? (re_ind == 0) : (1 <= re_ind <= GLPK.get_num_rows(lp) + GLPK.get_num_cols(lp))
    ae_max, ae_ind, re_max, re_ind = GLPK.check_kkt(lp, GLPK.IPT, GLPK.KKT_DE)
    @test ae_max == 0 ? (ae_ind == 0) : (1 <= ae_ind - GLPK.get_num_rows(lp) <= GLPK.get_num_cols(lp))
    @test re_max == 0 ? (re_ind == 0) : (1 <= re_ind - GLPK.get_num_rows(lp) <= GLPK.get_num_cols(lp))
    ae_max, ae_ind, re_max, re_ind = GLPK.check_kkt(lp, GLPK.IPT, GLPK.KKT_DB)
    @test ae_max == 0 ? (ae_ind == 0) : (1 <= ae_ind <= GLPK.get_num_rows(lp) + GLPK.get_num_cols(lp))
    @test re_max == 0 ? (re_ind == 0) : (1 <= re_ind <= GLPK.get_num_rows(lp) + GLPK.get_num_cols(lp))

    # verify results
    tol = 1e-10
    @test GLPK.get_row_stat(lp, 1) == GLPK.NU
    @test GLPK.get_row_prim(lp, 1) == 100.0
    @test abs(GLPK.get_row_dual(lp, 1) - 10.0 / 3) < tol

    @test GLPK.get_row_stat(lp, 2) == GLPK.NU
    @test GLPK.get_row_prim(lp, 2) == 600.0
    @test abs(GLPK.get_row_dual(lp, 2) - 2.0 / 3) < tol

    @test GLPK.get_row_stat(lp, 3) == GLPK.BS
    @test GLPK.get_row_prim(lp, 3) == 200.0
    @test abs(GLPK.get_row_dual(lp, 3) - 0.0) < tol

    @test GLPK.get_col_stat(lp, 1) == GLPK.BS
    @test abs(GLPK.get_col_prim(lp, 1) - 100.0 / 3) < tol
    @test abs(GLPK.get_col_dual(lp, 1) - 0.0) < tol

    @test GLPK.get_col_stat(lp, 2) == GLPK.BS
    @test abs(GLPK.get_col_prim(lp, 2) - 200.0 / 3) < tol
    @test abs(GLPK.get_col_dual(lp, 2) - 0.0) < tol

    @test GLPK.get_col_stat(lp, 3) == GLPK.NL
    @test abs(GLPK.get_col_prim(lp, 3) - 0.0) < tol
    @test abs(GLPK.get_col_dual(lp, 3) - (-8.0 / 3)) < tol

    @test abs(GLPK.get_obj_val(lp) - 2200.0 / 3) < tol

    # reset terminal output status
    GLPK.term_out(prev_term_out)
end

glpk_tst_1()
