import GLPK
import SparseArrays
using Test

# Basically, same example as in the GLPK manual
# but using many more functions just to test them out

function glpk_tst_1()

    # suppress output
    prev_term_out = GLPK.glp_term_out(GLPK.GLP_OFF)

    # define the problem as in the original GLPK manual example
    lp = GLPK.glp_create_prob()
    GLPK.glp_set_prob_name(lp, "sample")
    GLPK.glp_set_obj_name(lp, "OBJECTIVE")
    GLPK.glp_set_obj_dir(lp, GLPK.GLP_MAX)
    GLPK.glp_add_rows(lp, 3)
    GLPK.glp_set_row_name(lp, 1, "p")
    GLPK.glp_set_row_bnds(lp, 1, GLPK.GLP_UP, 0.0, 100.0)
    GLPK.glp_set_row_name(lp, 2, "q")
    GLPK.glp_set_row_bnds(lp, 2, GLPK.GLP_UP, 0.0, 600.0)
    GLPK.glp_set_row_name(lp, 3, "r")
    GLPK.glp_set_row_bnds(lp, 3, GLPK.GLP_UP, 0.0, 300.0)
    GLPK.glp_add_cols(lp, 3)
    GLPK.glp_set_col_name(lp, 1, "x1")
    GLPK.glp_set_col_bnds(lp, 1, GLPK.GLP_LO, 0.0, 0.0)
    GLPK.glp_set_obj_coef(lp, 1, 10.0)
    GLPK.glp_set_col_name(lp, 2, "x2")
    GLPK.glp_set_col_bnds(lp, 2, GLPK.GLP_LO, 0.0, 0.0)
    GLPK.glp_set_obj_coef(lp, 2, 6.0)
    GLPK.glp_set_col_name(lp, 3, "x3")
    GLPK.glp_set_col_bnds(lp, 3, GLPK.GLP_LO, 0.0, 0.0)
    GLPK.glp_set_obj_coef(lp, 3, 4.0)
    ia = zeros(Cint, 9)
    ja = zeros(Cint, 9)
    ar = zeros(Float64, 9)
    ia[1] = 1
    ja[1] = 1
    ar[1] = 1.0
    ia[2] = 1
    ja[2] = 2
    ar[2] = 1.0
    ia[3] = 1
    ja[3] = 3
    ar[3] = 1.0
    ia[4] = 2
    ja[4] = 1
    ar[4] = 10.0
    ia[5] = 2
    ja[5] = 2
    ar[5] = 4.0
    ia[6] = 2
    ja[6] = 3
    ar[6] = 5.0
    ia[7] = 3
    ja[7] = 1
    ar[7] = 2.0
    ia[8] = 3
    ja[8] = 2
    ar[8] = 2.0
    ia[9] = 3
    ja[9] = 3
    ar[9] = 6.0
    @test GLPK.glp_check_dup(
        3,
        3,
        length(ia),
        GLPK.offset(ia),
        GLPK.offset(ja),
    ) == 0
    GLPK.glp_load_matrix(
        lp,
        9,
        GLPK.offset(ia),
        GLPK.offset(ja),
        GLPK.offset(ar),
    )

    # mess up with the problem definition, test
    # alternative way to set rows and columns
    r = Cint[3]
    GLPK.glp_del_rows(lp, length(r), GLPK.offset(r))
    GLPK.glp_add_rows(lp, 1)
    GLPK.glp_set_row_name(lp, 3, "r")
    GLPK.glp_set_row_bnds(lp, 3, GLPK.GLP_UP, 0.0, 300.0)

    GLPK.glp_set_mat_row(lp, 1, 0, C_NULL, C_NULL)
    GLPK.glp_set_mat_col(lp, 1, 0, C_NULL, C_NULL)

    indices = Cint[1, 2, 3]
    indices_p = GLPK.offset(indices)
    row1 = [1.0, 1.0, 1.0]
    row2 = [10.0, 4.0, 5.0]
    row3 = [2.0, 2.0, 6.0]
    GLPK.glp_set_mat_row(lp, 1, 3, indices_p, GLPK.offset(row1))
    GLPK.glp_set_mat_row(lp, 2, 3, indices_p, GLPK.offset(row2))
    GLPK.glp_set_mat_row(lp, 3, 3, indices_p, GLPK.offset(row3))
    col1 = [1.0, 10.0, 2.0]
    col2 = [1.0, 4.0, 2.0]
    col3 = [1.0, 5.0, 6.0]
    GLPK.glp_set_mat_col(lp, 1, 3, indices_p, GLPK.offset(col1))
    GLPK.glp_set_mat_col(lp, 2, 3, indices_p, GLPK.offset(col2))
    GLPK.glp_set_mat_col(lp, 3, 3, indices_p, GLPK.offset(col3))

    # test random scaling of the coefficients and
    # letting the library figure out the scaling
    # and order by itself
    for i in 1:GLPK.glp_get_num_rows(lp)
        GLPK.glp_set_rii(lp, i, rand())
    end
    for i in 1:GLPK.glp_get_num_cols(lp)
        GLPK.glp_set_sjj(lp, i, rand())
    end
    GLPK.glp_scale_prob(lp, GLPK.GLP_SF_AUTO)
    GLPK.glp_unscale_prob(lp)
    GLPK.glp_sort_matrix(lp)

    # test general problem properties
    @test unsafe_string(GLPK.glp_get_prob_name(lp)) == "sample"
    @test unsafe_string(GLPK.glp_get_obj_name(lp)) == "OBJECTIVE"
    @test GLPK.glp_get_obj_dir(lp) == GLPK.GLP_MAX
    @test GLPK.glp_get_num_rows(lp) == 3
    @test GLPK.glp_get_num_cols(lp) == 3
    @test GLPK.glp_get_obj_coef(lp, 0) == 0.0
    @test GLPK.glp_get_num_nz(lp) == 9

    # test rows properties
    @test unsafe_string(GLPK.glp_get_row_name(lp, 1)) == "p"
    @test GLPK.glp_get_row_type(lp, 1) == GLPK.GLP_UP
    @test GLPK.glp_get_row_lb(lp, 1) == -floatmax(Float64)
    @test GLPK.glp_get_row_ub(lp, 1) == 100.0
    indices = zeros(Cint, 3)
    coefficients = zeros(Cdouble, 3)
    @test GLPK.glp_get_mat_row(
        lp,
        1,
        GLPK.offset(indices),
        GLPK.offset(coefficients),
    ) == 3
    @test indices == [1, 2, 3]
    @test coefficients == [1.0, 1.0, 1.0]
    @test GLPK.glp_get_rii(lp, 1) == 1.0

    @test unsafe_string(GLPK.glp_get_row_name(lp, 2)) == "q"
    @test GLPK.glp_get_row_type(lp, 2) == GLPK.GLP_UP
    @test GLPK.glp_get_row_lb(lp, 2) == -floatmax(Float64)
    @test GLPK.glp_get_row_ub(lp, 2) == 600.0
    @test GLPK.glp_get_mat_row(
        lp,
        2,
        GLPK.offset(indices),
        GLPK.offset(coefficients),
    ) == 3
    @test indices == [1, 2, 3]
    @test coefficients == [10.0, 4.0, 5.0]
    @test GLPK.glp_get_rii(lp, 2) == 1.0

    @test unsafe_string(GLPK.glp_get_row_name(lp, 3)) == "r"
    @test GLPK.glp_get_row_type(lp, 3) == GLPK.GLP_UP
    @test GLPK.glp_get_row_lb(lp, 3) == -floatmax(Float64)
    @test GLPK.glp_get_row_ub(lp, 3) == 300.0
    @test GLPK.glp_get_mat_row(
        lp,
        3,
        GLPK.offset(indices),
        GLPK.offset(coefficients),
    ) == 3
    @test indices == [1, 2, 3]
    @test coefficients == [2.0, 2.0, 6.0]
    @test GLPK.glp_get_rii(lp, 3) == 1.0

    # test columns properties
    @test unsafe_string(GLPK.glp_get_col_name(lp, 1)) == "x1"
    @test GLPK.glp_get_col_type(lp, 1) == GLPK.GLP_LO
    @test GLPK.glp_get_col_lb(lp, 1) == 0.0
    @test GLPK.glp_get_col_ub(lp, 1) == floatmax(Float64)
    @test GLPK.glp_get_obj_coef(lp, 1) == 10.0

    @test unsafe_string(GLPK.glp_get_col_name(lp, 2)) == "x2"
    @test GLPK.glp_get_col_type(lp, 2) == GLPK.GLP_LO
    @test GLPK.glp_get_col_lb(lp, 2) == 0.0
    @test GLPK.glp_get_col_ub(lp, 2) == floatmax(Float64)
    @test GLPK.glp_get_obj_coef(lp, 2) == 6.0

    @test unsafe_string(GLPK.glp_get_col_name(lp, 3)) == "x3"
    @test GLPK.glp_get_col_type(lp, 3) == GLPK.GLP_LO
    @test GLPK.glp_get_col_lb(lp, 3) == 0.0
    @test GLPK.glp_get_col_ub(lp, 3) == floatmax(Float64)
    @test GLPK.glp_get_obj_coef(lp, 3) == 4.0

    # test index creation for retreiving row/column indices
    # by name
    GLPK.glp_create_index(lp)
    for (ri, rn) in enumerate(["p", "q", "r"])
        @test GLPK.glp_find_row(lp, rn) == ri
    end
    @test GLPK.glp_find_row(lp, "?") == 0
    for (ci, cn) in enumerate(["x1", "x2", "x3"])
        @test GLPK.glp_find_col(lp, cn) == ci
    end
    @test GLPK.glp_find_col(lp, "?") == 0
    GLPK.glp_delete_index(lp)

    # mess up the status of rows and columns, and
    # test various basis-generation algorithms
    for i in 1:GLPK.glp_get_num_rows(lp)
        GLPK.glp_set_row_stat(lp, i, GLPK.GLP_BS)
    end
    for i in 1:GLPK.glp_get_num_cols(lp)
        GLPK.glp_set_col_stat(lp, i, GLPK.GLP_NF)
    end

    GLPK.glp_std_basis(lp)
    GLPK.glp_cpx_basis(lp)
    GLPK.glp_adv_basis(lp, 0)

    # # set parameters, reset them, set again
    param = GLPK.glp_smcp()
    GLPK.glp_init_smcp(param)
    param.msg_lev = GLPK.GLP_MSG_ERR
    param.presolve = GLPK.GLP_ON

    # verify status before solving
    @test GLPK.glp_get_status(lp) == GLPK.GLP_UNDEF
    @test GLPK.glp_get_prim_stat(lp) == GLPK.GLP_UNDEF
    @test GLPK.glp_get_dual_stat(lp) == GLPK.GLP_UNDEF

    # now solve, simplex algorithm
    flag = GLPK.glp_simplex(lp, param)

    # verify status after solving
    @test flag == 0
    @test GLPK.glp_get_status(lp) == GLPK.GLP_OPT
    @test GLPK.glp_get_prim_stat(lp) == GLPK.GLP_FEAS
    @test GLPK.glp_get_dual_stat(lp) == GLPK.GLP_FEAS
    @test GLPK.glp_get_unbnd_ray(lp) == 0

    ae_max = Ref{Cdouble}()
    ae_ind = Ref{Cint}()
    re_max = Ref{Cdouble}()
    re_ind = Ref{Cint}()
    GLPK.glp_check_kkt(
        lp,
        GLPK.GLP_SOL,
        GLPK.GLP_KKT_PE,
        ae_max,
        ae_ind,
        re_max,
        re_ind,
    )
    @test ae_max[] == 0 ? (ae_ind[] == 0) :
          (1 <= ae_ind[] <= GLPK.glp_get_num_rows(lp))
    @test re_max[] == 0 ? (re_ind[] == 0) :
          (1 <= re_ind[] <= GLPK.glp_get_num_rows(lp))
    GLPK.glp_check_kkt(
        lp,
        GLPK.GLP_SOL,
        GLPK.GLP_KKT_PB,
        ae_max,
        ae_ind,
        re_max,
        re_ind,
    )
    @test ae_max[] == 0 ? (ae_ind[] == 0) :
          (
        1 <= ae_ind[] <= GLPK.glp_get_num_rows(lp) + GLPK.glp_get_num_cols(lp)
    )
    @test re_max[] == 0 ? (re_ind[] == 0) :
          (
        1 <= re_ind[] <= GLPK.glp_get_num_rows(lp) + GLPK.glp_get_num_cols(lp)
    )

    # solve again, using simplex with exact arithmetics
    flag = GLPK.glp_exact(lp, param)

    # verify status after solving
    @test flag == 0
    @test GLPK.glp_get_status(lp) == GLPK.GLP_OPT
    @test GLPK.glp_get_prim_stat(lp) == GLPK.GLP_FEAS
    @test GLPK.glp_get_dual_stat(lp) == GLPK.GLP_FEAS
    @test GLPK.glp_get_unbnd_ray(lp) == 0

    ae_max = Ref{Cdouble}()
    ae_ind = Ref{Cint}()
    re_max = Ref{Cdouble}()
    re_ind = Ref{Cint}()
    GLPK.glp_check_kkt(
        lp,
        GLPK.GLP_SOL,
        GLPK.GLP_KKT_PE,
        ae_max,
        ae_ind,
        re_max,
        re_ind,
    )
    @test ae_max[] == 0 ? (ae_ind[] == 0) :
          (1 <= ae_ind[] <= GLPK.glp_get_num_rows(lp))
    @test re_max[] == 0 ? (re_ind[] == 0) :
          (1 <= re_ind[] <= GLPK.glp_get_num_rows(lp))
    GLPK.glp_check_kkt(
        lp,
        GLPK.GLP_SOL,
        GLPK.GLP_KKT_PB,
        ae_max,
        ae_ind,
        re_max,
        re_ind,
    )
    @test ae_max[] == 0 ? (ae_ind[] == 0) :
          (
        1 <= ae_ind[] <= GLPK.glp_get_num_rows(lp) + GLPK.glp_get_num_cols(lp)
    )
    @test re_max[] == 0 ? (re_ind[] == 0) :
          (
        1 <= re_ind[] <= GLPK.glp_get_num_rows(lp) + GLPK.glp_get_num_cols(lp)
    )

    # verify results
    tol = 1e-10
    @test GLPK.glp_get_row_stat(lp, 1) == GLPK.GLP_NU
    @test GLPK.glp_get_row_prim(lp, 1) == 100.0
    @test abs(GLPK.glp_get_row_dual(lp, 1) - 10.0 / 3) < tol

    @test GLPK.glp_get_row_stat(lp, 2) == GLPK.GLP_NU
    @test GLPK.glp_get_row_prim(lp, 2) == 600.0
    @test abs(GLPK.glp_get_row_dual(lp, 2) - 2.0 / 3) < tol

    @test GLPK.glp_get_row_stat(lp, 3) == GLPK.GLP_BS
    @test GLPK.glp_get_row_prim(lp, 3) == 200.0
    @test abs(GLPK.glp_get_row_dual(lp, 3) - 0.0) < tol

    @test GLPK.glp_get_col_stat(lp, 1) == GLPK.GLP_BS
    @test abs(GLPK.glp_get_col_prim(lp, 1) - 100.0 / 3) < tol
    @test abs(GLPK.glp_get_col_dual(lp, 1) - 0.0) < tol

    @test GLPK.glp_get_col_stat(lp, 2) == GLPK.GLP_BS
    @test abs(GLPK.glp_get_col_prim(lp, 2) - 200.0 / 3) < tol
    @test abs(GLPK.glp_get_col_dual(lp, 2) - 0.0) < tol

    @test GLPK.glp_get_col_stat(lp, 3) == GLPK.GLP_NL
    @test abs(GLPK.glp_get_col_prim(lp, 3) - 0.0) < tol
    @test abs(GLPK.glp_get_col_dual(lp, 3) - (-8.0 / 3)) < tol

    @test abs(GLPK.glp_get_obj_val(lp) - 2200.0 / 3) < tol

    # solve again, using interior point method
    param = GLPK.glp_iptcp()
    GLPK.glp_init_iptcp(param)
    param.msg_lev = GLPK.GLP_MSG_ERR
    flag = GLPK.glp_interior(lp, param)

    # verify status after solving
    @test flag == 0
    @test GLPK.glp_get_status(lp) == GLPK.GLP_OPT
    @test GLPK.glp_get_prim_stat(lp) == GLPK.GLP_FEAS
    @test GLPK.glp_get_dual_stat(lp) == GLPK.GLP_FEAS
    @test GLPK.glp_get_unbnd_ray(lp) == 0

    ae_max = Ref{Cdouble}()
    ae_ind = Ref{Cint}()
    re_max = Ref{Cdouble}()
    re_ind = Ref{Cint}()
    GLPK.glp_check_kkt(
        lp,
        GLPK.GLP_IPT,
        GLPK.GLP_KKT_PE,
        ae_max,
        ae_ind,
        re_max,
        re_ind,
    )
    @test ae_max[] == 0 ? (ae_ind[] == 0) :
          (1 <= ae_ind[] <= GLPK.glp_get_num_rows(lp))
    @test re_max[] == 0 ? (re_ind[] == 0) :
          (1 <= re_ind[] <= GLPK.glp_get_num_rows(lp))
    GLPK.glp_check_kkt(
        lp,
        GLPK.GLP_IPT,
        GLPK.GLP_KKT_PB,
        ae_max,
        ae_ind,
        re_max,
        re_ind,
    )
    @test ae_max[] == 0 ? (ae_ind[] == 0) :
          (
        1 <= ae_ind[] <= GLPK.glp_get_num_rows(lp) + GLPK.glp_get_num_cols(lp)
    )
    @test re_max[] == 0 ? (re_ind[] == 0) :
          (
        1 <= re_ind[] <= GLPK.glp_get_num_rows(lp) + GLPK.glp_get_num_cols(lp)
    )
    GLPK.glp_check_kkt(
        lp,
        GLPK.GLP_IPT,
        GLPK.GLP_KKT_DE,
        ae_max,
        ae_ind,
        re_max,
        re_ind,
    )
    @test ae_max[] == 0 ? (ae_ind[] == 0) :
          (
        1 <= ae_ind[] - GLPK.glp_get_num_rows(lp) <= GLPK.glp_get_num_cols(lp)
    )
    @test re_max[] == 0 ? (re_ind[] == 0) :
          (
        1 <= re_ind[] - GLPK.glp_get_num_rows(lp) <= GLPK.glp_get_num_cols(lp)
    )
    GLPK.glp_check_kkt(
        lp,
        GLPK.GLP_IPT,
        GLPK.GLP_KKT_DB,
        ae_max,
        ae_ind,
        re_max,
        re_ind,
    )
    @test ae_max[] == 0 ? (ae_ind[] == 0) :
          (
        1 <= ae_ind[] <= GLPK.glp_get_num_rows(lp) + GLPK.glp_get_num_cols(lp)
    )
    @test re_max[] == 0 ? (re_ind[] == 0) :
          (
        1 <= re_ind[] <= GLPK.glp_get_num_rows(lp) + GLPK.glp_get_num_cols(lp)
    )

    # verify results
    tol = 1e-10
    @test GLPK.glp_get_row_stat(lp, 1) == GLPK.GLP_NU
    @test GLPK.glp_get_row_prim(lp, 1) == 100.0
    @test abs(GLPK.glp_get_row_dual(lp, 1) - 10.0 / 3) < tol

    @test GLPK.glp_get_row_stat(lp, 2) == GLPK.GLP_NU
    @test GLPK.glp_get_row_prim(lp, 2) == 600.0
    @test abs(GLPK.glp_get_row_dual(lp, 2) - 2.0 / 3) < tol

    @test GLPK.glp_get_row_stat(lp, 3) == GLPK.GLP_BS
    @test GLPK.glp_get_row_prim(lp, 3) == 200.0
    @test abs(GLPK.glp_get_row_dual(lp, 3) - 0.0) < tol

    @test GLPK.glp_get_col_stat(lp, 1) == GLPK.GLP_BS
    @test abs(GLPK.glp_get_col_prim(lp, 1) - 100.0 / 3) < tol
    @test abs(GLPK.glp_get_col_dual(lp, 1) - 0.0) < tol

    @test GLPK.glp_get_col_stat(lp, 2) == GLPK.GLP_BS
    @test abs(GLPK.glp_get_col_prim(lp, 2) - 200.0 / 3) < tol
    @test abs(GLPK.glp_get_col_dual(lp, 2) - 0.0) < tol

    @test GLPK.glp_get_col_stat(lp, 3) == GLPK.GLP_NL
    @test abs(GLPK.glp_get_col_prim(lp, 3) - 0.0) < tol
    @test abs(GLPK.glp_get_col_dual(lp, 3) - (-8.0 / 3)) < tol

    @test abs(GLPK.glp_get_obj_val(lp) - 2200.0 / 3) < tol

    # reset terminal output status
    return GLPK.glp_term_out(prev_term_out)
end

glpk_tst_1()
