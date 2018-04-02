using Compat.Test
import GLPK

# Test Basis Factorization stuff, plus
# miscellaneous functions

function var_is_basic(lp::GLPK.Prob, k::Integer)
    rows = GLPK.get_num_rows(lp)
    if k <= rows
        return GLPK.get_row_stat(lp, k) == GLPK.BS
    else
        return GLPK.get_col_stat(lp, k-rows) == GLPK.BS
    end
end

function glpk_tst_4()
    prev_term_out = GLPK.term_out(GLPK.OFF)

    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    mip = GLPK.Prob()
    tran = GLPK.MathProgWorkspace()
    GLPK.mpl_read_model(tran, joinpath(datadir, "sudoku.mod"), 1)
    GLPK.mpl_read_data(tran, joinpath(datadir, "sudoku.dat"))
    GLPK.mpl_generate(tran, nothing)
    GLPK.mpl_build_prob(tran, mip)
    @test GLPK.simplex(mip, nothing) == 0
    @test GLPK.intopt(mip, nothing) == 0
    GLPK.mpl_postsolve(tran, mip, GLPK.MIP)

    @test GLPK.warm_up(mip) == 0

    @test GLPK.factorize(mip) == 0

    @test GLPK.bf_updated(mip) == 0

    bfp = GLPK.BasisFactParam()
    GLPK.get_bfcp(mip, bfp)
    bfp["piv_tol"] = 0.5
    #bfp.piv_tol = 0.5
    GLPK.set_bfcp(mip, bfp)

    rows = GLPK.get_num_rows(mip)
    cols = GLPK.get_num_cols(mip)

    for i = 1:rows
        bh = GLPK.get_bhead(mip, i)
        rb = GLPK.get_row_bind(mip, i)
    end
    if GLPK.jl_get_preemptive_check()
        @glpk_test_throws GLPK.get_bhead(mip, rows+1)
        @glpk_test_throws GLPK.get_row_bind(mip, -1)
    end
    for i = 1:cols
        cb = GLPK.get_col_bind(mip, i)
    end
    if GLPK.jl_get_preemptive_check()
        @glpk_test_throws GLPK.get_col_bind(mip, -1)
    end

    x = randn(Int(rows))
    GLPK.ftran(mip, x)
    GLPK.btran(mip, x)

    GLPK.adv_basis(mip)
    @test GLPK.factorize(mip) == 0

    for i = 1:rows+cols
        if !var_is_basic(mip, i)
            continue
        end

        ind = Int32[]
        val = Float64[]
        GLPK.eval_tab_row(mip, i, ind, val)

        ind, val = GLPK.eval_tab_row(mip, i)
    end

    for i = 1:rows+cols
        if var_is_basic(mip, i)
            continue
        end

        ind = Int32[]
        val = Float64[]
        GLPK.eval_tab_col(mip, i, ind, val)

        ind, val = GLPK.eval_tab_col(mip, i)
    end

    for i = 1:100
        ia, ja, val = findnz(sprand(Int(cols), 1, 0.5))
        # TODO: replace the map with round.(Int32, ia) when julia 0.5 support is dropped
        GLPK.transform_row(mip, map(x->round(Int32, x), ia), val)
    end

    for i = 1:100
        ia, ja, val = findnz(sprand(Int(rows), 1, 0.5))
        # TODO: replace the map with round.(Int32, ia) when julia 0.5 support is dropped
        GLPK.transform_col(mip, map(x->round(Int32, x), ia), val)
    end

    for i = 1:100
        ia, ja, val = findnz(sprand(Int(rows), 1, 0.5))
        dir = 2 * rand(Bool) - 1
        eps = 1e-9

        basic = map(x->var_is_basic(mip, x), ia)

        ia = ia[basic]
        val = val[basic]

        prt = GLPK.prim_rtest(mip, ia, val, dir, eps)
    end

    for i = 1:100
        ia, ja, val = findnz(sprand(Int(cols), 1, 0.5))
        dir = 2 * rand(Bool) - 1
        eps = 1e-9

        nonbasic = map(x->!var_is_basic(mip, x), ia)

        ia = ia[nonbasic]
        val = val[nonbasic]

        drt = GLPK.dual_rtest(mip, ia, val, dir, eps)
    end

    for i = 1:rows+cols
        if var_is_basic(mip, i)
            ac = GLPK.analyze_coef(mip, i)
        else
            ab = GLPK.analyze_bound(mip, i)
        end
    end

    GLPK.init_env()

    p = GLPK.malloc(100 * sizeof(Int))
    mu = GLPK.mem_usage()
    GLPK.free(p)
    mu = GLPK.mem_usage()
    GLPK.mem_limit(10000)

    p = GLPK.calloc(100, sizeof(Int))
    mu = GLPK.mem_usage()
    GLPK.free(p)
    mu = GLPK.mem_usage()

    #GLPK.free_env()

    sol_printable = tempname()
    sol_plain = tempname()

    try
        GLPK.print_sol(mip, sol_printable)
        GLPK.write_sol(mip, sol_plain)
        GLPK.read_sol(mip, sol_plain)
    finally
        isfile(sol_printable) && rm(sol_printable)
        isfile(sol_plain) && rm(sol_plain)
    end

    ipt_printable = tempname()
    ipt_plain = tempname()

    try
        GLPK.print_ipt(mip, ipt_printable)
        GLPK.write_ipt(mip, ipt_plain)
        GLPK.read_ipt(mip, ipt_plain)
    finally
        isfile(ipt_printable) && rm(ipt_printable)
        isfile(ipt_plain) && rm(ipt_plain)
    end

    mip_printable = tempname()
    mip_plain = tempname()

    try
        GLPK.print_mip(mip, mip_printable)
        GLPK.write_mip(mip, mip_plain)
        GLPK.read_mip(mip, mip_plain)
    finally
        isfile(mip_printable) && rm(mip_printable)
        isfile(mip_plain) && rm(mip_plain)
    end

    tee_out = tempname()
    try
        GLPK.open_tee(tee_out)
        GLPK.close_tee()
    finally
        isfile(tee_out) && rm(tee_out)
    end

    GLPK.term_out(prev_term_out)
end

glpk_tst_4()
