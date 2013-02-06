using Test
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

    datadir = joinpath(Pkg.dir(), "GLPK", "test", "data")

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
        @test_fails GLPK.get_bhead(mip, rows+1)
        @test_fails GLPK.get_row_bind(mip, -1)
    end
    for i = 1:cols
        cb = GLPK.get_col_bind(mip, i)
    end
    if GLPK.jl_get_preemptive_check()
        @test_fails GLPK.get_col_bind(mip, -1)
    end

    x = randn(int64(rows))
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
        ia, ja, val = findn_nzs(sprand(int64(cols), 1, 0.5))
        GLPK.transform_row(mip, int32(ia), val)
    end
    for i = 1:100
        ia, ja, val = findn_nzs(sprand(int64(rows), 1, 0.5))
        GLPK.transform_col(mip, int32(ia), val)
    end

    for i = 1:100
        ia, ja, val = findn_nzs(sprand(int64(rows), 1, 0.5))
        dir = 2 * randbool() - 1
        eps = 1e-9

        basic = map(x->var_is_basic(mip, x), ia)

        ia = ia[basic]
        val = val[basic]

        prt = GLPK.prim_rtest(mip, ia, val, dir, eps)
    end

    for i = 1:100
        ia, ja, val = findn_nzs(sprand(int64(cols), 1, 0.5))
        dir = 2 * randbool() - 1
        eps = 1e-9

        basic = map(x->var_is_basic(mip, x), ia)

        ia = ia[!basic]
        val = val[!basic]

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
        GLPK.print_sol(mip, "sol_printable.txt")
        GLPK.write_sol(mip, "sol.txt")
        GLPK.read_sol(mip, "sol.txt")
    finally
        if isfile(sol_printable)
            rm(sol_printable)
        end
        if isfile(sol_plain)
            rm(sol_plain)
        end
    end

    ipt_printable = tempname()
    ipt_plain = tempname()

    try
        GLPK.print_ipt(mip, "ipt_printable.txt")
        GLPK.write_ipt(mip, "ipt.txt")
        GLPK.read_ipt(mip, "ipt.txt")
    finally
        if isfile(ipt_printable)
            rm(ipt_printable)
        end
        if isfile(ipt_plain)
            rm(ipt_plain)
        end
    end

    mip_printable = tempname()
    mip_plain = tempname()

    try
        GLPK.print_mip(mip, "mip_printable.txt")
        GLPK.write_mip(mip, "mip.txt")
        GLPK.read_mip(mip, "mip.txt")
    finally
        if isfile(mip_printable)
            rm(mip_printable)
        end
        if isfile(mip_plain)
            rm(mip_plain)
        end
    end

    tee_out = tempname()
    try
        GLPK.open_tee(tee_out)
        GLPK.close_tee()
    finally
        if isfile(tee_out)
            rm(tee_out)
        end
    end

    GLPK.term_out(prev_term_out)
end

glpk_tst_4()
