# Copyright (c) 2012 Carlo Baldassi and GLPK.jl contributors
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

import GLPK
import SparseArrays
using Test

# Test Basis Factorization stuff, plus
# miscellaneous functions

function var_is_basic(lp::Ptr, k::Integer)
    rows = GLPK.glp_get_num_rows(lp)
    if k <= rows
        return GLPK.glp_get_row_stat(lp, k) == GLPK.GLP_BS
    else
        return GLPK.glp_get_col_stat(lp, k - rows) == GLPK.GLP_BS
    end
end

function glpk_tst_4()
    prev_term_out = GLPK.glp_term_out(GLPK.GLP_OFF)

    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    mip = GLPK.glp_create_prob()
    tran = GLPK.glp_mpl_alloc_wksp()
    GLPK.glp_mpl_read_model(tran, joinpath(datadir, "sudoku.mod"), 1)
    GLPK.glp_mpl_read_data(tran, joinpath(datadir, "sudoku.dat"))
    GLPK.glp_mpl_generate(tran, C_NULL)
    GLPK.glp_mpl_build_prob(tran, mip)
    @test GLPK.glp_simplex(mip, C_NULL) == 0
    @test GLPK.glp_intopt(mip, C_NULL) == 0
    GLPK.glp_mpl_postsolve(tran, mip, GLPK.GLP_MIP)

    @test GLPK.glp_warm_up(mip) == 0

    @test GLPK.glp_factorize(mip) == 0

    @test GLPK.glp_bf_updated(mip) == 0

    bfp = GLPK.glp_bfcp(
        0,
        0,
        0,
        0.0,
        0,
        0,
        0.0,
        0.0,
        0,
        0.0,
        0,
        0,
        tuple(fill(0.0, 38)...),
    )
    GLPK.glp_get_bfcp(mip, bfp)
    bfp.piv_tol = 0.5
    GLPK.glp_set_bfcp(mip, bfp)

    rows = GLPK.glp_get_num_rows(mip)
    cols = GLPK.glp_get_num_cols(mip)

    for i in 1:rows
        bh = GLPK.glp_get_bhead(mip, i)
        rb = GLPK.glp_get_row_bind(mip, i)
    end
    for i in 1:cols
        cb = GLPK.glp_get_col_bind(mip, i)
    end

    x = randn(rows)
    GLPK.glp_ftran(mip, x)
    GLPK.glp_btran(mip, x)

    GLPK.glp_adv_basis(mip, 0)
    @test GLPK.glp_factorize(mip) == 0

    for i in 1:rows+cols
        if !var_is_basic(mip, i)
            continue
        end

        ind = zeros(Cint, rows + cols)
        val = zeros(Cdouble, rows + cols)
        len = GLPK.glp_eval_tab_row(mip, i, GLPK.offset(ind), GLPK.offset(val))
        @test len <= length(ind)
    end

    for i in 1:rows+cols
        if var_is_basic(mip, i)
            continue
        end

        ind = zeros(Cint, rows + cols)
        val = zeros(Cdouble, rows + cols)
        len = GLPK.glp_eval_tab_col(mip, i, GLPK.offset(ind), GLPK.offset(val))
        @test len <= length(ind)
    end

    for i in 1:100
        ia, ja, val =
            SparseArrays.findnz(SparseArrays.sprand(Int(cols), 1, 0.5))
        indices = round.(Cint, ia)
        resize!(val, cols)
        resize!(indices, cols)
        GLPK.glp_transform_row(
            mip,
            length(ia),
            GLPK.offset(indices),
            GLPK.offset(val),
        )
    end

    for i in 1:100
        ia, ja, val =
            SparseArrays.findnz(SparseArrays.sprand(Int(rows), 1, 0.5))
        indices = round.(Cint, ia)
        resize!(val, rows)
        resize!(indices, rows)
        GLPK.glp_transform_col(
            mip,
            length(ia),
            GLPK.offset(indices),
            GLPK.offset(val),
        )
    end

    for i in 1:100
        ia, ja, val =
            SparseArrays.findnz(SparseArrays.sprand(Int(rows), 1, 0.5))
        dir = 2 * rand(Bool) - 1
        eps = 1e-9

        basic = map(x -> var_is_basic(mip, x), ia)

        ia = Cint.(ia[basic])
        val = val[basic]

        prt = GLPK.glp_prim_rtest(
            mip,
            length(ia),
            GLPK.offset(ia),
            GLPK.offset(val),
            dir,
            eps,
        )
    end

    for i in 1:100
        ia, ja, val =
            SparseArrays.findnz(SparseArrays.sprand(Int(cols), 1, 0.5))
        dir = 2 * rand(Bool) - 1
        eps = 1e-9

        nonbasic = map(x -> !var_is_basic(mip, x), ia)

        ia = Cint.(ia[nonbasic])
        val = val[nonbasic]

        drt = GLPK.glp_dual_rtest(
            mip,
            length(ia),
            GLPK.offset(ia),
            GLPK.offset(val),
            dir,
            eps,
        )
    end

    for i in 1:rows+cols
        if var_is_basic(mip, i)
            coef1 = Ref{Cdouble}()
            var1 = Ref{Cint}()
            value1 = Ref{Cdouble}()
            coef2 = Ref{Cdouble}()
            var2 = Ref{Cint}()
            value2 = Ref{Cdouble}()
            ac = GLPK.glp_analyze_coef(
                mip,
                i,
                coef1,
                var1,
                value1,
                coef2,
                var2,
                value2,
            )
        else
            limit1 = Ref{Cdouble}()
            var1 = Ref{Cint}()
            limit2 = Ref{Cdouble}()
            var2 = Ref{Cint}()
            ab = GLPK.glp_analyze_bound(mip, i, limit1, var1, limit2, var2)
        end
    end

    GLPK.glp_init_env()

    p = GLPK.glp_alloc(1, 100 * sizeof(Int))
    count, cpeak, total, tpeak =
        Ref{Cint}(), Ref{Cint}(), Ref{Cint}(), Ref{Cint}()
    GLPK.glp_mem_usage(count, cpeak, total, tpeak)
    GLPK.glp_free(p)
    mu = GLPK.glp_mem_usage(count, cpeak, total, tpeak)
    GLPK.glp_mem_limit(10000)

    p = GLPK.glp_alloc(100, sizeof(Int))
    mu = GLPK.glp_mem_usage(count, cpeak, total, tpeak)
    GLPK.glp_free(p)
    mu = GLPK.glp_mem_usage(count, cpeak, total, tpeak)

    # GLPK.glp_free_env()

    mktemp() do path, _
        return GLPK.glp_print_sol(mip, path)
    end

    mktemp() do path, _
        GLPK.glp_write_sol(mip, path)
        return GLPK.glp_read_sol(mip, path)
    end

    mktemp() do path, _
        return GLPK.glp_print_ipt(mip, path)
    end

    mktemp() do path, _
        GLPK.glp_write_ipt(mip, path)
        return GLPK.glp_read_ipt(mip, path)
    end

    mktemp() do path, _
        return GLPK.glp_print_mip(mip, path)
    end

    mktemp() do path, _
        GLPK.glp_write_mip(mip, path)
        return GLPK.glp_read_mip(mip, path)
    end

    mktemp() do path, _
        GLPK.glp_open_tee(path)
        return GLPK.glp_close_tee()
    end

    GLPK.glp_term_out(prev_term_out)
    GLPK.glp_mpl_free_wksp(tran)
    return GLPK.glp_delete_prob(mip)
end

glpk_tst_4()
