using Test
import GLPK

function cb_callback(tree::Ptr{Cvoid}, info::Ptr{Cvoid})
    reason = GLPK.glp_ios_reason(tree)
    @test reason in [
        GLPK.GLP_ISELECT,
        GLPK.GLP_IPREPRO,
        GLPK.GLP_IROWGEN,
        GLPK.GLP_IHEUR,
        GLPK.GLP_ICUTGEN,
        GLPK.GLP_IBRANCH,
        GLPK.GLP_IBINGO,
    ]

    prob = GLPK.glp_ios_get_prob(tree)

    a_cnt, n_cnt, t_cnt = Ref{Cint}(), Ref{Cint}(), Ref{Cint}()
    GLPK.glp_ios_tree_size(tree, a_cnt, n_cnt, t_cnt)
    cn = GLPK.glp_ios_curr_node(tree)
    bn = GLPK.glp_ios_best_node(tree)
    fan = GLPK.glp_ios_next_node(tree, 0)
    GLPK.glp_ios_next_node(tree, fan)
    lan = GLPK.glp_ios_prev_node(tree, 0)
    GLPK.glp_ios_prev_node(tree, lan)
    for i = 1:t_cnt[]
        try
            un = GLPK.glp_ios_up_node(tree, i)
            nl = GLPK.glp_ios_node_level(tree, i)
            nb = GLPK.glp_ios_node_bound(tree, i)
            nd = GLPK.glp_ios_node_data(tree, i)
        catch
        end
    end

    for i = 1:GLPK.glp_get_num_rows(prob)
        attr = GLPK.glp_attr(0, 0, 0, tuple(fill(0.0, 7)...))
        GLPK.glp_ios_row_attr(tree, i, attr)
    end
    for i = 1:GLPK.glp_get_num_cols(prob)
        cb = GLPK.glp_ios_can_branch(tree, i)
    end

    mg = GLPK.glp_ios_mip_gap(tree)
    if reason == GLPK.GLP_ISELECT
        @test cn == 0
        GLPK.glp_ios_select_node(tree, fan)
    elseif reason == GLPK.GLP_IPREPRO
    elseif reason == GLPK.GLP_IROWGEN
    elseif reason == GLPK.GLP_IHEUR
        #GLPK.ios_heur_sol(tree, zeros(Int, GLPK.get_num_cols(GLPK.ios_get_prob(tree))))
    elseif reason == GLPK.GLP_ICUTGEN
        ps = GLPK.glp_ios_pool_size(tree)
        indices, coefficients = Cint[], Cdouble[]
        GLPK.glp_ios_add_row(
            tree,
            "",
            0,
            0,
            Cint(length(indices)),
            pointer(indices) - sizeof(Cint),
            pointer(coefficients) - sizeof(Cdouble),
            GLPK.GLP_LO,
            0.0,
        )
        @test GLPK.glp_ios_pool_size(tree) == ps + 1
        indices, coefficients = Cint[1, 2, 3], Cdouble[7, 6, 5]
        GLPK.glp_ios_add_row(
            tree,
            "",
            0,
            0,
            Cint(length(indices)),
            pointer(indices) - sizeof(Cint),
            pointer(coefficients) - sizeof(Cdouble),
            GLPK.GLP_LO,
            0.0,
        )
        @test GLPK.glp_ios_pool_size(tree) == ps + 2
        GLPK.glp_ios_del_row(tree, 1)
        @test GLPK.glp_ios_pool_size(tree) == ps + 1
        GLPK.glp_ios_clear_pool(tree)
        @test GLPK.glp_ios_pool_size(tree) == 0
    elseif reason == GLPK.GLP_IBRANCH
        for i = 1:GLPK.glp_get_num_cols(prob)
            if GLPK.glp_ios_can_branch(tree, i) != 0
                GLPK.glp_ios_branch_upon(tree, i, GLPK.GLP_NO_BRNCH)
                break
            end
        end
    elseif reason == GLPK.GLP_IBINGO
    end

    return nothing
end

function glpk_tst_6()
    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    prev_term_out = GLPK.glp_term_out(GLPK.GLP_OFF)

    mip = GLPK.glp_create_prob()
    tran = GLPK.glp_mpl_alloc_wksp()
    GLPK.glp_mpl_read_model(tran, joinpath(datadir, "sudoku.mod"), 1)
    GLPK.glp_mpl_read_data(tran, joinpath(datadir, "sudoku.dat"))
    GLPK.glp_mpl_generate(tran, C_NULL)
    GLPK.glp_mpl_build_prob(tran, mip)
    #@test GLPK.glp_simplex(mip, nothing) == 0

    params = GLPK.glp_iocp()
    GLPK.glp_init_iocp(params)
    params.presolve = GLPK.GLP_ON
    params.cb_size = 1

    params.cb_func = @cfunction(cb_callback, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))

    @test GLPK.glp_intopt(mip, params) == 0

    GLPK.glp_mpl_postsolve(tran, mip, GLPK.GLP_MIP)

    GLPK.glp_term_out(prev_term_out)
    GLPK.glp_mpl_free_wksp(tran)
    GLPK.glp_delete_prob(mip)

end

glpk_tst_6()
