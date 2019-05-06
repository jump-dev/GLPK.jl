using Test
import GLPK

function cb_callback(tree::Ptr{Cvoid}, info::Ptr{Cvoid})
    reason = GLPK.ios_reason(tree)
    @test reason in [GLPK.ISELECT, GLPK.IPREPRO, GLPK.IROWGEN, GLPK.IHEUR,
                     GLPK.ICUTGEN, GLPK.IBRANCH, GLPK.IBINGO]

    prob = GLPK.ios_get_prob(tree)

    a_cnt, n_cnt, t_cnt = GLPK.ios_tree_size(tree)
    cn = GLPK.ios_curr_node(tree)
    bn = GLPK.ios_best_node(tree)
    fan = GLPK.ios_next_node(tree, 0)
    GLPK.ios_next_node(tree, fan)
    lan = GLPK.ios_prev_node(tree, 0)
    GLPK.ios_prev_node(tree, lan)
    for i = 1:t_cnt
        try
            un = GLPK.ios_up_node(tree, i)
            nl = GLPK.ios_node_level(tree, i)
            nb = GLPK.ios_node_bound(tree, i)
            nd = GLPK.ios_node_data(tree, i)
        catch
        end
    end

    for i = 1:GLPK.get_num_rows(prob)
        attr = GLPK.ios_row_attr(tree, i)
    end
    for i = 1:GLPK.get_num_cols(prob)
        cb = GLPK.ios_can_branch(tree, i)
    end

    mg = GLPK.ios_mip_gap(tree)
    if reason == GLPK.ISELECT
        @test cn == 0
        GLPK.ios_select_node(tree, fan)
    elseif reason == GLPK.IPREPRO
    elseif reason == GLPK.IROWGEN
    elseif reason == GLPK.IHEUR
        #GLPK.ios_heur_sol(tree, zeros(Int, GLPK.get_num_cols(GLPK.ios_get_prob(tree))))
    elseif reason == GLPK.ICUTGEN
        ps = GLPK.ios_pool_size(tree)
        GLPK.ios_add_row(tree, nothing, 0, 0, 0, Int32[], Float64[], GLPK.LO, 0.0)
        @test GLPK.ios_pool_size(tree) == ps + 1
        GLPK.ios_add_row(tree, nothing, 0, Int32[1,2,3], Float64[7.0,6.0,5.0], GLPK.LO, 0.0)
        @test GLPK.ios_pool_size(tree) == ps + 2
        GLPK.ios_del_row(tree, 1)
        @test GLPK.ios_pool_size(tree) == ps + 1
        GLPK.ios_clear_pool(tree)
        @test GLPK.ios_pool_size(tree) == 0
    elseif reason == GLPK.IBRANCH
        for i = 1:GLPK.get_num_cols(prob)
            if GLPK.ios_can_branch(tree, i) != 0
                GLPK.ios_branch_upon(tree, i, GLPK.NO_BRNCH)
                break
            end
        end
    elseif reason == GLPK.IBINGO
    end

    return nothing
end

function glpk_tst_6()
    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    prev_term_out = GLPK.term_out(GLPK.OFF)

    mip = GLPK.Prob()
    tran = GLPK.MathProgWorkspace()
    GLPK.mpl_read_model(tran, joinpath(datadir, "sudoku.mod"), 1)
    GLPK.mpl_read_data(tran, joinpath(datadir, "sudoku.dat"))
    GLPK.mpl_generate(tran, nothing)
    GLPK.mpl_build_prob(tran, mip)
    #@test GLPK.simplex(mip, nothing) == 0

    params = GLPK.IntoptParam()
    params.presolve = GLPK.ON
    params.cb_size = 1

    params.cb_func = @cfunction(cb_callback, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))

    @test GLPK.intopt(mip, params) == 0

    GLPK.mpl_postsolve(tran, mip, GLPK.MIP)

    GLPK.term_out(prev_term_out)
end

glpk_tst_6()
