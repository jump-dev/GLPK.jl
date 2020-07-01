using Test
import GLPK

# Test reading and writing MathProg files

function glpk_tst_3()
    prev_term_out = GLPK.glp_term_out(GLPK.GLP_OFF)

    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    lp = GLPK.glp_create_prob()
    tran = GLPK.glp_mpl_alloc_wksp()

    GLPK.glp_mpl_read_model(tran, joinpath(datadir, "queens.mod"), 0)

    outfile = tempname()
    try
        GLPK.glp_mpl_generate(tran, outfile)
        GLPK.glp_mpl_build_prob(tran, lp)
        GLPK.glp_mpl_free_wksp(tran)
    finally
        isfile(outfile) && rm(outfile)
    end

    @test GLPK.glp_simplex(lp, C_NULL) == 0

    GLPK.glp_term_out(prev_term_out)
    GLPK.glp_delete_prob(lp)
end

glpk_tst_3()
