using Test
import GLPK

function glpk_tst_5()
    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    prev_term_out = GLPK.glp_term_out(GLPK.GLP_OFF)

    lp = GLPK.glp_create_prob()

    cnffile = joinpath(datadir, "sat_tst.cnf")
    GLPK.glp_read_cnfsat(lp, cnffile)
    GLPK.glp_check_cnfsat(lp)
    cnfcopy = joinpath(datadir, "sat_tst-copy.cnf")
    try
        GLPK.glp_write_cnfsat(lp, cnfcopy)
    finally
        isfile(cnfcopy) && rm(cnfcopy)
    end
    if Sys.WORD_SIZE != sizeof(Csize_t) * 8
        # changed condition to:
        # if (sizeof(void *) != sizeof(size_t))
        expected_ret = GLPK.GLP_EFAIL
    else
        expected_ret = 0
    end
    @test GLPK.glp_minisat1(lp) == expected_ret

    @test GLPK.glp_intfeas1(lp, 0, 0) == expected_ret
    GLPK.glp_delete_prob(lp)
end

glpk_tst_5()
