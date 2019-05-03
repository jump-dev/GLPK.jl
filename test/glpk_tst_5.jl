using Test
import GLPK

function glpk_tst_5()
    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    prev_term_out = GLPK.term_out(GLPK.OFF)

    lp = GLPK.Prob()

    cnffile = joinpath(datadir, "sat_tst.cnf")
    GLPK.read_cnfsat(lp, cnffile)
    GLPK.check_cnfsat(lp)
    cnfcopy = joinpath(datadir, "sat_tst-copy.cnf")
    try
        GLPK.write_cnfsat(lp, cnfcopy)
    finally
        isfile(cnfcopy) && rm(cnfcopy)
    end
    if (4, 57) <= GLPK.version() <= (4, 61) && Sys.WORD_SIZE != sizeof(Cint) * 8
        # in api/minisat1.c, there is:
        #if (sizeof(void *) != sizeof(int))
        #{  xprintf("glp_minisat1: sorry, MiniSat solver is not supported "
        #      "on 64-bit platforms\n");
        #   ret = GLP_EFAIL;
        #   goto done;
        #}
        expected_ret = GLPK.EFAIL

    elseif (4, 62) <= GLPK.version() && Sys.WORD_SIZE != sizeof(Csize_t) * 8
        # changed condition to:
        # if (sizeof(void *) != sizeof(size_t))
        expected_ret = GLPK.EFAIL
    else
        expected_ret = 0
    end
    @test GLPK.minisat1(lp) == expected_ret

    @test GLPK.intfeas1(lp, 0, 0) == expected_ret
end

glpk_tst_5()
