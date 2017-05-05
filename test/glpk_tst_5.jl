using Base.Test
import GLPK

function glpk_tst_5()
    datadir = joinpath(dirname(@__FILE__), "data")
    isdir(datadir) || (datadir = joinpath(Pkg.dir(), "GLPK.jl", "test", "data"))

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
    if GLPK.version() >= (4, 57)
        # in api/minisat1.c, there is:
        #if (sizeof(void *) != sizeof(int))
        #{  xprintf("glp_minisat1: sorry, MiniSat solver is not supported "
        #      "on 64-bit platforms\n");
        #   ret = GLP_EFAIL;
        #   goto done;
        #}
        expected_ret = GLPK.EFAIL
    else
        expected_ret = 0
    end
    @test GLPK.minisat1(lp) == expected_ret

    @test GLPK.intfeas1(lp, 0, 0) == expected_ret
end

glpk_tst_5()
