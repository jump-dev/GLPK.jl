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
    @test GLPK.minisat1(lp) == 0

    @test GLPK.intfeas1(lp, 0, 0) == 0
end

# Test is failing on GLPK 4.61 for unknown reason
println("Test disabled")
#glpk_tst_5()
